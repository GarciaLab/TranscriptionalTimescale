% script to simulate burst time series with parameters from numerical and
% analytic results
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
FigurePath = ['../fig/bursting_supp/' project '/'];
mkdir(FigurePath)
DataPath = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'bursting_sim_struct.mat'])
load([DataPath 'bursting_chain_calc_struct.mat'])

% set basic plot parameters
t_max = 60;
ylimTrace = [-0.5 6.5];
n_bound_vec = 0:n_bcd_sites;

% sim name cell
sim_name_cell = {bursting_sim_struct.name};

% define colors
blue = [190 201 224]/255;
red = [246 141 100]/255;
green = [203 220 170]/255;
gray = [0.7020    0.7020    0.7020];
cmap1 = [green ; blue ;red];

% define resampling time res. Slower sampling time is neeeded to make plots
% intelligible
resamp_res = 0.5; % in seconds

% define time grid for resampling
time_rs = 0:resamp_res:3600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% koff-mediated cooperative binding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify appropriate index
coop_sim_index = find(contains(sim_name_cell,'koff-mediated'));

coop_plot_index = 1; 
sim_param_indices_coop = 176:4:201;
% plot results of stochastic simulations
trace_index = 3;%3

state_fig = figure;
% cmap2 = brewermap(9,'Set2');
hold on

% generate resampled trace (moving average, essentially)
time_raw = repelem(double(bursting_sim_struct(coop_sim_index).sim_time_cell{coop_plot_index,trace_index}),1);
trace_raw = repelem(double(bursting_sim_struct(coop_sim_index).sim_emission_cell{coop_plot_index,trace_index}),1);

trace_rs = interp1(time_raw,trace_raw,time_rs,'previous');

stairs(time_rs/60, trace_rs,'Color',[blue 0.0],'LineWidth',1.5);

ylim(ylimTrace)
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',n_bound_vec)
p = plot(0,0);
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';

saveas(state_fig,[FigurePath 'coop_trace.png'])
saveas(state_fig,[FigurePath 'coop_trace.pdf'])

% calculate fraction of time in eac state
dur_vec = diff([time_raw t_max*60]);
stateSums = accumarray(trace_raw'+1,dur_vec');
stateShares = stateSums/sum(stateSums);

hist_fig = figure('Position',[100 100 256 512]);
hold on
barh(n_bound_vec,stateShares,1,'FaceColor',blue);
  
xlabel('probability')
box on
p = plot(0,0);
% ax = gca;
% ax.YColor = 'black';
xlim([0 0.55])
% ax.XColor = 'black';
set(gca,'Fontsize',14,'xtick',0:.25:.5)
ylim([n_bound_vec(1)-0.5 n_bound_vec(end)+0.5])
StandardFigurePBoC(p,gca);

hist_fig.InvertHardcopy = 'off';
saveas(hist_fig,[FigurePath 'coop_hist.png'])
saveas(hist_fig,[FigurePath 'coop_hist.pdf'])

% %%
% % use mean dwell time inferred from exponential fit (see supp02)
% % mu = 0.23;
% upper_bound = 1;
% del_indices = find([diff(time_raw) 0]>=upper_bound & [diff(trace_raw) 0]==-1);
% time_raw_trunc = time_raw;
% for d = 1:length(del_indices)
%   deltaT = time_raw_trunc(del_indices(d)+1)-time_raw_trunc(del_indices(d))-upper_bound;
%   time_raw_trunc(del_indices(d)+1:end) = time_raw_trunc(del_indices(d)+1:end)-deltaT;
% end
% 
% state_fig2 = figure;
% cmap2 = brewermap(9,'Set2');
% hold on
% time_rs2 = 0:resamp_res:max(time_raw_trunc);
% trace_rs2 = interp1(time_raw_trunc,trace_raw,time_rs2,'previous');
% 
% stairs(time_rs2/60, trace_rs2,'Color',[blue 0.0],'LineWidth',1.5);
% 
% ylim(ylimTrace)
% xlim([0 max(time_rs2)/60])
% ylabel('transcription rate')
% xlabel('time (minutes)')
% box on
% set(gca,'Fontsize',14,'YTick',n_bound_vec)
% p = plot(0,0);
% StandardFigurePBoC(p,gca);
% state_fig2.InvertHardcopy = 'off';
% 
% saveas(state_fig2,[FigurePath 'coop_trace_trunc.png'])
% saveas(state_fig2,[FigurePath 'coop_trace_trunc.pdf'])
% 
% % calculate fraction of time in eac state
% dur_vec = diff([time_raw t_max*60]);
% stateSums = accumarray(trace_raw'+1,dur_vec');
% stateShares = stateSums/sum(stateSums);
% 
% hist_fig = figure('Position',[100 100 256 512]);
% hold on
% barh(n_bound_vec,stateShares,1,'FaceColor',blue);
%   
% xlabel('probability')
% box on
% p = plot(0,0);
% % ax = gca;
% % ax.YColor = 'black';
% xlim([0 0.55])
% % ax.XColor = 'black';
% set(gca,'Fontsize',14,'xtick',0:.25:.5)
% ylim([n_bound_vec(1)-0.5 n_bound_vec(end)+0.5])
% StandardFigurePBoC(p,gca);
% 
% hist_fig.InvertHardcopy = 'off';
% saveas(hist_fig,[FigurePath 'coop_hist.png'])
% saveas(hist_fig,[FigurePath 'coop_hist.pdf'])

%% make bar plots of effective on and off rates


% define matrices to use for
n = 0:n_bcd_sites;
n_states = length(n);
a = ones(n_states); m1 = tril(a,-1); m2 = tril(a,-2); m3 = triu(a,1); m4 = triu(a,2); m5 = ~~eye(n_states);
% extract transition rate matrix
Q = bursting_sim_struct(coop_sim_index).Q(:,:,coop_plot_index)';
k_minus_vec = [0 Q(m3&~m4)'];
k_plus_vec = [Q(m1&~m2)' 0];

% make on rate bar plot
kp_fig = figure;
bar(n,k_plus_vec,.3,'b')
xlabel('number of bound molecules')
ylabel('k_{+} (1/s)')
%ylabel('k_{+}')
set(gca,'YTick',[0.01 1 10^2])
ylim([0.01,10^2])
set(gcf,'position',[0, 0, 0.5600, 0.1560]*1E3)
set(gca,'YScale','log')
StandardFigure([],gca)
kp_fig.InvertHardcopy = 'off';
saveas(kp_fig,[FigurePath 'kp_bars.png'])
saveas(kp_fig,[FigurePath 'kp_bars.pdf'])


km_fig = figure;
bar(n,k_minus_vec,.3,'b')
xlabel('number of bound molecules')
ylabel('k_{-} (1/s)')
%ylabel('k_{+}')
set(gca,'YTick',[0.01 1 10^2])
ylim([0.01,10^2])
set(gcf,'position',[0,0, 0.5600, 0.1560]*1E3)
set(gca,'YScale','log')
StandardFigure([],gca)
km_fig.InvertHardcopy = 'off';
saveas(km_fig,[FigurePath 'km_bars.png'])
saveas(km_fig,[FigurePath 'km_bars.pdf'])



%%%%%%%%%%%%%%%%%%%
sim_param_indices_coop = 176:4:201;
kon = Q(7,6)
koff = Q(1,2)

load([DataPath 'bursting_chain_calc_struct.mat'])
omega_vec = exp(-bursting_chain_calc_struct(coop_sim_index).coopEnergies);
omega = omega_vec(sim_param_indices_coop(coop_plot_index))