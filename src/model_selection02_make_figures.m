% script to make figures examine waiting time distributions for 
% rate-limitiung step and and emergent cooperativity simulations
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
FigurePath = ['../fig/waiting_time_distributions/' project '/'];
mkdir(FigurePath)
DataPath = ['../out/waiting_time_distributions/' project '/'];
DataPathCalc = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'waiting_time_struct.mat'])
load([DataPathCalc 'bursting_chain_calc_struct.mat'])
% set basic plot parameters
t_max = 60;
ylimTrace = [-0.5 8];
n_bound_vec = 0:n_bcd_sites;

% sim name cell
sim_name_cell = {waiting_time_struct.name};

% define colors
purple = brighten([171 133 172]/255,.5);
blue = [190 201 224]/255;
red = [246 141 100]/255;
green = [203 220 170]/255;
gray = [0.7020    0.7020    0.7020];
% cmap1 = [green ; blue ;red];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Make figure illustrating passage time concept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify appropriate index
rateLim_sim_index = find(contains(sim_name_cell,'2rate-limiting'));
%
%  extract n bound vec
rateLim_sim_sub_index = 6;%length(waiting_time_struct(rateLim_sim_index).off_waiting_times_ideal);
% plot results of stochastic simulations
trace_index = 33;

state_fig = figure;%('Position',[100 100 1024 512]);
hold on

% extract trace data
% time_raw = double(bursting_sim_struct(rateLim_sim_index).sim_time_cell{rateLim_plot_index,trace_index});
% trace_raw = double(bursting_sim_struct(rateLim_sim_index).sim_emission_cell{rateLim_plot_index,trace_index});
trace_raw = waiting_time_struct(rateLim_sim_index).trace_array(trace_index,:,rateLim_sim_sub_index)-1;

% extract corresponding viterbi fit
viterbi_time = waiting_time_struct(rateLim_sim_index).time_vector;
viterbi_fit = waiting_time_struct(rateLim_sim_index).viterbi_traces(trace_index,:,rateLim_sim_sub_index)*n_bcd_sites;

stairs(viterbi_time/60, trace_raw,'Color',purple,'LineWidth',1);
stairs(viterbi_time/60, viterbi_fit,'Color','k','LineWidth',1.5);

ylim(ylimTrace)
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
% legend('raw trace','2 state fit')
set(gca,'Fontsize',14,'YTick',0:8)
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';
saveas(state_fig,[FigurePath 'rateLim_trace.png'])
saveas(state_fig,[FigurePath 'rateLim_trace.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) Make figures showing passage times for rate-limiting step mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_step_vec = [1 2 5 15];
rateLim_sim_indices = find(ismember(sim_name_cell,{'1rate-limiting steps','2rate-limiting steps','5rate-limiting steps','15rate-limiting steps'}));
close all

plot_index = 6;%length(waiting_time_struct(rateLim_sim_indices(1)).off_waiting_times_ideal);
pt_mean = mean(waiting_time_struct(rateLim_sim_indices(1)).off_waiting_times_ideal{plot_index})/60;
% define bins for grouping waiting time measurements
wt_bins = linspace(0,4,50);
wt_centers = (wt_bins(1:end-1)+wt_bins(2:end))/2;

% define gamma function
gamma_fun = @(x) x(1) * x(2)^x(3) .* wt_centers.^(x(3)-1).*exp(-x(2)*wt_centers) ./ gamma(x(3));
exp_fun = @(x) x(1) * exp(-x(2) .* wt_centers);
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',5e4);

% perform fits
rl_fit_struct = struct;
for i = 1:length(rateLim_sim_indices)
  wt_vec_raw = waiting_time_struct(rateLim_sim_indices(i)).off_waiting_times_ideal{plot_index};
  wt_vec = wt_vec_raw  / pt_mean / 60;% mean(wt_vec_raw);
  p_vec = histcounts(wt_vec,wt_bins);
  if i == 1
    fit_obj = @(x) exp_fun(x)-p_vec;
    rl_fit_struct(i).params = lsqnonlin(fit_obj,[1 1], [0 0],[Inf Inf],options);
    pd_vec = exp_fun(rl_fit_struct(i).params);
  else
    fit_obj = @(x) gamma_fun(x)-p_vec;
    rl_fit_struct(i).params = lsqnonlin(fit_obj,[1 1 1], [0 0 0],[Inf Inf Inf],options);
    pd_vec = gamma_fun(rl_fit_struct(i).params);
  end  
%   rl_fit_struct(i).wt_axis = wt_centers * mean(wt_vec_raw);
  rl_fit_struct(i).pd_vec = pd_vec /sum(pd_vec);
  rl_fit_struct(i).p_vec = p_vec / sum(p_vec);
  rl_fit_struct(i).var = var(wt_vec_raw);
  rl_fit_struct(i).mean = mean(wt_vec_raw);
end


offset = 2;
rl_hist_fig = figure;
cmap1_raw = brewermap([],'Pastel2');
cmap2_raw = brewermap([],'Set2');
start1 = brighten(cmap1_raw(2,:),-0.8);
end1 = 0.3*brighten(cmap1_raw(2,:),0.8)+[.7 .7 .7];
start2 = brighten(cmap2_raw(2,:),-.8);
end2 = .3*brighten(cmap2_raw(2,:),.8)+[.7 .7 .7];


cmap1 = interp1([1,length(rateLim_sim_indices)],[start1 ; end1],linspace(1,length(rateLim_sim_indices),length(rateLim_sim_indices)+offset));
cmap2 = interp1([1,length(rateLim_sim_indices)],[start2 ; end2],linspace(1,length(rateLim_sim_indices),length(rateLim_sim_indices)+offset));

hold on
iter = length(rateLim_sim_indices);
pl = [];
lgd_str = {};
for i = fliplr(1:length(rateLim_sim_indices))
  p_vec = rl_fit_struct(i).p_vec;
%   wt_axis = rl_fit_struct(i).wt_axis/60;
  if i == length(rateLim_sim_indices)
    bar(wt_centers*pt_mean,p_vec,1,'FaceColor',cmap1(i+offset/2,:),'FaceAlpha',1,'EdgeAlpha',0.3);   
  elseif i > 1
    bar(wt_centers*pt_mean,p_vec,1,'FaceColor',cmap1(i+offset/2,:),'FaceAlpha',.5,'EdgeAlpha',0.3);  
  else
    bar(wt_centers*pt_mean,p_vec,1,'FaceColor',cmap1(i+offset/2,:),'FaceAlpha',0.5,'EdgeAlpha',0.3);  
  end
  pd_vec = rl_fit_struct(i).pd_vec;
  pl = [pl plot(wt_centers*pt_mean,pd_vec,'Color',cmap2(i+offset/2,:),'LineWidth',2)]; 
  if i == 1
    lgd_str = [lgd_str{:} {[num2str(n_step_vec(i)) ' step']}];
  else
    lgd_str = [lgd_str{:} {[num2str(n_step_vec(i)) ' steps']}];
  end
end
% for i = 1:length(rateLim_sim_indices)
%   pd_vec = rl_fit_struct(i).pd_vec;
%   plot(wt_centers,pd_vec,'Color',cmap(i,:),'LineWidth',2);  
% end

xlim([0 wt_bins(end)*pt_mean])
ylim([0 .13])
set(gca,'Fontsize',14)
ylabel('probability')
xlabel('first-pasage time (minutes)')
p = plot(0,0);
legend(pl,lgd_str{:});
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
rl_hist_fig.InvertHardcopy = 'off';
saveas(rl_hist_fig,[FigurePath 'rateLim_wt_hist.png'])
saveas(rl_hist_fig,[FigurePath 'rateLim_wt_hist.pdf'])


coop_sim_index = find(ismember(sim_name_cell,{'kon-mediated cooperativity'}));
sim_vec = 1:length(waiting_time_struct(i).off_waiting_times_ideal);
% plot_index = length(waiting_time_struct(i).off_waiting_times_ideal);

blueDark = brighten(blue,-0.5);

% perform fits
coop_fit_struct = struct;
for i = sim_vec
  wt_vec_raw = waiting_time_struct(coop_sim_index).off_waiting_times_ideal{i};
  wt_vec = wt_vec_raw / pt_mean / 60;    
  p_vec = histcounts(wt_vec,wt_bins);
  
  fit_obj = @(x) exp_fun(x)-p_vec;
  coop_fit_struct(i).params = lsqnonlin(fit_obj,[1 1], [0 0],[Inf Inf],options);
  pd_vec = exp_fun(coop_fit_struct(i).params);

  coop_fit_struct(i).pd_vec = pd_vec/sum(pd_vec);
  coop_fit_struct(i).p_vec = p_vec / sum(p_vec);
  coop_fit_struct(i).var = var(wt_vec_raw);
  coop_fit_struct(i).mean = mean(wt_vec_raw);
end


coop_hist_fig = figure;
offset = length(rateLim_sim_indices);
hold on

p_vec = coop_fit_struct(plot_index).p_vec;  
bar(wt_centers*pt_mean,p_vec,1,'FaceColor',blue,'FaceAlpha',.7,'EdgeAlpha',0.3);   
pd_vec = coop_fit_struct(plot_index).pd_vec;
plot(wt_centers*pt_mean,pd_vec,'Color',blueDark,'LineWidth',2);  

xlim([0 wt_bins(end)*pt_mean])
set(gca,'Fontsize',14)
ylabel('probability')
xlabel('first-pasage time (minutes)')
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
coop_hist_fig.InvertHardcopy = 'off';
saveas(coop_hist_fig,[FigurePath 'coop_wt_hist.png'])
saveas(coop_hist_fig,[FigurePath 'coop_wt_hist.pdf'])


xmax = 11;
ymax = 11;
n_boots = 10;
% make fano factor figure
% marker_cell = {'o','o','^','d','s'};
marker_cell = {'o','o','o','o','o'};
plot_indices = [1 rateLim_sim_indices];
plot_colors = [blueDark ; cmap2(2:end-1,:)];

fano_struct = struct;

fano_fig = figure;
hold on
% plot CV=1
plot(linspace(0,max([xmax,ymax])),linspace(0,max([xmax,ymax])),'--','Color','k','LineWidth',1);
% plot CV=2
plot(linspace(0,max([xmax,ymax])),2*linspace(0,max([xmax,ymax])),'--','Color','k','LineWidth',1);
% plot CV=1/2
plot(linspace(0,2*max([xmax,ymax])),linspace(0,max([xmax,ymax])),'--','Color','k','LineWidth',1);

iter = 1;
for p = plot_indices
  wt_cell = waiting_time_struct(p).off_waiting_times_ideal;
  mean_array = NaN(n_boots,size(wt_cell,2));
  std_array = NaN(n_boots,size(wt_cell,2));
  for w = 1:length(wt_cell)
    wt_times = wt_cell{w};
    for n = 1:n_boots
      boot_indices = randsample(1:length(wt_times),length(wt_times),true);
      mean_array(n,w) = mean(wt_times(boot_indices));
      std_array(n,w) = std(wt_times(boot_indices));
    end
  end
  fano_struct(iter).mean_vec = nanmean(mean_array)/60;
  fano_struct(iter).mean_se_vec = nanstd(mean_array)/60;
  fano_struct(iter).std_vec = nanmean(std_array)/60;
  fano_struct(iter).std_se_vec = nanstd(std_array)/60;
    
  iter = iter + 1;
end
% for f = 1:length(fano_struct)
%   yse = fano_struct(f).std_se_vec;
%   xse = fano_struct(f).mean_se_vec;
%   errorbar(fano_struct(f).mean_vec,fano_struct(f).std_vec,yse,yse,xse,xse,'.','CapSize',0,'Color','k')
% end
for f = 1:length(fano_struct)
    % plot
  if f~=1
    scatter(fano_struct(f).mean_vec,fano_struct(f).std_vec,50,'s',marker_cell{f},'MarkerFaceColor',...
      plot_colors(f,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1)
  else
    scatter(fano_struct(f).mean_vec,fano_struct(f).std_vec,75,'o',marker_cell{f},'MarkerFaceColor',...
      plot_colors(f,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1)
  end
end

xlim([0 xmax])
ylim([0 ymax])
% grid on
set(gca,'Fontsize',14)
ylabel('standard deviation')
xlabel('mean first passage time (minutes)')
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
fano_fig.InvertHardcopy = 'off';
saveas(fano_fig,[FigurePath 'fano_scatter.png'])
saveas(fano_fig,[FigurePath 'fano_scatter.pdf'])