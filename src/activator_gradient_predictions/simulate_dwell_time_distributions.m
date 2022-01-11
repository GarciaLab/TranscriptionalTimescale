% script to simulate full binding/unbinding model to check average
% unbinding rate for on-rate-mediated cooperativity model
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('../utilities')

% set paths
FigurePath = ['.'];
mkdir(FigurePath)
DataPath = ['../../out/emergent_bursting/' project '/'];

% load data
% load([DataPath 'bursting_sim_struct.mat'])
load([DataPath 'bursting_chain_calc_struct.mat'])


% sim name cell
sim_name_cell = {bursting_chain_calc_struct.name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% koff-mediated cooperative binding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify appropriate index
coop_on_sim_index = find(contains(sim_name_cell,'kon-mediated'));
coop_off_sim_index = find(contains(sim_name_cell,'koff-mediated'));

coop_plot_index = 1; 
sim_param_indices_coop = 176:4:201;

% extract transition rate matrix
Q_on = bursting_chain_calc_struct(coop_on_sim_index).Q(:,:,sim_param_indices_coop(coop_plot_index))';
SS_on = bursting_chain_calc_struct(coop_on_sim_index).SS(:,sim_param_indices_coop(coop_plot_index))';

Q_off = bursting_chain_calc_struct(coop_off_sim_index).Q(:,:,sim_param_indices_coop(coop_plot_index))';
SS_off = bursting_chain_calc_struct(coop_off_sim_index).SS(:,sim_param_indices_coop(coop_plot_index))';

%% %%%%%%%%%% Conduct simulations
n_sim = 100; % number of simulations
T = 1e4; % total time to simulate in seconds

% %%%%%%%%%%%%%%%%%%%% build model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(123);
tic
[unbinding_vec_off, binding_vec_off] = microscopic_binding_sim(Q_off,SS_off,100,T);
toc
[unbinding_vec_on, binding_vec_on] = microscopic_binding_sim(Q_on,SS_on,25,T);

%% %%%%%%%%%%%% fit exponential to simulated dwell times and plot %%%%%%%%%
time_step = 0.5;
ub_bins = 0:time_step:1e3;
ub_centers = ub_bins(1:end-1)+diff(ub_bins)/2;

c_vec_off = histcounts(unbinding_vec_off,ub_bins);
c_vec_on = histcounts(unbinding_vec_on,ub_bins);
p_vec_off = c_vec_off/sum(c_vec_off);
p_vec_on = c_vec_on/sum(c_vec_on);

% fit raw counts
% exp_fun = @(x) x(1) * exp(-ub_centers./x(2));
% fit_obj = @(x) exp_fun(x)-c_vec_off;
% 
% options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',5e4);
f_off = fit(ub_centers',p_vec_off','smoothingspline','SmoothingParam',0.5);
p_trend_off = feval(f_off,ub_centers);
f_on = fit(ub_centers',p_vec_on','smoothingspline','SmoothingParam',0.5);
p_trend_on = feval(f_on,ub_centers);

% Make plots
close all

%% make exponential fit figure first
exp_fig_on = figure;
hold on
% bar(ub_centers,p_vec_off,1,'FaceColor',cmap(3,:),'FaceAlpha',.9,'EdgeAlpha',1,'EdgeColor','k');   

% pd_vec = params(1)*exp(-[0.01 ub_centers]./params(2))/sum(c_vec_off);
plot(ub_centers,p_vec_on,'Color','k','LineWidth',3);
xlim([0 15])
% legend('independent binding','cooperative binding')
xlabel('activator dwell times (s)')
ylabel('probability')
box on
% StandardFigurePBoC([],gca);
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'Fontsize',14)

% set(gca,'Color','w')
exp_fig_on.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(exp_fig_on,['activator_dwell_exp_on_coop_6bs.png'])
saveas(exp_fig_on,['activator_dwell_exp_on_coop_6bs.pdf'])

%%
exp_fig_off = figure;
hold on
% bar(ub_centers,p_vec_off,1,'FaceColor',cmap(3,:),'FaceAlpha',.9,'EdgeAlpha',1,'EdgeColor','k');   

% pd_vec = params(1)*exp(-[0.01 ub_centers]./params(2))/sum(c_vec_off);
plot(ub_centers,p_trend_off,'Color','k','LineWidth',3);
xlim([0 15])
% legend('independent binding','cooperative binding')
xlabel('activator dwell times (s)')
ylabel('probability')
box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
% set(gca,'Color','w')
exp_fig_off.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(exp_fig_off,['activator_dwell_exp_off_coop_6bs.png'])
saveas(exp_fig_off,['activator_dwell_exp_off_coop_6bs.pdf'])


%% now try to capture non-exponential behavior
close all
log_bins = logspace(log10(min(ub_centers)/2),log10(1e3),1e2);
log_centers = log_bins(1:end-1)+diff(log_bins)/2;

c_on_log = histcounts(unbinding_vec_on,log_bins);
c_on_log = c_on_log./diff(log_bins) * 2;
p_on_log = c_on_log/sum(c_on_log);

c_off_log = histcounts(unbinding_vec_off,log_bins);
c_off_log = c_off_log./diff(log_bins) * 2;
p_off_log = c_off_log/sum(c_off_log);

% p_vec_off_interp = interp1(ub_centers,p_vec_off,log_centers);
% ft_off = ~isnan(p_vec_off_interp)&p_vec_off_interp>0;
% p_vec_off_interp = p_vec_off_interp(ft_off);
% p_vec_on_interp = interp1(ub_centers,p_vec_on,log_centers);
% ft_on = ~isnan(p_vec_on_interp)&p_vec_on_interp>0;
% p_vec_on_interp = p_vec_on_interp(ft_on);

f_off_log = fit(log_centers',p_off_log','smoothingspline','SmoothingParam',0.3);
p_trend_off_log = feval(f_off_log,log_centers');
% f_off_log2 = fit(log_centers',p_trend_off_log,'smoothingspline','SmoothingParam',0.2);
% p_trend_off_log2 = feval(f_off_log2,log_centers');
f_on_log = fit(log_centers',p_on_log','smoothingspline','SmoothingParam',0.5);
p_trend_on_log = feval(f_on_log,log_centers');

% Make plots
close all

%% make exponential fit figure first
pwr_fig_on = figure;
hold on
% bar(ub_centers,p_vec_off,1,'FaceColor',cmap(3,:),'FaceAlpha',.9,'EdgeAlpha',1,'EdgeColor','k');   

% pd_vec = params(1)*exp(-[0.01 ub_centers]./params(2))/sum(c_vec_off);
plot(log_centers,imgaussfilt(p_trend_on_log,1),'Color','k','LineWidth',3);

set(gca,'Yscale','log')
set(gca,'Xscale','log')

% set(gca,'Xtick',xTickPoints,'xTickLabels',xTickLabelsStr)

ylim([5e-6 1e-1])
ax = gca;
ax.XLim = [min(ub_centers),1e2];

xlabel('activator dwell times (s)')
ylabel('probability')
% legend('independent binding','cooperative binding','Location','southwest')
set(gcf,'Color','w')
box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
pwr_fig.InvertHardcopy = 'off';

% legend('independent binding','cooperative binding')
xlabel('activator dwell times (s)')
ylabel('probability')
box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
% set(gca,'Color','w')
pwr_fig_on.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(pwr_fig_on,['activator_dwell_log_on_coop_6bs.png'])
saveas(pwr_fig_on,['activator_dwell_log_on_coop_6bs.pdf'])
%%
pwr_fig_off = figure;
hold on
% bar(ub_centers,p_vec_off,1,'FaceColor',cmap(3,:),'FaceAlpha',.9,'EdgeAlpha',1,'EdgeColor','k');   

% pd_vec = params(1)*exp(-[0.01 ub_centers]./params(2))/sum(c_vec_off);
plot(log_centers,imgaussfilt(p_trend_off_log,1),'Color','k','LineWidth',3);

set(gca,'Yscale','log')
set(gca,'Xscale','log')

% set(gca,'Xtick',xTickPoints,'xTickLabels',xTickLabelsStr)


ax = gca;
ax.XLim = [min(ub_centers),1e2];
ylim([5e-6 1e-1])
xlabel('activator dwell times (s)')
ylabel('probability')
% legend('independent binding','cooperative binding','Location','southwest')
set(gcf,'Color','w')
box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
pwr_fig.InvertHardcopy = 'off';

% legend('independent binding','cooperative binding')
xlabel('activator dwell times (s)')
ylabel('probability')
box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
% set(gca,'Color','w')
pwr_fig_off.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(pwr_fig_off,['activator_dwell_log_off_coop_6bs.png'])
saveas(pwr_fig_off,['activator_dwell_log_off_coop_6bs.pdf'])

% save data
results_struct = struct;
results_struct.dwell_time_axis = ub_centers;
results_struct.prob_vec_on_coop_raw = p_vec_on;
results_struct.prob_vec_off_coop_raw = p_vec_off;
results_struct.prob_vec_on_coop = p_trend_on;
results_struct.prob_vec_off_coop = p_trend_off;
results_struct.dwell_time_axis_log = log_centers;
results_struct.prob_vec_on_coop_log = p_trend_on_log;
results_struct.prob_vec_off_coop_log = p_trend_off_log;
save('results_struct.mat','results_struct')