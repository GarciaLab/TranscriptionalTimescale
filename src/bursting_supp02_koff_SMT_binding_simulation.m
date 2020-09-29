% script to simulate full binding/unbinding model to check average
% unbinding rate
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


% sim name cell
sim_name_cell = {bursting_sim_struct.name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% koff-mediated cooperative binding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify appropriate index
coop_sim_index = find(contains(sim_name_cell,'koff-mediated'));

coop_plot_index = 1; 
sim_param_indices_coop = 176:4:201;

% extract transition rate matrix
Q = bursting_sim_struct(coop_sim_index).Q(:,:,coop_plot_index)';
SS = bursting_sim_struct(coop_sim_index).SS(:,coop_plot_index)';

% extract microscopic rates
n_vec = 0:n_bcd_sites;
n_states = length(n_vec);
a = ones(n_states); m1 = tril(a,-1); m2 = tril(a,-2); m3 = triu(a,1); m4 = triu(a,2); m5 = ~~eye(n_states);

% effective rates
k_minus_vec = [0 Q(m3&~m4)'];
k_plus_vec = [Q(m1&~m2)' 0];

% extract microscopic rates
k_unbind_vec = k_minus_vec ./n_vec;
k_unbind_vec(1) = 0;
k_bind_vec = k_plus_vec ./(n_bcd_sites-n_vec);
k_bind_vec(end) = 0;

%% %%%%%%%%%%%%%%%%%%%% build model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_sim = 10; % number of simulations
T = 1e4; % total time to simulate in seconds
rng(123);
% cell arrays to store results
event_time_cell = cell(1,n_sim);
bound_state_cell = cell(1,n_sim);

for n = 1:n_sim
  % initialize
  event_time_vec = [0];
  
  nBound = randsample(n_vec,1,true,SS);
  s_bound = randsample(1:n_bcd_sites,nBound,false);
  bound_state_array = zeros(1,n_bcd_sites);
  bound_state_array(s_bound) = 1;
  
  simTime = 0; 
  
  while simTime < T
    
    currState = bound_state_array(end,:);
    nBound = sum(currState);
    
    % draw expected jump time
    tau = 1 / (k_minus_vec(nBound+1) + k_plus_vec(nBound+1));
    dt = exprnd(tau);
    
    simTime = simTime + dt;
    
    % determine whether we bind or unbind
    eventType = randsample([-1 1],1,true,[k_minus_vec(nBound+1) k_plus_vec(nBound+1)]);
    
    if eventType == 1      
      options = find(currState==0);      
    else
      options = find(currState==1);
    end  
    % select binding site to update    
    bsSwitch = randsample(repelem(options,2),1,false);% repelem prevents undesirable behavior when only 1 option
    currState(bsSwitch) = currState(bsSwitch) + eventType;    
    bound_state_array(end+1,:) = currState;
    event_time_vec(end+1) = simTime;    
  end
  
  bound_state_cell{n} = bound_state_array;
  event_time_cell{n} = event_time_vec;
  
end

% extract single-molecules unbinding times
unbinding_vec = [];
binding_vec = [];
for n = 1:n_sim
  bound_state_array = bound_state_cell{n};
  event_time_vec = event_time_cell{n};
  
  for i = 1:n_bcd_sites
    binding_state_vec = bound_state_array(:,i);
    
    change_vec = [0 diff(binding_state_vec')];
    
    binding_times_raw = event_time_vec(change_vec==1);
    unbinding_times_raw = event_time_vec(change_vec==-1);
     
    if ~isempty(binding_times_raw) && ~isempty(unbinding_times_raw)            
      binding_times1 = binding_times_raw(binding_times_raw<unbinding_times_raw(end));
      unbinding_times1 = unbinding_times_raw(unbinding_times_raw>binding_times_raw(1));
      
      unbinding_vec = [unbinding_vec unbinding_times1-binding_times1];
      
      binding_times2 = binding_times_raw(binding_times_raw>unbinding_times_raw(1));
      unbinding_times2 = unbinding_times_raw(unbinding_times_raw<binding_times_raw(end));
      
      binding_vec = [binding_vec binding_times2-unbinding_times2];
    end
  end   
    
end

%% %%%%%%%%%%%% fit exponential to simulated dwell times and plot %%%%%%%%%
ub_bins = 0:0.1:max(unbinding_vec);
ub_centers = ub_bins(1:end-1)+diff(ub_bins)/2;

c_vec = histcounts(unbinding_vec,ub_bins);
p_vec = c_vec/sum(c_vec);

% fit raw counts
exp_fun = @(x) x(1) * exp(-ub_centers./x(2));
fit_obj = @(x) exp_fun(x)-c_vec;

options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',5e4);

params = lsqnonlin(fit_obj,[1 1], [0 0],[Inf Inf],options);

% Make plots
close all

% make exponential fit figure first
bar_fig = figure;
cmap = brewermap(9,'Set2');

hold on
bar(ub_centers,p_vec,1,'FaceColor',cmap(3,:),'FaceAlpha',.9,'EdgeAlpha',1,'EdgeColor','k');   

pd_vec = params(1)*exp(-[0.01 ub_centers]./params(2))/sum(c_vec);
plot([0.01 ub_centers],pd_vec,'Color','k','LineWidth',2);

xlim([0 2])
legend('activator dwell times',['exponential fit (\mu = ' num2str(round(params(2),2)) 's)'])
xlabel('activator dwell times (s)')
ylabel('probability')
% box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
bar_fig.InvertHardcopy = 'off';
saveas(bar_fig,[FigurePath 'activator_dwell_exp.png'])
saveas(bar_fig,[FigurePath 'activator_dwell_exp.pdf'],'pdf')


% now try to capture non-exponential behavior
close all

log_bins = logspace(-4,log10(max(unbinding_vec)),3e2);
log_centers = log_bins(1:end-1)+diff(log_bins)/2;

p_vec_interp = interp1(ub_centers,p_vec,log_centers);

pwr_fig = figure;
cmap = brewermap(9,'Set2');
hold on
% scatter(log_centers,log_counts)
% b = bar(log_prob,1,'FaceColor',cmap(3,:),'FaceAlpha',1,'EdgeAlpha',1); 
scatter(log_centers,p_vec_interp,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
plot([.01 ub_centers],pd_vec,'Color','k','LineWidth',2);

set(gca,'Yscale','log')
set(gca,'Xscale','log')

% set(gca,'Xtick',xTickPoints,'xTickLabels',xTickLabelsStr)

ylim([ 1e-10 1e0])
ax = gca;
ax.XLim = [min(ub_centers),1e3];
xlabel('activator dwell times (sec)')
ylabel('probability')
lgd = legend('dwell time pdf',['exponential fit (\mu = ' num2str(round(params(2),2)) 's)']);
box off
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
pwr_fig.InvertHardcopy = 'off';

lgd.FontSize=12;

saveas(pwr_fig,[FigurePath 'activator_dwell_pwr.png'])
saveas(pwr_fig,[FigurePath 'activator_dwell_pwr.pdf'])

%% %%%%%%% Now examine what happens if we forbid long binding events %%%%%5
max_dwell = 10;
event_time_trunc_cell = cell(1,n_sim);
bound_state_trunc_cell = cell(1,n_sim);

for n = 1:n_sim
  % initialize
  event_time_vec = [0];
  
  nBound = randsample(n_vec,1,true,SS);
  s_bound = randsample(1:n_bcd_sites,nBound,false);
  bound_state_array = zeros(1,n_bcd_sites);
  bound_state_array(s_bound) = 1;
  
  % keep track of duration of current binding state
  bound_time_vec = rand(1,n_bcd_sites)*max_dwell.*bound_state_array(1,:);
  
  simTime = 0; 
  
  while simTime < T
    
    currState = bound_state_array(end,:);
    nBound = sum(currState);
        
    % draw expected jump time
    tau = 1 / (k_minus_vec(nBound+1) + k_plus_vec(nBound+1));
    dt = exprnd(tau);
    
    % check whether any site will reach max bound time before this
    % transition
    bound_indices = find(currState==1);
    [dtDwell, mi] = min(max_dwell-bound_time_vec(bound_indices));
    
    bsDwell = bound_indices(mi);
    
    trunc_flag = 0;
    if dtDwell <= dt
      dt = dtDwell;
      trunc_flag = 1;
      eventType = -1;
      bsSwitch = bsDwell;
    else
            
      % determine whether we bind or unbind
      eventType = randsample([-1 1],1,true,[k_minus_vec(nBound+1) k_plus_vec(nBound+1)]);

      if eventType == 1      
        options = find(currState==0);      
      else
        options = find(currState==1);
      end  
      % select binding site to update    
      bsSwitch = randsample(repelem(options,2),1,false);% repelem prevents undesirable behavior when only 1 option
    end
    
    simTime = simTime + dt;
    
    currState(bsSwitch) = currState(bsSwitch) + eventType;  
    
    bound_time_vec = (bound_time_vec + dt).*currState;
              
    bound_state_array(end+1,:) = currState;
    event_time_vec(end+1) = simTime;  
        
    if eventType == -1
      bound_time_vec(bsSwitch) = 0;
    end
  end
  
  bound_state_trunc_cell{n} = bound_state_array;
  event_time_trunc_cell{n} = event_time_vec;
  
end

%% %%%%%%%% plot the trace and histogram for truncated binding case %%%%%%%
plot_index = 2;
time_raw = event_time_trunc_cell{plot_index};
trace_raw = sum(bound_state_trunc_cell{plot_index},2);
time_rs = 0:resamp_res:3600;
t_max = max(time_rs)/60;

resamp_res = 0.5; % in seconds
blue = [190 201 224]/255;
ylimTrace = [-0.5 6.5];
n_bound_vec = 0:6;

state_fig_truncated = figure;
cmap2 = brewermap(9,'Set2');
hold on

trace_rs = interp1(time_raw, trace_raw,time_rs,'previous');

stairs(time_rs/60, trace_rs,'Color',[blue 0.0],'LineWidth',1.5);

ylim(ylimTrace)
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',n_bound_vec)
p = plot(0,0);
StandardFigurePBoC(p,gca);
state_fig_truncated.InvertHardcopy = 'off';

saveas(state_fig_truncated,[FigurePath 'coop_trace_trunc.png'])
saveas(state_fig_truncated,[FigurePath 'coop_trace_trunc.pdf'])

%% calculate fraction of time in eac state
dur_vec = diff(time_raw);
stateSums = accumarray(trace_raw(1:end-1)+1,dur_vec');
stateShares = stateSums/sum(stateSums);

hist_fig = figure('Position',[100 100 256 512]);
hold on
barh(n_bound_vec,stateShares,1,'FaceColor',blue);
  
xlabel('probability')
box on
p = plot(0,0);
% ax = gca;
% ax.YColor = 'black';
% xlim([0 0.55])
% ax.XColor = 'black';
set(gca,'Fontsize',14,'xtick',0:.25:1)
ylim([n_bound_vec(1)-0.5 n_bound_vec(end)+0.5])
StandardFigurePBoC(p,gca);

hist_fig.InvertHardcopy = 'off';
saveas(hist_fig,[FigurePath 'coop_hist_trunc.png'])
saveas(hist_fig,[FigurePath 'coop_hist_trunc.pdf'])
