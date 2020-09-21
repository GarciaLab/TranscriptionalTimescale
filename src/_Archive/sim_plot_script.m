% script to simulate burst time series with parameters from numerical and
% analytic results
clear
close all
% load numeric results
project = 'n11_on6_koff33';
addpath('utilities')
% set paths
FigPath = ['../fig/emergent_bursting/' project '/'];
DataPath = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'results_struct.mat'])

% extract relevant parameters
kon_vec_fun = matlabFunction(results_struct.kon_vec);
koff_vec_fun = matlabFunction(results_struct.koff_vec);
tau_vec_fun = matlabFunction(results_struct.tau_vec);
on_states = results_struct.on_states;
off_states = results_struct.off_states;
n_bound_vec = results_struct.n_bound_vec;
n_states = numel(n_bound_vec);
% pick result to use
max_rate = 1; % in minutes
min_rate = .5;
rng(123)
options = find(results_struct.kon_eff_vec <= max_rate & results_struct.koff_eff_vec <= max_rate &...
                results_struct.kon_eff_vec >= min_rate & results_struct.koff_eff_vec >= min_rate);
sim_index = randsample(options,1);

% extract specific simulation parameters
kon_sim = results_struct.kon_samples(sim_index)*60;
b_sim = results_struct.b_samples(sim_index);
koff_sim = results_struct.koff_val * 60;
kon_val_vec = kon_vec_fun(b_sim,kon_sim);
koff_val_vec = koff_vec_fun(koff_sim);
tau_val_vec = tau_vec_fun(b_sim,koff_sim,kon_sim);

kon_eff_val = results_struct.kon_eff_fun(koff_sim,kon_sim);
koff_eff_val = results_struct.koff_eff_fun(b_sim,koff_sim,kon_sim);
% other simulation parameteres
t_sim = 60; % minutes
n_sim = 10; % number of simulations

%% simulate emeregent bursting
sim_results_emergent = struct;

% jump
for n = 1:n_sim
    state_val_vec = [1];
    jump_time_vec = [0];
    t_curr = 0;
    while t_curr < t_sim    
        % randomly select time til next state jump
        state_curr = state_val_vec(end);    
        next_jump = exprnd(tau_val_vec(state_curr));       
        t_next = min([t_curr + next_jump,t_sim]);
        % randomly choose next state
        up_thresh = kon_val_vec(state_curr) * tau_val_vec(state_curr);
        if rand() < up_thresh
            state_next = state_curr + 1;
        else
            state_next = state_curr - 1;
        end
        if t_next < t_sim
            state_val_vec = [state_val_vec state_next];
            jump_time_vec = [jump_time_vec t_next];
        end
        t_curr = t_next;
    end
    sim_results_emergent(n).state_val_vec = state_val_vec;
    sim_results_emergent(n).state_val_vec_eff = (state_val_vec-1)>=min(on_states);
    sim_results_emergent(n).jump_time_vec = jump_time_vec;
end    

% save results
save([DataPath 'sim_results_emergent.mat'],'sim_results_emergent')


% Plot results
close all
t_max = 30;
plot_index = 2;
blue = .7*[27 117 188]/256 + .3*[1 1 1];
green = [203 220 170]/256;

state_fig_discrete = figure;
cmap = brewermap([],'Set2');
hold on
s1 = patchline(sim_results_emergent(plot_index).jump_time_vec,(sim_results_emergent(plot_index).state_val_vec-1)/(n_states-1),...
    'EdgeColor',blue,'LineWidth',0.5,'edgealpha',1);

stairs(sim_results_emergent(plot_index).jump_time_vec,sim_results_emergent(plot_index).state_val_vec_eff,...
    'Color',[.3 .3 .3],'LineWidth',2)
% s1.FaceAlpha = 0;
ylim([0 1.05])
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',0:0.2:1)
p = plot(0,0);
StandardFigurePBoC(p,gca);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
state_fig_discrete.InvertHardcopy = 'off';

saveas(state_fig_discrete,[FigPath 'emergent_binding_plot_discrete.png'])
saveas(state_fig_discrete,[FigPath 'emergent_binding_plot_discrete.pdf'])


state_fig = figure;
hold on
s1 = patchline(sim_results_emergent(plot_index).jump_time_vec,(sim_results_emergent(plot_index).state_val_vec-1)/(n_states-1),...
    'EdgeColor',blue,'LineWidth',0.75,'edgealpha',.75);
ylim([0 1.05])
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',0:0.2:1)
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';

saveas(state_fig,[FigPath 'emergent_binding_plot.png'])
saveas(state_fig,[FigPath 'emergent_binding_plot.pdf'])


%% simulate simple 2 state on/off kinetics
state_options_simp = [1 2];
rate_vec_simp = [kon_eff_val, koff_eff_val];
sim_results_2state = struct;
% jump
for n = 1:n_sim
    state_val_vec = [1];
    jump_time_vec = [0];
    t_curr = 0;
    while t_curr < t_sim    
        % randomly select time til next state jump
        state_curr = state_val_vec(end);    
        next_jump = exprnd(1/rate_vec_simp(state_curr));       
        t_next = min([t_curr + next_jump,t_sim]);
        
        % next state
        state_next = state_options_simp(state_options_simp~=state_curr);
        
        if t_next < t_sim
            state_val_vec = [state_val_vec state_next];
            jump_time_vec = [jump_time_vec t_next];
        end
        t_curr = t_next;
    end
    sim_results_2state(n).state_val_vec = state_val_vec-1;    
    sim_results_2state(n).jump_time_vec = jump_time_vec;
end    

% save results
save([DataPath 'sim_results_2state.mat'],'sim_results_2state')

%% plot results
close all
plot_index_2s = 5;


state_fig_discrete = figure;
hold on
cmap = brewermap([],'Set2');

stairs(sim_results_2state(plot_index_2s).jump_time_vec,sim_results_2state(plot_index_2s).state_val_vec,...
    'Color',blue,'LineWidth',2)
ylim([0 1.05])
xlim([0 t_max])

ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',[0 1])

p = plot(0,0);
StandardFigurePBoC(p,gca);
state_fig_discrete.InvertHardcopy = 'off';

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';

saveas(state_fig_discrete,[FigPath 'two_state_plot.png'])
saveas(state_fig_discrete,[FigPath 'two_state_plot.pdf'])