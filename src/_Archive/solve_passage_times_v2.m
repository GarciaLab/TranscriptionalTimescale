% script to solve for effective off/on rates in linear chain binding model
% script also performs parameter sweeps to identify "bursty" regimes that
% can be simulated to generate example traces for figures and further
% analysis

% Created by: Nick Lammers
% Date: 4/17/2020

% clear workspace
clear
close all
addpath('utilities')

% specify basic model characteristics
n_bound_vec = 0:6; % number of different binding tates
n_states = numel(n_bound_vec);

% ON and OFF blocks (used for calculating effective bursting rates)
on_states = 5:7; 
off_states = 1:3;
midpoint = 4; % state where nonlinear rate is at half-max

% specify binding dynamics
syms kon koff amp n_hill
nonlin_rate = 'kon';
on_nonlin_flag = strcmp(nonlin_rate,'kon');
off_nonlin_flag = strcmp(nonlin_rate,'koff');
% define null rates without nonlinearity
koff_vec = koff*n_bound_vec; % assume Bcd unbinds independently
kon_vec = kon*(n_states-n_bound_vec-1); % assume Bcd unbinds independently

if on_nonlin_flag    
    kon_vec = kon_vec.*(amp .* n_bound_vec.^n_hill ./ ((midpoint-1).^n_hill + n_bound_vec.^n_hill));
end
if off_nonlin_flag   
    koff_vec = koff_vec.*(amp .* (n_states-n_bound_vec).^n_hill ./ ((midpoint-1).^n_hill + (n_states-n_bound_vec-1).^n_hill));
end

% state dwell time vec
tau_vec = 1./(kon_vec + koff_vec);
% koff_val (from Mir et al 2018)
koff_val = 1/3;
% specify simulation type
project = ['n' num2str(max(n_bound_vec)) '_' nonlin_rate];
% set figure path
FigPath = ['../fig/emergent_bursting/' project '/'];
mkdir(FigPath)
DataPath = ['../out/emergent_bursting/' project '/'];
mkdir(DataPath)

%% generate expressions for first passage times between ON and OFF blocks

kon_states = [off_states];
% OFF to ON
kon_eq_list = {};
kon_sym_list = [];
kon_str_list = {};

% generate list of symbols to use for calculating effective on rate
for i = kon_states
    % initialize symbol
    sym_str = ['ET' num2str(i-1) num2str(midpoint-1)];
    kon_str_list{i} = sym_str;
    eval(['syms ' sym_str])
    kon_sym_list = [kon_sym_list eval(sym_str)];
end

% generate equations
for i = (kon_states)
    next_on = 0;
    if i < max(kon_states)
        next_on = kon_sym_list(i+1);
    end
    next_off = 0;
    if i > min(kon_states)
        next_off = kon_sym_list(i-1);
    end
    kon_eq_list{i} = kon_sym_list(i) == tau_vec(i)*(1 + kon_vec(i)*next_on + koff_vec(i)*next_off);
end

% solve 
solKon = solve(kon_eq_list,kon_sym_list);
eval(['off_dwell = solKon.' kon_str_list{end} ';']);
kon_eff_fun = matlabFunction(1/off_dwell);


% ON to OFF
koff_states = [on_states];

koff_eq_list = {};
koff_sym_list = [];
koff_str_list = {};
% generate list of symbols to use
for i = koff_states
    % initialize symbol
    sym_str = ['ET' num2str(i-1) num2str(midpoint-1)];
    koff_str_list = [koff_str_list{:} {sym_str}];
    eval(['syms ' sym_str])
    koff_sym_list = [koff_sym_list eval(sym_str)];
end

% generate equation
os = numel(off_states)+1;
iter = 1;
for i = koff_states 
    next_on = 0;
    if iter < numel(koff_sym_list)
        next_on = koff_sym_list(iter+1);
    end
    next_off = 0;
    if iter > 1
        next_off = koff_sym_list(iter-1);
    end
    koff_eq_list = [koff_eq_list {koff_sym_list(iter) == tau_vec(i)*(1 + kon_vec(i)*next_on + koff_vec(i)*next_off)}];
    iter = iter + 1;
end
% solve 
solKoff = solve(koff_eq_list,koff_sym_list);
eval(['on_dwell = solKoff.' koff_str_list{1} ';']);
koff_eff_fun = matlabFunction(1/on_dwell);


%% conduct simple parameter sweep of possible solutions'
n_samp = 1e4;
amp_samples = 10.^(rand(1,n_samp)*3);
kon_samples = 10.^(rand(1,n_samp)*2-2);
n_hill = 100;

kon_eff_vec = NaN(1,n_samp);
koff_eff_vec = NaN(1,n_samp);
tic
for i = 1:n_samp
    kon_eff_vec(i) = 60*kon_eff_fun(amp_samples(i),koff_val,kon_samples(i),n_hill);%(kon_samples(i),1,b_samples(i)));
    koff_eff_vec(i) = 60*koff_eff_fun(amp_samples(i),koff_val,kon_samples(i),n_hill);%double(koff_eff(kon_samples(i),1,b_samples(i)));
end
toc

slow_ft = kon_eff_vec <= 1 & koff_eff_vec <= 1;
cmap = brewermap([],'Set2');
close all

scatter_fig = figure;
hold on
scatter(kon_eff_vec,koff_eff_vec,20,cmap(3,:))
scatter(kon_eff_vec(slow_ft),koff_eff_vec(slow_ft),20,cmap(2,:))
set(gca,'YScale','log','XScale','log')    
set(gca,'Fontsize',14)
xlabel('effective on rate (events/min)')
ylabel('effective off rate (events/min)')
grid on  

%%
saveas(scatter_fig,[FigPath 'param_sweep.png'])    
    
% save analytic results
results_struct = struct;
results_struct.kon_eff_fun = kon_eff_fun;
results_struct.koff_eff_fun = koff_eff_fun;
results_struct.koff_val = koff_val;

results_struct.kon_eff_vec = kon_eff_vec;
results_struct.koff_eff_vec = koff_eff_vec;
results_struct.kon_samples = kon_samples;
results_struct.b_samples = amp_samples;

results_struct.n_bound_vec = n_bound_vec;
results_struct.on_states = on_states;
results_struct.off_states = off_states;
results_struct.kon_vec = kon_vec;
results_struct.koff_vec = koff_vec;
results_struct.tau_vec = tau_vec;

save([DataPath 'results_struct.mat'],'results_struct')