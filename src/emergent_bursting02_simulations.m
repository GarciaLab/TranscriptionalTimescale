% conduct small-scale stochastic simulations to generate illustrative 
% traces for Figure 3
clear
close all
addpath('utilities')

% load numeric results
project = 'n6';

% set paths
DataPath = ['../out/emergent_bursting/' project '/'];

load([DataPath 'bursting_chain_calc_struct.mat'])
load([DataPath 'bursting_step_calc_struct.mat'])

% set stochastic simulation parameters
n_sim = 25; % number of indivudal traces to simulate
t_sim = 60*60; % duration of each simulation (in seconds)

p = gcp('nocreate');
if isempty(p)
  parpool;
end

% get list of model names
model_names = {bursting_chain_calc_struct.name};

%%%%%%%%%%%%%%%%%%%%%%%
% independent burst chain  first
%%%%%%%%%%%%%%%%%%%%%%%
ind_index = find(strcmp(model_names,'independent')); 
% pick intermediate parameter values where K_d=1, such that <r>=N/2
sim_param_index_ind = ceil(size(bursting_chain_calc_struct(ind_index).Q,3)/2);  
% seed for consistency (this doesn't do much)
r_seed = 123;

bursting_sim_struct = struct;

disp('Simulating independent binding condition...')
tic  
bursting_temp = stochastic_simulation_wrapper(bursting_chain_calc_struct(ind_index),sim_param_index_ind,n_sim,t_sim,r_seed);
fnames = fieldnames(bursting_temp);
for f = 1:length(fnames)
  bursting_sim_struct(1).(fnames{f}) = bursting_temp.(fnames{f});
end
toc

%%%%%%%%%%%%%%%%%%%%%%%
% cooperative binding
%%%%%%%%%%%%%%%%%%%%%%%

% now cooperative chain
coop_indices = find(contains(model_names,'cooperativity')); 
sim_param_indices_coop = 176:4:201;

% set seed for consistency
r_seed_vec = round(rand(1,length(bursting_step_calc_struct)+2)*1000);
disp('Simulating cooperative binding condition...')
tic
iter = 1;
for i = coop_indices
  bursting_temp = stochastic_simulation_wrapper(bursting_chain_calc_struct(i),sim_param_indices_coop,n_sim,t_sim,r_seed_vec(iter));
  fnames = fieldnames(bursting_temp);
  for f = 1:length(fnames)
    bursting_sim_struct(i).(fnames{f}) = bursting_temp.(fnames{f});
  end
  iter = iter + 1;
end
toc

% now iterate through compound chains for rate limiting step
rng(122);
r_seed_rl_vec = round(rand(1,length(bursting_step_calc_struct))*1000);
rl_index = 2; % only simulate 2 rate-limiting step case
sim_param_indices_RL = sim_param_indices_coop; % simulate systems with same dynamics as cooperative case 
offset = length(bursting_sim_struct); % for book-keeping

disp('Simulating rate-limiting step codition...')
tic
for i = rl_index
  bursting_temp = stochastic_simulation_wrapper(bursting_step_calc_struct(i),sim_param_indices_RL,n_sim,t_sim,r_seed_rl_vec(iter));
  fnames = fieldnames(bursting_temp);
  for f = 1:length(fnames)
    bursting_sim_struct(iter+1).(fnames{f}) = bursting_temp.(fnames{f});
  end
  iter = iter + 1;
end
toc
% save
save([DataPath, 'bursting_sim_struct.mat'],'bursting_sim_struct');