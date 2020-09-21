% script to generate data for wiating time distributions of different
% models
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
ReadPath = ['../out/emergent_bursting/' project '/'];
WritePath = ['../out/waiting_time_distributions/' project '/'];
mkdir(WritePath);

% load data
load([ReadPath 'bursting_step_calc_struct.mat'])
load([ReadPath 'bursting_chain_calc_struct.mat'])

% simulation parameters
n_sim = 150; % number of indivudal traces to simulate
t_sim = 20*60*60; % duration of each simulation (in seconds)

%% %%%%%%%%%%%%%%%%%% run stochastic simulation script %%%%%%%%%%%%%%%%%%%%
p = gcp('nocreate');
if isempty(p)
  p = parpool;
end

%%%%%%%%%%%%%%%%%%%%%%%
% cooperative binding
%%%%%%%%%%%%%%%%%%%%%%%

simIndices = 176:4:201;

% set seed for consistency
r_seed_vec = round(rand(1,2+length(bursting_step_calc_struct))*1000);

disp('Simulating cooperative binding condition...')
tic
iter = 1;
for i = [2 3]
  bursting_temp = stochastic_simulation_wrapper(bursting_chain_calc_struct(i),simIndices,n_sim,t_sim,r_seed_vec(iter));
  fnames = fieldnames(bursting_temp);
  for f = 1:length(fnames)
    bursting_sim_struct(iter).(fnames{f}) = bursting_temp.(fnames{f});
  end
  iter = iter + 1;
end
toc

% now iterate through compound chains for rate limiting step
rng(122);

% simIndices = simIndicesCoop;%[simIndicesCoop(end-2) simIndicesCoop(end)];
offset = length(bursting_sim_struct);

%%%%%%%%%%%%%%%%%%%%%%%
% rate-limiting step binding
%%%%%%%%%%%%%%%%%%%%%%%

disp('Simulating rate-limiting step codition...')
tic
for i = 1:length(bursting_step_calc_struct)
  bursting_temp = stochastic_simulation_wrapper(bursting_step_calc_struct(i),simIndices,n_sim,t_sim,r_seed_vec(i));
  fnames = fieldnames(bursting_temp);
  for f = 1:length(fnames)
    bursting_sim_struct(i+offset).(fnames{f}) = bursting_temp.(fnames{f});
  end
  iter = iter + 1;
end
toc


%% %%%%%%%%%%%%%%%%%%% Extract waiting times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% set basic parameters
emergent_indices = 1:2;
rate_lim_indices = 3:length(bursting_sim_struct);
sim_indices = [emergent_indices rate_lim_indices];

dT = 5; % time res for interpolated data in seconds (effectively our experimental time res)
n_bound_vec = 0:6;
n_bs = 6;
waiting_time_struct = struct;
iter = 1;

for s = sim_indices
  
  sub_index_vec = 1:size(bursting_sim_struct(s).SS,2);    
  wt_off_cell_ideal = cell(1,length(sim_indices));

  for i = 1:length(sub_index_vec)    
    wt_off_vec_ideal = [];
    for n = 1:n_sim    
      
      % extract idealized waiting times
      raw_trace = bursting_sim_struct(s).sim_emission_cell{i,n};
      raw_times = bursting_sim_struct(s).sim_time_cell{i,n};
      % select for high and low periods
      hl_indices = find(raw_trace==0|raw_trace==n_bcd_sites);     
      condensed_trace = raw_trace(hl_indices);
      % find rise and fall events
      diff_trace = diff([0 condensed_trace]);
      drop_indices_raw = find(diff_trace<0);
      rise_indices_raw = find(diff_trace>0);
    
      if ~isempty(rise_indices_raw) && ~isempty(drop_indices_raw)            
        rise_indices = rise_indices_raw(rise_indices_raw>drop_indices_raw(1));
        drop_indices = drop_indices_raw(drop_indices_raw<rise_indices_raw(end));
        % convert back to real time space
        drop_times = raw_times(hl_indices(drop_indices));
        rise_times = raw_times(hl_indices(rise_indices));
        % record
        wt_off_vec_ideal = [wt_off_vec_ideal rise_times-drop_times];
      end
    end    
    wt_off_cell_ideal{i} = wt_off_vec_ideal;
  end

  % store results in a data structure    
  waiting_time_struct(iter).trace_array = trace_array; 
  waiting_time_struct(iter).off_waiting_times_ideal = wt_off_cell_ideal;
  waiting_time_struct(iter).name = bursting_sim_struct(s).name;
  iter = iter + 1;
end

% save
save([WritePath 'waiting_time_struct.mat'],'waiting_time_struct')