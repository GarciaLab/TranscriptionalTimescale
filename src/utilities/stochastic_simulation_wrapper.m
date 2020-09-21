% script to simulate burst time series with parameters from analytic
% calculations
function bursting_sim_struct = stochastic_simulation_wrapper(bursting_calc_struct,simSubIndices,n_sim,t_sim,r_seed)


% independent burst chain  first
bursting_sim_struct = struct;
% set seed for consistency
rng(r_seed)

% record key info
bursting_sim_struct.name = bursting_calc_struct.name;
bursting_sim_struct.E = bursting_calc_struct.E;
for p = 1:length(simSubIndices)
  % record network characteristics
  bursting_sim_struct.Q(:,:,p) = bursting_calc_struct.Q(:,:,simSubIndices(p));
  bursting_sim_struct.SS(:,p) = bursting_calc_struct.SS(:,simSubIndices(p));   
  if isfield(bursting_calc_struct,'SSFull')
    bursting_sim_struct.SSFull(:,p) = bursting_calc_struct.SSFull(:,p);
  else
    bursting_sim_struct.SSFull(:,p) = bursting_sim_struct.SS(:,p);
  end
  
  % call core stochastic simulation function
  [bursting_sim_struct.sim_emission_cell(p,:), bursting_sim_struct.sim_emission_cell(p,:),bursting_sim_struct.sim_time_cell(p,:)] = ...
        stochastic_sim_fun(bursting_sim_struct.Q(:,:,p), bursting_sim_struct.SSFull(:,p),bursting_sim_struct.E,n_sim,t_sim);

end