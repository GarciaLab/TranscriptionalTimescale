% this script calculates transition rate matrices and state probabilities
% for various activator binding models discussed in the text. These
% quantities are used as the basis for stochastic simulations conducted in
% emergent_bursting02
clear
close all
addpath('utilities')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define basic hyperparameters and symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_bs = 6; % number of Bcd binding sites
n_states = n_bs + 1; % total number of states
n_bound_vec = 0:n_bs; % vector encoding # bound in each state
off_rate_basal = 1/2; % in seconds. This sets overall system timescales (from Mir et al, 2018)
n_calc_points = 201; % number of distinct binding and cooperative energy values to explore
rate_lim_steps_vec = [1:2 5 10 15]; % numbers of rate-limiting steps to simulate (all in OFF pathway)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Calculations for simple chain model with and without cooperativity
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since we are assuming that our system is a linear markov chain, it follows 
% that the chain must be in thermodynamic equilbrium. Thus we will take 
% the approach of first calculating how the relative occupancies of each 
% state in the chain change as a function of binding/unbinding cooparativity 
% terms. Kinetics can then be specified by invoking measured timescales for
% Bcd unbinding (Mir et al 2018). Note that all energies are in kT units

% set save paths
project = ['n' num2str(n_bs) ]; % project identifier

DataPath = ['../out/emergent_bursting/' project '/'];
mkdir(DataPath)

% accounting for state state multplicities

% Since our states correspond to a certain # bound and NOT the
% bound/unbound state of each site, there are multiple ways to realize each
% state. For instance, there are 6! / (4! x 2!) = 15 ways to get 2 Bcd moleCoopules
% bound to the enhancer. These multiplicities must be accounted for in our
% effective state energies

mult_vec = NaN(size(n_bound_vec));
for w = 1:length(mult_vec)
  mult_vec(w) = factorial(n_bs) ./ (factorial(n_bs-n_bound_vec(w)) .*factorial(n_bound_vec(w)));
end
mu_vec = -log(mult_vec); % energy contribution from state multiplicities

% specify different binding energies to explore
ecMax = 2; % ec corresponds to negative log of omega term in text (ec=-log(w))
coopEnergies = fliplr(linspace(-ecMax,ecMax,n_calc_points)); 
activatorEnergyVec = -0.5*coopEnergies*(n_bs-1); % this definition ensures symmetric state probabilities (simply for convenience)

% define a function to calculate state probabilities for a given binding
% energy and cooperativity. For simplicity, we consider only simple pairwise 
% cooperativity between, which leads to a linear energy landscape. 

stateEnergy_fun = @(eb,ec) eb*n_bound_vec + ec*n_bound_vec.*(n_bound_vec-1)/2 + mu_vec;
stateProb_fun = @(eb,ec) exp(-stateEnergy_fun(eb,ec)) / sum(exp(-stateEnergy_fun(eb,ec)));

% for each energy landscape, we can calculate the implied on rate for a 
% single Bcd molecule (absent cooperative effects) unsing known koff and
% detailed balance
on_rate_basal_vec = NaN(1,n_calc_points); 

% initialize 3D arrays to store transition rate (generator) matrices for each
% condition
Q_ind_array = zeros(n_states,n_states,n_calc_points); % independent binding
Q_coop_off_array = zeros(n_states,n_states,n_calc_points); % off rate-mediated
Q_coop_on_array = zeros(n_states,n_states,n_calc_points); % on rate-mediated

P_coop_array = NaN(n_states,n_calc_points);
P_ind_array = NaN(n_states,n_calc_points);
% define matrices used for indexing
a = ones(n_states); m1 = tril(a,-1); m2 = tril(a,-2); m3 = triu(a,1); m4 = triu(a,2); m5 = ~~eye(n_states);

% calculate metrics for (a) and (b) first
for eb = 1:length(activatorEnergyVec)  
  eBind = activatorEnergyVec(eb);
  eCoop = coopEnergies(eb);
  % independent binding
  state_probs_ind = stateProb_fun(eBind,0);  
  P_ind_array(:,eb) = state_probs_ind;
  % cooperativity
  state_probs_coop = stateProb_fun(eBind,eCoop);    
  P_coop_array(:,eb) = state_probs_coop;
  
  % calculate basal on rate using detailed balance
  on_rate_basal_vec(eb) = off_rate_basal * state_probs_ind(2) / state_probs_ind(1) / n_bs;

  % independent   
  rate_ind_temp(:,1) = off_rate_basal * n_bound_vec;
  rate_ind_temp(1:end-1,2) = rate_ind_temp(2:end,1)' .* state_probs_ind(2:end) ./ state_probs_ind(1:end-1);

  % generator matrix for ind
  c_ind_slice = zeros(n_states);
  c_ind_slice(m1&~m2) = rate_ind_temp(1:end-1,2);
  c_ind_slice(m3&~m4) = rate_ind_temp(2:end,1);
  c_ind_slice(m5) = -sum(c_ind_slice);
  Q_ind_array(:,:,eb) = c_ind_slice';

  % on rate-mediated
  rate_coop_on_temp(:,1) = off_rate_basal * n_bound_vec;
  rate_coop_on_temp(1:end-1,2) = rate_coop_on_temp(2:end,1)' .* state_probs_coop(2:end) ./ state_probs_coop(1:end-1);

  % generator matrix for ind
  c_on_slice = zeros(n_states);
  c_on_slice(m1&~m2) = rate_coop_on_temp(1:end-1,2);
  c_on_slice(m3&~m4) = rate_coop_on_temp(2:end,1);
  c_on_slice(m5) = -sum(c_on_slice);
  Q_coop_on_array(:,:,eb) = c_on_slice';

  % off rate-mediated
  rate_coop_off_temp(:,2) = on_rate_basal_vec(1,eb) * (n_bs-n_bound_vec);
  rate_coop_off_temp(2:end,1) = rate_coop_off_temp(1:end-1,2)' .* state_probs_coop(1:end-1) ./ state_probs_coop(2:end); 
  
  % re-adjust rates to ensure consistency with experimental measurements
  % this WILL NOT yield average microscopic rate (taken across all events)
  % will be precisely equal to Mir2018 cvalue but it will ensure that 
  % micrsocopic rates out of each state are on a reasonably scale
  propensity_vec = rate_coop_off_temp(2:end,1)';
  basal_vec = rate_coop_off_temp(2:end,1)'./n_bound_vec(2:end);
  af = off_rate_basal/(sum(state_probs_coop(2:end).*basal_vec)/sum(state_probs_coop(2:end)));
  rate_coop_off_temp = rate_coop_off_temp*af;
  
  % generator matrix for ind
  c_off_slice = zeros(n_states);
  c_off_slice(m1&~m2) = rate_coop_off_temp(1:end-1,2);
  c_off_slice(m3&~m4) = rate_coop_off_temp(2:end,1);
  c_off_slice(m5) = -sum(c_off_slice);
  Q_coop_off_array(:,:,eb) = c_off_slice';

end

% solve for effective on and off rates (0<->6). This is not strictly necessary but
% will be used as a consistency check that stochastic simulations are
% behaving as expected

% initialize arrays
eff_ton_ind_vec = NaN(1,n_calc_points);
eff_toff_ind_vec = NaN(1,n_calc_points);

eff_ton_on_coop_vec = NaN(1,n_calc_points);
eff_toff_on_coop_vec = NaN(1,n_calc_points);

eff_ton_off_coop_vec = NaN(1,n_calc_points);
eff_toff_off_coop_vec = NaN(1,n_calc_points);

% specify states from which to calculate on/off rates. This choice is
% somewhat arbitrary, but we know kon must be from a low state to a high
% state and koff must by high->low
calc_vec = [1 n_states];

for eb = 1:length(activatorEnergyVec) 
  % ind
  [eff_ton_ind_vec(eb), eff_toff_ind_vec(eb)] = pt_solve(Q_ind_array(:,:,eb),calc_vec(1),calc_vec(2));
  
  % koff-mediated cooperativity
  [eff_ton_on_coop_vec(eb), eff_toff_on_coop_vec(eb)] = pt_solve(Q_coop_on_array(:,:,eb),calc_vec(1),calc_vec(2));

  % kon-mediated cooperativity
  [eff_ton_off_coop_vec(eb), eff_toff_off_coop_vec(eb)] = pt_solve(Q_coop_off_array(:,:,eb),calc_vec(1),calc_vec(2)); 
end

% initialize data structure
bursting_chain_calc_struct = struct;

% independent binding
bursting_chain_calc_struct(1).name = 'independent';
bursting_chain_calc_struct(1).Q = Q_ind_array;
bursting_chain_calc_struct(1).eff_on_states = calc_vec;
bursting_chain_calc_struct(1).eff_rates = [1./eff_ton_ind_vec' 1./eff_toff_ind_vec'];
bursting_chain_calc_struct(1).SS = P_ind_array;
% kon-mediated cooperative binding/unbinding 
bursting_chain_calc_struct(2).name = 'kon-mediated cooperativity';
bursting_chain_calc_struct(2).Q = Q_coop_on_array;
bursting_chain_calc_struct(2).eff_on_states = calc_vec;
bursting_chain_calc_struct(2).eff_rates = [1./eff_ton_on_coop_vec' 1./eff_toff_on_coop_vec'];
bursting_chain_calc_struct(2).SS = P_coop_array;
% koff-mediated cooperative binding/unbinding 
bursting_chain_calc_struct(3).name = 'koff-mediated cooperativity';
bursting_chain_calc_struct(3).Q = Q_coop_off_array;
bursting_chain_calc_struct(3).eff_on_states = calc_vec;
bursting_chain_calc_struct(3).eff_rates = [1./eff_ton_off_coop_vec' 1./eff_toff_off_coop_vec'];
bursting_chain_calc_struct(3).SS = P_coop_array;
% generic fields
for i = 1:3
  bursting_chain_calc_struct(i).edge_score = bursting_chain_calc_struct(i).SS(1,:)+bursting_chain_calc_struct(i).SS(end,:);
  bursting_chain_calc_struct(i).E = n_bound_vec;
  bursting_chain_calc_struct(i).off_rate_basal = off_rate_basal;
  bursting_chain_calc_struct(i).on_rate_basal_vec = on_rate_basal_vec;
  bursting_chain_calc_struct(i).activatorEnergies = activatorEnergyVec;
  bursting_chain_calc_struct(i).coopEnergies = coopEnergies;
end
save([DataPath, 'bursting_chain_calc_struct.mat'],'bursting_chain_calc_struct');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) generate transition rate matrices for rate-limiting step models
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now calculate metrics for (c)
ind_state_vec = [26 176]; % define the two states to lie near the extrema

% set range of switching kinetics to explore (use kon coop model as ref
% point)
slow_kinetics_array = [1./eff_ton_on_coop_vec' 1./eff_toff_on_coop_vec']; 

% iterate through different numbers of rate-limiting steps
bursting_step_calc_struct = struct;

for n =1:length(rate_lim_steps_vec)
  
  % generate transcription rate vec
  emission_vec = repmat(n_bound_vec,1,rate_lim_steps_vec(n)+1);
  
  % initialize array to store rate matrics
  Q_rate_lim_array = zeros((rate_lim_steps_vec(n)+1)*n_states,(rate_lim_steps_vec(n)+1)*n_states,n_calc_points); % on rate-mediated
  P_lim_array = NaN(n_states,n_calc_points);
  P_lim_array_full = NaN((rate_lim_steps_vec(n)+1)*n_states,n_calc_points);
  % indexing vec for convenience
  sub_vec = 1:n_states;

  for eb = 1:n_calc_points
    
    % extract component rate arrays
    QON = Q_ind_array(:,:,ind_state_vec(1));
    QOFF = Q_ind_array(:,:,ind_state_vec(2));
    
    % estimate expected effective occupancies (assuming strong timescale
    % separation)
    SSON = P_ind_array(:,ind_state_vec(1));
    SSOFF = P_ind_array(:,ind_state_vec(2));
    
    % effective ss
    pon = slow_kinetics_array(eb,1)/(slow_kinetics_array(eb,2)+slow_kinetics_array(eb,1));
    P_lim_array(:,eb) = SSON*pon + SSOFF*(1-pon);
    
    % full ss (same dims as Q)
    P_lim_array_full(:,eb) = [repmat(SSOFF/rate_lim_steps_vec(n)*(1-pon),rate_lim_steps_vec(n),1) ; SSON*pon];
    
    % initialize combined rate array
    Q_slice = Q_rate_lim_array(:,:,eb);
    
    % add "OFF" blocks
    for m = 1:rate_lim_steps_vec(n)
      % fill in binding rate block
      Q_slice(sub_vec+(m-1)*n_states,sub_vec+(m-1)*n_states) = QOFF';
      
      % calculate indices for slow exchange terms
      linear_indices = sub2ind(size(Q_slice),sub_vec+m*n_states,sub_vec+(m-1)*n_states);
      
      % add exchange terms
      Q_slice(linear_indices) = slow_kinetics_array(eb,1)*rate_lim_steps_vec(n); % on rate
            
    end
    % add "ON" block
    Q_slice(sub_vec+rate_lim_steps_vec(n)*n_states,sub_vec+rate_lim_steps_vec(n)*n_states) = QON';
    linear_indices = sub2ind(size(Q_slice),sub_vec,sub_vec+rate_lim_steps_vec(n)*n_states);
    Q_slice(linear_indices) = slow_kinetics_array(eb,2);
    
    % enforce proper form
    Q_slice(eye(size(Q_slice))==1) = 0;
    Q_slice(eye(size(Q_slice))==1) = -sum(Q_slice);
    Q_rate_lim_array(:,:,eb) = Q_slice';   
  end
  
  % record
  bursting_step_calc_struct(n).name = [num2str(rate_lim_steps_vec(n)) 'rate-limiting steps'];
  bursting_step_calc_struct(n).Q = Q_rate_lim_array; 
  bursting_step_calc_struct(n).SS = P_lim_array; 
  bursting_step_calc_struct(n).SSFull =P_lim_array_full;
  bursting_step_calc_struct(n).edge = P_lim_array(1,:)+P_lim_array(end,:); 
  bursting_step_calc_struct(n).eff_on_states = calc_vec;
  bursting_step_calc_struct(n).eff_rates = [1./eff_ton_off_coop_vec' 1./eff_toff_off_coop_vec'];
  bursting_step_calc_struct(n).E = emission_vec;
  bursting_step_calc_struct(n).off_rate_basal = off_rate_basal;
  bursting_step_calc_struct(n).on_rate_basal_vec = on_rate_basal_vec;
  bursting_step_calc_struct(n).activatorEnergies = activatorEnergyVec;
  bursting_step_calc_struct(n).coopEnergies = coopEnergies;
  bursting_step_calc_struct(n).ind_state_vec = ind_state_vec;
end

save([DataPath, 'bursting_step_calc_struct.mat'],'bursting_step_calc_struct');
