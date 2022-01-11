% this script calculates transition rate matrices and state probabilities
% for various activator binding models discussed in the text. These
% quantities are used as the basis for stochastic simulations conducted in
% emergent_bursting02
clear
close all
addpath('../utilities')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define basic hyperparameters and symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_bs = 6; % number of Bcd binding sites
n_states = n_bs + 1; % total number of states
n_bound_vec = 0:n_bs; % vector encoding # bound in each state
off_rate_basal = 0.5; % in seconds. This sets overall system timescales 
n_calc_points = 201; % number of distinct activator concentrations
activator_grad_vec = logspace(-log10(10),log10(10),n_calc_points);

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

DataPath = ['../../out/activator_gradient/' project '/'];
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
% ecMax = 2; % ec corresponds to negative log of omega term in text (ec=-log(w))
% 6.6859
coopEnergies = repelem(-1.5,length(activator_grad_vec));%fliplr(linspace(-ecMax,ecMax,n_calc_points)); 
activatorEnergyVec = -0.5*coopEnergies*(n_bs-1) - log(activator_grad_vec); 

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
%   propensity_vec = rate_coop_off_temp(2:end,1)';
%   basal_vec = rate_coop_off_temp(2:end,1)'./n_bound_vec(2:end);
%   af = off_rate_basal/(sum(state_probs_coop(2:end).*basal_vec)/sum(state_probs_coop(2:end)));
%   rate_coop_off_temp = rate_coop_off_temp*af;
  
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

%%

close all
koff_vec = 1./eff_toff_off_coop_vec;
kon_vec = 1./eff_ton_off_coop_vec;

k_fig = figure;
hold on
plot(activator_grad_vec,60*kon_vec)
plot(activator_grad_vec,60*koff_vec)
legend('k_{on}','k_{off}','Location','southeast')
grid on
set(gca,'yscale','log')
set(gca,'xscale','log')

xlabel('[Activator]')
ylabel('events per minute')

set(gca,'FontSize',14)
ylim([1e-6 1e1])

saveas(k_fig,'kon_koff_vs_activator.png')
saveas(k_fig,'kon_koff_vs_activator.pdf')

r_fig = figure;
hold on
plot(activator_grad_vec,sum((0:6)' .* P_coop_array))
% legend('k_{on}','k_{off}','Location','southeast')
grid on
set(gca,'xscale','log')
xlabel('[Activator]')
ylabel('mean transcription rate')

set(gca,'FontSize',14)

saveas(r_fig,'mean rate_vs_activator_unbinding_coop.png')
saveas(r_fig,'mean_rate_vs_activator_unbinding_coop.pdf')

%%
koff_vec = 1./eff_toff_on_coop_vec;
kon_vec = 1./eff_ton_on_coop_vec;

k_fig = figure;
hold on
plot(activator_grad_vec,60*kon_vec)
plot(activator_grad_vec,60*koff_vec)
legend('k_{on}','k_{off}','Location','southeast')
grid on
set(gca,'yscale','log')
set(gca,'xscale','log')

xlabel('[Activator]')
ylabel('events per minute')

set(gca,'FontSize',14)
ylim([1e-5 1e2])

saveas(k_fig,'kon_koff_vs_activator_binding_coop.png')
saveas(k_fig,'kon_koff_vs_activator_binding_coop.pdf')

r_fig = figure;
hold on
plot(activator_grad_vec,sum((0:6)' .* P_coop_array))
% legend('k_{on}','k_{off}','Location','southeast')
grid on
set(gca,'xscale','log')
xlabel('[Activator]')
ylabel('mean transcription rate')

set(gca,'FontSize',14)

saveas(r_fig,'mean rate_vs_activator_binding_coop.png')
saveas(r_fig,'mean_rate_vs_activator_binding_coop.pdf')