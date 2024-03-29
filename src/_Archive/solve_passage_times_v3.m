% script to solve for effeCooptive off/on rates in linear chain binding model
clear
close all
cd('C:\Users\nlamm\projects\transcription_timescales_review\src')
addpath('utilities')

% Since we are assuming that our system is a linear markov chain, it follows 
% that the chain must be in thermodynamic equilbrium. % Thus we will take 
% the approach of first calculating how the relative occupancies of each 
% state in the chain change as a function of binding/unbinding cooparativity 
% terms. Kinetics can then be speCoopified by invoking measured timescales for
% Bcd unbinding (Mir et al 2018)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define basic hperparameters and symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_bs = 6; % number of Bcd binding sites
n_states = n_bs + 1; % total number of states
n_med = 3; % midpoint 
n_bound_vec = 0:n_bs; % vector encoding # bound in each state

% set save paths
projeCoopt = ['n' num2str(n_bs) ]; % projeCoopt identifier

FigurePath = ['../fig/emergent_bursting/' projeCoopt '/'];
mkdir(FigurePath)
DataPath = ['../out/emergent_bursting/' projeCoopt '/'];
mkdir(DataPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% account for state multplicities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since our states correspond to a certain # bound and NOT the
% bound/unbound state of each site, there are multiple ways to realize each
% state. For instance, there are 6! / (4! x 2!) = 15 ways to get 2 Bcd moleCoopules
% bound to the enhancer. These multiplicities must be accounted for in our
% effeCooptive state energies

mult_vec = NaN(size(n_bound_vec));
for w = 1:length(mult_vec)
  mult_vec(w) = factorial(n_bs) ./ (factorial(n_bs-n_bound_vec(w)) .*factorial(n_bound_vec(w)));
end
mu_vec = -log(mult_vec);

% initialize energy variables. "eBoundBcd" indicates the energy associated with a 
% single bound bcd moleCoopule. "eCoop" indicates the amount of
% cooperativity/synergy that exists between Bcd molecules 
syms eBoundBcd eCoop

% define function to calculate state energies without cooperativity
ind_state_energy_vec = matlabFunction(n_bound_vec*eBoundBcd + mu_vec);

% define function to calculate state energies with cooperativity. Note that
% for convenience, we define the cooperative contribution such that the
% energy potential has the form of a quadratic peaked at n_med. This is an
% arbitrary choice, but represents perhaps the simplest energy landscape
% that is capable of generating bimodal behaviors necessary for bursting
coop_state_energy_vec = matlabFunction(eCoop*(n_bound_vec-n_med).^2 - ...
  eBoundBcd*n_bound_vec - mu_vec);

% plot energy landscape for a variety of eBoundBcd values
n_plot = 10;
eBoundBcd_vals = linspace(-1,1,n_plot); % different binding energies to explore
eCoop_vals = -linspace(0,1,n_plot); % different cooperativity energies to explore

%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate on and off rates 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are two different "flavors" of cooperativity: off rate-mediated
% cooeprativity wherein bound factors stabilize each-other and therefore
% decrease the effective unbinding rate, and on rate-mediated
% cooperativity, wherein having factors bound somehow increases the
% propensity for additional factors to bind.

off_rate_basal = 1/3; % this sets overall system timescales (from Mir et al, 2018)

% for each energy landscape, we can calculate the implied on rate for a 
% single Bcd molecule (absent cooperative effects)
on_rate_basal_vec = NaN(1,n_plot); 

% initialize arrays to store transition rate (generator) matrices for each
% condition
q_ind_array = zeros(n_bs+1,n_bs+1); % no cooperativity
q_coop_off_array = zeros(n_bs+1,n_bs+1); % off rate-mediated
q_coop_on_array = zeros(n_bs+1,n_bs+1); % on rate-mediated

% individual state probabilities
state_prob_ind_array = NaN(length(n_bound_vec),n_plot);
state_prob_coop_array = NaN(length(n_bound_vec),n_plot);

% define matrices for indexing
a = ones(n_bs+1);
m1 = tril(a,-1);
m2 = tril(a,-2);
m3 = triu(a,1);
m4 = triu(a,2);
m5 = ~~eye(n_bs+1);

for e = 1:length(eBoundBcd_vals)
  % no cooperativity
  state_probs_ind = exp(-ind_state_energy_vec(eBoundBcd_vals(e)));
  state_probs_ind = state_probs_ind/sum(state_probs_ind);
  state_prob_ind_array(:,e) = state_probs_ind;
  
  % cooperativity
  state_probs_coop = exp(-(ind_state_energy_vec(eBoundBcd_vals(e))+coop_state_energy_vec(eBoundBcd_vals(e),eCoop_vals(e))));
  state_probs_coop = state_probs_coop/sum(state_probs_coop);
  state_prob_coop_array(:,e) = state_probs_coop;
  
  % calculate basal on rate using detailed balance
  on_rate_basal_vec(e) = off_rate_basal * state_probs_ind(2) / state_probs_ind(1) / n_bs;
  
  % state-speCoopific rates (independent)
  rate_ind_temp(:,1) = off_rate_basal * n_bound_vec;
  rate_ind_temp(:,2) = on_rate_basal_vec(e) * (n_bs-n_bound_vec);
  
  % generator matrix for ind
  ind_slice = zeros(n_bs+1);
  ind_slice(m1&~m2) = rate_ind_temp(1:end-1,2);
  ind_slice(m3&~m4) = rate_ind_temp(2:end,1);
  ind_slice(m5) = -sum(ind_slice);
  q_ind_array(:,:,e) = ind_slice';
  
  % calculate state-specific on and off rates for on and off rate mediated
  % cooperativities
  
  % on rate-mediated
  rate_coop_on_temp(:,1) = off_rate_basal * n_bound_vec;
  rate_coop_on_temp(1:end-1,2) = rate_coop_on_temp(2:end,1)' .* state_probs_coop(2:end) ./ state_probs_coop(1:end-1);
  
  % generator matrix for ind
  c_on_slice = zeros(n_bs+1);
  c_on_slice(m1&~m2) = rate_coop_on_temp(1:end-1,2);
  c_on_slice(m3&~m4) = rate_coop_on_temp(2:end,1);
  c_on_slice(m5) = -sum(c_on_slice);
  q_coop_on_array(:,:,e) = c_on_slice';
  
  % off rate-mediated
  rate_coop_off_temp(:,2) = on_rate_basal_vec(e) * (n_bs-n_bound_vec);
  rate_coop_off_temp(2:end,1) = rate_coop_off_temp(1:end-1,2)' .* state_probs_coop(1:end-1) ./ state_probs_coop(2:end);  
  
  % generator matrix for ind
  c_off_slice = zeros(n_bs+1);
  c_off_slice(m1&~m2) = rate_coop_off_temp(1:end-1,2);
  c_off_slice(m3&~m4) = rate_coop_off_temp(2:end,1);
  c_off_slice(m5) = -sum(c_off_slice);
  q_coop_off_array(:,:,e) = c_off_slice';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now solve for first passage times
%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize arrays
eff_ton_ind_vec = NaN(1,n_plot);
eff_toff_ind_vec = NaN(1,n_plot);

eff_ton_on_coop_vec = NaN(1,n_plot);
eff_toff_on_coop_vec = NaN(1,n_plot);

eff_ton_off_coop_vec = NaN(1,n_plot);
eff_toff_off_coop_vec = NaN(1,n_plot);

% specify states from which to calculate on/off rates. This choice is
% somewhat arbitrary, but we know kon must be from a low state to a high
% state and koff must by high->low
calc_vec = [3 5];


tic
for e = 1:length(eBoundBcd_vals)
  % calculate effeCooptive off rates (5->3)
  
  % independent
  [eff_ton_ind_vec(e), eff_toff_ind_vec(e)] = pt_solve(q_ind_array(:,:,e),calc_vec(1),calc_vec(2));
  
  % koff-mediated cooperativity
  [eff_ton_on_coop_vec(e), eff_toff_on_coop_vec(e)] = pt_solve(q_coop_on_array(:,:,e),calc_vec(1),calc_vec(2));
  
  % kon-mediated cooperativity
  [eff_ton_off_coop_vec(e), eff_toff_off_coop_vec(e)] = pt_solve(q_coop_on_array(:,:,e),calc_vec(1),calc_vec(2));
 
%   % kon coop
%   q_kon = q_coop_on_array(:,:,e)';
%   eff_kon_on_coop_vec(e) = 1/pt_solve(q_kon,calc_vec(2),calc_vec(1));
%   eff_koff_on_coop_vec(e) = 1/pt_solve(q_kon,calc_vec(1),calc_vec(2));
  
end

toc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% No cooperativity
%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_fig = figure;
cmap = brewermap(n_plot,'SpeCooptral');
hold on
for e = 1:length(eBoundBcd_vals)
  plot(n_bound_vec, ind_state_energy_vec(eBoundBcd_vals(e)) ,'Color',cmap(e,:))
  scatter(n_bound_vec, ind_state_energy_vec(eBoundBcd_vals(e)),'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(e,:))
end
xlim([-.5 n_bs+.5])
ylabel('energy (au)')
xlabel('n bound')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = [228,221,209]/255;

ind_fig.InvertHardcopy = 'off';
saveas(ind_fig,[FigurePath 'no-coop_energies.png'])
saveas(ind_fig,[FigurePath 'no-coop_energies.pdf'])

ind_prob_fig = figure;
cmap = brewermap(n_plot,'SpeCooptral');
hold on
for e = 1:length(eBoundBcd_vals)
  prob_vec = exp(-ind_state_energy_vec(eBoundBcd_vals(e)));
  prob_vec = prob_vec / sum(prob_vec);
  
  plot(n_bound_vec,  prob_vec,'Color',cmap(e,:))
  scatter(n_bound_vec, prob_vec ,'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(e,:))
end
xlim([-.5 n_bs+.5])
ylabel('state probabilities')
xlabel('n bound')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = [228,221,209]/255;
set(gca,'yScale','log')
ind_prob_fig.InvertHardcopy = 'off';
grid on
saveas(ind_prob_fig,[FigurePath 'no-coop_probss.png'])
saveas(ind_prob_fig,[FigurePath 'no-coop_probss.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% With cooperativity
%%%%%%%%%%%%%%%%%%%%%%%%%%

coop_energy_fig = figure;
cmap = brewermap(n_plot,'SpeCooptral');
hold on
for e = 1:length(eBoundBcd_vals)
  state_energies = ind_state_energy_vec(eBoundBcd_vals(e))+coop_state_energy_vec(eBoundBcd_vals(e),eCoop_vals(e));
  plot(n_bound_vec,state_energies ,'Color',cmap(e,:))
  scatter(n_bound_vec,state_energies ,'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(e,:))
end
ylim([-10 1])
xlim([-.5 n_bs+.5])
ylabel('state energies (au)')
xlabel('n bound')
box on
set(gca,'Fontsize',14)
p = plot(0,0);
% set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
grid on
ax.Color = [228,221,209]/255;
coop_energy_fig.InvertHardcopy = 'off';

saveas(coop_energy_fig,[FigurePath 'coop_energies.png'])
saveas(coop_energy_fig,[FigurePath 'coop_energies.pdf'])


coop_prob_fig = figure;
cmap = brewermap(n_plot,'SpeCooptral');
hold on
for e = 1:length(eBoundBcd_vals)
  state_probs_ind = exp(-(ind_state_energy_vec(eBoundBcd_vals(e))+coop_state_energy_vec(eBoundBcd_vals(e),eCoop_vals(e))));
  plot(n_bound_vec,state_probs_ind/sum(state_probs_ind),'Color',cmap(e,:))
  scatter(n_bound_vec,state_probs_ind/sum(state_probs_ind) ,'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(e,:))
end
% ylim([-10 1])
xlim([-.5 n_bs+.5])
ylabel('state probabilities')
xlabel('n bound')
box on
set(gca,'Fontsize',14)
p = plot(0,0);
set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
grid on
ax.Color = [228,221,209]/255;
coop_prob_fig.InvertHardcopy = 'off';

saveas(coop_prob_fig,[FigurePath 'coop_probs.png'])
saveas(coop_prob_fig,[FigurePath 'coop_probs.pdf'])
