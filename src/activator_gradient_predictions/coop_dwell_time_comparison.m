% this script calculates transition rate matrices and state probabilities
% for various activator binding models and runs simple simulations to get
% dwell time distributions
clear
close all
addpath('../utilities')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define basic hyperparameters and symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_bs = 3; % number of Bcd binding sites
n_states = n_bs + 1; % total number of states
n_bound_vec = 0:n_bs; % vector encoding # bound in each state
off_rate_basal = 1/1.2640; % in seconds. This sets overall system timescales (from Mir et al, 2018)
n_calc_points = 1;%201; % number of distinct binding and cooperative energy values to explore

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
% ecMax = 2; % ec corresponds to negative log of omega term in text (ec=-log(w))
coopEnergies = -2;%fliplr(linspace(-ecMax,ecMax,n_calc_points)); 
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
  af = 1;%off_rate_basal/(sum(state_probs_coop(2:end).*basal_vec)/sum(state_probs_coop(2:end)));
  rate_coop_off_temp = rate_coop_off_temp*af;
  
  % generator matrix for ind
  c_off_slice = zeros(n_states);
  c_off_slice(m1&~m2) = rate_coop_off_temp(1:end-1,2);
  c_off_slice(m3&~m4) = rate_coop_off_temp(2:end,1);
  c_off_slice(m5) = -sum(c_off_slice);
  Q_coop_off_array(:,:,eb) = c_off_slice';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Conduct simulations
n_sim = 25; % number of simulations
T = 1e4; % total time to simulate in seconds

% %%%%%%%%%%%%%%%%%%%% build model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(123);

[unbinding_vec_off, binding_vec_off] = microscopic_binding_sim(Q_coop_off_array',P_coop_array,n_sim,T);

[unbinding_vec_ind, binding_vec_ind] = microscopic_binding_sim(Q_ind_array',P_ind_array,n_sim,T);

%%
ub_bins = 0:0.1:100;
ub_centers = ub_bins(1:end-1)+diff(ub_bins)/2;

c_vec_off = histcounts(unbinding_vec_off,ub_bins);
c_vec_on = histcounts(unbinding_vec_ind,ub_bins);
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

% make exponential fit figure first
exp_fig = figure;
cmap = brewermap(9,'Set2');

hold on
% bar(ub_centers,p_vec_off,1,'FaceColor',cmap(3,:),'FaceAlpha',.9,'EdgeAlpha',1,'EdgeColor','k');   

% pd_vec = params(1)*exp(-[0.01 ub_centers]./params(2))/sum(c_vec_off);
plot(ub_centers,p_trend_on,'Color',cmap(3,:),'LineWidth',3);
plot(ub_centers,p_trend_off,'Color',cmap(5,:),'LineWidth',3);


xlim([0 15])
legend('independent binding','cooperative binding')
xlabel('activator dwell times (s)')
ylabel('probability')
box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
% set(gca,'Color','w')
exp_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(exp_fig,['activator_dwell_exp_coop_vs_ind.png'])
saveas(exp_fig,['activator_dwell_exp_coop_vs_ind.pdf'])


%% now try to capture non-exponential behavior
close all

% log_bins = logspace(-4,log10(max(unbinding_vec)),3e2);
% log_centers = log_bins(1:end-1)+diff(log_bins)/2;

f_off_log = fit(log10(ub_centers'),p_vec_off','smoothingspline','SmoothingParam',0.9);
p_trend_off_log = feval(f_off_log,log10(ub_centers));
f_on_log = fit(log10(ub_centers'),p_vec_on','smoothingspline','SmoothingParam',0.9);
p_trend_on_log = feval(f_on_log,log10(ub_centers));


pwr_fig = figure;
cmap = brewermap(9,'Set2');
hold on
% scatter(log_centers,log_counts)
% b = bar(log_prob,1,'FaceColor',cmap(3,:),'FaceAlpha',1,'EdgeAlpha',1); 

plot(ub_centers,p_trend_on_log,'Color',cmap(3,:),'LineWidth',3);
plot(ub_centers,p_trend_off_log,'Color',cmap(5,:),'LineWidth',3);


set(gca,'Yscale','log')
set(gca,'Xscale','log')

% set(gca,'Xtick',xTickPoints,'xTickLabels',xTickLabelsStr)

ylim([ 1e-5 1e-1])

ax = gca;
ax.XLim = [min(ub_centers),1e2];
xlabel('activator dwell times (s)')
ylabel('probability')
legend('independent binding','cooperative binding','Location','southwest')
set(gcf,'Color','w')
box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
pwr_fig.InvertHardcopy = 'off';

% lgd.FontSize=12;

saveas(pwr_fig,['activator_dwell_pwr_coop_vs_ind.png'])
saveas(pwr_fig,['activator_dwell_pwr_coop_vs_ind.pdf'])

