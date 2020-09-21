% script to simulate burst time series with parameters from numerical and
% analytic results
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
n_states = n_bcd_sites+1;
n_bound_vec = 0:n_bcd_sites;
n_calc_points = 201;
off_rate_basal = 1/2; % in seconds. This sets overall system timescales (from Mir et al, 2018)

project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
FigurePath = ['../fig/emergent_bursting/' project '/'];
mkdir(FigurePath)
DataPath = ['../out/emergent_bursting/' project '/'];

%% calculate transition matrices and state probabilities from scratch
mult_vec = NaN(size(n_bound_vec));
for w = 1:length(mult_vec)
  mult_vec(w) = factorial(n_bcd_sites) ./ (factorial(n_bcd_sites-n_bound_vec(w)) .*factorial(n_bound_vec(w)));
end
mu_vec = -log(mult_vec); % energy contribution from state multiplicities

coopEnergies = -linspace(0,1.9,201); 
activatorEnergyVec = repelem(4.75,n_calc_points); % fix to value used for coop plot %-0.5*coopEnergies*(n_bcd_sites-1);

stateEnergy_fun = @(eb,ec) eb*n_bound_vec + ec*n_bound_vec.*(n_bound_vec-1)/2 + mu_vec;
stateProb_fun = @(eb,ec) exp(-stateEnergy_fun(eb,ec)) / sum(exp(-stateEnergy_fun(eb,ec)));

% for each energy landscape, we can calculate the implied on rate for a 
% single Bcd molecule (absent cooperative effects) unsing known koff and
% detailed balance
on_rate_basal_vec = NaN(1,n_calc_points); 

% initialize 3D array to store transition rate (generator) matrices for each
% condition
Q_coop_on_array = zeros(n_states,n_states,n_calc_points); % on rate-mediated

P_coop_array = NaN(n_states,n_calc_points);
% define matrices used for indexing
a = ones(n_states); m1 = tril(a,-1); m2 = tril(a,-2); m3 = triu(a,1); m4 = triu(a,2); m5 = ~~eye(n_states);

% calculate metrics for (a) and (b) first
for ec = 1:length(activatorEnergyVec)  
  eBind = activatorEnergyVec(ec);
  eCoop = coopEnergies(ec);
  
  % state probabilities
  state_probs_coop = stateProb_fun(eBind,eCoop);    
  P_coop_array(:,ec) = state_probs_coop;
  % calculate basal on rate using detailed balance  
  on_rate_basal_vec(ec) = off_rate_basal * state_probs_coop(2) / state_probs_coop(1) / n_bcd_sites;
  % calculate implied rates
  rate_coop_on_temp(:,1) = off_rate_basal * n_bound_vec;
  rate_coop_on_temp(1:end-1,2) = rate_coop_on_temp(2:end,1)' .* state_probs_coop(2:end) ./ state_probs_coop(1:end-1);

  % generator matrix 
  c_on_slice = zeros(n_states);
  c_on_slice(m1&~m2) = rate_coop_on_temp(1:end-1,2);
  c_on_slice(m3&~m4) = rate_coop_on_temp(2:end,1);
  c_on_slice(m5) = -sum(c_on_slice);
  Q_coop_on_array(:,:,ec) = c_on_slice';

end


%% Generate figures 
% define colors
blue = [190 201 224]/255;
pboc = [228,221,209]/255;

% extract cooperativity value vector
omega_vec = exp(-coopEnergies);

% initialize arrays 
eff_64_coop_vec = NaN(1,size(Q_coop_on_array,3));
eff_63_coop_vec = NaN(1,size(Q_coop_on_array,3));
eff_02_coop_vec = NaN(1,size(Q_coop_on_array,3));
calc_vec = [1 n_states];

for ec = 1:length(coopEnergies)     
  [~, eff_64_coop_vec(ec)] = pt_solve(Q_coop_on_array(:,:,ec),calc_vec(1),calc_vec(2)); 
%   [~, eff_63_coop_vec(ec)] = pt_solve(Q_coop_on_array(:,:,ec),calc_vec(1)-1,calc_vec(2)); 
%   [~, eff_64_coop_vec(ec)] = pt_solve(Q_coop_on_array(:,:,ec),calc_vec(1),calc_vec(2)); 
end


et64_fig = figure;
cmap = brewermap([],'Set2');
hold on

yyaxis left
% plot(omega_vec,eff_64_coop_vec,'-','Color',cmap(2,:),'LineWidth',2.5)
plot(omega_vec,eff_64_coop_vec,'-','Color',cmap(2,:),'LineWidth',2.5)
set(gca,'Yscale','log')
ylabel('passage time from 6 to 0 (seconds)')
ax = gca;
ax.YColor = cmap(2,:);

yyaxis right
plot(omega_vec,P_coop_array(end,:),'-','Color',cmap(3,:),'LineWidth',2.5)
ylabel('probability of state 6')
ax = gca;
ax.YColor = cmap(3,:);
% ylim(ylimTrace)
xlim([1 7])
xlabel('cooperativity strength (\omega)')
box on
set(gca,'Fontsize',14)
% set(gca,'Yscale','log')
ax = gca;

ax.XColor = 'black';
% legend('
ax.Color = pboc;
% StandardFigurePBoC(p,gca);
et64_fig.InvertHardcopy = 'off';
saveas(et64_fig,[FigurePath 'et64_vs_omega.png'])
saveas(et64_fig,[FigurePath 'et64_vs_omega.pdf'])