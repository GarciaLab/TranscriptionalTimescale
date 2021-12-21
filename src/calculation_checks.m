%Script to check consistency between matrix-based state probability 
%calculations and formulas cited in appendix
clear
close all

addpath('utilities')

% load numeric results
project = 'n6';

% set paths
DataPath = ['../out/emergent_bursting/' project '/'];

load([DataPath 'bursting_chain_calc_struct.mat'])


N = 6;
n_vec = 0:N;

%%%%%%%%%%%%%%%%%%%
% independent binding model checks
%%%%%%%%%%%%%%%%%%%
ind_index = 1;
Q_mat_ind = bursting_chain_calc_struct(ind_index).Q;

ss_array_ind_matrix = NaN(size(Q_mat_ind,1),size(Q_mat_ind,3));
ss_array_ind_manuscript = NaN(size(Q_mat_ind,1),size(Q_mat_ind,3));

% precalculate N weights
WN = [];
for i = 1:length(n_vec)
  WN(i) = factorial(N)/(factorial(N-n_vec(i))*factorial(n_vec(i)));
end

for i = 1:size(Q_mat_ind,3)
  Q = Q_mat_ind(:,:,i)';
  
  % matrix-based approach
  [V,D] = eig(Q);
  [~,mi] = max(diag(D));
  ss_array_ind_matrix(:,i) = V(:,mi)/sum(V(:,mi));
  
  % manuscript approach
  kon = Q(2,1)/N;
  koff = Q(1,2);
  K_d = koff/kon;
  state_weights = WN.*K_d.^-n_vec;
  ss_array_ind_manuscript(:,i) = state_weights / sum(state_weights);
end
all_ind_match = all(round(ss_array_ind_manuscript(:),3)==round(ss_array_ind_matrix(:),3))


%%%%%%%%%%%%%%%%%%%
% cooperative binding model checks
%%%%%%%%%%%%%%%%%%%
coop_index = 2;
Q_mat_coop = bursting_chain_calc_struct(coop_index).Q;
omega_vec = exp(-bursting_chain_calc_struct(coop_index).coopEnergies);

choose_2 = n_vec.*(n_vec-1)/2;

ss_array_coop_matrix = NaN(size(Q_mat_coop,1),size(Q_mat_coop,3));
ss_array_coop_manuscript = NaN(size(Q_mat_coop,1),size(Q_mat_coop,3));

for i = 197%1:size(Q_mat_coop,3)
  Q = Q_mat_coop(:,:,i)';
  
  % matrix-based approach
  [V,D] = eig(Q);
  [~,mi] = max(diag(D));
  ss_array_coop_matrix(:,i) = V(:,mi)/sum(V(:,mi));
  
  % manuscript approach
  kon = Q(2,1)/N;
  koff = Q(1,2);
  K_d = koff/kon;
  state_weights = WN.*K_d.^-n_vec.*omega_vec(i).^choose_2;
  ss_array_coop_manuscript(:,i) = state_weights / sum(state_weights);
end
all_coop_match = all(round(ss_array_coop_manuscript(:),3)==round(ss_array_coop_matrix(:),3))