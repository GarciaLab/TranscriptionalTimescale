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

% load data
load([DataPath 'bursting_sim_struct.mat'])
load([DataPath 'bursting_chain_calc_struct.mat'])
coop_index = 2;
% define colors
blue = [190 201 224]/255;
pboc = [228,221,209]/255;
% extract cooperativity value vector
ec_vec = bursting_chain_calc_struct(coop_index).coopEnergies;
omega_vec = exp(-ec_vec);
ec_filter = ec_vec <=0; % select for only synergistic interactions

% calculate amount of time neaded to go from 6 to 4 as a function of
% cooperativity strength
Q_coop_on_array = bursting_chain_calc_struct(coop_index).Q;

% initialize arrays 
eff_64_coop_vec = NaN(1,size(Q_coop_on_array,3));
calc_vec = [1 n_states];

for ec = 1:length(ec_vec)     
  [~, eff_64_coop_vec(ec)] = pt_solve(Q_coop_on_array(:,:,ec),calc_vec(1),calc_vec(2)); 
end


et64_fig = figure;
cmap = brewermap([],'Set2');
hold on

yyaxis left
plot(omega_vec(ec_filter),eff_64_coop_vec(ec_filter)/60,'-','Color',cmap(2,:),'LineWidth',2.5)
ylabel('time to escape from 6 to 0 (minutes)')
ax = gca;
ax.YColor = cmap(2,:);

yyaxis right
plot(omega_vec(ec_filter),bursting_chain_calc_struct(coop_index).SS(7,ec_filter),'-','Color',...
  cmap(3,:),'LineWidth',2.5)
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