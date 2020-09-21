% script to make burst model calculations appendix figure
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
FigurePath = ['../fig/bursting_appendix/' project '/'];
mkdir(FigurePath)
DataPath = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'bursting_chain_calc_struct.mat'])
load([DataPath 'bursting_step_calc_struct.mat'])

% define colors
pboc = [228,221,209]/255;
purple = brighten([171 133 172]/255,.5);
blue = [190 201 224]/255;
red = [246 141 100]/255;
green = [203 220 170]/255;
cmap0 = [green ; blue ;red];

% key calc quantities
n_bound_vec = 0:6;
mult_vec = NaN(size(n_bound_vec));
for w = 1:length(mult_vec)
  mult_vec(w) = factorial(n_bcd_sites) ./ (factorial(n_bcd_sites-n_bound_vec(w)) .*factorial(n_bound_vec(w)));
end
mu_vec = -log(mult_vec); % energy contribution from state multiplicities

% define functions
stateEnergy_fun = @(eb,ec) eb*n_bound_vec + ec*n_bound_vec.*(n_bound_vec-1)/2 + mu_vec;
stateProb_fun = @(eb,ec) exp(-stateEnergy_fun(eb,ec)) / sum(exp(-stateEnergy_fun(eb,ec)));

% state multiplicity plots

multiplicity_fig = figure('Position',[100 100 512 256]);
hold on
area(n_bound_vec,mult_vec,'FaceColor',purple,'LineWidth',1)
p = plot(0,0);

xlabel('binding state')
ylabel('state multiplicity')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';

StandardFigurePBoC(p,gca);
multiplicity_fig.InvertHardcopy = 'off';
saveas(multiplicity_fig,[FigurePath 'multiplicity_fig.png'])
saveas(multiplicity_fig,[FigurePath 'multiplicity_fig.pdf'])

% independent binding figs
ind_name = 'independent';
chain_calc_names = {bursting_chain_calc_struct.name};
ind_index = find(contains(chain_calc_names,ind_name));

% get energy vectors
eb_energies = bursting_chain_calc_struct(ind_index).activatorEnergies;
plot_index = 1:10:length(eb_energies);

% ind_energy_fig = figure;
% hold on
cmap1 = brewermap(length(plot_index),'RdYlBu');
% 
% for e = 1:length(plot_index)
%   e_vec = stateEnergy_fun(eb_energies(plot_index(e)),0);
%   plot(n_bound_vec,e_vec,'Color',cmap1(e,:))
% end
% 
% p = plot(0,0);
% % ylim(ylimTrace)
% % xlim([0 t_max])
% xlabel('binding state')
% ylabel('state energy (k_bT)')
% box on
% set(gca,'Fontsize',14)
% % set(gca,'yScale','log')
% ax = gca;
% ax.YColor = 'black';
% ax.XColor = 'black';
% % ax.Color = pboc;
% StandardFigurePBoC(p,gca);
% ind_energy_fig.InvertHardcopy = 'off';
% saveas(ind_energy_fig,[FigurePath 'ind_energies.png'])
% saveas(ind_energy_fig,[FigurePath 'ind_energies.pdf'])
% 
% 
% ind_prob_fig = figure;
% hold on
% cmap1 = brewermap(length(plot_index),'RdYlBu');
% for e = 1:length(plot_index)
%   p_vec = stateProb_fun(eb_energies(plot_index(e)),0);
%   plot(n_bound_vec,p_vec,'Color',cmap1(e,:))
% end
% % p = plot(0,0);
% % ylim(ylimTrace)
% % xlim([0 t_max])
% xlabel('binding state')
% ylabel('probability')
% box on
% set(gca,'Fontsize',14)
% % set(gca,'yScale','log')
% ax = gca;
% ax.YColor = 'black';
% ax.XColor = 'black';
% % ax.Color = pboc;
% StandardFigurePBoC(p,gca);
% ind_prob_fig.InvertHardcopy = 'off';
% saveas(ind_prob_fig,[FigurePath 'ind_probs.png'])
% saveas(ind_prob_fig,[FigurePath 'ind_probs.pdf'])
% 
% 
exi = [1 ceil(length(plot_index)/2) length(plot_index)];

ind_ex_fig1 = figure;
hold on
for e = 1:length(exi)
  e_vec = stateEnergy_fun(eb_energies(plot_index(exi(e))),0);
  plot(n_bound_vec,e_vec,'Color',cmap0(e,:),'LineWidth',3)
end
p = plot(0,0);
% ylim(ylimTrace)
% xlim([0 t_max])
legend('positive \epsilon_b', 'negligible \epsilon_b','negative \epsilon_b','Location','northwest')
xlabel('binding state')
ylabel('energy (k_bT)')
box on
set(gca,'Fontsize',14)
% set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
% ax.Color = pboc;
StandardFigurePBoC(p,gca);
ind_ex_fig1.InvertHardcopy = 'off';
saveas(ind_ex_fig1,[FigurePath 'ind_example_energies.png'])
saveas(ind_ex_fig1,[FigurePath 'ind_example_energies.pdf'])


ind_ex_fig2 = figure;
hold on
for e = 1:length(exi)
  e_vec = stateProb_fun(eb_energies(plot_index(exi(e))),0);
  plot(n_bound_vec,e_vec,'Color',cmap0(e,:),'LineWidth',3)
end
p = plot(0,0);
% ylim(ylimTrace)
% xlim([0 t_max])
xlabel('binding state')
ylabel('probability')
box on
set(gca,'Fontsize',14)
legend('positive \epsilon_b', 'negligible \epsilon_b','negative \epsilon_b','Location','northwest')
% set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
% ax.Color = pboc;
StandardFigurePBoC(p,gca);
ind_ex_fig2.InvertHardcopy = 'off';
saveas(ind_ex_fig2,[FigurePath 'ind_example_probs.png'])
saveas(ind_ex_fig2,[FigurePath 'ind_example_probs.pdf'])



% get index
coop_name = 'kon-mediated';
coop_index = find(contains(chain_calc_names,coop_name));

% get energy vectors
eb_energies = bursting_chain_calc_struct(coop_index).activatorEnergies;
ec_energies = bursting_chain_calc_struct(coop_index).coopEnergies;


coop_energy_fig = figure;
hold on
colormap(cmap1);
cmap1 = brewermap(length(plot_index),'RdYlBu');
for e = 1:length(plot_index)
  e_vec = stateEnergy_fun(eb_energies(plot_index(e)),ec_energies(plot_index(e)));
  plot(n_bound_vec,e_vec,'Color',cmap1(e,:))
end
p = plot(0,0);
% ylim(ylimTrace)
% xlim([0 t_max])
xlabel('binding state')
ylabel('state energy (k_bT)')
box on
set(gca,'Fontsize',14)
% set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
% ax.Color = pboc;
StandardFigurePBoC(p,gca);
coop_energy_fig.InvertHardcopy = 'off';
saveas(coop_energy_fig,[FigurePath 'coop_energies.png'])
saveas(coop_energy_fig,[FigurePath 'coop_energies.pdf'])


coop_prob_fig = figure;
hold on
cmap1 = brewermap(length(plot_index),'RdYlBu');
for e = 1:length(plot_index)
  p_vec = stateProb_fun(eb_energies(plot_index(e)),ec_energies(plot_index(e)));
  plot(n_bound_vec,p_vec,'Color',cmap1(e,:))
end
% p = plot(0,0);
% ylim(ylimTrace)
% xlim([0 t_max])
xlabel('binding state')
ylabel('state probability')
box on
set(gca,'Fontsize',14)
% set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
% ax.Color = pboc;
StandardFigurePBoC(p,gca);
coop_prob_fig.InvertHardcopy = 'off';
saveas(coop_prob_fig,[FigurePath 'coop_probs.png'])
saveas(coop_prob_fig,[FigurePath 'coop_probs.pdf'])


% Plot binding affinity versus passage times

et_fig = figure;
hold on
cmap2 = brewermap(9,'Set2');
scatter(bursting_chain_calc_struct(coop_index).activatorEnergies,...
  1/60./bursting_chain_calc_struct(coop_index).eff_rates(:,1),'MarkerFaceColor',cmap2(3,:),'MarkerEdgeColor','k')
% ylim(ylimTrace)
xlim([0 5])
xlabel('binding energy (k_bT)')
ylabel('first-passage time (minutes)')
box on
% grid on
set(gca,'Fontsize',14)
set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = pboc;
grid on
% StandardFigurePBoC(p,gca);
et_fig.InvertHardcopy = 'off';
saveas(et_fig,[FigurePath 'et_vs_eb.png'])
saveas(et_fig,[FigurePath 'et_vs_eb.pdf'])