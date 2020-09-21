% script to make burst model calculations appendix figure
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 2:6;
project = 'n_binding_site_comparisons';


% set paths
FigurePath = ['../fig/emergent_bursting/' project '/'];
mkdir(FigurePath)
DataPath = ['../out/emergent_bursting/' project '/'];

% load data
master_struct = struct;
for i = 1:length(n_bcd_sites)
  loadPath = ['../out/emergent_bursting/n' num2str(n_bcd_sites(i)) '/'];
  load([loadPath 'bursting_chain_calc_struct.mat'])
  master_struct(i).bursting_chain_calc_struct = bursting_chain_calc_struct;
  master_struct(i).n_sites = n_bcd_sites(i);
end

% define colors
pboc = [228,221,209]/255;
purple = brighten([171 133 172]/255,.5);
blue = [190 201 224]/255;
red = [246 141 100]/255;
green = [203 220 170]/255;
cmap0 = [green ; blue ;red];

% specify index to show
coop_index = 2;
n_index = find(n_bcd_sites==6);

% cooperativity energy
ec_vec = master_struct(n_index).bursting_chain_calc_struct(coop_index).coopEnergies;
omega_vec = exp(-ec_vec);
ec_filter = ec_vec <=0; % select for only synergistic interactions

ss6_fig = figure;
hold on
cmap = brewermap(14,'Blues');
color_ind = [101 81 31 11 31 81 101];

plot(omega_vec(ec_filter),sum(master_struct(n_index).bursting_chain_calc_struct(coop_index)....
  .SS([1 7],ec_filter)),'-','Color',cmap(end,:),'LineWidth',2.5)
plot(omega_vec(ec_filter),sum(master_struct(n_index).bursting_chain_calc_struct(coop_index)....
  .SS(2:6,ec_filter)),'-','Color',cmap(3,:),'LineWidth',2.5)
 
ylabel('state probability')
% ylim(ylimTrace)
xlim([1 7])
xlabel('cooperativity strength (\omega)')
legend('end states (0,6)', 'intermediate states (1-5)','Location','southwest')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = pboc;
% StandardFigurePBoC(p,gca);
ss6_fig.InvertHardcopy = 'off';
saveas(ss6_fig,[FigurePath 'ss6_vs_omega.png'])
saveas(ss6_fig,[FigurePath 'ss6_vs_omega.pdf'])

et_fig = figure;
hold on
cmap2 = brewermap(12,'YlOrRd');
lgd_str = {};
% define quantities to plot
for n = 1:length(n_bcd_sites)
  konEff = master_struct(n).bursting_chain_calc_struct(coop_index).eff_rates(:,1);
  koffEff = master_struct(n).bursting_chain_calc_struct(coop_index).eff_rates(:,2);
  kappa_vec =  (konEff+koffEff)./(konEff.*koffEff) / 60;

  plot(omega_vec(ec_filter),kappa_vec(ec_filter),'Color',cmap2(2*n,:),'LineWidth',2)
  
  lgd_str = [lgd_str{:} {[num2str(n_bcd_sites(n)) ' binding sites']}];
end
set(gca,'yScale','log')
set(gca,'YColor','k')
ylabel('bursting timescale (minutes)')
% ylim([0.2 20])
xlim([1 7])
xlabel('cooperativity strength (\omega)')
legend(lgd_str{:},'Location','northwest')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = pboc;
% StandardFigurePBoC(p,gca);
et_fig.InvertHardcopy = 'off';
saveas(et_fig,[FigurePath 'et_vs_omega.png'])
saveas(et_fig,[FigurePath 'et_vs_omega.pdf'])


ss_fig = figure;
hold on
lgd_str = {};
% define quantities to plot
for n = 1:length(n_bcd_sites)
  plot(omega_vec(ec_filter),sum(master_struct(n).bursting_chain_calc_struct(coop_index)....
  .SS([1 end],ec_filter)),'-','Color',cmap2(2*n,:),'LineWidth',2.5)
  
  lgd_str = [lgd_str{:} {[num2str(n_bcd_sites(n)) ' binding sites']}];
end
% set(gca,'yScale','log')
set(gca,'YColor','k')
ylabel('probability of end states (0,N)')
% ylim([0.2 20])
xlim([1 7])
xlabel('cooperativity strength (\omega)')
legend(lgd_str{:},'Location','southeast')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = pboc;
% StandardFigurePBoC(p,gca);
ss_fig.InvertHardcopy = 'off';
saveas(ss_fig,[FigurePath 'ss_vs_omega.png'])
saveas(ss_fig,[FigurePath 'ss_vs_omega.pdf'])


