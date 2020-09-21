% Script to pull and plot illustrative Hunchback P2P traces
clear
% add utilties folder
addpath('C:\Users\nlamm\projects\transcription_timescales_review\src\utilities')
% set path to data
DataPath = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\Dl-Ven_hbP2P-mCh_v2\';
% set write path
FigPath = 'C:\Users\nlamm\projects\transcription_timescales_review\fig\';
mkdir(FigPath);
% load data
load([DataPath 'hmm_input_output_w6_K3.mat'])

%% Plot illustrative trace
close all

% specify trace to use
trace_id = 679;
% extract vectors to plot
time_vec = hmm_input_output(trace_id).time/60;
fluo_vec = hmm_input_output(trace_id).fluo;
hmm_vec = hmm_input_output(trace_id).z_vec;

fluo_fig = figure('Position',[0 0 512 256]);
cmap = brewermap([],'Set2');
hold on
p = plot(0,0);
plot(time_vec,fluo_vec,'-','Color','black','LineWidth',1.5');    
% e.CapSize = 0;
scatter(time_vec,fluo_vec,'MarkerFaceColor','black','MarkerEdgeAlpha',0)
xlabel('time (minutes)')
ylabel('{\it hbP2P} spot fluorescence (AU)')
box on
xlim([10,max(time_vec)])
StandardFigure(p,gca);
saveas(fluo_fig,[FigPath 'fluo_trend.pdf'])
saveas(fluo_fig,[FigPath 'fluo_trend.png'])

% now hmm fit
fluo_hmm_fig = figure('Position',[0 0 512 256]);
hold on
p = plot(0,0);
sh = stairs(time_vec,hmm_vec>1,'-','Color',cmap(2,:),'LineWidth',1.5);  
% fill in bursts
bottom = min(sh.YData);
x = [sh.XData(1),repelem(sh.XData(2:end),2)];
y = [repelem(sh.YData(1:end-1),2),sh.YData(end)];
fill([x,fliplr(x)],[y,bottom*ones(size(y))], cmap(2,:))

xlabel('time (minutes)')
ylabel('promoter state')
box on
xlim([10,max(time_vec)])
ylim([0 1.1])
StandardFigure(p,gca);
saveas(fluo_hmm_fig,[FigPath 'hmm_trend.pdf'])
saveas(fluo_hmm_fig,[FigPath 'hmm_trend.png'])

% overlay
overlay_fig = figure('Position',[0 0 512 256]);
hold on
p = plot(0,0);
plt = plot(time_vec,fluo_vec,'-','Color','black','LineWidth',1.5');    
y_lim = get(gca,'Ylim');
top = y_lim(2);
bottom = y_lim(1);
fill([x,fliplr(x)],[top*y,bottom*ones(size(x))], cmap(2,:),'FaceAlpha',0.4,'EdgeAlpha',0)
%     scatter(time,fluo,'MarkerFaceColor','black','MarkerEdgeColor','black')
xlabel('time (minutes)')
ylabel('{\it hbP2P} spot fluorescence (AU)')
box on
xlim([10,max(time_vec)])
ylim(y_lim)
suffix = '_standard';
StandardFigure(p,gca);

saveas(overlay_fig,[FigPath 'fluo_w_hmm_trend.pdf'])
saveas(overlay_fig,[FigPath 'fluo_w_hmm_trend.png'])