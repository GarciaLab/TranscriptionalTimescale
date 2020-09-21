% script to simulate simple models for connection between bcd binding and
% transcriptional bursting
clear
close all
FigPath = ['../fig/bcd_hb_timscale_sims/'];
mkdir(FigPath)
% specify basic parameters
t_elong = 140; 
rnap_init_rate = [realmin 1/3];
bcd_off_rate = 1/3;
bcd_on_rate = 1/3;
t_sim = 40*60;
state_options = [1 2];

t_vec = 0:t_sim;
MS2_kernel = ones(1,t_elong);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model 1: 1 binding site all or nothing
kon1 = bcd_on_rate;
koff1 = bcd_off_rate;
R1 = [-kon1, koff1; kon1, -koff1];
% simulate 
state_val_vec1 = [1];
jump_time_vec1 = [];
init_time_vec1 = [];
t_curr = 0;
while t_curr < t_sim    
    % randomly select time til next state jump
    state_curr = state_val_vec1(end);    
    next_jump = exprnd(-1/R1(state_curr,state_curr));    
    % sample initiation events    
    t_init = t_curr;
    t_next = min([t_curr + next_jump,t_sim]);
    while t_init < t_next
        next_init = exprnd(1/rnap_init_rate(state_curr));
        t_init = t_init + next_init;
        if t_init < t_next            
            init_time_vec1 = [init_time_vec1 t_init];
        end
    end    
    if t_next < t_sim
        state_val_vec1 = [state_val_vec1 state_options(state_options~=state_curr)];
        jump_time_vec1 = [jump_time_vec1 t_next];
    end
    t_curr = t_next;
end
% calculate observed MS2 signal 
init_vec_reg = histcounts(init_time_vec1,t_vec);
MS2_vec_reg = conv(init_vec_reg,MS2_kernel,'full');
MS2_vec_reg = MS2_vec_reg(1:numel(t_vec));
% make plot
fig1 = figure;
hold on
plot(t_vec/60,MS2_vec_reg,'Color','black','LineWidth',1.5)
y_lim = ylim;
plot([init_time_vec1; init_time_vec1]/60,[repelem(y_lim(1),numel(init_time_vec1)); repelem(y_lim(2),numel(init_time_vec1))],...
    'Color',[0.8500, 0.3250, 0.0980,.1])

xlabel('time (minutes)')
ylabel('number of elongating RNAP')
set(gca,'Fontsize',14)
saveas(fig1,[FigPath 'model01_sim_1site.png'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model 2: 6 binding sites any
state_options2 = 1:7;
kon_vec2 = fliplr(state_options2)*kon1 - kon1;
koff_vec2 = state_options2*koff1 - koff1;
kappa_vec2 = kon_vec2+koff_vec2;
init_vec2 = [realmin repelem(rnap_init_rate(2),6)];
% simulate 
state_val_vec2 = [2];
jump_time_vec2 = [];
init_time_vec2 = [];
t_curr = 0;
while t_curr < t_sim        
    % randomly select time til next state jump
    state_curr = state_val_vec1(end);    
    next_jump = exprnd(1/kappa_vec2(state_curr));  
    state_options2 = [state_curr-1 state_curr+1];
    next_state = randsample(state_options2,1,true,[koff_vec2(state_curr) kon_vec2(state_curr)]);
    % sample initiation events    
    t_init = t_curr;
    t_next = min([t_curr + next_jump,t_sim]);
    while t_init < t_next
        next_init = exprnd(1/init_vec2(state_curr));
        t_init = t_init + next_init;
        if t_init < t_next            
            init_time_vec2 = [init_time_vec2 t_init];
        end
    end    
    if t_next < t_sim
        state_val_vec2 = [state_val_vec2 next_state];
        jump_time_vec2 = [jump_time_vec2 t_next];
    end
    t_curr = t_next;
end
% calculate observed MS2 signal 
init_vec_reg2 = histcounts(init_time_vec2,t_vec);
MS2_vec_reg2 = conv(init_vec_reg2,MS2_kernel,'full');
MS2_vec_reg2 = MS2_vec_reg2(1:numel(t_vec));
% make plot
fig2 = figure;
hold on
plot(t_vec/60,MS2_vec_reg2,'Color','black','LineWidth',1.5)
y_lim = ylim;
plot([init_time_vec2; init_time_vec2]/60,[repelem(y_lim(1),numel(init_time_vec2)); repelem(y_lim(2),numel(init_time_vec2))],...
    'Color',[0.8500, 0.3250, 0.0980,.1])

xlabel('time (minutes)')
ylabel('number of elongating RNAP')
set(gca,'Fontsize',14)
saveas(fig2,[FigPath 'model02_sim_6sites_any.png'])
    

%% model 3: 6 binding sites all
state_options2 = 1:7;
kon_vec2 = fliplr(state_options2)*kon1 - kon1;
koff_vec2 = state_options2*koff1 - koff1;
kappa_vec2 = kon_vec2+koff_vec2;
init_vec2 = [repelem(realmin,6) rnap_init_rate(2)];
% simulate 
state_val_vec2 = [2];
jump_time_vec2 = [];
init_time_vec2 = [];
t_curr = 0;
while t_curr < t_sim        
    % randomly select time til next state jump
    state_curr = state_val_vec1(end);    
    next_jump = exprnd(1/kappa_vec2(state_curr));  
    state_options2 = [state_curr-1 state_curr+1];
    next_state = randsample(state_options2,1,true,[koff_vec2(state_curr) kon_vec2(state_curr)]);
    % sample initiation events    
    t_init = t_curr;
    t_next = min([t_curr + next_jump,t_sim]);
    while t_init < t_next
        next_init = exprnd(1/init_vec2(state_curr));
        t_init = t_init + next_init;
        if t_init < t_next            
            init_time_vec2 = [init_time_vec2 t_init];
        end
    end    
    if t_next < t_sim
        state_val_vec2 = [state_val_vec2 next_state];
        jump_time_vec2 = [jump_time_vec2 t_next];
    end
    t_curr = t_next;
end
% calculate observed MS2 signal 
init_vec_reg2 = histcounts(init_time_vec2,t_vec);
MS2_vec_reg2 = conv(init_vec_reg2,MS2_kernel,'full');
MS2_vec_reg2 = MS2_vec_reg2(1:numel(t_vec));
% make plot
fig2 = figure;
hold on
plot(t_vec/60,MS2_vec_reg2,'Color','black','LineWidth',1.5)
y_lim = ylim;
plot([init_time_vec2; init_time_vec2]/60,[repelem(y_lim(1),numel(init_time_vec2)); repelem(y_lim(2),numel(init_time_vec2))],...
    'Color',[0.8500, 0.3250, 0.0980,.1])

xlabel('time (minutes)')
ylabel('number of elongating RNAP')
set(gca,'Fontsize',14)
saveas(fig2,[FigPath 'model03_sim_6sites_all.png'])
    