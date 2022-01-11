function [unbinding_vec, binding_vec] = microscopic_binding_sim(Q,SS,n_sim,T)

n_bs = size(Q,1)-1;

% extract microscopic rates
n_vec = 0:n_bs;
n_states = length(n_vec);
a = ones(n_states); m1 = tril(a,-1); m2 = tril(a,-2); m3 = triu(a,1); m4 = triu(a,2); m5 = ~~eye(n_states);

% effective rates
k_minus_vec = [0 Q(m3&~m4)'];
k_plus_vec = [Q(m1&~m2)' 0];

% % extract microscopic rates
% k_unbind_vec = k_minus_vec ./ n_vec;
% k_unbind_vec(1) = 0;
% k_bind_vec = k_plus_vec ./(n_bs-n_vec);
% k_bind_vec(end) = 0;

% cell arrays to store results
event_time_cell = cell(1,n_sim);
bound_state_cell = cell(1,n_sim);

for n = 1:n_sim
  % initialize
  event_time_vec = [0];
  
  nBound = randsample(n_vec,1,true,SS);
  s_bound = randsample(1:n_bs,nBound,false);
  bound_state_array = zeros(1,n_bs);
  bound_state_array(s_bound) = 1;
  
  simTime = 0; 
  
  while simTime < T
    
    currState = bound_state_array(end,:);
    nBound = sum(currState);
    
    % draw expected jump time
    tau = 1 / (k_minus_vec(nBound+1) + k_plus_vec(nBound+1));
    dt = exprnd(tau);
    
    simTime = simTime + dt;
    
    % determine whether we bind or unbind
    eventType = randsample([-1 1],1,true,[k_minus_vec(nBound+1) k_plus_vec(nBound+1)]);
    
    if eventType == 1      
      options = find(currState==0);      
    else
      options = find(currState==1);
    end  
    % select binding site to update    
    bsSwitch = randsample(repelem(options,2),1,false);% repelem prevents undesirable behavior when only 1 option
    currState(bsSwitch) = currState(bsSwitch) + eventType;    
    bound_state_array(end+1,:) = currState;
    event_time_vec(end+1) = simTime;    
  end
  
  bound_state_cell{n} = bound_state_array;
  event_time_cell{n} = event_time_vec;
  
end

% extract single-molecules unbinding times
unbinding_vec = [];
binding_vec = [];
for n = 1:n_sim
  bound_state_array = bound_state_cell{n};
  event_time_vec = event_time_cell{n};
  
  for i = 1:n_bs
    binding_state_vec = bound_state_array(:,i);
    
    change_vec = [0 diff(binding_state_vec')];
    
    binding_times_raw = event_time_vec(change_vec==1);
    unbinding_times_raw = event_time_vec(change_vec==-1);
     
    if ~isempty(binding_times_raw) && ~isempty(unbinding_times_raw)            
      binding_times1 = binding_times_raw(binding_times_raw<unbinding_times_raw(end));
      unbinding_times1 = unbinding_times_raw(unbinding_times_raw>binding_times_raw(1));
      
      unbinding_vec = [unbinding_vec unbinding_times1-binding_times1];
      
      binding_times2 = binding_times_raw(binding_times_raw>unbinding_times_raw(1));
      unbinding_times2 = unbinding_times_raw(unbinding_times_raw<binding_times_raw(end));
      
      binding_vec = [binding_vec binding_times2-unbinding_times2];
    end
  end   
    
end