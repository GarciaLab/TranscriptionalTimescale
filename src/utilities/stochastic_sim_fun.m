function [sim_state_cell, sim_emission_cell, sim_time_cell] = stochastic_sim_fun(Q,ss_vec,emission_vec,n_sim,t_sim)
  % ARGUMENTS
% Q: N+1xN+1 Transition rate matrix for the system 
% ss_vec: State probability vec. Used to initialize simulations
% emission_vec: 1 X N+1 vector giving transcription rates for each state
% n_sim: number of independent simulations to run
% t_sim: duration to simulate


% OUTPUT
% sim_state_cell : 1 X n_sim cell array with each entry containing state
%                 sequence for corresponding simulation
% sim_emission_cell : 1 X n_sim cell array with each entry containing
%                 transcription rate sequence 
% sim_time_cell : 1 X n_sim cell array with each entry containing
%                 time of jump from one state to the next


  % state options
  state_option_vec = 1:size(Q,1);
  
  % initialize cell arrays
  sim_state_cell = cell(1,n_sim);
  sim_emission_cell = cell(1,n_sim);
  sim_time_cell = cell(1,n_sim);
 
  % iterate
  parfor n = 1:n_sim
      state_val_vec = [randsample(state_option_vec,1,true,ss_vec)];
      jump_time_vec = [0];
      t_curr = 0;
      while t_curr < t_sim            
          state_curr = state_val_vec(end);   % current state  
          next_jump = exprnd(-1/Q(state_curr,state_curr)); % time til next jump                        
          t_next = t_curr + next_jump;

          % randomly choose next state
          weight_vec = Q(state_curr,:);
          weight_vec(state_curr) = 0;
          weight_vec = weight_vec / sum(weight_vec);
          state_next = randsample(state_option_vec,1,true,weight_vec);
          
          % if within alotted time, add state to state vector
          if t_next < t_sim
              state_val_vec = [state_val_vec state_next];
              jump_time_vec = [jump_time_vec t_next];
          end
          t_curr = t_next;
      end
      sim_state_cell{n} = int8(state_val_vec);    
      sim_emission_cell{n} = int8(emission_vec(state_val_vec));         
      sim_time_cell{n} = double(jump_time_vec);
  end  