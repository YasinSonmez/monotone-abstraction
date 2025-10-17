% get abstraction
[con, veh_turn_abs] = get_3D_monotone_abs_v4();
% make u_res more coarse
% veh_turn_abs.set_u_res(300);

% ego position
s_min = veh_turn_abs.x_range(1,1);
s_max = veh_turn_abs.x_range(1,2);
s_res = veh_turn_abs.x_res(1);
% ego velocity
v_min = veh_turn_abs.x_range(2,1);
v_max = veh_turn_abs.x_range(2,2);
v_res = veh_turn_abs.x_res(2);
% oncoming vehicle position
r_min = veh_turn_abs.x_range(3,1);
r_max = veh_turn_abs.x_range(3,2);
r_res = veh_turn_abs.x_res(3);

% size of state-space grid
x_num_cells = veh_turn_abs.x_num_cells;

% initialize feasible_inputs cell
tic;
feasible_inputs = cell(veh_turn_abs.x_num_cells');
num_feasible_inputs = prod(x_num_cells) * veh_turn_abs.u_num_cells;
for s_idx = 1:x_num_cells(1)
    for v_idx = 1:x_num_cells(2)
        for r_idx = 1:x_num_cells(3)
            
            disp(['Initializing inputs for (s_idx, v_idx, r_idx) = (', ...
                    num2str(s_idx), ', ', ...
                    num2str(v_idx), ', ', ...
                    num2str(r_idx), ')']);
            
            % get state index
            x_idx = [s_idx; v_idx; r_idx];
            cell_idx = num2cell(x_idx);
            
            % add all inputs to be tested
            feasible_inputs{cell_idx{:}} = 1:veh_turn_abs.u_num_cells;
            
        end
    end
end

%%% INVARIANT SET ALGORITHM %%%
% iterate over all states and do the following:
% 1) check if any inputs lead to unsafe states
% 2) remove inputs that lead to unsafe states from feasible input list
iter = 0;
num_unsafe_states = 0;
while(iter == 0 || removed_inputs)
    
    % tracking fixed-point iterations
    removed_inputs = 0;
    iter = iter + 1;
    
    for s_idx = x_num_cells(1):-1:1
        for v_idx = x_num_cells(2):-1:1
            for r_idx = x_num_cells(3):-1:1
                
                % get state index
                x_max_idx = [s_idx; v_idx; r_idx];
                cell_idx = num2cell(x_max_idx);
                % get max / min state values
                x_max_val = ...
                    veh_turn_abs.get_priority_state_at_idx(x_max_idx);
                x_min_val = [max(x_max_val(1) - s_res, s_min);
                             max(x_max_val(2) - v_res, v_min);
                             min(x_max_val(3) + r_res, r_max)];
                
                % create list of unsafe inputs
                unsafe_inputs = [];
                
                % check all feasible inputs
                for i = 1:length(feasible_inputs{cell_idx{:}})
                    
                    % get input value
                    u_idx = feasible_inputs{cell_idx{:}}(i);
                    u_val = veh_turn_abs.get_priority_input_at_idx(u_idx);
                    
                    % print current state & input index
                    clc;
                    disp(['Fixed-point algorithm iteration: ', ...
                        num2str(iter)]);
                    disp(['Test state: (h_idx, v_idx, v_L_idx) = (', ...
                        num2str(s_idx), ', ', ...
                        num2str(v_idx), ', ', ...
                        num2str(r_idx), ')']);
                    disp(['Test input: (u_idx) = ', num2str(u_idx)]);
                    disp(['Number of feasible inputs: ', ...
                        num2str(num_feasible_inputs - removed_inputs)]);
                    disp(['Unsafe inputs removed: ', ...
                        num2str(removed_inputs)]);
                    disp(['Unsafe states removed: ', ...
                        num2str(num_unsafe_states)]);
                    
                    % compute max / min transitions
                    x_next_max_val = ...
                        veh_turn_expr_v4_no_loops(x_max_val, u_val, zeros(0,1));
                    x_next_min_val = ...
                        veh_turn_expr_v4_no_loops(x_min_val, u_val, zeros(0,1));
                    % if both states reached the goal set,
                    % we are safe no matter what
                    if(x_next_max_val(3) > 10 && x_next_min_val(3) > 10)
                        continue;
                    end
                    % do thresholding if only the min state reached the
                    % goal set
                    x_next_min_val(3) = min(x_next_min_val(3), 10);
                    % get max / min state indices
                    x_next_max_idx = ...
                        veh_turn_abs.get_state_idx(x_next_max_val);
                    x_next_min_idx = ...
                        veh_turn_abs.get_state_idx(x_next_min_val);
                    
                    if(sum(x_next_min_idx == -1) || sum(x_next_max_idx == -1))
                        % we can transition out of the domain
                        % remember to remove current input from list
                        unsafe_inputs = [unsafe_inputs i];
                        removed_inputs = removed_inputs + 1;
                        continue;
                    end
                    
                    % this should always hold
                    assert(sum(x_next_min_idx <= x_next_max_idx) == 3);
                    
                    % look for any unsafe states in the set of successors
                    found_unsafe = false;
                    for s_next_idx = x_next_min_idx(1):x_next_max_idx(1)
                        for v_next_idx = x_next_min_idx(2):x_next_max_idx(2)
                            for r_next_idx = x_next_min_idx(3):x_next_max_idx(3)
                                
                                % get index for the successor cell
                                cell_next_idx = ...
                                    num2cell([s_next_idx; v_next_idx; r_next_idx]);
                                
                                if(isempty(feasible_inputs{cell_next_idx{:}}))
                                    % successor is an unsafe state
                                    found_unsafe = true;
                                    break;
                                end
                                
                            end
                        end
                    end
                    
                    if(found_unsafe)
                        % there is an unsafe successor
                        % remember to remove current input from list
                        unsafe_inputs = [unsafe_inputs i];
                        removed_inputs = removed_inputs + 1;
                        continue;
                    end
                
                end
                
                % remove unsafe inputs
                if(~isempty(unsafe_inputs))
                    feasible_inputs{cell_idx{:}}(unsafe_inputs) = [];
                end
                
                % if no feasible inputs remain, the state is unsafe
                if(~isempty(unsafe_inputs) && ...
                    isempty(feasible_inputs{cell_idx{:}}))
                    num_unsafe_states = num_unsafe_states + 1;
                end
                
            end
        end
    end
    
    % update number of remaining feasible inputs
    num_feasible_inputs = num_feasible_inputs - removed_inputs;
    
end
toc;
