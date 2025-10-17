% get abstraction
[con, veh_follow_abs] = get_3D_box_abs_v2();

% size of state-space grid
x_num_cells = veh_follow_abs.x_num_cells;

% initialize feasible_inputs cell
feasible_inputs = cell(veh_follow_abs.x_num_cells');
num_feasible_inputs = prod(x_num_cells) * veh_follow_abs.u_num_cells;
for h_idx = 1:x_num_cells(1)
    for v_idx = 1:x_num_cells(2)
        for v_L_idx = 1:x_num_cells(3)
            
            disp(['Initializing inputs for (h_idx, v_idx, v_L_idx) = (', ...
                    num2str(h_idx), ', ', ...
                    num2str(v_idx), ', ', ...
                    num2str(v_L_idx), ')']);
            
            % get state index
            x_idx = [h_idx; v_idx; v_L_idx];
            cell_idx = num2cell(x_idx);
            
            % add all inputs to be tested
            feasible_inputs{cell_idx{:}} = 1:veh_follow_abs.u_num_cells;
            
        end
    end
end

%%% PRECOMPUTATION STEP %%%
% iterate over discretized headways
for h_point = con.h_min:con.h_res:con.h_max
    % iterate over discretized lead vehicle velocities
    for v_lead_point = con.v_min:con.v_res:con.v_max
        clc;
        disp('Precomputation step - removing unsafe states.');
        disp(['Test point: (h, v_L) = (', ...
              num2str(h_point), ', ', ...
              num2str(v_lead_point), ')']);
        v_follow_point = outer_approx_boundary(h_point, ...
                                               v_lead_point, ...
                                               abs(con.a_min_brake), ...
                                               abs(con.a_max_brake), ...
                                               con.h_min);
        % threshold at maximum velocity
        v_follow_point = min(v_follow_point, con.v_max);
        % quantize v_follow to the priority state
        threshold = 1e-5;
        v_follow_point = con.v_min + ceil((v_follow_point - con.v_min) ...
                                   / con.v_res - threshold) * con.v_res;
        % expand the v_follow point upwards until an unsafe collision
        % is possible
        while(v_follow_point <= con.v_max + threshold && ...
              ~unsafe_collision_check(h_point, v_follow_point, ...
                                      v_lead_point, con.a_min_brake, ...
                                      con.a_max_brake, con.v_allow, 1))
            v_follow_point = v_follow_point + con.v_res;
        end
        v_follow_point = v_follow_point - con.v_res;
        % quantize to the last safe state we saw (within range)
        x_bound = [h_point; v_follow_point; v_lead_point];
        x_bound_idx = veh_follow_abs.get_state_idx(x_bound);
        % now remove all unsafe states above this one
        if((x_bound_idx(1) == 20 || x_bound_idx(2) == 21) ...
                && x_bound_idx(3) == 21)
            1;
        end
        while(x_bound_idx(2) < x_num_cells(2))
            x_bound_idx = x_bound_idx + [0; 1; 0];
            cell_idx = num2cell(x_bound_idx);
            feasible_inputs{cell_idx{:}} = zeros(1,0);
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
    
    for h_idx = x_num_cells(1):-1:1
        for v_idx = x_num_cells(2):-1:1
            for v_L_idx = x_num_cells(3):-1:1
                
                % get state index
                x_idx = [h_idx; v_idx; v_L_idx];
                cell_idx = num2cell(x_idx);
                
                % create list of unsafe inputs
                unsafe_inputs = [];
                
                % check all feasible inputs
                for i = 1:length(feasible_inputs{cell_idx{:}})
                    
                    % get input and disturbance index
                    u_idx = feasible_inputs{cell_idx{:}}(i);
                    w_min_idx = 1;
                    w_max_idx = veh_follow_abs.w_num_cells;
                    
                    % print current state & input index
                    clc;
                    disp(['Fixed-point algorithm iteration: ', ...
                        num2str(iter)]);
                    disp(['Test state: (h_idx, v_idx, v_L_idx) = (', ...
                        num2str(h_idx), ', ', ...
                        num2str(v_idx), ', ', ...
                        num2str(v_L_idx), ')']);
                    disp(['Test input: (u_idx) = ', num2str(u_idx)]);
                    disp(['Number of feasible inputs: ', ...
                        num2str(num_feasible_inputs - removed_inputs)]);
                    disp(['Unsafe inputs removed: ', ...
                        num2str(removed_inputs)]);
                    disp(['Unsafe states removed: ', ...
                        num2str(num_unsafe_states)]);
                    
                    % compute min and max transitions
                    x_next_min_idx = ...
                        veh_follow_abs.get_min_transition_state(x_idx, u_idx, w_min_idx);
                    x_next_max_idx = ...
                        veh_follow_abs.get_max_transition_state(x_idx, u_idx, w_max_idx);
                    
                    if(sum(x_next_min_idx == -1) || sum(x_next_max_idx == -1))
                        % we can transition out of the domain
                        % remember to remove current input from list
                        unsafe_inputs = [unsafe_inputs i];
                        removed_inputs = removed_inputs + 1;
                        continue;
                    end
                    
                    % look for any unsafe states in the set of successors
                    found_unsafe = false;
                    for h_next_idx = x_next_min_idx(1):x_next_max_idx(1)
                        for v_next_idx = x_next_min_idx(2):x_next_max_idx(2)
                            for v_L_next_idx = x_next_min_idx(3):x_next_max_idx(3)
                                
                                % get index for the successor cell
                                cell_next_idx = ...
                                    num2cell([h_next_idx; v_next_idx; v_L_next_idx]);
                                
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
