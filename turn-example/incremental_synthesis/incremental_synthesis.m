%%% This script is used to check whether we can safely cut in between %%%
%%% two oncoming vehicles while executing an unprotected left turn. %%%

% set the initial position of the second oncoming vehicle, as well as the
% initial distance between the two oncoming vehicles
s2_init = -70;
s12_dist_init = 70;

% get precomputed abstractions and 
mdl = run_setup_scenario1;

% iterate over the state cells for the ego first abstraction
% (since it uses a finer resolution for velocity values)

% ego position
s_min = mdl.abs_ego_first.x_range(1,1);
s_max = mdl.abs_ego_first.x_range(1,2);
s_res = mdl.abs_ego_first.x_res(1);
% ego velocity
v_min = mdl.abs_ego_first.x_range(2,1);
v_max = mdl.abs_ego_first.x_range(2,2);
v_res = mdl.abs_ego_first.x_res(2);
% oncoming vehicle 2 position
s2_min = mdl.abs_ego_first.x_range(3,1);
s2_max = mdl.abs_ego_first.x_range(3,2);
s2_res = mdl.abs_ego_first.x_res(3);
% oncoming vehicle 1 position
s1_min = mdl.abs_oncoming_first.x_range(3,1);
s1_max = mdl.abs_oncoming_first.x_range(3,2);
% ego wheel torque
T_min = mdl.abs_ego_first.u_range(1);
T_max = mdl.abs_ego_first.u_range(2);
T_res = mdl.abs_ego_first.u_res;

% for checking if we are out of range
s_max_oncoming_first = mdl.abs_oncoming_first.x_range(1,2);

% list of feasible control inputs
% feasible inputs will be indexed with the following convention
% (we will refer to this as a FINE index, since v resolution is finer):
numCells_s = mdl.abs_ego_first.x_numCells(1); % s lower idx
numCells_v = mdl.abs_ego_first.x_numCells(2); % v lower idx
numCells_s2 = mdl.abs_ego_first.x_numCells(3); % s2 lower idx
feasible_inputs = cell(numCells_s, numCells_v, numCells_s2);
num_feasible_inputs = 0;
num_removed_inputs = 0;

% iterate over all states and do the following:
% 1) check for a feasible input
% 2) if some feasible inputs exist, add them to list for current state
num_unsafe_states = 0;
for s_val = s_min:s_res:s_max
    for v_val = v_min:v_res:v_max
        for s2_val = s2_min:s2_res:s2_max
            
            % domain of test points:          s \in [-70m, 10m]
            %                                 v \in [0m/s, 12m/s]
            %                                 s1 \in [0m, 60m]
            %                                 s2 \in [-70m, -10m]
            %
            % print current test point
            s1_val = convert_s2_to_s1(s2_val, s2_init, s12_dist_init, ...
                                        mdl.v0_min, mdl.v0_max);
            clc;
            disp('Feasible inputs initialization step.');
            disp(['Test point: (s, v, s1, s2) = (', ...
              num2str(s_val), ', ', ...
              num2str(v_val), ', ', ...
              num2str(s1_val), ', ', ...
              num2str(s2_val), ')']);
            disp(['Unsafe states identified: ', ...
              num2str(num_unsafe_states)]);
            disp(['Initial inputs identified: ', ...
              num2str(num_feasible_inputs)]);
            
            % compute input bound ensuring oncoming vehicle 1
            % goes safely before ego vehicle
            %
            % domain of upper safe set:       s \in [-70m, -10m]
            %                                 v \in [0m/s, 12m/s]
            %                                 s1 \in [-70m, 10m]
            %
            % check if we are in the domain:
            if((s_val <= -10) && (s1_val <= 10))
                % still need to avoid collision
                x_upper_val = [s_val; v_val; s1_val];
                x_upper_idx = mdl.abs_oncoming_first.get_state_idx(x_upper_val);
                s_idx = x_upper_idx(1);
                v_idx = x_upper_idx(2);
                s1_idx = x_upper_idx(3);
                torque_upper_idx = ...
                    mdl.oncoming_first_controller(s_idx, v_idx, s1_idx);
                if(torque_upper_idx == 0)
                    continue; % this state is not in the upper safe set
                end
                torque_upper_val = ...
                    mdl.abs_oncoming_first.get_priority_input_at_idx(torque_upper_idx);
            elseif(s1_val > 10)
                % oncoming vehicle cleared intersection - we are good!
                torque_upper_val = T_max;
            else
                % we are in an unsafe configuration
                torque_upper_val = T_min;
            end
            
            % compute input bound ensuring ego vehicle
            % goes safely before oncoming vehicle 2
            %
            % domain of lower safe set:       s \in [-70m, 10m]
            %                                 v \in [0m/s, 12m/s]
            %                                 s2 \in [-70m, -10m]
            %
            % (always in the domain here)
            x_lower_val = [s_val; v_val; s2_val];
            x_lower_idx = mdl.abs_ego_first.get_state_idx(x_lower_val);
            s_idx = x_lower_idx(1);
            v_idx = x_lower_idx(2);
            s2_idx = x_lower_idx(3);
            torque_lower_idx = ...
                mdl.ego_first_controller(s_idx, v_idx, s2_idx);
            if(torque_lower_idx == 0)
                continue; % this state is not in the lower safe set
            end
            torque_lower_val = ...
                mdl.abs_ego_first.get_priority_input_at_idx(torque_lower_idx);
            
            % check if the state is safe
            if(torque_lower_val <= torque_upper_val)
                % we have feasible inputs, add to feasible input list
                for u_val = torque_lower_val:T_res:torque_upper_val
                    u_idx = mdl.abs_ego_first.get_input_idx(u_val);
                    cell_idx = num2cell(x_lower_idx);
                    feasible_inputs{cell_idx{:}}(end + 1) = u_idx;
                    num_feasible_inputs = num_feasible_inputs + 1;
                end
            else
                num_unsafe_states = num_unsafe_states + 1;
            end
            
        end
    end
end

% iterate over all states and do the following:
% 1) check if any inputs lead to unsafe states
% 2) remove inputs that lead to unsafe states from feasible input list
% (currently converges after 3 iterations)
iter = 0;
while(iter == 0 || removed_inputs)
    
    % tracking fixed-point iterations
    removed_inputs = 0;
    iter = iter + 1;
    
    for s_val = s_min:s_res:s_max
        for v_val = v_min:v_res:v_max
            for s2_val = s2_min:s2_res:s2_max

                % print current test point
                s1_val = convert_s2_to_s1(s2_val, s2_init, s12_dist_init, ...
                                            mdl.v0_min, mdl.v0_max);
                clc;
                disp(['Fixed-point algorithm iteration: ', ...
                    num2str(iter)]);
                disp(['Test point: (s, v, s1, s2) = (', ...
                    num2str(s_val), ', ', ...
                    num2str(v_val), ', ', ...
                    num2str(s1_val), ', ', ...
                    num2str(s2_val), ')']);
                disp(['Number of feasible inputs: ', ...
                    num2str(num_feasible_inputs)]);
                disp(['Unsafe inputs removed: ', ...
                    num2str(removed_inputs)]);
                disp(['Unsafe states identified: ', ...
                num2str(num_unsafe_states)]);
                
                % need to do box abstraction techniques here
                % 
                %
                % max / min initial states
                x_max_val = [s_val; v_val; s2_val];
                x_max_idx = ...
                    mdl.abs_ego_first.get_state_idx(x_max_val);
                cell_idx = num2cell(x_max_idx);
                x_min_val = [min(s_val + s_res, s_max); ...
                             min(v_val + v_res, v_max); ...
                             max(s1_val - s_res, s1_min)];
                % check if the min state is safe
                ego_occupancy = ...
                    vehicle_occupancy_scenario2(x_min_val(1));
                oncoming_occupancy = ...
                    vehicle_occupancy_scenario2(x_min_val(3));
                if(~ego_occupancy.is_ego_safe(oncoming_occupancy))
                    feasible_inputs{cell_idx{:}} = [];
                    removed_inputs = ...
                        removed_inputs + length(feasible_inputs{cell_idx{:}});
                    continue;
                end
                
                % create list of unsafe inputs
                unsafe_inputs = [];

                % check all feasible inputs
                for i = 1:length(feasible_inputs{cell_idx{:}})

                    % get input and disturbance index
                    %
                    % input value / index
                    u_idx = feasible_inputs{cell_idx{:}}(i);
                    u_val = mdl.abs_ego_first.get_priority_input_at_idx(u_idx);
                    % no disturbance
                    w_idx = zeros(1,0);
                    
                    % the abstraction self-loops we used for computing the
                    % safe sets previously are problematic here. so, here:
                    %
                    % 1) we use transition functions which still check for
                    %    unsafe states, but have no self-loops when a goal
                    %    state is reached
                    % 2) when a goal state is reached, we threshold the
                    %    state value so it stays in the domain
                    %
                    % compute max / min transitions
                    x_next_max_val = ...
                        veh_turn_expr_v3_no_loops(x_max_val, u_val, zeros(0,1));
                    x_next_min_val = ...
                        veh_turn_expr_v4_no_loops(x_min_val, u_val, zeros(0,1));
                    % if both states reached the goal set,
                    % we are safe no matter what
                    if(x_next_max_val(1) > 10 && x_next_min_val(1) > 10)
                        continue;
                    end
                    % do thresholding if only the min state reached the
                    % goal set
                    x_next_min_val(1) = min(x_next_min_val(1), 10);
                    % convert to ego first index
                    x_next_max_idx = ...
                        mdl.abs_ego_first.get_state_idx(x_next_max_val);
                    x_next_min_val_ego_first = ...
                        [x_next_min_val(1); ...
                         x_next_min_val(2); ...
                         s2_val - s_res + mdl.v_max * mdl.dt];
                    x_next_min_idx = ...
                        mdl.abs_ego_first.get_state_idx( ...
                        x_next_min_val_ego_first);
                    
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
                            for s2_next_idx = x_next_min_idx(3):x_next_max_idx(3)
                                
                                % get index for the successor cell
                                cell_next_idx = ...
                                    num2cell([s_next_idx; v_next_idx; s2_next_idx]);
                                
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

% function for converting the position of oncoming vehicle 2 
%                      to the position of oncoming vehicle 1
function s1_val = convert_s2_to_s1(s2_val, s2_init, s12_dist_init, v0_min, v0_max)
    t0 = (s2_val - s2_init) / v0_max;
    s1_val = t0 * (v0_min) + s12_dist_init + s2_init;
end
