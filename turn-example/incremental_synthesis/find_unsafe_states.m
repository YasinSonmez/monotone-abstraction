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
s1_max = mdl.abs_oncoming_first.x_range(3,2);

% for checking if we are out of range
s_max_oncoming_first = mdl.abs_oncoming_first.x_range(1,2);

% list of unsafe states
unsafe_states = {};

% iterate over all states and check for a feasible input
for s_val = s_min:s_res:s_max
    for v_val = v_min:v_res:v_max
        for s2_val = s2_min:s2_res:s2_max
            
            % check if we are out of range
            if(s_val > s_max_oncoming_first)
                continue;
            end
            
            % print current test point
            t0 = (s2_val - s2_init) / mdl.v0_max;
            s1_val = t0 * (mdl.v0_min) + s12_dist_init + s2_init;
            s1_val = min(s1_val, s1_max);
            % s1_val = min(s2_val + s12_dist, s1_max);
            disp(['Test point: (s, v, s1, s2) = (', ...
              num2str(s_val), ', ', ...
              num2str(v_val), ', ', ...
              num2str(s1_val), ', ', ...
              num2str(s2_val), ')']);
            
            % compute input bound ensuring oncoming vehicle 1
            % goes safely before ego vehicle
            x_upper_val = [s_val; v_val; s1_val];
            x_upper_idx = mdl.abs_oncoming_first.get_state_idx(x_upper_val);
            s_idx = x_upper_idx(1);
            v_idx = x_upper_idx(2);
            s1_idx = x_upper_idx(3);
            accel_upper_idx = ...
                mdl.oncoming_first_controller(s_idx, v_idx, s1_idx);
            if(accel_upper_idx == 0)
                continue; % this state is not in the upper safe set
            end
            accel_upper_val = ...
                mdl.abs_oncoming_first.get_priority_input_at_idx(accel_upper_idx);
            
            % compute input bound ensuring ego vehicle 
            % goes safely before oncoming vehicle 2
            x_lower_val = [s_val; v_val; s2_val];
            x_lower_idx = mdl.abs_ego_first.get_state_idx(x_lower_val);
            s_idx = x_lower_idx(1);
            v_idx = x_lower_idx(2);
            s2_idx = x_lower_idx(3);
            accel_lower_idx = ...
                mdl.ego_first_controller(s_idx, v_idx, s2_idx);
            if(accel_lower_idx == 0)
                continue; % this state is not in the lower safe set
            end
            accel_lower_val = ...
                mdl.abs_ego_first.get_priority_input_at_idx(accel_lower_idx);
            
            % check if the state is safe
            if(accel_lower_val > accel_upper_val)
                % This state is in the intersection of the upper and
                % lower safe set, but there is no feasible input.
                unsafe_states{end + 1} = [s_val; v_val; s1_val; s2_val];
            end
            
        end
    end
end
