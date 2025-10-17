function [x_safe_traj, x_unsafe_traj] = ...
                test_safe_set_point(x_test_point, abs_file, bound_file)

% load safe set for the scenario
load(abs_file);

% find the unsafe behavior
x_unsafe_idx = veh_turn_abs.get_state_idx(x_test_point);
x_unsafe_idx = x_unsafe_idx + [0; 1; 0];
% store first state
x_unsafe_traj = {};
x_unsafe_traj{1} = ...
            veh_turn_abs.get_priority_state_at_idx(x_unsafe_idx);
% roll out the dynamics
while(x_unsafe_idx ~= -1)
    % compute transition
    x_unsafe_idx = ...
        veh_turn_abs.get_transition_state(x_unsafe_idx, 1, zeros(1, 0));
    % store state
    if(x_unsafe_idx ~= -1)
        x_val = veh_turn_abs.get_priority_state_at_idx(x_unsafe_idx);
        x_unsafe_traj{end + 1} = x_val;
        % corner case: check if inter-sample behavior is unsafe
        eval(['v_bound = ', bound_file, '(x_val(1), x_val(3));']);
        if(veh_turn_abs.x_priority(2) == 1 && v_bound < x_val(2))
            break; % unsafe inter-sample behavior!
        end
        if(veh_turn_abs.x_priority(2) == 0 && v_bound > x_val(2))
            break; % unsafe inter-sample behavior!
        end
    end
end

% verify the safe behavior
x_safe_idx = veh_turn_abs.get_state_idx(x_test_point);
x_safe_plus_idx = ...
        veh_turn_abs.get_transition_state(x_safe_idx, 1, zeros(1, 0));
% store first two states
x_safe_traj = {};
x_safe_traj{1} = ...
        veh_turn_abs.get_priority_state_at_idx(x_safe_idx);
x_safe_traj{2} = ...
        veh_turn_abs.get_priority_state_at_idx(x_safe_plus_idx);
% roll out the dynamics
while(sum(x_safe_plus_idx == x_safe_idx) ~= 3)
    % compute transition
    x_safe_idx = x_safe_plus_idx;
    x_safe_plus_idx = ...
        veh_turn_abs.get_transition_state(x_safe_idx, 1, zeros(1, 0));
    % store state
    x_safe_traj{end + 1} = ...
        veh_turn_abs.get_priority_state_at_idx(x_safe_idx);
end

end
