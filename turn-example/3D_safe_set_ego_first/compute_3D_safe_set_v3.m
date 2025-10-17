%%% Get motion constraint constants and 3D monotone abstraction %%%
%%% for vehicle following scenario %%%
[con, veh_turn_abs] = get_3D_monotone_abs_v3();
veh_turn_abs.initialize_safe_set();

%%% PRECOMPUTATION STEP %%%
% iterate over ego vehicle positions (along its pathway)
% for s_point = con.s_min:con.s_res:con.s_max
%     % iterate over oncoming vehicle positions (along its pathway)
%     for s0_point = con.s0_min:con.s0_res:con.s0_max
%         disp(['Test point: (s, s0) = (', ...
%               num2str(s_point), ', ', num2str(s0_point), ')']);
%         v_ego_point = outer_approx_boundary_scenario1(s_point, s0_point);
%         % quantize to discrete state
%         x_bound = [s_point; v_ego_point; s0_point];
%         x_bound_idx = veh_turn_abs.get_state_idx(x_bound);
%         % ensure we don't quantize to an unsafe state
%         if(mod(v_ego_point, con.v_res) > 1e-5 && x_bound_idx(2) > 1)
%             x_bound_idx = x_bound_idx - [0; 1; 0];
%         end
%         veh_turn_abs.add_basis_idx(x_bound_idx);
%     end
% end

%%% INVARIANT SET ALGORITHM %%%
% currently converges after 75 iterations %
tic;
iter = 0;
while(true)
    unsafe_basis_indices = {};
    iter = iter + 1;
    for i = 1:veh_turn_abs.basis_size
        % print status message
        clc;
        disp(['Iteration: ', num2str(iter)]);
        disp(['Checking basis index: ', num2str(i)]);
        disp(['Unsafe basis elements found: ', num2str(length(unsafe_basis_indices))]);
        disp(['Number of basis elements: ', num2str(veh_turn_abs.basis_size)]);
        % test points on the boundary of the safe set
        x_idx = veh_turn_abs.get_basis_idx(i);
        % use max acceleration force for the ego vehicle
        % (veh_follow_expr_3.m ensures max velocity is not exceeded)
        u_idx = 1;
        % no disturbance in this scenario
        w_idx = zeros(1,0);
        x_plus_idx = veh_turn_abs.get_transition_state(x_idx, u_idx, w_idx);
        % check if the point is safe
        if(~veh_turn_abs.x_in_safe_set(x_plus_idx))
            unsafe_basis_indices{end + 1} = x_idx;
        end
    end
    if(isempty(unsafe_basis_indices))
        % all basis elements are safe and we are done!
        break;
    else
        for i = 1:length(unsafe_basis_indices)
            % print status message
            clc;
            disp(['Iteration: ', num2str(iter)]);
            disp(['Unsafe basis elements: ', num2str(length(unsafe_basis_indices))]);
            disp(['Removing unsafe basis index: ', num2str(i)]);
            veh_turn_abs.remove_basis_idx(unsafe_basis_indices{i});
        end
    end
end
toc;
