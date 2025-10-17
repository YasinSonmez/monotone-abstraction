%%% Get motion constraint constants and 3D monotone abstraction %%%
%%% for vehicle following scenario %%%
[con, veh_follow_abs] = get_3D_monotone_abs_v2();
veh_follow_abs.initialize_safe_set();

%%% PRECOMPUTATION STEP %%%
% iterate over discretized headways
% for h_point = con.h_min:con.h_res:con.h_max
%     % iterate over discretized lead vehicle velocities
%     for v_lead_point = con.v_min:con.v_res:con.v_max
%         disp(['Test point: (h, v_L) = (', ...
%               num2str(h_point), ', ', ...
%               num2str(v_lead_point), ')']);
%         v_follow_point = outer_approx_boundary(h_point, ...
%                                                v_lead_point, ...
%                                                abs(con.a_min_brake), ...
%                                                abs(con.a_max_brake), ...
%                                                con.h_min);
%         % threshold at maximum velocity
%         v_follow_point = min(v_follow_point, con.v_max);
%         % quantize v_follow to the priority state
%         threshold = 1e-5;
%         v_follow_point = con.v_min + ceil((v_follow_point - con.v_min) ...
%                                    / con.v_res - threshold) * con.v_res;
%         % expand the v_follow point upwards until an unsafe collision
%         % is possible
%         while(v_follow_point <= con.v_max + threshold && ...
%               ~unsafe_collision_check(h_point, v_follow_point, ...
%                                       v_lead_point, con.a_min_brake, ...
%                                       con.a_max_brake, con.v_allow, 1))
%             v_follow_point = v_follow_point + con.v_res;
%         end
%         v_follow_point = v_follow_point - con.v_res;
%         % quantize to the last safe state we saw (within range)
%         x_bound = [h_point; v_follow_point; v_lead_point];
%         x_bound_idx = veh_follow_abs.get_state_idx(x_bound);
%         veh_follow_abs.add_basis_idx(x_bound_idx);
%     end
% end

%%% INVARIANT SET ALGORITHM %%%
% currently converges after 56 iterations %
iter = 0;
while(true)
    unsafe_basis_indices = {};
    iter = iter + 1;
    for i = 1:veh_follow_abs.basis_size
        % print status message
        clc;
        disp(['Iteration: ', num2str(iter)]);
        disp(['Checking basis index: ', num2str(i)]);
        disp(['Unsafe basis elements found: ', num2str(length(unsafe_basis_indices))]);
        disp(['Number of basis elements: ', num2str(veh_follow_abs.basis_size)]);
        % test points on the boundary of the safe set
        x_idx = veh_follow_abs.get_basis_idx(i);
        % use max possible braking force for both vehicles
        % (veh_follow_expr.m ensures vehicles do not reverse)
        u_idx = 1;
        w_idx = veh_follow_abs.w_num_cells;
        x_plus_idx = veh_follow_abs.get_transition_state(x_idx, u_idx, w_idx);
        % check if the point is safe
        if(~veh_follow_abs.x_in_safe_set(x_plus_idx))
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
            veh_follow_abs.remove_basis_idx(unsafe_basis_indices{i});
        end
    end
end
