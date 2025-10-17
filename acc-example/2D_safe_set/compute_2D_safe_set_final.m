% create vehicle following dynamics

% h = headway
% v = ego velocity
% v_L = lead velocity
% a = ego acceleration

n_x = 2; % states are [h; v]
n_u = 1; % input is a
n_w = 1; % disturbance is v_L

% vehicle following dynamics, exact discretization with dt = 0.1s
% x(t + 0.1) = x(1) + w(1)*0.1 - x(2)*0.1 - 0.5*u(1)*(0.1)^2
% x(t + 0.1) = x(2) + u(1)*0.1

veh_follow_dyn = monotone_dyn(n_x, n_u, n_w, '[x(1) + w(1)*0.1 - x(2)*0.1 - 0.5*u(1)*(0.1)^2; x(2) + u(1)*0.1]');

% create abstraction
veh_follow_abs = monotone_abstraction(veh_follow_dyn);

% motion constraints
h_min = 5; % m
h_max = 80; % m
v_min = 0; % m/s
v_max = 20; % m/s
a_min = -3; % m/s^2
a_max = 3; % m/s^2

% state constraints
x_range = [h_min h_max;
           v_min v_max];
veh_follow_abs.set_x_range(x_range);

% input constraints
u_range = [a_min a_max];
veh_follow_abs.set_u_range(u_range);

% disturbance constraints
w_range = [v_min v_max];
veh_follow_abs.set_w_range(w_range);

% abstraction resolution
h_res = 0.2;
v_res = 0.1;
a_res = 0.1;

% state resolution
x_res = [h_res;
         v_res];
veh_follow_abs.set_x_res(x_res);

% input resolution
u_res = a_res;
veh_follow_abs.set_u_res(u_res);

% disturbance resolution
w_res = v_res;
veh_follow_abs.set_w_res(w_res);

% state priorities
x_priority = [0;  % smaller headway more unsafe
              1]; % larger velocity more unsafe
veh_follow_abs.set_x_priority(x_priority);

% input priorities
u_priority = 1; % larger acceleration more unsafe
veh_follow_abs.set_u_priority(u_priority);

% disturbance priorities
w_priority = 0; % smaller velocity more unsafe
veh_follow_abs.set_w_priority(w_priority);

%%% INVARIANT SET ALGORITHM (WORKAROUND VERSION) %%%
% We'll implement our own version that works around the library bugs

fprintf('Starting invariant set computation...\n');
tic; % Start timing
veh_follow_abs.initialize_safe_set();
iter = 1;
while(true)
    unsafe_basis_indices = {};
    iter = iter + 1;
    disp(iter);
    for i = 1:veh_follow_abs.basis_size
        % test points on the boundary of the safe set
        x_idx = veh_follow_abs.get_basis_idx(i);
        x_val = veh_follow_abs.get_priority_state_at_idx(x_idx);
        u_val = max(a_min, -10*x_val(2)); % use max possible braking force
        u_idx = veh_follow_abs.get_input_idx(u_val);
        w_val = v_min; % assume front vehicle not moving
        w_idx = veh_follow_abs.get_disturbance_idx(w_val);
        
        % COMPUTE TRANSITION DIRECTLY (bypassing all problematic methods)
        % Implement the dynamics directly to avoid library bugs
        dt = 0.1;
        x_plus_val = [x_val(1) + w_val*dt - x_val(2)*dt - 0.5*u_val*dt^2;
                      x_val(2) + u_val*dt];
        
        % convert to next state index
        x_plus_idx = veh_follow_abs.get_state_idx(x_plus_val);
        
        % check if the point should be safe (analytically)
        headway = x_val(1);
        velocity = x_val(2);
        stop_dist = headway - velocity^2/(2*abs(u_val));
        % check if the point is safe
        if(~veh_follow_abs.x_in_safe_set(x_plus_idx))
            unsafe_basis_indices{end + 1} = x_idx;
            if(stop_dist >= 5)
                1; % should be safe!
            end
        end
    end
    if(isempty(unsafe_basis_indices))
        % all basis elements are safe and we are done!
        break;
    else
        for i = 1:length(unsafe_basis_indices)
            veh_follow_abs.remove_basis_idx(unsafe_basis_indices{i});
        end
    end
end

% Stop timing and display results
computation_time = toc;
fprintf('Converged! All basis elements are safe.\n');
fprintf('Safe set computation completed in %d iterations.\n', iter);
fprintf('Final safe set has %d basis elements.\n', veh_follow_abs.basis_size);
fprintf('Computation completed in %.2f seconds (%.0f ms)\n', computation_time, computation_time*1000);

% Print some basis elements
fprintf('\nFirst 10 Safe Set Basis Elements:\n');
for i = 1:min(10000, veh_follow_abs.basis_size)
    x_idx = veh_follow_abs.get_basis_idx(i);
    x_val = veh_follow_abs.get_priority_state_at_idx(x_idx);
    fprintf('Basis %d: [%.1f, %.1f]\n', i, x_val(1), x_val(2));
end
