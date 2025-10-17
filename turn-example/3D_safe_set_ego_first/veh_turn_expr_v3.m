function x_plus = veh_turn_expr_v3(x_val, u_val, w_val)

% s = ego vehicle position along its path
% v = ego velocity along its path
% s_0 = oncoming vehicle position along its path
% a = ego vehicle acceleration along its path

n_x = 3; % states are [h; v; v_l]
n_u = 1; % input is a
n_w = 0; % disturbance is a_L

% ensure state dimensions match
if(size(x_val,1) ~= n_x || size(x_val,2) ~= 1)
    error(['Invalid state dimension, expected ', num2str(dyn.n_x), ' by 1 column vector.']);
end
if(size(u_val,1) ~= n_u || size(u_val,2) ~= 1)
    error(['Invalid input dimension, expected ', num2str(dyn.n_u), ' by 1 column vector.']);
end
if(size(w_val,1) ~= n_w || size(w_val,2) ~= 1)
    error(['Invalid disturbance dimension, expected ', num2str(dyn.n_w), ' by 1 column vector.']);
end

% use modelling parameters defined in constants file
con = veh_turn_const_v3();

% from the modelling parameters, compute a, b, and c for each vehicle,
% where we represent the dynamics for each vehicle as follows:
%
% \dot{p} = v,
% \dot{v} = a + b * v + c * v^2,
%
% using our model, we have:
%
% a = (1 / con.M) * (T / con.R_w - con.alpha)
%
% b = - (1 / con.M) * con.beta
%
% c = - (1 / con.M) * con.gamma

% choose worst-case parameter values for ego vehicle
% here, we are trying to MINIMIZE a, b and c
R_w = ternary(u_val > 0, con.R_w_max, con.R_w_min);
M = ternary((u_val / R_w - con.alpha_max) > 0, con.M_max, con.M_min);
a = (1 / M) * (u_val / R_w - con.alpha_max);
b = - (1 / con.M_min) * con.beta_max;
c = - (1 / con.M_min) * con.gamma_max;

% using numerical solution to vehicle dynamics eqn, compute next state
x_plus = zeros(n_x, 1);
[x_plus(1), x_plus(2)] = veh_dyn(x_val(1), x_val(2), a, b, c, con.dt);
x_plus(3) = x_val(3) + con.v0_max * con.dt;

% do thresholding on the ego vehicle velocity
% (this is necessary because sometimes the ode45 solver lets
% the ego vehicle go slightly outside of its velocity bounds)
x_plus(2) = min(max(x_plus(2), con.v_min), con.v_max);

% determine which zone the two vehicles are in
ego_occupancy = vehicle_occupancy_scenario1(x_plus(1));
oncoming_occupancy = vehicle_occupancy_scenario1(x_plus(3));

% check which situation we are in
if(~ego_occupancy.is_ego_safe(oncoming_occupancy))
    % this state is unsafe => map to an unsafe state
    x_plus = con.unsafe_state;
elseif(ego_occupancy.has_ego_reached_goal(oncoming_occupancy))
    % the ego vehicle reached the goal set => use a self-loop
    x_plus = x_val;
end

% compute next state
function [pos_plus, vel_plus] = veh_dyn(p0, v0, a, b, c, dt)
    
    % use ode45 solver to get state at next time step
    [t, x] = ode45(@f, [0 dt], [p0; v0; a; b; c]);
    
    % make sure we computed the state up until the sample time
    threshold = 1e-5;
    assert(abs(dt - t(end)) < threshold);
    
    % extract the position and velocity at time t = dt
    pos_plus = x(end, 1);
    vel_plus = x(end, 2);
    
end

function dxdt = f(t, x)

    % state convention is as follows:
    % x(1) = position,
    % x(2) = velocity,
    % x(3) = a,
    % x(4) = b,
    % x(5) = c
    
    % equations of motion, with state defined as above:
    dxdt = [x(2);
            x(3) + x(4) * x(2) + x(5) * x(2)^2;
            0;
            0;
            0];
    
    % ensure the vehicle does not reverse
    if(x(2) <= 0 && dxdt(2) < 0)
        dxdt(2) = 0;
    end
    
    % ensure the vehicle does not exceed its maximum velocity
    if(x(2) >= 12 && dxdt(2) > 0)
        dxdt(2) = 0;
    end
    
end

% helper function implementing ternary operator for convenience
function output = ternary(condition, input1, input2)
    if(condition)
        output = input1;
    else
        output = input2;
    end
end

end
