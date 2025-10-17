function x_plus = veh_follow_expr(x_val, u_val, w_val)

% h = headway
% v = ego velocity
% v_L = lead velocity
% a = ego acceleration
% a_L = lead acceleration

n_x = 3; % states are [h; v; v_l]
n_u = 1; % input is a
n_w = 1; % disturbance is a_L

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
con = veh_follow_const();

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
% here, we are trying to MAXIMIZE a, b and c
R_w_ego = ternary(u_val > 0, con.R_w_min, con.R_w_max);
M_ego = ternary((u_val / R_w_ego - con.alpha_min) > 0, con.M_min, con.M_max);
a = (1 / M_ego) * (u_val / R_w_ego - con.alpha_min);
b = - (1 / con.M_max) * con.beta_min;
c = - (1 / con.M_max) * con.gamma_min;

% choose worst-case parameter values for lead vehicle
% here, we are trying to MINIMIZE a, b and c
R_w_L = ternary(w_val > 0, con.R_w_max, con.R_w_min);
M_L = ternary((w_val / R_w_L - con.alpha_max) > 0, con.M_max, con.M_min);
a_L = (1 / M_L) * (w_val / R_w_L - con.alpha_max);
b_L = - (1 / con.M_min) * con.beta_max;
c_L = - (1 / con.M_min) * con.gamma_max;

% using numerical solution to vehicle dynamics eqn, compute next state
[pos_plus, vel_plus] = veh_dyn(0, x_val(2), a, b, c, con.dt, 0);
[posL_plus, velL_plus] = veh_dyn(0, x_val(3), a_L, b_L, c_L, con.dt, 1);
x_plus = [x_val(1) + posL_plus - pos_plus; vel_plus; velL_plus];

% do thresholding on the ego / lead vehicle velocities
% (this is necessary because sometimes the ode45 solver lets
% the vehicles go slightly outside of their velocity bounds)
x_plus(2) = min(max(x_plus(2), con.v_min), con.v_max);
x_plus(3) = min(max(x_plus(3), con.v_min), con.v_max);

% compute next state
function [pos_plus, vel_plus] = veh_dyn(p0, v0, a, b, c, dt, lead)
    
    % use ode45 solver to get state at next time step
    [t, x] = ode45(@f, [0 dt], [p0; v0; a; b; c; lead]);
    
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
    % x(5) = c,
    % x(6) = 1 if this is the lead vehicle, 0 otherwise
    
    % equations of motion, with state defined as above:
    dxdt = [x(2);
            x(3) + x(4) * x(2) + x(5) * x(2)^2;
            0;
            0;
            0;
            0];
    
    % ensure the vehicle does not reverse
    if(x(2) <= 0 && dxdt(2) < 0)
        dxdt(2) = 0;
    end
    
    % ensure the LEAD vehicle does not exceed v_max
    if(x(6) == 1 && x(2) >= 20 && dxdt(2) > 0)
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
