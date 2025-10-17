% test vehicle parameters (from NSF-TTP paper)
M = 2044;        % vehicle weight (kg)
R_w = 0.3074;    % wheel radius (m)
beta = 339.1329; % friction force modelling parameter
gamma = 15;      % friction force modelling parameter (larger for testing)
Tb = 500;        % nominal braking force (Nm)

% put dynamics into the following form:
% vdot = - b * v^2 - a
a = (1 / M) * (beta + Tb / R_w);
b = (1 / M) * gamma;

% initial conditions
p0 = 0;  % (m)
v0 = 12; % (m / s)

% solve w ode45
[t, x] = ode45(@veh_dyn, [0 20], [p0; v0]);

% estimate stopping distance
stop_dist_est = estimate_stop_pos(t, x);

% stopping distance (exact solution)
stop_dist_exact = 1 / b * log(1 + b / a * v0^2);

% get x_exact from the exact solution
x_exact = [1 / b * log(cos(atan(sqrt(b / a) * v0) - sqrt(a * b) * t) ...
                 / cos(atan(sqrt(b / a) * v0))), ...
           sqrt(a / b)* tan(atan(sqrt(b / a) * v0) - sqrt(a * b) * t)];

% plot results
plot(t, x);
title('x simulated');
figure;
plot(t, x_exact);
title('x exact');

function dxdt = veh_dyn(t, x)

    % test vehicle parameters (from NSF-TTP paper)
    M = 2044;        % vehicle weight (kg)
    R_w = 0.3074;    % wheel radius (m)
    beta = 339.1329; % friction force modelling parameter
    gamma = 15;      % friction force modelling parameter (larger for testing)
    Tb = 500;        % nominal braking force (Nm)

    % put dynamics into the following form:
    % vdot = - b * v^2 - a
    a = (1 / M) * (beta + Tb / R_w);
    b = (1 / M) * gamma;
    
    % equations of motion:
    % vdot = 1 / M * (Tb / R_w - beta - gamma * v^2)
    dxdt = [x(2); - b * x(2)^2 - a];
    
    % ensure the vehicle does not reverse
    if(x(2) <= 0 && dxdt(2) < 0)
        dxdt(2) = 0;
    end
    
end

function x_stop = estimate_stop_pos(t, x)
    
    % find first negative velocity
    idx = find(x(:,2) < 0, 1);
    t0 = t(idx - 1);
    t1 = t(idx);
    x0 = x(idx - 1, 1);
    x1 = x(idx, 1);
    v0 = x(idx - 1, 2);
    v1 = x(idx, 2);
    
    % estimate first time instant when v = 0
    t_stop = t0 + v0 / (v0 - v1) * (t1 - t0);
    % use t_stop to estimate x_stop (assuming constant deceleration)
    x_stop = x0 + (t_stop - t0) / (t1 - t0) * x1;
    
end
