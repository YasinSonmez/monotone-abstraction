% test vehicle parameters (from NSF-TTP paper)
M = 2044;        % vehicle weight (kg)
R_w = 0.3074;    % wheel radius (m)
Tb = 500;        % nominal braking force (Nm)

% frictional force modelling paramters:
% Ff = alpha + beta*v + gamma*v^2
% friction force modelling parameter (constant term)
alpha = 339.1329;
% friction force modelling parameter (linear term, larger for test)
beta = 5;
% friction force modelling parameter (squared term, larger for test)
gamma = 15;

% initial conditions
p0 = 0;  % (m)
v0 = 12; % (m / s)

% time step
dt = 0.2;

% solve w ode45
[t_step, x_step] = ode45(@veh_dyn, [0 dt], [p0; v0]);
[t_full, x_full] = ode45(@veh_dyn, [0 20], [p0; v0]);

% plot results
plot(t_step, x_step);
title('x step');
figure;
plot(t_full, x_full);
title('x full');

function dxdt = veh_dyn(t, x)

    % test vehicle parameters (from NSF-TTP paper)
    M = 2044;        % vehicle weight (kg)
    R_w = 0.3074;    % wheel radius (m)
    Tb = 500;        % nominal braking force (Nm)

    % frictional force modelling paramters:
    % Ff = alpha + beta*v + gamma*v^2
    % friction force modelling parameter (constant term)
    alpha = 339.1329;
    % friction force modelling parameter (linear term, larger for test)
    beta = 5;
    % friction force modelling parameter (squared term, larger for test)
    gamma = 15;

    % put dynamics into the following form:
    % vdot = - a - b * v - c * v^2
    a = (1 / M) * (alpha + Tb / R_w);
    b = (1 / M) * beta;
    c = (1 / M) * gamma;

    % equations of motion:
    % vdot = 1 / M * (Tb / R_w - beta - gamma * v^2)
    dxdt = [x(2); - a - b * x(2) - c * x(2)^2];
    
    % ensure the vehicle does not reverse
    if(x(2) <= 0 && dxdt(2) < 0)
        dxdt(2) = 0;
    end
    
end
