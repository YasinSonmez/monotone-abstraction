function con = veh_follow_const()
% provides all modelling parameters for the vehicle following example

% motion constraints
con.h_min = 0; % m
con.h_max = 80; % m
con.v_min = 0; % m/s
con.v_max = 20; % m/s

% wheel torque constraints
con.T_max = 1200; % Nm
con.T_brake_min = -1800; % Nm
con.T_brake_max = -2400; % Nm

% modelling parameters
con.M_min = 2000; % kg
con.M_max = 2500; % kg
con.R_w_min = 0.30; % m
con.R_w_max = 0.35; % m

% frictional force modelling parameters
con.alpha_min = 300;
con.alpha_max = 350;
con.beta_min = 0.10;
con.beta_max = 0.25;
con.gamma_min = 0.30;
con.gamma_max = 0.65;

% using model with air drag, compute upper / lower bounds on the vehicle
% deceleration while braking (to be used in precomputation step)
%
% to overapproximate the safe set, be:
% - pessimistic about how fast the lead vehicle brakes
% - optimistic about how fast the ego vehicle brakes
% (still with the worst-case parameter values in each case, though)
%
% lead vehicle:
Ff_lead = con.alpha_max;
con.a_lead_brake = (1 / con.M_min) * ...
                    (abs(con.T_brake_max) / con.R_w_min + Ff_lead);
% ego vehicle:
Ff_ego = con.alpha_min + con.beta_min * con.v_max + ...
            con.gamma_min * con.v_max^2;
con.a_ego_brake = (1 / con.M_max) * ...
                    (abs(con.T_brake_min) / con.R_w_max + Ff_ego);

% abstraction resolution
con.h_res = 0.8; % 100 pts
con.v_res = 0.4; % 50 pts
con.T_res = 150; % 20 pts for u, 24 pts for w

% vehicle dynamics sampling time
con.dt = 0.4;

end
