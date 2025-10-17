function con = veh_follow_const_v2()
% provides all modelling parameters for the vehicle following example

% motion constraints
con.h_min = 0; % m
con.h_max = 80; % m
con.v_min = 0; % m/s
con.v_max = 20; % m/s

% wheel torque constraints
con.T_max = 1200; % Nm
con.T_brake_min = -1800; % Nm
con.T_brake_max = -2500; % Nm

% modelling parameters
con.M_min = 2000; % kg
con.M_max = 2250; % kg
con.R_w_min = 0.30; % m
con.R_w_max = 0.35; % m

% frictional force modelling parameters
con.alpha_min = 300;
con.alpha_max = 350;
con.beta_min = 0.10;
con.beta_max = 0.25;
con.gamma_min = 0.30;
con.gamma_max = 0.65;

% using model with air drag, compute upper and lower bounds for
% deceleration rate while braking (to be used in precomputation step)
% Ff_max = con.alpha_max + con.beta_max * con.v_max + ...
%             con.gamma_max * con.v_max^2;
% Ff_min = con.alpha_min;
% con.a_max_brake = (1 / con.M_min) * ...
%                     (abs(con.T_brake_max) / con.R_w_min + Ff_max);
% con.a_min_brake = (1 / con.M_max) * ...
%                     (abs(con.T_brake_min) / con.R_w_max + Ff_min);

% compute lower / upper bounds on deceleration rate while braking
% want to ensure a(t) - a_L(t) <= 0 while the vehicles are braking
% in this case, v(t) - v_L(t) will DECREASE over time
Ff_min = con.alpha_min;
Ff_max = con.alpha_max + con.beta_max * con.v_max + ...
            con.gamma_max * con.v_max^2;
con.a_min_brake_ego = (1 / con.M_max) * ...
                        (abs(con.T_brake_max) / con.R_w_max + Ff_min);
con.a_max_brake_lead = (1 / con.M_min) * ...
                        (abs(con.T_brake_min) / con.R_w_min + Ff_max);
assert(con.a_min_brake_ego >= con.a_max_brake_lead);

% abstraction resolution
con.h_res = 2; % 40 cells
con.v_res = 1; % 20 cells
con.T_res = 100; % 37 cells for u, 30 cells for w

% vehicle dynamics sampling time
con.dt = 0.5;

% allowable relative velocity at time of impact
% (for computing expanded safe set)
con.v_allow = 0; % m/s, from "AHS Safe Control Laws for Platoon Leaders"

end
