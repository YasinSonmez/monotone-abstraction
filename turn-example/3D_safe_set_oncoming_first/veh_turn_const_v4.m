function con = veh_turn_const_v4()
% provides all modelling parameters for the vehicle following example

% motion constraints
con.s_min = -70; % m
con.s_max = -10; % m
con.s0_min = -70; % m
con.s0_max = 10; % m
con.v_min = 0; % m/s
con.v_max = 12; % m/s (for the ego vehicle, this is ~25 mph which seems
                %      realistic for executing a turn)

% wheel torque constraints
con.T_max = 1200; % Nm
con.T_min = -1800; % Nm

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

% using model with air drag, compute upper / lower bounds on the vehicle
% acceleration (to be used in precomputation step)
%
% for this scenario, compute the min acceleration
% (still with the worst-case parameter values, though)
Ff_max = con.alpha_min + con.beta_min * con.v_max + ...
            con.gamma_min * con.v_max^2;
con.a_min_ego = (1 / con.M_max) * (con.T_min / con.R_w_max + Ff_max);

% abstraction resolution
con.s_res = 2; % 30 cells
con.s0_res = 2; % 40 cells
con.v_res = 0.5; % 24 cells
con.T_res = 100; % 30 cells

% sampling time for vehicle dynamics discretization
con.dt = 0.5;

% constants for unprotected left turn scenario
con.v0_min = 8; % m/s
con.v0_max = 12; % m/s
con.collision_zone_width = 10; % m (collision zone width, 
                               %    length of Hyundai Ioniq = ~4.5m)
                               
% unsafe used for re-mapping the transition function when needed
% (any unsafe state will do)
con.unsafe_state = [con.s_min - 15; con.v_max + 5; con.s0_min - 15];

end
