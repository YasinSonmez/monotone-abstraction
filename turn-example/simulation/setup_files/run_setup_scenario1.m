function mdl = run_setup_scenario1
%%% Description %%%
% Set up model struct for unprotected left turn simulation

% load the precomputed abstractions
load('veh_turn_abs_ego_first_5_7_2021.mat', 'veh_turn_abs');
mdl.abs_ego_first = veh_turn_abs;
load('veh_turn_abs_oncoming_first_5_7_2021.mat', 'veh_turn_abs');
mdl.abs_oncoming_first = veh_turn_abs;

% load the precomputed controllers
load('ego_first_controller_5_7_2021.mat', 'ego_first_controller');
mdl.ego_first_controller = ego_first_controller;
load('oncoming_first_controller_5_7_2021.mat', 'oncoming_first_controller');
mdl.oncoming_first_controller = oncoming_first_controller;

% load the controller
load('feasible_inputs_5_7_2021.mat', 'feasible_inputs');
mdl.feasible_inputs = feasible_inputs;

% set simulation parameters
mdl.x0 = [-38; 9.5; 0; -70];
mdl.v0 = 12;
mdl.v_min = 0;
mdl.v_max = 12;
mdl.v0_min = 8;
mdl.v0_max = 12;
mdl.T_min = -1800;
mdl.T_max = 1200;

% choose the parameters for our vehicle dynamics model
mdl.M = 2250;
mdl.R_w = 0.325;
mdl.alpha = 325;
mdl.beta = 0.175;
mdl.gamma = 0.475;

% sampling time
mdl.dt = 0.5;     % (for abstraction computations)
mdl.dt_sim = 0.1; % (for vehicle dynamics block in simulation)

end
