% load safe set for ego first scenario
abs_file = 'veh_turn_abs_oncoming_first_11_2.mat';
bound_file = 'outer_approx_boundary_scenario2';

% test point 1 (left area)
x_test_point1 = [-36.4; 10.2; -66];
% test point 2 (middle area)
x_test_point2 = [-26.8; 8.4; -22.8];
% test point 3 (right area)
x_test_point3 = [-15.6; 6; -1.2];

% make a distinction between the safe and unsafe behavior
[x_safe_traj1, x_unsafe_traj1] = ...
    test_safe_set_point(x_test_point1, abs_file, bound_file);
[x_safe_traj2, x_unsafe_traj2] = ...
    test_safe_set_point(x_test_point2, abs_file, bound_file);
[x_safe_traj3, x_unsafe_traj3] = ...
    test_safe_set_point(x_test_point3, abs_file, bound_file);
