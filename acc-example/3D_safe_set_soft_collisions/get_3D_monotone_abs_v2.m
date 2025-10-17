function [con, veh_follow_abs] = get_3D_monotone_abs_v2()
% create 3D monotone abstraction for vehicle following scenario:

% h = headway
% v = ego velocity
% v_L = lead velocity
% a = ego acceleration
% a_L = lead acceleration

n_x = 3; % states are [h; v; v_l]
n_u = 1; % input is a
n_w = 1; % disturbance is a_L

% vehicle following dynamics, exact discretization with dt = 0.1s
% [x(1)]            [ x(1) + x(3)*0.1 + 0.5*w(1)*(0.1)^2 - x(2)*0.1 - 0.5*u(1)*(0.1)^2 ]
% [x(2)](t + 0.1) = [ x(2) + u(1)*0.1                                                  ]
% [x(3)]            [ x(3) + w(1)*0.1                                                  ]

% (see file veh_follow_expr.m)

veh_follow_dyn = ...
    monotone_dyn(n_x, n_u, n_w, 'veh_follow_expr_v2');

% create abstraction
veh_follow_abs = monotone_abstraction(veh_follow_dyn);

% get modelling parameters
con = veh_follow_const_v2();

% state constraints
x_range = [con.h_min con.h_max;
           con.v_min con.v_max;
           con.v_min con.v_max];
veh_follow_abs.set_x_range(x_range);

% input constraints
u_range = [con.T_brake_max con.T_max];
veh_follow_abs.set_u_range(u_range);

% disturbance constraints
w_range = [con.T_brake_min con.T_max];
veh_follow_abs.set_w_range(w_range);

% state resolution
x_res = [con.h_res;
         con.v_res;
         con.v_res];
veh_follow_abs.set_x_res(x_res);

% input resolution
u_res = con.T_res;
veh_follow_abs.set_u_res(u_res);

% disturbance resolution
w_res = con.T_res;
veh_follow_abs.set_w_res(w_res);

% state priorities
x_priority = [0;  % smaller headway more unsafe
              1;  % larger ego velocity more unsafe
              0]; % smaller lead velocity more unsafe
veh_follow_abs.set_x_priority(x_priority);

% input priorities
u_priority = 1; % larger ego vehicle wheel torque more unsafe
veh_follow_abs.set_u_priority(u_priority);

% disturbance priorities
w_priority = 0; % smaller lead vehicle wheel torque more unsafe
veh_follow_abs.set_w_priority(w_priority);

end
