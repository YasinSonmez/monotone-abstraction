function [con, veh_follow_abs] = get_3D_monotone_abs_v3()
% create 3D monotone abstraction for vehicle following scenario:

% s = ego vehicle position along path
% v = ego velocity
% s_0 = oncoming vehicle position along path
% a = ego vehicle acceleration

n_x = 3; % states are [s; v; s_0]
n_u = 1; % input is a
n_w = 0; % no disturbance

% vehicle following dynamics, exact discretization with dt = 0.1s
% [x(1)]            [ x(1) + x(2)*0.1 + 0.5*u(1)*(0.1)^2 ]
% [x(2)](t + 0.1) = [ x(2) + u(1)*0.1                    ]
% [x(3)]            [ x(3) + v_0*0.1                     ]

% (see file veh_follow_expr_v3.m)

veh_follow_dyn = ...
    monotone_dyn(n_x, n_u, n_w, 'veh_turn_expr_v3');

% create abstraction
veh_follow_abs = monotone_abstraction(veh_follow_dyn);

% get modelling parameters
con = veh_turn_const_v3();

% state constraints
x_range = [con.s_min con.s_max;
           con.v_min con.v_max;
           con.s0_min con.s0_max];
veh_follow_abs.set_x_range(x_range);

% input constraints
u_range = [con.T_min con.T_max];
veh_follow_abs.set_u_range(u_range);

% disturbance constraints
% w_range = [con.a_max_brake con.a_max];
% veh_follow_abs.set_w_range(w_range);

% state resolution
x_res = [con.s_res;
         con.v_res;
         con.s0_res];
veh_follow_abs.set_x_res(x_res);

% input resolution
u_res = con.T_res;
veh_follow_abs.set_u_res(u_res);

% disturbance resolution
% w_res = con.a_res;
% veh_follow_abs.set_w_res(w_res);

% state priorities
% x_priority = [0;  % smaller headway more unsafe
%               1;  % larger ego velocity more unsafe
%               0]; % smaller lead velocity more unsafe
x_priority = [0;  % smaller ego position more unsafe
              0;  % smaller ego velocity more unsafe
              1]; % larger oncoming position more unsafe
veh_follow_abs.set_x_priority(x_priority);

% input priorities
% u_priority = 1; % larger ego acceleration more unsafe
u_priority = 0; % smaller ego acceleration more unsafe
veh_follow_abs.set_u_priority(u_priority);

% shouldn't need this?
% disturbance priorities
% w_priority = 0; % smaller lead acceleration more unsafe
% veh_follow_abs.set_w_priority(w_priority);

end
