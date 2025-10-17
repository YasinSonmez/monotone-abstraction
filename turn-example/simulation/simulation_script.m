%%%  script for running the simulation  %%%
% (need to do this because simulink is having an issue with the hybrid
%  vehicle dynamics for some reason...)
mdl = run_setup_scenario1;
out.x_t = mdl.x0';
out.u_t = zeros(0, 1);
out.u_bounds = zeros(0, 2);

for t = 0:0.1:6
  
  % extract simulation state
  x_sim = out.x_t(end, :)';
    
  if(mod(t, mdl.dt) == 0)
      
      % use the MPC on this time-step!
      
      % extract state
      s_val = x_sim(1);
      v_val = x_sim(2);
      s2_val = x_sim(4);
      
      % velocity bounds not violated
      assert(0 <= v_val && v_val <= 12);
      
      % get rid of weird initialization bug
      if(norm([s_val; v_val; s2_val]) < 1e-5)
          % set to initial conditions in this case
          s_val = mdl.x0(1);
          v_val = mdl.x0(2);
          s2_val = mdl.x0(4);
      end
      
      if(s_val > 10)
          % already passed through the intersection - no constraints needed
          torque_lower_val = mdl.T_min;
          torque_upper_val = mdl.T_max;
          u_bounds = [mdl.T_min; mdl.T_max];
      else
          % get correct index for feasible inputs structure
          x_val = [s_val; v_val; s2_val];
          x_idx = mdl.abs_ego_first.get_state_idx(x_val);
          cell_idx = num2cell(x_idx);
          
          % get input bounds (for plotting)
          torque_lower_bound_idx = max(mdl.feasible_inputs{cell_idx{:}});
          torque_upper_bound_idx = min(mdl.feasible_inputs{cell_idx{:}});
          torque_lower_bound = ...
              mdl.abs_ego_first.get_priority_input_at_idx(torque_lower_bound_idx);
          torque_upper_bound = ...
              mdl.abs_ego_first.get_priority_input_at_idx(torque_upper_bound_idx);
          u_bounds = [torque_lower_bound; torque_upper_bound];
          
          % get a range of feasible inputs (from incremental synthesis)
          lower_input_idx = mdl.feasible_inputs{cell_idx{:}}(1);
          upper_input_idx = lower_input_idx;
          i = 1;
          while((i < length(mdl.feasible_inputs{cell_idx{:}})) && ...
                  (upper_input_idx - 1 == mdl.feasible_inputs{cell_idx{:}}(i + 1)))
              upper_input_idx = upper_input_idx - 1;
              i = i + 1;
          end
          torque_lower_val = ...
              mdl.abs_ego_first.get_priority_input_at_idx(lower_input_idx);
          torque_upper_val = ...
              mdl.abs_ego_first.get_priority_input_at_idx(upper_input_idx);
      end
      
      % set input, using quadratic programming (1-step MPC here)
      
      % the MPC problem below uses a kinematic model with vehicle acceleration
      % as the control input for simplicity
      
      % first, compute the lower / upper bounds on acceleration from the
      % lower / upper bounds on wheel torque
      
      % a = (1 / M) * (T / R_w - Ff)
      Ff = mdl.alpha + mdl.beta * v_val + mdl.gamma * v_val^2;
      accel_lower_val = (1 / mdl.M) * (torque_lower_val / mdl.R_w - Ff);
      accel_upper_val = (1 / mdl.M) * (torque_upper_val / mdl.R_w - Ff);
      
      % state cost
      x0 = x_sim(1:2);
      v_des = 11.5;
      Q = [0  0;
          0  1];
      Q = 1.0 * Q;
      q = [0; -v_des];
      q = 1.0 * q;
      Qf = Q;
      qf = q;
      % input cost
      R = 1.0;
      r = 0;
      
      % dynamics constraint matrices
      dt = 0.5;
      A = [1  dt;
          0  1];
      B = [0.5 * dt^2; dt];
      theta = [0; 0];
      
      %%% SYSTEM_CONSTRAINTS %%%
      constraints = [];
      n = size(x0,1); % state dimension
      m = size(B,2); % input dimension
      x_bar = []; % x = [x(t+1); x(t+2); ... x(t+T)]
      u_bar = []; % u = [u(t); u(t+1); ... u(t+T-1)] (both col vectors)
      T = 5; % 5 step time horizon
      
      % create state decision variables for each time index
      for i = 1:T
          x{i} = sdpvar(n,1);
          u{i} = sdpvar(m,1);
          x_bar = [x_bar; x{i}];
          u_bar = [u_bar; u{i}];
      end
      
      % require that system updates satisfy x(t+1) = Ax(t) + Bu(t) + theta
      [G, L] = make_sys_constr(T, A, B, theta, x0);
      constraints = [constraints, x_bar <= G*u_bar + L];
      constraints = [constraints, x_bar >= G*u_bar + L];
      
      %%% CONSTRAINTS %%%
      % input constraints
      Hu = [1; -1];
      hu = [accel_upper_val; -accel_lower_val];
      % require that inputs satisfy Hu*u(t) <= hu for all t
      [Hu_bar, hu_bar] = make_input_constr(T, Hu, hu);
      constraints = [constraints, Hu_bar*u_bar <= hu_bar];
      
      %%% OBJECTIVE_FUNCTION %%%
      [Q_bar, q_bar, R_bar, r_bar] = make_QP_costs(T,Q,Qf,q,qf,R,r);
      obj_fun = 1/2*(x_bar'*Q_bar*x_bar + u_bar'*R_bar*u_bar) + ...
          q_bar'*x_bar + r_bar'*u_bar;
      
      %%% CALL SOLVER %%%
      optimize(constraints, obj_fun, sdpsettings('solver','quadprog'));
      accel_opt = value(u_bar(1));
      
      % now, convert the desired acceleration from the above MPC to a
      % desired torque to be implemented
      
      %     a = (1 / M) * (T / R_w - Ff)
      % ==> T = (a * M + Ff) * R_w
      u_opt = (accel_opt * mdl.M + Ff) * mdl.R_w;
      
  else
      
      % don't use MPC on this time-step!
      u_opt = out.u_t(end);
      u_bounds = out.u_bounds(end, :)';
      
  end
  
  % record input, simulate vehicle dynamics
  out.u_t = [out.u_t; u_opt];
  out.u_bounds = [out.u_bounds; u_bounds'];
  if(t < 6)
    out.x_t = [out.x_t; sys_dyn(x_sim, u_opt)'];
  end
  
end

function x_plus = sys_dyn(x_val, u_val)

    % get the modelling parameters
    mdl = run_setup_scenario1;

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
    a = (1 / mdl.M) * (u_val / mdl.R_w - mdl.alpha);
    b = - (1 / mdl.M) * mdl.beta;
    c = - (1 / mdl.M) * mdl.gamma;

    % using numerical solution to vehicle dynamics eqn, compute next state
    x_plus = zeros(4, 1);
    [x_plus(1), x_plus(2)] = veh_dyn(x_val(1), x_val(2), a, b, c, mdl.dt_sim);
    x_plus(3) = x_val(3) + mdl.v0_min * mdl.dt_sim;
    x_plus(4) = x_val(4) + mdl.v0_max * mdl.dt_sim;

    % do thresholding on the ego vehicle velocity
    % (this is necessary because sometimes the ode45 solver lets
    % the ego vehicle go slightly outside of its velocity bounds)
    x_plus(2) = min(max(x_plus(2), mdl.v_min), mdl.v_max);
    
end

% compute next state
function [pos_plus, vel_plus] = veh_dyn(p0, v0, a, b, c, dt)
    
    % use ode45 solver to get state at next time step
    [t, x] = ode45(@f, [0 dt], [p0; v0; a; b; c]);
    
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
    % x(5) = c
    
    % equations of motion, with state defined as above:
    dxdt = [x(2);
            x(3) + x(4) * x(2) + x(5) * x(2)^2;
            0;
            0;
            0];
    
    % ensure the vehicle does not reverse
    if(x(2) <= 0 && dxdt(2) < 0)
        dxdt(2) = 0;
    end
    
    % ensure the vehicle does not exceed its maximum velocity
    if(x(2) >= 12 && dxdt(2) > 0)
        dxdt(2) = 0;
    end
    
end
