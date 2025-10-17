function ego_controller(block)

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 2;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  block.InputPort(1).Dimensions        = 4;
  block.InputPort(1).DirectFeedthrough = false;
  
  block.OutputPort(1).Dimensions       = 1;
  block.OutputPort(2).Dimensions       = 2;
  
  %% Set block sample time to inherited
  block.SampleTimes = [-1 0];
  
  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Run accelerator on TLC
  block.SetAccelRunOnTLC(true);
  
  %% Register methods
  block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
  block.RegBlockMethod('InitializeConditions',     @InitConditions);  
  block.RegBlockMethod('Outputs',                  @Output);  
  
%endfunction

function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = fd;
  block.OutputPort(1).SamplingMode  = fd;
  block.OutputPort(2).SamplingMode  = fd;
  
%endfunction

function InitConditions(block) 
  block.OutputPort(1).Data = 0;
  block.OutputPort(2).Data = [0; 0];
  
%endfunction

function Output(block)

  % input ports
  global mdl;
  x_sim = block.InputPort(1).Data;
  
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
      T_bounds = [mdl.T_min; mdl.T_max];
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
      T_bounds = [torque_lower_bound; torque_upper_bound];
      
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
  v_des = 9;
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
  T_opt = (accel_opt * mdl.M + Ff) * mdl.R_w;
  
  block.OutputPort(1).Data = T_opt;
  block.OutputPort(2).Data = T_bounds;
  
%endfunction
