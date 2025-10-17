function vehicle_dynamics(block)

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register number of dialog parameters   
  block.NumDialogPrms = 0;

  %% Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = 4;
  block.InputPort(1).DirectFeedthrough = false;
  
  block.InputPort(2).Dimensions        = 1;
  block.InputPort(2).DirectFeedthrough = false;
  
  block.OutputPort(1).Dimensions       = 4;
  
  %% Set block sample time to continuous
  block.SampleTimes = [0 0];
  
  %% Setup Dwork
  block.NumContStates = 4;

  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivative);  
  
%endfunction

function InitConditions(block)

  %% Initialize Dwork
  global mdl;
  
  block.ContStates.Data = mdl.x0;
  
%endfunction

function Output(block)

  block.OutputPort(1).Data = block.ContStates.Data;
  
%endfunction

function Derivative(block)

  % input ports
  global mdl;
  x = block.InputPort(1).Data;
  u = block.InputPort(2).Data;

  % vehicle dynamics
  dx = zeros(4,1);
  % s (position of ego vehicle)
  dx(1) = x(2);
  % v (velocity of ego vehicle)
  Ff = mdl.alpha + mdl.beta * x(2) + mdl.gamma * x(2)^2;
  dx(2) = (1 / mdl.M) * (u / mdl.R_w - Ff);
  % s_1 (position of oncoming vehicle 1)
  dx(3) = mdl.v0_min;
  % s_2 (position of oncombin vehicle 2)
  dx(4) = mdl.v0_max;
  
  % ensure the vehicle does not reverse
  if(x(2) <= 0 && dx(2) < 0)
      dx(2) = 0;
  end
  
  % ensure the vehicle does not exceed its maximum velocity
  if(x(2) >= 12 && dx(2) > 0)
      dx(2) = 0;
  end

  block.Derivatives.Data = dx;
  
%endfunction