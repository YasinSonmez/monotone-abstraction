classdef monotone_dyn < handle
    
    % MONOTONE_ABSTRACTION: Creates an monotone dynamical system object.
    % ==================================================================
    % 
    % SYNTAX
    % ------
    %
    %   dyn = monotone_dynamics(n_x, n_u, n_w, f_expr)
    %
    % DESCRIPTION
    % -----------
    %   An instance of this class represents a monotone dynamical system.
    % 
    % INPUT
    % -----
    %   n_x     Integer representing system state dimension
    %
    %   n_u     Integer representing system input dimension
    %
    %   n_w     Integer representing system disturbance dimension
    %
    %   f_expr  String representing the system dynamics
    %           e.g. for an integrator system: '[x(2); u(1)]'
    %           to represent: x'(1) = x(2),
    %                         x'(2) = u(1)
    %
    % METHODS
    % -------
    %   dyn.next_state(x, u, w)  get the state transitioned to from state
    %                            x, with input u and disturbance w
    %
    %   dyn.get_state_dim()  get the state dimension of the system
    %
    %   dyn.get_input_dim()  get the input dimension of the system
    %
    %   dyn.get_disturbance_dim()  get the disturbance dimension for the
    %                              system
    %
    
    properties (SetAccess=protected)
        n_x;            % state dimension
        n_u;            % input dimension
        n_w;            % disturbance dimension
        f_expr;         % transition function
    end
    
    methods
        
        % Constructor
        function dyn = monotone_dyn(n_x, n_u, n_w, f_expr)
            % system dimensions
            dyn.n_x = n_x;
            dyn.n_u = n_u;
            dyn.n_w = n_w;
            % f_expr is the name of the user-defined dynamics function
            dyn.f_expr = f_expr;
        end
        
        % transition function
        function x_plus = next_state(dyn, x, u, w)
            % check to ensure state, input, disturbance
            % dimensions are respected
            if(size(x,1) ~= dyn.n_x || size(x,2) ~= 1)
                error(['Invalid state dimension, expected ', num2str(dyn.n_x), ' by 1 column vector.']);
            end
            if(size(u,1) ~= dyn.n_u || size(u,2) ~= 1)
                error(['Invalid input dimension, expected ', num2str(dyn.n_u), ' by 1 column vector.']);
            end
            if(size(w,1) ~= dyn.n_w || size(w,2) ~= 1)
                error(['Invalid disturbance dimension, expected ', num2str(dyn.n_w), ' by 1 column vector.']);
            end
            % evaluate the user-defined dynamics function
            x_plus = zeros(3,1);
            try
                eval(['x_plus = ', dyn.f_expr, '(x, u, w);']);
            catch
                error('Error using dynamics expression - are you sure your system dimensions match?');
            end
        end
        
        % get state dimension
        function ret = get_state_dim(dyn)
            ret = dyn.n_x;
        end
        
        % get input dimension
        function ret = get_input_dim(dyn)
            ret = dyn.n_u;
        end
        
        % get state dimension
        function ret = get_disturbance_dim(dyn)
            ret = dyn.n_w;
        end
        
    end
    
end
