classdef monotone_abstraction < handle
    
    % MONOTONE_ABSTRACTION: Creates an abstraction of a dynamical system
    %                       with monotone dynamics.
    % ==================================================================
    % 
    % SYNTAX
    % ------
    %
    %   ma = monotone_abstraction()
    %
    % Description
    % -----------
    %   An instance of this class represents a finite abstraction of a
    %   system with monotone dynamics
    %
    
    properties (SetAccess=protected)
        % monotone dynamical system
        dyn;
        
        % range of admissible values for the state, input, and disturbance
        x_range;
        u_range;
        w_range;
        
        % resolution for state, input, and disturbance
        x_res;
        u_res;
        w_res;
        
        % number of cells across one dimension of the abstraction
        % for example, the state abstraction is a
        %        x_cellRange(1) by x_cellRange(2) by x_cellRange(3)
        % boolean matrix for a system with 3D state.
        % also, the resolution for each dimension of the abstraction
        x_numCells;
        u_numCells;
        w_numCells;
        
        % priority direction for state, input, and disturbance
        % 1 indicates that larger numbers have priority
        % 0 indicates that smaller numbers have priority
        x_priority;
        u_priority;
        w_priority;
        
        % controllable set - set of states where we can enforce safety
        controllable_set_basis = [];
        
        % current set - to be used for invariant set computation
        safe_set_basis = [];
        
        % store previously computed safe sets here
        safe_set_list = {};
        
        % transitions from states using min input and max disturbance
        % values
        min_transitions = [];
        
    end
    
    methods
        %% constructor
        function abs = monotone_abstraction(dyn)
            abs.dyn = dyn;
            % initialize matrices for controllable set basis,
            % safe set basis, min transitions using the correct state dim.
            abs.controllable_set_basis = zeros(dyn.n_x, 0);
            abs.safe_set_basis = zeros(dyn.n_x, 0);
            abs.min_transitions = zeros(dyn.n_x, 0);
        end
        
        function ret = x_dim(abs)
            ret = abs.dyn.get_state_dim();
        end
        
        function ret = u_dim(abs)
            ret = abs.dyn.get_input_dim();
        end
        
        function ret = w_dim(abs)
            ret = abs.dyn.get_disturbance_dim();
        end
        
        function ret = x_num_cells(abs)
            ret = abs.x_numCells;
        end
        
        function ret = u_num_cells(abs)
            ret = abs.u_numCells;
        end
        
        function ret = w_num_cells(abs)
            ret = abs.w_numCells;
        end
        
        %% helper function for setting abstraction resolution
        function z_numCells = set_res_helper(~, z_dim, z_range, dz)
            
            if(~(size(dz, 2) == 1 && (size(dz, 1) == 1 || size(dz, 1) == z_dim)))
                error(['Input should either be a scalar (each state dimension has same resolution) or a ', num2str(z_dim), ' by 1 vector otherwise.']);
            end
            
            if(isempty(z_range))
                error('The range of this variable must be set before setting its resolution.');
            end
            
            z_numCells = zeros(z_dim,1);
            % adding threshold here to prevent floating point errors
            threshold = 1e-5;
            for i = 1:z_dim
                if(size(dz, 1) == 1)
                    % each dimension has same resolution
                    idx = 1;
                else
                    % each dimension has different resolution
                    idx = i;
                end
                z_numCells(i) = ...
                    ceil((z_range(i,2) - z_range(i,1)) / dz(idx)) + 1;
                if(mod(z_range(i,2) - z_range(i,1), dz(idx)) > threshold)
                    error('The range of this variable cannot be spanned by an integer number of cells with the given resolution.');
                end
            end
            
        end
        
        % set resolution of state abstraction
        function abs = set_x_res(abs, dx)
            abs.x_numCells = abs.set_res_helper(abs.x_dim, abs.x_range, dx);
            if(size(dx,1) == 1)
                abs.x_res = dx*ones(abs.x_dim,1);
            else
                abs.x_res = dx;
            end
            % create min transitions cell for caching transitions from
            % states with min input and disturbance values
            abs.min_transitions = zeros(abs.dyn.n_x, prod(abs.x_numCells));
        end
        
        % set resolution of input abstraction
        function abs = set_u_res(abs, du)
            u_dim = abs.u_dim;
            abs.u_numCells = abs.set_res_helper(u_dim, abs.u_range, du);
            if(size(du,1) == 1)
                abs.u_res = du*ones(u_dim,1);
            else
                abs.u_res = du;
            end
        end
        
        % set resolution of disturbance abstraction
        function abs = set_w_res(abs, dw)
            w_dim = abs.w_dim;
            abs.w_numCells = abs.set_res_helper(w_dim, abs.w_range, dw);
            if(size(dw,1) == 1)
                abs.w_res = dw*ones(w_dim,1);
            else
                abs.w_res = dw;
            end
        end
        
        %% check if range input is valid
        function check_z_range(~, z_dim, z_range)
            ok_rows = sum(z_range(:,1) <= z_range(:,2));
            if(~(size(z_range,1) == z_dim && size(z_range,2) == 2 && ok_rows == z_dim))
                error(['Input should be a ', num2str(z_dim), ' by 2 vector, with min (max) values in left (right) column).']);
            end
        end
        
        % set the range of values for each state dimension
        function abs = set_x_range(abs, x_range)
            abs.check_z_range(abs.x_dim, x_range);
            abs.x_range = x_range;
        end
        
        % set the range of values for each input dimension
        function abs = set_u_range(abs, u_range)
            abs.check_z_range(abs.u_dim, u_range);
            abs.u_range = u_range;
        end
        
        % set the range of values for each disturbance dimension
        function abs = set_w_range(abs, w_range)
            abs.check_z_range(abs.w_dim, w_range);
            abs.w_range = w_range;
        end
        
        %% helper for setting priorities (whether larger / smaller values are prioritized)
        function check_z_priority(~, z_dim, z_priority)
            if(~(size(z_priority,1) == z_dim && size(z_priority,2) == 1) || ...
               ~(sum(z_priority == 0 | z_priority == 1) == z_dim))
                error(['Input should be a ', num2str(z_dim), ' by 1 vector of 1s and 0s']);
            end
        end
        
        % set priority for each state dimension
        function abs = set_x_priority(abs, x_priority)
            abs.check_z_priority(abs.x_dim, x_priority);
            abs.x_priority = x_priority;
        end
        
        % set priority for each input dimension
        function abs = set_u_priority(abs, u_priority)
            abs.check_z_priority(abs.u_dim, u_priority);
            abs.u_priority = u_priority;
        end
        
        % set priority for each disturbance dimension
        function abs = set_w_priority(abs, w_priority)
            abs.check_z_priority(abs.w_dim, w_priority);
            abs.w_priority = w_priority;
        end
        
        %% helper function for getting the priority abstraction value
        function z = get_priority_value_at_idx(~, z_dim, z_num_cells, z_range, z_res, z_priority, z_idx)
            z = zeros(z_dim, 1);
            if(sum((ones(z_dim, 1) <= z_idx) & (z_idx <= z_num_cells)) ~= z_dim)
                error('Cannot get priority value at index - index out of range.');
            end
            for i = 1:z_dim
                if(z_priority(i) == 1)
                    % larger values have priority
                    z(i) = z_range(i,1) + (z_idx(i) - 1)*z_res(i);
                else
                    % smaller values have priority
                    z(i) = z_range(i,2) - (z_idx(i) - 1)*z_res(i);
                end
            end
        end
        
        % get the state value at a particular index
        function x = get_priority_state_at_idx(abs, x_idx)
            % assuming state index is valid
            x = abs.get_priority_value_at_idx(abs.x_dim, abs.x_num_cells, abs.x_range, abs.x_res, abs.x_priority, x_idx);
        end
        
        % get the input value at a particular index
        function u = get_priority_input_at_idx(abs, u_idx)
            % assuming input index is valid
            u = abs.get_priority_value_at_idx(abs.u_dim, abs.u_num_cells, abs.u_range, abs.u_res, abs.u_priority, u_idx);
        end
        
        % get the disturbance value at a particular index
        function w = get_priority_disturbance_at_idx(abs, w_idx)
            % assuming disturbance index is valid
            w = abs.get_priority_value_at_idx(abs.w_dim, abs.w_num_cells, abs.w_range, abs.w_res, abs.w_priority, w_idx);
        end
        
        %% helper function for converting a state / input / disturbance value to an index
        function z_idx = get_z_idx(~, z_dim, z_range, z_res, z_priority, z_val)
            % using threshold here to prevent floating point errors
            threshold = 1e-5;
            for i = 1:z_dim
                if(~(z_range(i,1) <= z_val(i) + threshold && ...
                     z_val(i) <= z_range(i,2) + threshold))
                    z_idx = -1;
                    return;
                end
            end
            z_idx = zeros(z_dim, 1);
            for i = 1:z_dim
                if(z_priority(i) == 1)
                    % larger values have priority
                    z_idx(i) = ceil((z_val(i) - z_range(i,1)) ...
                                        / z_res(i) - threshold) + 1;
                                    % subtract threshold to prevent
                                    % rounding errors. using 1 indexing.
                else
                    % smaller values have priority
                    z_idx(i) = ceil((z_range(i,2) - z_val(i)) ...
                                        / z_res(i) - threshold) + 1;
                                    % subtract threshold to prevent
                                    % rounding errors. using 1 indexing.
                end
            end
        end
        
        % convert state value to state index (in the abstraction)
        function x_idx = get_state_idx(abs, x_val)
            x_idx = abs.get_z_idx(abs.x_dim, abs.x_range, abs.x_res, abs.x_priority, x_val);
        end
        
        % convert input value to input index (in the abstraction)
        function u_idx = get_input_idx(abs, u_val)
            u_idx = abs.get_z_idx(abs.u_dim, abs.u_range, abs.u_res, abs.u_priority, u_val);
        end
        
        % convert input value to input index (in the abstraction)
        function w_idx = get_disturbance_idx(abs, w_val)
            w_idx = abs.get_z_idx(abs.w_dim, abs.w_range, abs.w_res, abs.w_priority, w_val);
        end
        
        %% functions for invariant set computation
        
        % initialize the safe set to be the entire state domain
        function abs = initialize_safe_set(abs)
            abs.safe_set_basis(:, 1) = abs.x_numCells;
        end
        
        % check if state x2 has higher priority than state x1
        function res = x1_leq_x2(abs, x1_idx, x2_idx)
            res = (sum(x1_idx <= x2_idx) == abs.x_dim);
        end
        
        % check if state x1 and x2 are equal
        function res = x1_eq_x2(abs, x1_idx, x2_idx)
            res = (sum(x1_idx == x2_idx) == abs.x_dim);
        end
        
        % check if a state is on the safe set boundary
        function res = x_on_safe_set_boundary(abs, x_idx)
            res = false;
            if(x_idx == -1)
                return;
            end
            for i = 1:size(abs.safe_set_basis, 2)
                if(abs.x1_eq_x2(x_idx, abs.safe_set_basis(:, i)))
                    res = true;
                    break;
                end
            end
        end
        
        % check if a state is in the controllable set
        function res = x_in_controllable_set(abs, x_idx)
            res = false;
            if(x_idx == -1)
                return;
            end
            for i = 1:length(abs.controllable_set_basis)
                if(abs.x1_leq_x2(x_idx, abs.controllable_set_basis(:, i)))
                    res = true;
                    break;
                end
            end
        end
        
        % check if a state is in the safe set
        function res = x_in_safe_set(abs, x_idx, varargin)
            res = false;
            if(x_idx == -1)
                return;
            end
            if(isempty(varargin))
            % check if the state is in the current safe set
                for i = 1:size(abs.safe_set_basis, 2)
                    if(abs.x1_leq_x2(x_idx, abs.safe_set_basis(:, i)))
                        res = true;
                        break;
                    end
                end
            elseif(length(varargin) == 1)
            % check if the state is in a safe set in our stored list
            %
            % note to self: think there was actually a bug here before.
            % we should be iterating from:
            % i = 1:size(abs.safe_set_list{set_idx}, 2)
            % as below. fixing this now - Stan 5/4/2021
                set_idx = varargin{1};
                for i = 1:size(abs.safe_set_list{set_idx}, 2)
                    if(abs.x1_leq_x2(x_idx, abs.safe_set_list{set_idx}(:, i)))
                        res = true;
                        break;
                    end
                end
            else
                error('Unexpected number of input arguments.');
            end
        end
        
        % remove an element from the safe set basis
        function abs = remove_basis_idx(abs, basis_idx)
            % find the basis element
            found = false;
            for i = 1:size(abs.safe_set_basis, 2)
                if(sum(abs.safe_set_basis(:, i) == basis_idx) == abs.x_dim)
                    removed = abs.safe_set_basis(:, i);
                    abs.safe_set_basis(:, i) = [];
                    found = true;
                    break;
                end
            end
            if(~found)
                disp('Element to be removed not found! Basis left unchanged.');
                return;
            end
            % add new element to the basis if needed
            for i = 1:abs.x_dim
                if(removed(i) == 1)
                    continue; % out of bounds
                else
                    translation = zeros(abs.x_dim,1);
                    translation(i) = 1;
                    candidate = removed - translation;
                    if(~abs.x_in_safe_set(candidate))
                        abs.safe_set_basis(:, end + 1) = candidate;
                    end
                end
            end
        end
        
        % get the number of elements in the basis
        function num = basis_size(abs)
            num = size(abs.safe_set_basis, 2);
        end
        
        % get the state index for a particular basis element
        function x_idx = get_basis_idx(abs, basis_idx)
            x_idx = abs.safe_set_basis(:, basis_idx);
        end
        
        % helper function for using the min transitions object
        function flat_x_idx = flatten_x_idx(abs, x_idx)
            % TODO: implement this function for any state dim.
            % (currently implemented for special case n = 3)
            if(~(size(x_idx, 1) == 3 && size(x_idx, 2) == 1))
                error('Only implemented for n_x = 3!');
            end
            flat_x_idx = ...
                  x_idx(1) + ...
                 (x_idx(2) - 1) * abs.x_numCells(1) + ...
                 (x_idx(3) - 1) * abs.x_numCells(1) * abs.x_numCells(2);
        end
        
        % get the index of the state transitioned to
        function x_plus_idx = get_transition_state(abs, x_idx, u_idx, w_idx)
            % check if we already computed this transition before
            flat_x_idx = abs.flatten_x_idx(x_idx);
            threshold = 1e-5;
            if(norm(abs.min_transitions(:, flat_x_idx)) > threshold)
                % we already computed this transition before
                % NOTE: we assume u_idx and w_idx will always be the same 
                % everytime this function is called, therefore we can
                % cache transitions
                % TODO: add a check to make sure the cache is being
                % used correctly
                x_plus_idx = abs.min_transitions(:, flat_x_idx);
            else
                % we haven't computed this transition before
                % get next state value
                x_val = abs.get_priority_state_at_idx(x_idx);
                u_val = abs.get_priority_input_at_idx(u_idx);
                w_val = abs.get_priority_disturbance_at_idx(w_idx);
                x_plus_val = abs.dyn.next_state(x_val, u_val, w_val);
                % convert to next state index
                x_plus_idx = abs.get_state_idx(x_plus_val);
                % store this transition so we don't need to recompute it
                abs.min_transitions(:, flat_x_idx) = x_plus_idx;
            end
        end
        
        % add a new element to the safe set basis
        function abs = add_basis_idx(abs, x_idx)
            % check to make sure the element isn't redundant
            if(~abs.x_in_safe_set(x_idx))
                abs.safe_set_basis(:, end + 1) = x_idx;
            end
        end
        
        % stores the current safe set in the safe set list
        function abs = store_safe_set(abs)
            abs.safe_set_list{end + 1} = abs.safe_set_basis;
        end
        
        % save the current safe set in the controllable set variable
        function abs = set_controllable_set(abs)
            abs.controllable_set_basis = abs.safe_set_basis;
        end
        
        % clear cached transition values
        function abs = clear_cached_transitions(abs)
            % create a new min transitions cell for caching transitions
            abs.min_transitions = zeros(abs.dyn.n_x, prod(abs.x_numCells));
        end
        
        % invariant set algorithm (implemented for special case)
        function abs = compute_safe_set(abs, u_idx, use_controllable_set)
            %  Make sure we are in the special case
            if(~(abs.u_dim == 1 && (abs.w_dim == 0 || abs.w_dim == 1)))
                error('Only implemented for u_dim = 1, w_dim = 0 or 1.');
            end
            % Invariant set algorithm (special case)
            while(true)
                unsafe_basis_indices = zeros(abs.x_dim, 0);
                for i = 1:abs.basis_size
                    % test points on the boundary of the safe set
                    x_idx = abs.get_basis_idx(i);
                    % get disturbance
                    if(abs.w_dim == 1)
                        % disturbance is lead vehicle acceleration force
                        w_idx = abs.w_num_cells;
                    else
                        % no disturbance in this scenario
                        w_idx = zeros(1,0);
                    end
                    x_plus_idx = abs.get_transition_state(x_idx, u_idx, w_idx);
                    % check if the point is safe
                    unsafe = ...
                        (use_controllable_set && ...
                         ~abs.x_in_controllable_set(x_plus_idx)) || ...
                        (~use_controllable_set && ...
                         ~abs.x_in_safe_set(x_plus_idx));
                    if(unsafe)
                        unsafe_basis_indices(:, end + 1) = x_idx;
                    end
                end
                if(isempty(unsafe_basis_indices))
                    % all basis elements are safe and we are done!
                    break;
                else
                    for i = 1:size(unsafe_basis_indices, 2)
                        abs.remove_basis_idx(unsafe_basis_indices(:, i));
                    end
                end
            end
        end
        
    end
    
end
