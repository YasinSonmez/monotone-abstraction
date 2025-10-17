classdef box_abstraction < handle
    
    % BOX_ABSTRACTION: Creates a box abstraction of a dynamical system.
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
        
        % whether the lower / upper bound of the state, input, and
        % and disturbance domain intervals are unbounded or not
        % 1 indicates that the lower / upper bound is unbounded
        % 0 indicates that the lower / upper bound is not unbounded
        x_unbounded;
        u_unbounded;
        w_unbounded;
        
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
        
        % cached transitions using min / max state, input, and
        % disturbance values
        min_transitions = [];
        max_transitions = [];
        
    end
    
    methods
        %% constructor
        function abs = box_abstraction(dyn)
            abs.dyn = dyn;
            % initialize matrices min / max transitions using the correct
            % state dim.
            abs.min_transitions = zeros(dyn.n_x, 0);
            abs.max_transitions = zeros(dyn.n_x, 0);
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
            threshold = 1e-6;
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
            % create min / max transitions cell for caching transitions 
            % from states with min / max input and disturbance values
            if(~isempty(abs.u_numCells))
                abs.min_transitions = zeros(abs.dyn.n_x, prod(abs.x_numCells) * abs.u_numCells);
                abs.max_transitions = zeros(abs.dyn.n_x, prod(abs.x_numCells) * abs.u_numCells);
            end

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
            % create min / max transitions cell for caching transitions 
            % from states with min / max input and disturbance values
            if(~isempty(abs.x_numCells))
                abs.min_transitions = zeros(abs.dyn.n_x, prod(abs.x_numCells) * abs.u_numCells);
                abs.max_transitions = zeros(abs.dyn.n_x, prod(abs.x_numCells) * abs.u_numCells);
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
        
        %% check if range input is valid
        function check_z_unbounded(~, z_dim, z_unbounded)
            ok_rows = sum(z_unbounded(:,1) <= z_unbounded(:,2));
            if(~(size(z_unbounded,1) == z_dim && size(z_unbounded,2) == 2 && ok_rows == z_dim))
                error(['Input should be a ', num2str(z_dim), ' by 2 vector, with min (max) values in left (right) column).']);
            end
        end
        
        % set the range of values for each state dimension
        function abs = set_x_unbounded(abs, x_unbounded)
            abs.check_z_range(abs.x_dim, x_unbounded);
            abs.x_unbounded = x_unbounded;
        end
        
        % set the range of values for each input dimension
        function abs = set_u_unbounded(abs, u_unbounded)
            abs.check_z_range(abs.u_dim, u_unbounded);
            abs.u_unbounded = u_unbounded;
        end
        
        % set the range of values for each disturbance dimension
        function abs = set_w_unbounded(abs, w_unbounded)
            abs.check_z_range(abs.w_dim, w_unbounded);
            abs.w_unbounded = w_unbounded;
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
        
        %% helper function for getting the maximum abstraction value for the current cell
        function z = get_max_value_at_idx(~, z_dim, z_num_cells, z_range, z_res, z_priority, z_idx)
            % get the MAXIMUM value (w.r.t the given priority)
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
        
        % get the max state value at a particular index
        function x = get_max_state_at_idx(abs, x_idx)
            % assuming state index is valid
            x = abs.get_max_value_at_idx(abs.x_dim, abs.x_num_cells, abs.x_range, abs.x_res, abs.x_priority, x_idx);
        end
        
        % get the max input value at a particular index
        function u = get_max_input_at_idx(abs, u_idx)
            % assuming input index is valid
            u = abs.get_max_value_at_idx(abs.u_dim, abs.u_num_cells, abs.u_range, abs.u_res, abs.u_priority, u_idx);
        end
        
        % get the max disturbance value at a particular index
        function w = get_max_disturbance_at_idx(abs, w_idx)
            % assuming disturbance index is valid
            w = abs.get_max_value_at_idx(abs.w_dim, abs.w_num_cells, abs.w_range, abs.w_res, abs.w_priority, w_idx);
        end
        
        %% helper function for getting the minimum abstraction value for the current cell
        function z = convert_max_to_min_value(~, z, z_dim, z_range, z_res, z_priority)
            % convert to the MINIMUM value (w.r.t. the given priority),
            % but ensure we don't go outside the domain
            for i = 1:z_dim
                if(z_priority(i) == 1)
                    % larger values have priority
                    z(i) = max(z(i) - z_res(i), z_range(i,1));
                else
                    % smaller values have priority
                    z(i) = min(z(i) + z_res(i), z_range(i,2));
                end
            end
        end
        
        % get the min state value at a particular index
        function x = get_min_state_at_idx(abs, x_idx)
            % assuming state index is valid
            x = abs.get_max_value_at_idx(abs.x_dim, abs.x_num_cells, abs.x_range, abs.x_res, abs.x_priority, x_idx);
            x = abs.convert_max_to_min_value(x, abs.x_dim, abs.x_range, abs.x_res, abs.x_priority);
        end
        
        % get the min input value at a particular index
        function u = get_min_input_at_idx(abs, u_idx)
            % assuming input index is valid
            u = abs.get_max_value_at_idx(abs.u_dim, abs.u_num_cells, abs.u_range, abs.u_res, abs.u_priority, u_idx);
            u = abs.convert_max_to_min_value(u, abs.u_dim, abs.u_range, abs.u_res, abs.u_priority);
        end
        
        % get the min disturbance value at a particular index
        function w = get_min_disturbance_at_idx(abs, w_idx)
            % assuming disturbance index is valid
            w = abs.get_max_value_at_idx(abs.w_dim, abs.w_num_cells, abs.w_range, abs.w_res, abs.w_priority, w_idx);
            w = abs.convert_max_to_min_value(w, abs.w_dim, abs.w_range, abs.w_res, abs.w_priority);
        end
        
        %% helper function for converting a state / input / disturbance value to an index
        function z_idx = get_z_idx(~, z_dim, z_range, z_res, z_priority, z_unbounded, z_val)
            % for each state, check if its lower / upper bound is actually
            % unbounded (for example, in the vehicle following example, 
            % the domain for the headway h is [0, \inf). if so, then we
            % have to do thresholding to ensure we quantize correctly (see
            % [15, Section V-B-3] for details, as cited in our paper)
            for i = 1:z_dim
                if(z_unbounded(i,1))
                    z_val(i) = max(z_val(i), z_range(i,1));
                elseif(z_unbounded(i,2))
                    z_val(i) = min(z_val(i), z_range(i,2));
                end
            end
            % using threshold here to prevent floating point errors
            threshold = 1e-6;
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
            x_idx = abs.get_z_idx(abs.x_dim, abs.x_range, abs.x_res, abs.x_priority, abs.x_unbounded, x_val);
        end
        
        % convert input value to input index (in the abstraction)
        function u_idx = get_input_idx(abs, u_val)
            u_idx = abs.get_z_idx(abs.u_dim, abs.u_range, abs.u_res, abs.u_priority, abs.u_unbounded, u_val);
        end
        
        % convert input value to input index (in the abstraction)
        function w_idx = get_disturbance_idx(abs, w_val)
            w_idx = abs.get_z_idx(abs.w_dim, abs.w_range, abs.w_res, abs.w_priority, abs.w_unbounded, w_val);
        end
        
        %% functions for invariant set computation
        
        % helper function for using the min / max transitions objects
        function flat_x_u_idx = flatten_x_u_idx(abs, x_idx, u_idx)
            % TODO: implement this function for any state / input dim.
            % (currently implemented for special case n = 3, m = 1)
            if(~(size(x_idx, 1) == 3 && size(x_idx, 2) == 1))
                error('Only implemented for n_x = 3!');
            end
            if(~(size(u_idx, 1) == 1 && size(u_idx, 1) == 1))
                error('Only implemented for n_u = 1!');
            end
            flat_x_u_idx = ...
                  x_idx(1) + ...
                 (x_idx(2) - 1) * abs.x_numCells(1) + ...
                 (x_idx(3) - 1) * abs.x_numCells(1) * abs.x_numCells(2) + ...
                 (u_idx - 1) * prod(abs.x_numCells);
        end
        
        % get the index of the min state transitioned to
        function x_plus_idx = get_min_transition_state(abs, x_idx, u_idx, w_idx)
            % check if we already computed this transition before
            flat_x_u_idx = abs.flatten_x_u_idx(x_idx, u_idx);
            threshold = 1e-6;
            if(norm(abs.min_transitions(:, flat_x_u_idx)) > threshold)
                % we already computed this transition before
                % NOTE: we assume u_idx and w_idx will always be the same 
                % everytime this function is called, therefore we can
                % cache transitions
                % TODO: add a check to make sure the cache is being
                % used correctly
                x_plus_idx = abs.min_transitions(:, flat_x_u_idx);
            else
                % we haven't computed this transition before
                % get next state value
                x_val = abs.get_min_state_at_idx(x_idx);
                u_val = abs.get_min_input_at_idx(u_idx);
                w_val = abs.get_min_disturbance_at_idx(w_idx);
                x_plus_val = abs.dyn.next_state(x_val, u_val, w_val);
                % convert to next state index
                x_plus_idx = abs.get_state_idx(x_plus_val);
                % store this transition so we don't need to recompute it
                abs.min_transitions(:, flat_x_u_idx) = x_plus_idx;
            end
        end
        
        % get the index of the max state transitioned to
        function x_plus_idx = get_max_transition_state(abs, x_idx, u_idx, w_idx)
            % check if we already computed this transition before
            flat_x_u_idx = abs.flatten_x_u_idx(x_idx, u_idx);
            threshold = 1e-6;
            if(norm(abs.max_transitions(:, flat_x_u_idx)) > threshold)
                % we already computed this transition before
                % NOTE: we assume u_idx and w_idx will always be the same 
                % everytime this function is called, therefore we can
                % cache transitions
                % TODO: add a check to make sure the cache is being
                % used correctly
                x_plus_idx = abs.max_transitions(:, flat_x_u_idx);
            else
                % we haven't computed this transition before
                % get next state value
                x_val = abs.get_max_state_at_idx(x_idx);
                u_val = abs.get_max_input_at_idx(u_idx);
                w_val = abs.get_max_disturbance_at_idx(w_idx);
                x_plus_val = abs.dyn.next_state(x_val, u_val, w_val);
                % convert to next state index
                x_plus_idx = abs.get_state_idx(x_plus_val);
                % store this transition so we don't need to recompute it
                abs.max_transitions(:, flat_x_u_idx) = x_plus_idx;
            end
        end
        
        % clear cached transition values
        function abs = clear_cached_transitions(abs)
            % create a new min transitions cell for caching transitions
            abs.min_transitions = zeros(abs.dyn.n_x, prod(abs.x_numCells));
        end
        
    end
    
end
