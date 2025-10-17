% get abstraction
[~, veh_follow_box_abs] = get_3D_box_abs_v2();

% size of state-space grid
x_num_cells = veh_follow_box_abs.x_num_cells;

% load the precomputed feasible_inputs cell
% load('feasible_inputs_5_3_2021.mat')

% first ensure that the set is well-structured, essentially meaning that
% it is a convex set (we expect this due to the monotonicity property)
for h_idx = 1:x_num_cells(1)
    for v_idx = 1:x_num_cells(2)
        for v_L_idx = 1:x_num_cells(3)
            % check all points adjacent to the current cell which have a 
            % smaller index
            h_check_idx = max(h_idx - 1, 1);
            v_check_idx = max(v_idx - 1, 1);
            v_L_check_idx = max(v_L_idx - 1, 1);
            cell_idx = num2cell([h_idx; v_idx; v_L_idx]);
            cell_check1_idx = num2cell([h_check_idx; v_idx; v_L_idx]);
            cell_check2_idx = num2cell([h_idx; v_check_idx; v_L_idx]);
            cell_check3_idx = num2cell([h_idx; v_idx; v_L_check_idx]);
            % check if any of the adjacent cells are not in the safe set,
            % but the current cell is
            if(~isempty(feasible_inputs{cell_idx{:}}) && ...
                    (isempty(feasible_inputs{cell_check1_idx{:}}) || ...
                     isempty(feasible_inputs{cell_check2_idx{:}}) || ...
                     isempty(feasible_inputs{cell_check3_idx{:}})))
                 % we can't plot this set - it is not convex!
                error('Safe set is not well-structured!');
            end
        end
    end
end

% we are good!
disp('Safe set is well-structured!');

% get monotone abstraction
[~, veh_follow_abs] = get_3D_monotone_abs_v2();

for h_idx = x_num_cells(1):-1:1
    for v_L_idx = x_num_cells(3):-1:1
        % first check if the safe set lies over the states corresponding
        % to (h_idx, v_L_idx)
        v_idx = 1;
        cell_idx = num2cell([h_idx; v_idx; v_L_idx]);
        if(isempty(feasible_inputs{cell_idx{:}}))
            continue;
        end
        % find upper boundary of safe set
        while(v_idx < x_num_cells(2) && ...
                ~isempty(feasible_inputs{cell_idx{:}}))
            v_idx = v_idx + 1;
            cell_idx = num2cell([h_idx; v_idx; v_L_idx]);
        end
        % now, add the boundary point to the safe set basis
        veh_follow_abs.add_basis_idx([h_idx; v_idx; v_L_idx]);
    end
end
