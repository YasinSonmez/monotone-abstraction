%%% Maximal safety controller algorithm %%%
load('veh_turn_abs_ego_first_5_7_2021.mat');
disp('Loaded precomputed controllable set and monotone abstraction.');

% check if transitions stay in the controllable set
use_controllable_set = true;
% store current safe set and set controllable safe set variable
veh_turn_abs.store_safe_set();
veh_turn_abs.set_controllable_set();
% compute the remaining safe sets
for u_idx = 2:veh_turn_abs.u_num_cells
    disp(['Computing safe set for input index ', num2str(u_idx), ...
          ' of ', num2str(veh_turn_abs.u_num_cells), '.']);
    veh_turn_abs.clear_cached_transitions();
    veh_turn_abs.compute_safe_set(u_idx, use_controllable_set);
    veh_turn_abs.store_safe_set();
    clc;
end

% set state space dimensions
x_numCells = veh_turn_abs.x_num_cells;
s_numCells = x_numCells(1);
v_numCells = x_numCells(2);
s0_numCells = x_numCells(3);
% construct maximal safety controller from safe sets computed above
ego_first_controller = zeros(s_numCells, v_numCells, s0_numCells);
total_idx = prod(x_numCells);
idx = 0;
for s_idx = 1:s_numCells
    for v_idx = 1:v_numCells
        for s0_idx = 1:s0_numCells
            clc;
            idx = idx + 1;
            percent_done = round(idx / total_idx * 100, 2);
            disp('Constructing maximal safety controller...');
            disp(['Current state index: (s, v, s0) = (', ...
                   num2str(s_idx), ',', ...
                   num2str(v_idx), ',', ...
                   num2str(s0_idx), ')']);
            disp(['Percent done: ', num2str(percent_done), '%']);
            x_idx = [s_idx; v_idx; s0_idx];
            set_idx = 1;
            while(set_idx <= veh_turn_abs.u_num_cells && ...
                  veh_turn_abs.x_in_safe_set(x_idx, set_idx))
                set_idx = set_idx + 1;
            end
            ego_first_controller(s_idx, v_idx, s0_idx) = set_idx - 1;
        end
    end
end
