%%% Plot 3D invariant set after computing it %%%
load('veh_turn_abs_ego_first_5_7_2021.mat');
con = veh_turn_const_v3();

% for plotting
S = con.s_min:con.s_res:con.s_max;
V = con.v_min:con.v_res:con.v_max;
S0 = con.s0_min:con.s0_res:con.s0_max;
[X,Y] = meshgrid(S, S0);

Z = con.v_max*ones(size(X));

% plot 
for i = 1:veh_turn_abs.basis_size
    x_idx = veh_turn_abs.get_basis_idx(i);
    x_val = veh_turn_abs.get_priority_state_at_idx(x_idx);
    % since states with small values have higher priority, indices in the
    % abstraction are reversed compared to the indices used for plotting
    S_idx = size(Z,2) - x_idx(1) + 1;
    % no need to do reversal for S0 state (small values have low priority)
    S0_idx = x_idx(3);
    Z(S0_idx, S_idx) = x_val(2);
end

% traverse plot along the S direction
% fill in missing values, since we know the safe set boundary (i.e. the
% value of v) will monotonically decrease
for i = 1:length(S0)
    min_v = con.v_max;
    for j = 1:length(S)
        % v decreases as s increases
        Z(i, j) = min(Z(i, j), min_v);
        min_v = Z(i, j);
    end
end

% same, traversing along the S0 direction
for i = 1:length(S)
    min_v = con.v_max;
    for j = 1:length(S0)
        % v decreases as s0 decreases => need to reverse s0 index
        S0_idx = size(Z,1) - j + 1;
        Z(S0_idx, i) = min(Z(S0_idx, i), min_v);
        min_v = Z(S0_idx, i);
    end
end

% plot final set
figure();
h = surf(S, S0, Z);
hold on;
set(h,'LineStyle','none');
set(gca,'FontSize',32);
% specify axis range, tick values
xlim([con.s_min, con.s_max]);
ylim([con.s0_min, con.s0_max]);
zlim([con.v_min, con.v_max]);
xticks([-60 -40 -20 0]);
yticks([-70 -60 -50 -40 -30 -20 -10]);
zticks([0 2 4 6 8 10 12]);
% title('Safe set for unprotected left turn scenario (ego first)');
xlabel('Ego Position (m)','Interpreter','tex');
ylabel('Oncoming Position (m)','Interpreter','tex');
zlabel('Ego Velocity (m/s)','Interpreter','tex');
view(55,35);
