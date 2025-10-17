%%% Plot 3D invariant set after computing it %%%

% load abstraction (with safe set) and modelling constants
% load('veh_turn_abs_oncoming_first_5_7_2021.mat');
con = veh_turn_const_v4();

% for plotting
S = con.s_min:con.s_res:con.s_max;
V = con.v_min:con.v_res:con.v_max;
S0 = con.s0_min:con.s0_res:con.s0_max;
[X,Y] = meshgrid(S, S0);

Z = con.v_min*ones(size(X));

% plot 
for i = 1:veh_turn_abs.basis_size
    x_idx = veh_turn_abs.get_basis_idx(i);
    x_val = veh_turn_abs.get_priority_state_at_idx(x_idx);
    % since S0 states with small values have higher priority, abstraction
    % indices are reversed compared to the indices used for plotting
    S_idx = x_idx(1);
    % no need to do reversal for S state (small values have low priority)
    S0_idx = size(Z,1) - x_idx(3) + 1;
    Z(S0_idx, S_idx) = x_val(2);
end

% traverse plot along the S direction
% fill in missing values, since we know the safe set boundary (i.e. the
% value of v) will monotonically increases
for i = 1:length(S0)
    max_v = con.v_min;
    for j = 1:length(S)
        % v increases as s decreases
        S_idx = size(Z,2) - j + 1;
        Z(i, S_idx) = max(Z(i, S_idx), max_v);
        max_v = Z(i, S_idx);
    end
end

% same, traversing along the S0 direction
for i = 1:length(S)
    max_v = con.v_min;
    for j = 1:length(S0)
        % v increases as s0 increases => need to reverse s0 index
        Z(j, i) = max(Z(j, i), max_v);
        max_v = Z(j, i);
    end
end

% plot final set
h = surf(S, S0, Z);
hold on;
set(h,'LineStyle','none');
set(gca,'FontSize',32);
% specify axis range, tick values
xlim([con.s_min, con.s_max]);
ylim([con.s0_min, con.s0_max]);
zlim([con.v_min, con.v_max]);
xticks([-60 -40 -20]);
yticks([-60 -40 -20 0]);
zticks([0 2 4 6 8 10 12]);
% title('Safe set for unprotected left turn scenario (oncoming first)');
xlabel('Ego Position (m)','Interpreter','tex');
ylabel('Oncoming Position (m)','Interpreter','tex');
zlabel('Ego Velocity (m/s)','Interpreter','tex');
view(55,35);