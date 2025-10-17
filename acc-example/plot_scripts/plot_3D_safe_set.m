%%% Plot 3D invariant set after computing it %%%
%%% (using script compute_3D_invariant_set.m) %%%

con = veh_follow_const_v2();

% for plotting
H = con.h_min:con.h_res:con.h_max;
V_L = con.v_min:con.v_res:con.v_max;
[X,Y] = meshgrid(H,V_L);

Z = zeros(size(X));

% plot 
for i = 1:veh_follow_abs.basis_size
    x_idx = veh_follow_abs.get_basis_idx(i);
    x_val = veh_follow_abs.get_priority_state_at_idx(x_idx);
    % since states with small values have higher priority, indices in the
    % abstraction are reversed compared to the indices used for plotting
    V_L_idx = size(Z,1) - x_idx(3) + 1;
    H_idx = size(Z,2) - x_idx(1) + 1;
    Z(V_L_idx,H_idx) = x_val(2);
end

% traverse plot along the H direction
% fill in missing values, since we know the safe set boundary (i.e. the
% value of v) will monotonically increase
for i = 1:length(V_L)
    max_v = 0;
    for j = 1:length(H)
        Z(i,j) = max(Z(i,j),max_v);
        max_v = Z(i,j);
    end
end

% same, traversing along the V_L direction
for i = 1:length(H)
    max_v = 0;
    for j = 1:length(V_L)
        Z(j,i) = max(Z(j,i),max_v);
        max_v = Z(j,i);
    end
end

% adjustment used to prevent 'z-stitching' for first example as discussed here:
% https://www.mathworks.com/matlabcentral/answers/6623-plotting-overlapping-surfaces
if(exist('adjustment','var'))
    Z = Z + adjustment;
end

% plot final set
h = surf(H,V_L,Z);
hold on;
set(h,'LineStyle','none');
set(gca,'FontSize',32);
% specify axis range, tick values
xlim([con.h_min, con.h_max]);
ylim([con.v_min, con.v_max]);
zlim([con.v_min, con.v_max + 0.1]);
xticks([0 20 40 60 80]);
yticks([0 5 10 15 20]);
zticks([0 5 10 15 20]);
% title('Safe set for vehicle-following scenario');
xlabel('Headway (m)','Interpreter','tex');
ylabel('Lead Velocity (m/s)','Interpreter','tex');
zlabel('Ego Velocity (m/s)','Interpreter','tex');