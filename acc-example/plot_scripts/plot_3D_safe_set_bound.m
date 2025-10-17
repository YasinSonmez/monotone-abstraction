%%% Plot invariant set boundary %%%
[con, veh_follow_abs] = get_3D_monotone_abs();

% for plotting
H = con.h_min:con.h_res:con.h_max;
V_L = con.v_min:con.v_res:con.v_max;
[X,Y] = meshgrid(H,V_L);

Z = zeros(size(X));

for i = 1:length(V_L)
    for j = 1:length(H)
        h_val = H(j);
        v_L_val = V_L(i);
        % plot directly from safe-set definition
        v_val = outer_approx_boundary(h_val, ...
                                      v_L_val, ...
                                      abs(con.a_lead_brake), ...
                                      abs(con.a_ego_brake), ...
                                      con.h_min);
        v_val = min(v_val, con.v_max);
        Z(i,j) = v_val;
    end
end

% plot final set
h = surf(H,V_L,Z);
hold on;
set(h,'LineStyle','none');
caxis([con.v_min con.v_max]);
title('Safe set bound for vehicle-following scenario');
xlabel('Headway (m)');
ylabel('Lead Velocity (m/s)');
zlabel('Ego Velocity (m/s)');
