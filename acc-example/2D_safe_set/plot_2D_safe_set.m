% Plot the boundary of the invariant set

% motion constraints
h_min = 0; % m
h_max = 80; % m
v_min = 0; % m/s
v_max = 20; % m/s

% resolution
v_res = 0.1;
h_res = 0.2;

% for plotting
H = h_min:h_res:h_max;
V = v_min:v_res:v_max;
[X,Y] = meshgrid(H,V);

Z = zeros(size(X));

for i = 1:length(V)
    i
    for j = 1:length(H)
        h_val = H(j);
        v_val = V(i);
        x_val = [h_val; v_val];
        x_idx = veh_follow_abs.get_state_idx(x_val);
        Z(i,j) = veh_follow_abs.x_in_safe_set(x_idx);
    end
end

h = surf(H,V,Z);
hold on;
set(h,'LineStyle','none');
caxis([0 1]);

% add closed form solution plot boundary

V_closed = v_min:v_res:v_max;
H_closed = V_closed.^2/(2*3) + 5;

plot(H_closed, V_closed, 'LineWidth', 1.35, 'color', 'red');
title('Vehicle-following invariant set');
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
