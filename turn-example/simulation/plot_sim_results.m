% plot results from unprotected left turn example
t = 0:0.1:6;
spacing = 0.025;

% shift inputs by one sampling interval due to simulink timing weirdness.
% otherwise the plot is going to be very confusing due to the low 
% sampling frequency
% out.u_t = [out.u_t(6:end); (mdl.T_min + mdl.T_max) / 2 .* ones(5,1)];
% out.u_bounds = [out.u_bounds(6:end, :); ...
%                 mdl.T_min .* ones(5,1), mdl.T_max .* ones(5,1)];

% position
figure();
subaxis(3,1,1,'Spacing', spacing);
plot(t, out.x_t(:,1), '-', ...
     t, out.x_t(:,3), '--', ...
     t, out.x_t(:,4), '-.', 'LineWidth', 1.35);
hold on;
plot(t, -10*ones(size(t)), ':', ...
     t, 10*ones(size(t)), ':', 'Color', [0.65, 0, 1], 'LineWidth', 1.35);
ylabel('Pos. (m)');
% title('Unprotected Left Turn Scenario');
legend('Ego', 'Oncoming 1', '2', ...
       'Location', 'northoutside', 'Orientation', 'horizontal');
grid on;
ylim([-75, 25]);
yticks([-60 -40 -20 0 20]);
set(gca,'XTickLabel',[]);

% velocity
subaxis(3,1,2,'Spacing', spacing);
plot(t, out.x_t(:,2), '-', ...
     t, 8*ones(size(t)), '--', ...
     t, 12*ones(size(t)), '-.', 'LineWidth', 1.35);
grid on;
ylim([4, 16]);
yticks([5 10 15]);
ylabel('Vel. (m / s)');
set(gca,'XTickLabel',[]);

% torque
subaxis(3,1,3,'Spacing', spacing);
plot(t, out.u_t ./ 1000, '-', 'LineWidth', 1.35);
hold on;
plot(t, out.u_bounds(:,1) ./ 1000, ':', ...
     t, out.u_bounds(:,2) ./ 1000, ':', 'Color', [1 0 0], 'LineWidth', 1.35);
grid on;
ylim([-2, 2]);
yticks([-2 -1 0 1 2]);
ylabel('Torque (kNm)');
xlim([0 6]);
xticks([0 1 2 3 4 5 6]);
xlabel('Time (s)');

% align ylabels
align_Ylabels(gcf);

% make tikz file
matlab2tikz('unprotected_left_turn_sim_new.tex', ...
            'height', '8cm', 'width', '12cm');