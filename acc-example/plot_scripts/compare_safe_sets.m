% set default figure renderer
set(0, 'DefaultFigureRenderer', 'opengl');

% specify figure size
figure('Renderer', 'painters', 'Position', [10 10 1500 1200]);

% load baseline safe set
clear;
load('veh_follow_abs_3D_5_9_2021.mat');
plot_3D_safe_set;
colormap winter;
freezeColors;

% load expanded safe set
clear;
load('veh_follow_abs_3D_expanded_5_9_2021.mat', 'veh_follow_abs');
% adjustment used to prevent 'z-stitching' as discussed here:
% https://www.mathworks.com/matlabcentral/answers/6623-plotting-overlapping-surfaces
adjustment = 0.075;
plot_3D_safe_set;
colormap default;
view(-35,35);
