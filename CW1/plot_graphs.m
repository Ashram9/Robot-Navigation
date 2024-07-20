% plot trajectory of lawnmower
addpath('Coursework 1 Data/');

% load .csv files of 3 methods
data_method1 = readmatrix('results_GNSS_KF_with_OD.csv');
latitude_GNSS = data_method1(:, 2);
longitude_GNSS = data_method1(:, 3);

data_method2 = readmatrix('results_DR.csv');
latitude_DR = data_method2(:, 2);
longitude_DR = data_method2(:, 3);

data_method3 = readmatrix('motion_profile.csv');
latitude_integrated = data_method3(:, 2);
longitude_integrated = data_method3(:, 3);
north_velocity_data = data_method3(:, 4);
east_velocity_data = data_method3(:, 5);
heading_deg = data_method3(:, 6);


% Plot Trajectory Comparison
figure;
plot(longitude_GNSS, latitude_GNSS, '-.b', 'LineWidth', 1.5, 'DisplayName', 'GNSS');
hold on;
plot(longitude_DR, latitude_DR, '--g', 'LineWidth', 1.5, 'DisplayName', 'DR');
plot(longitude_integrated, latitude_integrated, '-r', 'LineWidth', 2, 'DisplayName', 'DR/GNSS Integration');
legend;
title('Lawnmower Trajectory of 3 Different Methods');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
grid on;
saveas(gcf, 'trajectory_comparison.png'); 

% Plot Velocity Vector (NE)
figure;
% Use quiver function
quiver(longitude_integrated, latitude_integrated, east_velocity_data, north_velocity_data, 'b');
title('Velocity Vector Plot (m/s)');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
grid on;
saveas(gcf, 'velocity_arrow_plot.png');

% Plot Heading
figure;
heading_rad = deg2rad(heading_deg);
north_component = cos(heading_rad);
east_component = sin(heading_rad);
% Use the quiver function to plot arrows for heading angles,
% which determined by the corresponding east-bound and north-bound components
quiver(longitude_integrated, latitude_integrated, east_component, north_component, 'b');
title('Heading Quiver Plot');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
saveas(gcf, 'heading_quiver_plot.png');

