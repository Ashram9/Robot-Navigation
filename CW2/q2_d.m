% This script runs Q2(d)

% Create the configuration object.
configuration = drivebot.SimulatorConfiguration();

% Enable the laser to support pure SLAM
configuration.enableLaser = true;

% For this part of the coursework, this should be set to false.
configuration.perturbWithNoise = false;

% Set this value to truncate the run at a specified timestep rather than
% run through the whole simulation to its end.
configuration.maximumStepNumber = 1200;

% Set up the simulator
simulator = drivebot.DriveBotSimulator(configuration, 'q2_d');

% Create the localization system
drivebotSLAMSystem = drivebot.DriveBotSLAMSystem(configuration);
drivebotSLAMSystem.setRecommendOptimizationPeriod(inf);

% Q2d:
% Explore the  timestep where the loop closure occurs, and get
% results just before and after the loop closure event
% warning('q2_d:unimplemented', ...
%         'Analyse loop closure behaviour for Q2d.')

% This tells the SLAM system to do a very detailed check that the input
% appears to be correct but can make the code run slowly. Once you are
% confident your code is working safely, you can set this to false.
drivebotSLAMSystem.setValidateGraph(false);

% Run the main loop and correct results
results = minislam.mainLoop(simulator, drivebotSLAMSystem);

% Minimal output plots. For your answers, please provide titles and label
% the axes.
title('Simulation');
xlabel('x');
ylabel('y');
saveas(gcf,'Figures/q2_d_simulation_beforeLC.png');

% Plot optimisation times.
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(results{1}.optimizationTimes, '*')
hold on
saveas(gcf,'Figures/q2_d_Optimization_times_beforeLC.png');

% Plot the error curves.
minislam.graphics.FigureManager.getFigure('Errors');
clf
% plot(results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory')
errors = results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory';
% wrap phi in [-pi, pi]
errors(:,3) = g2o.stuff.normalize_thetas(errors(:,3));
plot(errors)
legend('x', 'y', 'phi')
title('Errors')
saveas(gcf,'Figures/q2_d_Errors_beforeLC.png');

% Plot covariance.
minislam.graphics.FigureManager.getFigure('Vehicle Covariances');
clf
plot(results{1}.vehicleCovarianceHistory')
hold on
legend('x', 'y', 'phi')
title('Vehicle Covariances')
saveas(gcf,'Figures/q2_d_Vehicle_Covariances_beforeLC.png');

% Plot chi2 values.
minislam.graphics.FigureManager.getFigure('chi2 values');
clf
plot(results{1}.chi2Time, results{1}.chi2History)
hold on
saveas(gcf,'Figures/q2_d_chi2_beforeLC.png');

%% after loop closing

% Create the configuration object.
configuration = drivebot.SimulatorConfiguration();

% Enable the laser to support pure SLAM
configuration.enableLaser = true;

% For this part of the coursework, this should be set to false.
configuration.perturbWithNoise = false;

% Set this value to truncate the run at a specified timestep rather than
% run through the whole simulation to its end.
configuration.maximumStepNumber = 2000;

% Set up the simulator
simulator = drivebot.DriveBotSimulator(configuration, 'q2_d');

% Create the localization system
drivebotSLAMSystem = drivebot.DriveBotSLAMSystem(configuration);
drivebotSLAMSystem.setRecommendOptimizationPeriod(inf);

% Q2d:
% Explore the  timestep where the loop closure occurs, and get
% results just before and after the loop closure event
% warning('q2_d:unimplemented', ...
%         'Analyse loop closure behaviour for Q2d.')

% This tells the SLAM system to do a very detailed check that the input
% appears to be correct but can make the code run slowly. Once you are
% confident your code is working safely, you can set this to false.
drivebotSLAMSystem.setValidateGraph(false);

% Run the main loop and correct results
results = minislam.mainLoop(simulator, drivebotSLAMSystem);

% Minimal output plots. For your answers, please provide titles and label
% the axes.
title('Simulation');
xlabel('x');
ylabel('y');
saveas(gcf,'Figures/q2_d_simulation_afterLC.png');

% Plot optimisation times.
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(results{1}.optimizationTimes, '*')
hold on
saveas(gcf,'Figures/q2_d_Optimization_times_afterLC.png');

% Plot the error curves.
minislam.graphics.FigureManager.getFigure('Errors');
clf
% plot(results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory')
errors = results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory';
% wrap phi in [-pi, pi]
errors(:,3) = g2o.stuff.normalize_thetas(errors(:,3));
plot(errors)
legend('x', 'y', 'phi')
title('Errors')
saveas(gcf,'Figures/q2_d_Errors_afterLC.png');

% Plot covariance.
minislam.graphics.FigureManager.getFigure('Vehicle Covariances');
clf
plot(results{1}.vehicleCovarianceHistory')
hold on
legend('x', 'y', 'phi')
title('Vehicle Covariances')
saveas(gcf,'Figures/q2_d_Vehicle_Covariances_afterLC.png');

% Plot chi2 values.
minislam.graphics.FigureManager.getFigure('chi2 values');
clf
plot(results{1}.chi2Time, results{1}.chi2History)
hold on
saveas(gcf,'Figures/q2_d_chi2_afterLC.png');

