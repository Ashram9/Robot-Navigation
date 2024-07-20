% This script runs Q1(e)

% Create the configuration object.
configuration = drivebot.SimulatorConfiguration();

% Since we are doing prediction and GPS, disable the SLAM sensor
configuration.enableGPS = true;

% Set to true for part ii
configuration.enableCompass = true; %false;

% Set up the simulator
simulator = drivebot.DriveBotSimulator(configuration, 'q1_e');

% Create the localization system
drivebotSLAMSystem = drivebot.DriveBotSLAMSystem(configuration);

% Q1(e)i:
% Use the method "setRecommendOptimizationPeriod" in DriveBotSLAMSystem
% to control the rate at which the optimizer runs
drivebotSLAMSystem.setRecommendOptimizationPeriod(1);

% This tells the SLAM system to do a very detailed check that the input
% appears to be correct but can make the code run slowly. Once you are
% confident your code is working safely, you can set this to false.
drivebotSLAMSystem.setValidateGraph(false);

% Run the main loop and correct results
results = minislam.mainLoop(simulator, drivebotSLAMSystem);

% Minimal output plots. For your answers, please provide titles and label
% the axes.
title('Simulation')
xlabel('x')
ylabel('y')
saveas(gcf,'Figures/q1_e_simulation.png')

% Plot optimisation times
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(results{1}.vehicleStateTime, results{1}.optimizationTimes, '*')
hold on
saveas(gcf,'Figures/q1_e_Optimization_times.png');

% Plot the error curves
minislam.graphics.FigureManager.getFigure('Errors');
clf
% plot(results{1}.vehicleStateTime, results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory')
errors = results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory';
% wrap phi in [-pi, pi]
errors(:,3) = g2o.stuff.normalize_thetas(errors(:,3));
plot(errors);
legend('x', 'y', 'phi')
title('Errors')
saveas(gcf,'Figures/q1_e_Errors.png');

% Plot covariance
minislam.graphics.FigureManager.getFigure('Vehicle Covariances');
clf
plot(results{1}.vehicleStateTime, results{1}.vehicleCovarianceHistory')
hold on
legend('x', 'y', 'phi')
title('Vehicle Covariances')
saveas(gcf,'Figures/q1_e_Vehicle_Covariances.png');

% Plot chi2 values
minislam.graphics.FigureManager.getFigure('chi2 values');
clf
plot(results{1}.chi2Time, results{1}.chi2History)
hold on
saveas(gcf,'Figures/q1_e_chi2.png');
