% This script runs Q1(b)

% Create the configuration object.
configuration = drivebot.SimulatorConfiguration();

% Since we are doing just prediction, all the other sensors are disabled.
% This is the default setting.

% If you set this parameter to false, the simulator generates measurements
% with no noise in them. You might find this useful for debugging.
% However, unless specified otherwise, any submitted results must have this
% value set to true.
configuration.perturbWithNoise = true;

% Set up the simulator
simulator = drivebot.DriveBotSimulator(configuration, 'q1_b');

% Create the localization system
drivebotSLAMSystem = drivebot.DriveBotSLAMSystem(configuration);

% This tells the SLAM system to do a very detailed check that the input
% appears to be correct but can make the code run slowly. Once you are
% confident your code is working safely, you can set this to false.
drivebotSLAMSystem.setValidateGraph(false);

% Run the main loop and correct results
results = minislam.mainLoop(simulator, drivebotSLAMSystem);

title('Simulation')
xlabel('x')
ylabel('y')
saveas(gcf,'Figures/q1_b_simulation.png')

% Minimal output plots. These just show you an example of how to plot the
% results. For your report, you need to look at improving these figures
% including labelling axes, etc.

% Plot optimisation times
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(results{1}.vehicleStateTime, results{1}.optimizationTimes, '*')
hold on
saveas(gcf,'Figures/q1_b_Optimization_times.png');

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
saveas(gcf,'Figures/q1_b_Errors.png');

% Plot covariance
minislam.graphics.FigureManager.getFigure('Vehicle Covariances');
clf
plot(results{1}.vehicleStateTime, results{1}.vehicleCovarianceHistory')
hold on
legend('x', 'y', 'phi')
title('Vehicle Covariances')
saveas(gcf,'Figures/q1_b_Vehicle_Covariances.png');