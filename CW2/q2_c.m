% This script runs Q2(c)

% Create the configuration object.
configuration = drivebot.SimulatorConfiguration();

% Enable the laser to support pure SLAM
configuration.enableLaser = true;

% If you set this parameter to false, the simulator generates measurements
% with no noise in them. You might find this useful for debugging.
% However, unless specified otherwise, any submitted results must have this
% value set to true.
configuration.perturbWithNoise = true;

% Set up the simulator
simulator = drivebot.DriveBotSimulator(configuration, 'q2_c');

% Create the localization system
drivebotSLAMSystem = drivebot.DriveBotSLAMSystem(configuration);
drivebotSLAMSystem.setRecommendOptimizationPeriod(inf);

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
saveas(gcf,'Figures/q2_c_simulation.png');

% Plot optimisation times
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(results{1}.optimizationTimes, '*')
hold on
saveas(gcf,'Figures/q2_c_Optimization_times.png');

% Plot the error curves
minislam.graphics.FigureManager.getFigure('Errors');
clf
% plot(results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory')
errors = results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory';
% wrap phi in [-pi, pi]
errors(:,3) = g2o.stuff.normalize_thetas(errors(:,3));
plot(errors)
legend('x', 'y', 'phi')
title('Errors')
saveas(gcf,'Figures/q2_c_Errors.png');

% Plot covariance
minislam.graphics.FigureManager.getFigure('Vehicle Covariances');
clf
plot(results{1}.vehicleCovarianceHistory')
hold on
legend('x', 'y', 'phi')
title('Vehicle Covariances')
saveas(gcf,'Figures/q2_c_Vehicle_Covariances.png');

% Plot chi2 values
minislam.graphics.FigureManager.getFigure('chi2 values');
clf
plot(results{1}.chi2Time, results{1}.chi2History)
hold on
saveas(gcf,'Figures/q2_c_chi2.png');

% This is how to extract the graph from the optimizer
graph = drivebotSLAMSystem.optimizer();

% This is how to extract cell arrays of the vertices and edges from the
% graph
allVertices = graph.vertices();
allEdges = graph.edges();

% Work out the number of vehicle poses and landmarks. 
numVehicleVertices = 0;
numLandmarks = 0;

landmarkObservationsPerVehicleVertex = 0;
observationsPerLandmarkVertex = 0;

% Q2c:
% Finish implementing the code to capture information about the graph
% structure.
% warning('q2_c:unimplemented', ...
%         'Implement the rest of the graph query code for Q2c.')

% numbers of all vertices and edges
numVertices = length(allVertices);
numEdges = length(allEdges);

numVehicleVertices = length(results{1}.vehicleStateHistory');
fprintf('the number of vehicle poses stored: %d.\n', numVehicleVertices);
numLandmarks = numVertices - numVehicleVertices;
fprintf('the number of landmarks initialized: %d.\n', numLandmarks);

% numbers of observation edges
numObsEdge = numEdges - numVehicleVertices + 1;

landmarkObservationsPerVehicleVertex = numObsEdge / numVehicleVertices;
fprintf('the average number of observations made by a robot at each timestep: %f.\n', landmarkObservationsPerVehicleVertex);
observationsPerLandmarkVertex = numObsEdge / numLandmarks;
fprintf('the average number of observations received by each landmark: %f.\n', observationsPerLandmarkVertex);


