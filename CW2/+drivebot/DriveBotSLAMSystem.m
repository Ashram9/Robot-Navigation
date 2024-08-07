% This class implements an event-based estimation system using g2o and
% the barebones for building up a minimal, ideal SLAM system. The system is
% event-based and responds to a sequence of events which are time stamped
% and served in order. To implement your SLAM system, you will need to
% implement various methods which mostly involve working with the graph.
% These methods are initially stubbed out and will generate exceptions if
% you try to call them.

classdef DriveBotSLAMSystem < minislam.slam.SLAMSystem
    
    properties(Access = public, Constant)
        % Platform state dimension
        NP = 3;
        
        % Landmark dimension
        NL = 2;
        
        % Initial cache size; might help a bit with performance
        INITIAL_CACHE_SIZE = 10000;
    end
    
    properties(Access = protected)
        
        % The most recently created vehicle vertex.
        currentVehicleVertex;
        
        % The set of all vertices associated with the vehicle state over
        % time.
        vehicleVertices;
        vehicleVertexId;
        
        % The set of all prediction edges. These are removed from the graph
        % afterwards if we don't use prediction
        processModelEdges;
        numProcessModelEdges;
        
        % The landmark vertices. Confusingly enough, "Map" here refers to
        % the data structure which is used to store the landmarks. (It
        % allows random access of landmarkID to landmark object.)
        landmarkIDStateVectorMap;
        
        % How often we recommend running the optimization
        recommendOptimizationPeriod;
        
        % Flag to show if we should prune the edges. This is needed for
        % question Q3a
        removePredictionEdgesFromGraph;
        keepFirstPredictionEdge;
        
        % add a flag for graph pruning
        graphPrune = false;
        
    end
    
    methods(Access = public)
        
        % Create the localization system and start it up.
        function this = DriveBotSLAMSystem(configuration)
            
            % Call the base class constructor
            this = this@minislam.slam.SLAMSystem(configuration);
            
            % Preallocate for convenience
            this.vehicleVertices = cell(1, this.INITIAL_CACHE_SIZE);
            
            % No vehicle vertices initally set
            this.vehicleVertexId = 0;
            
            % The set of prediction edges, initially empty
            this.processModelEdges = cell(1, this.INITIAL_CACHE_SIZE);
            this.numProcessModelEdges = 0;
            
            % Allocate the landmark map
            this.landmarkIDStateVectorMap = containers.Map('KeyType', 'int64', 'ValueType', 'any');
            
            % By default, run very infrequently
            this.recommendOptimizationPeriod = inf;
            
            this.removePredictionEdgesFromGraph = false;
            this.keepFirstPredictionEdge = false;
        end
        
        % Destroy the graph when we destroy the SLAM system.
        % Without this, MATLAB will crash whenever this object is destroyed.

        function delete(this)
            vertices = this.graph.vertices();

            for v = 1 : length(vertices)
                this.graph.removeVertex(vertices{v});
            end
        end
        
        % Recommend if an optimization is a good idea. Based on an event,
        % some activities (e.g., such as loop closing) can have a very big
        % impact on the estimates. The logic we have here just recommends
        % an optimization if a fixed number of steps have been completed.
        
        function recommendation = recommendOptimization(this)
            
            % This is how to do it after every 100 steps
            recommendation = rem(this.stepNumber, ...
                this.recommendOptimizationPeriod) == 0;
        end
        
        % Set the value of how often recommend optimization should return
        % true
        function setRecommendOptimizationPeriod(this, newRecommendOptimizationPeriod)
            this.recommendOptimizationPeriod = newRecommendOptimizationPeriod;
        end
        
        % Return the current mean and covariance estimate of the robot.
        % This is only valid after optimization has been called.
        function [x, P] = platformEstimate(this)
            [xS, PS] = this.graph.computeMarginals(this.currentVehicleVertex);
            x=full(xS);
            P=full(PS);
        end
        
        % Returns the entire history of the platform estimates. Suppose
        % there are n vehicle vertices. T is a 1 by N dimensional vector of
        % timesteps. X is a 3 by N dimensional vector of vehicle state (x,
        % y, theta). P is a 3 by N dimensional vector where the nth column
        % are the diagonals from the covariance matrix.
        function [T, X, P] = platformEstimateHistory(this)
            
            % Extract the graph
            [xS, PS] = this.graph.computeMarginals();
            
            % Create the output array
            X = zeros(this.NP, this.vehicleVertexId);
            P = zeros(this.NP, this.vehicleVertexId);
            T = zeros(1, this.vehicleVertexId);
            
            % Copy the outputs over
            for v = 1 : this.vehicleVertexId
                idx = this.vehicleVertices{v}.hessianIndex();
                
                T(v) = this.vehicleVertices{v}.time();
                
                % Copy the estimate into the array. If the vertices is
                % fixed (conditioned), its estimate is okay. The covariance
                % is not explicitly defined, but has a value of zero.
                % Therefore we fill this manually.
                if (isempty(idx) == true)
                    X(:, v) = this.vehicleVertices{v}.estimate();
                    P(:, v) = zeros(3, 1);
                else
                    X(:, v) = full(xS(idx));
                    P(:, v) = full(diag(PS(idx, idx)));
                end
            end
        end
        
        % Return the means and covariances of the landmark estimates. These
        % are only valid after optimization has been called.
        function [x, P, landmarkIds] = landmarkEstimates(this)
            
            landmarkVertices = values(this.landmarkIDStateVectorMap);
            
            numberOfLandmarks = length(landmarkVertices);
            
            landmarkIds = NaN(1, numberOfLandmarks);
            x = NaN(this.NL, numberOfLandmarks);
            P = NaN(this.NL, this.NL, numberOfLandmarks);
            
            [xS, PS] = this.graph.computeMarginals();
            
            for l = 1 : numberOfLandmarks
                landmarkIds(l) = landmarkVertices{l}.landmarkId();
                idx = landmarkVertices{l}.hessianIndex();
                x(:, l) = full(xS(idx));
                if (isempty(idx == true))
                    P(:, :, l) = zeros(3, 3);
                else
                    P(:, :, l) = full(PS(idx, idx));
                end
            end
        end
        
        % We overload the optimize method so that you can add additional
        % logic here
        function chi2 = optimize(this, maximumNumberOfOptimizationSteps)
            
            % Remove the prediction edges if requested.
            if (this.removePredictionEdgesFromGraph == true)
                this.deleteVehiclePredictionEdges();
            end
            
            %Q3b.
            % Pruning graph if requested.
            if (this.graphPrune == true)
                this.graphPruning();
            end
            
            % Now call the actual optimizer. Let it handle the default if
            % no steps are specified.
            if (nargin > 1)
                chi2 = optimize@minislam.slam.SLAMSystem(this, ...
                    maximumNumberOfOptimizationSteps);
            else
                chi2 = optimize@minislam.slam.SLAMSystem(this);
            end
        end
        
        function setRemovePredictionEdges(this, removeEdges, keepFirst)
            this.removePredictionEdgesFromGraph = removeEdges;
            this.keepFirstPredictionEdge = keepFirst;
            
        end

        function setGraphPruning(this, toPruneOrNot)
            this.graphPrune = toPruneOrNot;
        end 
    end
    
    % These are the methods you will need to overload
    methods(Access = protected)
        
        % Handle the initial condition
        
        function handleInitialConditionEvent(this, event)
            
            % Create the first vertex, set its estimate to the initial
            % value and add it to the graph.
            this.currentVehicleVertex = drivebot.graph.VehicleStateVertex(this.currentTime);
            this.currentVehicleVertex.setEstimate(event.data);
            this.graph.addVertex(this.currentVehicleVertex);
            
            % Set the book keeping for this initial vertex.
            this.vehicleVertexId = 1;
            this.vehicleVertices{this.vehicleVertexId} = this.currentVehicleVertex;
            
            % If the covariance is 0, the vertex is known perfectly and so
            % we set it as fixed. If the covariance is non-zero, add a
            % unary initial prior condition edge instead. This adds a soft
            % constraint on where the state can be.
            if (det(event.covariance) < 1e-6)
                this.currentVehicleVertex.setFixed(true);
            else
                initialPriorEdge = drivebot.graph.InitialPriorEdge();
                initialPriorEdge.setMeasurement(event.data);
                initialPriorEdge.setInformation(inv(event.covariance));
                initialPriorEdge.setVertex(this.currentVehicleVertex);
                this.graph.addEdge(initialPriorEdge);
            end
        end
        
        function handleNoPrediction(~)
            % Nothing to do
        end
        
        function handleHeartbeatEvent(this, ~)
            % Nothing to do
        end
        
        function handlePredictToTime(this, time, dT)

            % Create the next vehicle vertex and add it to the graph
            
            this.currentVehicleVertex = drivebot.graph.VehicleStateVertex(time);
            
            % Q1b:
            % Implement prediction code here
            % warning('drivebotslam:handlepredicttotime:unimplemented', ...
            %     'Implement the rest of this method for Q1b.');

            % Create a process edge and add it to the graph
            processModelEdge = drivebot.graph.VehicleKinematicsEdge(dT);
            processModelEdge.setVertex(1, this.vehicleVertices{this.vehicleVertexId});
            processModelEdge.setVertex(2, this.currentVehicleVertex); 
            processModelEdge.setMeasurement(this.u);
            processModelEdge.setInformation(inv(this.uCov));
            processModelEdge.initialize();
            this.graph.addEdge(processModelEdge);

            % Add the prediciton to the current vertex and add it to the
            % graph
            this.graph.addVertex(this.currentVehicleVertex);

            this.numProcessModelEdges = this.numProcessModelEdges + 1;
            this.processModelEdges{this.numProcessModelEdges} = processModelEdge;
            
            % Bump the indices
            this.vehicleVertexId = this.vehicleVertexId + 1;
            this.vehicleVertices{this.vehicleVertexId} = this.currentVehicleVertex;
        end
        
        function handleGPSObservationEvent(this, event)

            % Q1d:
            % Create a GPS measurement edge and add it to the graph
            % warning('drivebotslam:handlegpsobservationevent:unimplemented', ...
            %     'Implement the rest of this method for Q1d.');
            
            % Create a GPS measurement edge and add it to the graph
            gpsMeasurementEdge = drivebot.graph.GPSMeasurementEdge(this.configuration.gpsPositionOffset);
            gpsMeasurementEdge.setVertex(1, this.currentVehicleVertex);
            gpsMeasurementEdge.setMeasurement(event.data);
            gpsMeasurementEdge.setInformation(inv(event.covariance));
            this.graph.addEdge(gpsMeasurementEdge);
        end
        
        function handleCompassObservationEvent(this, event)
            
            % Q1c
            % Create a compass measurement edge and add it to the graph
            compassMeasurementEdge = drivebot.graph.CompassMeasurementEdge(this.configuration.compassAngularOffset);
            compassMeasurementEdge.setVertex(1, this.currentVehicleVertex);
            compassMeasurementEdge.setMeasurement(event.data);
            compassMeasurementEdge.setInformation(inv(event.covariance));
            this.graph.addEdge(compassMeasurementEdge);
        end
        
        function handleLandmarkObservationEvent(this, event)
            
            % Iterate over all the landmark measurements
            for l = 1 : length(event.landmarkIds)
                
                % Get the landmark vertex associated with this measurement.
                % If necessary, a new landmark vertex is created and added
                % to the graph.
                [landmarkVertex, newVertexCreated] = this.createOrGetLandmark(event.landmarkIds(l));
                z = event.data(:, l);

                % Q2b:
                % Complete the implementation
                % warning('drivebotslamsystem:handlelandmarkobservationevent:unimplemented', ...
                %     'Implement the rest of this method for Q2b.');
                % Basically, make a new edge
                landmarkRangeBearingEdge = drivebot.graph.LandmarkRangeBearingEdge();
                landmarkRangeBearingEdge.setVertex(1, this.vehicleVertices{this.vehicleVertexId});
                landmarkRangeBearingEdge.setVertex(2, landmarkVertex);
                landmarkRangeBearingEdge.setMeasurement(z); 
                landmarkRangeBearingEdge.setInformation(inv(event.covariance));
                landmarkRangeBearingEdge.initialize();
                this.graph.addEdge(landmarkRangeBearingEdge);
            end
        end
        
        function deleteVehiclePredictionEdges(this)

            % Q3a:            
            % warning('drivebotslam:deletevehiclepredictionedges:unimplemented', ...
            %     'Implement the rest of this method for Q3a.');

            % extract cell arrays of the all edges from the graph
            edges = this.graph.edges();
            numOfEdges = length(edges);

            % initialise the edge count 
            count = 0;

            % traverse all the edges
            for i = 1:numOfEdges
                % check if the edge is a prediction edge 
                if isa(edges{i}, 'drivebot.graph.VehicleKinematicsEdge')
                    count = count + 1;

                    % check if we need to keep the first prediction edge 
                    if (this.keepFirstPredictionEdge == true && count == 1)                       
                        continue
                    else 
                       % otherwise, remove the current edges 
                        this.graph.removeEdge(edges{i});
                    end
                end
            end

        end
        
        
        % This method returns a landmark associated with landmarkId. If a
        % landmark exists already, it is returned. If it does not exist, a
        % vertex is created and is added to the graph.
        function [landmarkVertex, newVertexCreated] = createOrGetLandmark(this, landmarkId)

            
            % If the landmark exists already, return it
            if (isKey(this.landmarkIDStateVectorMap, landmarkId) == true)
                landmarkVertex = this.landmarkIDStateVectorMap(landmarkId);
                newVertexCreated = false;
                return
            end
            
            fprintf('Creating landmark %d\n', landmarkId);
            
            % Create the new landmark add it to the graph
            landmarkVertex = drivebot.graph.LandmarkStateVertex(landmarkId);
            this.landmarkIDStateVectorMap(landmarkId) = landmarkVertex;
            
            this.graph.addVertex(landmarkVertex);
            
            newVertexCreated = true;
        end
        
        function storeStepResults(this)
            % Nothing
        end

        %Q3b:
        % Implement the graph pruning function
        function graphPruning(this)
            vertices = this.graph.vertices();
            numVertices = length(vertices);
        
            previousVertexLandmarks = []; % Store the set of landmark IDs observed at the previous vertex
            verticesToRemove = []; % Stores the index of the vertex to be removed
        
            % Traverse all the vertices
            for i = 1:numVertices
                vertex = vertices{i};
        
                if isa(vertex, 'drivebot.VehicleStateVertex')
                    vertexEdges = vertex.edges();
                    currentVertexLandmarks = [];
        
                    % Collects all landmark IDs observed at the current vertex
                    for j = 1:length(vertexEdges)
                        edge = vertexEdges{j};
                        if isa(edge, 'drivebot.LandmarkRangeBearingEdge')
                            lm = edge.vertex(2);
                            lmID = lm.landmarkId();
                            currentVertexLandmarks = [currentVertexLandmarks, lmID];
                        end
                    end
        
                    % Check whether the landmarks observed by the current
                    % vertex are same as the previous one
                    if ~isempty(previousVertexLandmarks) && all(sort(previousVertexLandmarks) == sort(currentVertexLandmarks)) && (i ~= 1) && (i ~= 2)
                        % If same, store the previous vertex for further
                        % removal (redundant)
                        verticesToRemove = [verticesToRemove, i-1];
                    end
        
                    % Update the landmark observation set for the previous vertex
                    previousVertexLandmarks = currentVertexLandmarks;
                end
            end
        
            % Remove marked vertices and associated edges
            for i = 1:length(verticesToRemove)
                vertexIndex = verticesToRemove(i);
                vertex = vertices{vertexIndex};
        
                % Remove edges associated with the vertex
                vertexEdges = vertex.edges();
                for j = 1:length(vertexEdges)
                    this.graph.removeEdge(vertexEdges{j});
                end
        
                % Remove the vertex
                this.graph.removeVertex(vertex);
            end
        end


        
    end
end
