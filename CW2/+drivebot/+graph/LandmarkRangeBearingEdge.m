% This edge encodes a 3D range bearing measurement.
%
% The measurement is in spherical polar coordinates

% Jacobian from https://github.com/petercorke/robotics-toolbox-matlab/blob/master/RangeBearingSensor.m

classdef LandmarkRangeBearingEdge < g2o.core.BaseBinaryEdge
    
    methods(Access = public)
    
        function this = LandmarkRangeBearingEdge()
            this = this@g2o.core.BaseBinaryEdge(2);
        end
        
        function initialize(this)
            % Q2b:
            % Complete implementation
            % warning('landmarkrangebearingedge:initialize:unimplemented', ...
            %     'Implement the rest of this method for Q2b.');
            x = this.edgeVertices{1}.x;
            landmark = zeros(2, 1);
            landmark(1) = x(1) + this.z(1) * cos(this.z(2) + x(3));
            landmark(2) = x(2) + this.z(1) * sin(this.z(2) + x(3));
            this.edgeVertices{2}.setEstimate(landmark);
        end
        
        function computeError(this)

            % Q2b:
            % Complete implementation
            % warning('landmarkrangebearingedge:computeerror:unimplemented', ...
            %     'Implement the rest of this method for Q2b.');
            % Vertex states
            x = this.edgeVertices{1}.estimate();
            landmark = this.edgeVertices{2}.estimate();
            dx = landmark(1:2) - x(1:2);

            this.errorZ(1) = norm(dx) - this.z(1); % Difference between expected distance and measured distance
            this.errorZ(2) = g2o.stuff.normalize_theta(atan2(dx(2), dx(1)) - x(3) - this.z(2)); % Difference between expected relative angle and actual measured angle
        end
        
        function linearizeOplus(this)
            % Q2b:
            % Complete implementation
            % warning('landmarkrangebearingedge:linearizeoplus:unimplemented', ...
            %     'Implement the rest of this method for Q2b.');
            x = this.edgeVertices{1}.estimate();
            landmark = this.edgeVertices{2}.estimate();
            dx = landmark(1:2) - x(1:2);
            r = norm(dx);
            
            this.J{1} = ...
                [-dx(1)/r -dx(2)/r 0; % How does the robot position affect the range error
                dx(2)/r^2 -dx(1)/r^2 -1]; % How does the robot position affect the bearing error
            this.J{2} = - this.J{1}(1:2, 1:2); % The same thing as the vehicle edge
        end        
    end
end