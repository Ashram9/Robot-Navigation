% This class uses a slightly simpler model for the vehicle kinematics used
% in the lectures. This is the more standard built in type for estimate.
%
% The model assumes that the vehicle speed is specified in the vehicle
% frame and is then projected into the world frame. Specifically,
%
% M = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0;0 0 1];
%
% The process model has the form:
%
% x = x + M * [vx;vy;theta]
%
% where vx, vy and vtheta are the velocities.
%
% The error model 
% eTheta = 

classdef VehicleKinematicsEdge < g2o.core.BaseBinaryEdge
    
    properties(Access = protected)
        % The length of the time step
        dT;
    end
    
    methods(Access = public)
        function this = VehicleKinematicsEdge(dT)
            assert(dT >= 0);
            this = this@g2o.core.BaseBinaryEdge(3);            
            this.dT = dT;
        end
       
        function initialize(this)
            
                        
            priorX = this.edgeVertices{1}.x;

            c = cos(priorX(3));
            s = sin(priorX(3));
            
            M = this.dT * [c -s 0;
                s c 0;
                0 0 1];
            
            % Compute the posterior assming no noise
            this.edgeVertices{2}.x = this.edgeVertices{1}.x + M * this.z;

            % Wrap the heading to -pi to pi
            this.edgeVertices{2}.x(3) = g2o.stuff.normalize_theta(this.edgeVertices{2}.x(3));

        end
        
        function computeError(this)
    
            % Q1b:
            % Complete implementation
            % warning('vehiclekinematicsedge:computeerror:unimplemented', ...
            %         'Implement the rest of this method for Q1b.');
            
            % Rotation matrix from prior state
            priorX = this.edgeVertices{1}.x;

            c = cos(priorX(3));
            s = sin(priorX(3));
            
            Mi = [c s 0;
                -s c 0;
                0 0 1];

            % Compute the error.
            dx = this.edgeVertices{2}.x - priorX;
            dx(3) = g2o.stuff.normalize_theta(dx(3));
            this.errorZ = Mi * (dx) / this.dT - this.z;
            
            % Wrap the heading error to -pi to pi
            this.errorZ(3) = g2o.stuff.normalize_theta(this.errorZ(3));
        end
        
        % Compute the Jacobians
        function linearizeOplus(this)

            % Q1b:
            % Complete implementation
            % warning('vehiclekinematicsedge:linearizeoplus:unimplemented', ...
            %     'Implement the rest of this method for Q1b.');
            priorX = this.edgeVertices{1}.x; % Vertex 1
            c = cos(priorX(3));
            s = sin(priorX(3));
            dx = this.edgeVertices{2}.x - priorX; %Extract change
            Mi = [c s 0;
                -s c 0;
                0 0 1]; % The same inverse transformation as above
            this.J{2} = Mi / this.dT;
            this.J{1}(1, 1) = - c; % Partial derivative of the x direction, wrt the x direction
            this.J{1}(1, 2) = - s; % Partial derivative of the x direction wrt the y direction
            this.J{1}(1, 3) = -dx(1) * s + dx(2) * c; % How orientation changes x
            this.J{1}(2, 1) = s; % How x changes y
            this.J{1}(2, 2) = - c; % How y changes y
            this.J{1}(2, 3) = -dx(1) * c - dx(2) * s; % How orientation changes y
            this.J{1}(3, 3) = -1; % How orientation changes orientation
            this.J{1} = this.J{1} / this.dT; % Normalise by time step

        end
    end    
end