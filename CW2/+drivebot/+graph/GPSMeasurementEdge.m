classdef GPSMeasurementEdge < g2o.core.BaseUnaryEdge
   
    properties(Access = protected)
        
        xyOffset;
        
    end
    
    methods(Access = public)
    
        function this = GPSMeasurementEdge(xyOffset)
            this = this@g2o.core.BaseUnaryEdge(2);
            this.xyOffset = xyOffset;
        end
        
        function computeError(this)

	    % Q1d:
        % Implement the code
        % warning('gpsmeasurementedge:computeerror:unimplemented', ...
        %         'Implement the rest of this method for Q1d.');
            x = this.edgeVertices{1}.estimate();
            % Transformation matrix for GPS offset
            c = cos(x(3));
            s = sin(x(3));
            
            M = [c -s;
                s c];

            this.errorZ = (x(1:2) + M * this.xyOffset) - this.z;
        end
        
        function linearizeOplus(this)

	    % Q1d:
        % Implement the code
        % warning('gpsmeasurementedge:lineareizeoplus:unimplemented', ...
        %         'Implement the rest of this method for Q1d.');
        % Basically, one to one correspondance  
            x = this.edgeVertices{1}.estimate();
            % Transformation matrix for GPS offset
            c = cos(x(3));
            s = sin(x(3));
            
            % That crazy thing on orientation just corrects for the GPS
            % offset. It doesn't change that much.
            this.J{1} = ...
                [1 0 -c*this.xyOffset(2)-s*this.xyOffset(1); % How each of the vehicle state affects x error
                0 1 c*this.xyOffset(1)-s*this.xyOffset(2)]; % How each of the vehicle state affects y error


        end
    end
end
