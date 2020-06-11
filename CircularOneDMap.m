classdef circularOneDMap < MapGeometry
    %OneDMap is descendant for piecewise linear map
    %   Constructor contains one argument which is number of nodes
    
    methods
        function map = circularOneDMap(N)
            % Create map
            map@MapGeometry(1);
            % Store size
            map.sizes = N;
            % Calculate internal coordinates of nodes
	    x = (1:N)'/N;
            map.internal = [cos(2*pi*x), sin(2*pi*x)];
            % Form array of links
            map.links = [1:N-1;2:N]';
	    map.links = [map.links;[N,1]]
            % Form array of ribs
            map.ribs = [1:N-2;2:N-1;3:N]';
	    map.ribs = [map.ribs;[N-1, N, 1]];
	    map.ribs = [map.ribs;[N, 1, 2]]	
            % Set mapped coordinates to empty set
            map.mapped = [];
        end

        function newMap = extendPrim(map)
	    newMap = map;
        end

                
        function res = getBorder(map)
            res = [];
        end
    end
    
end

