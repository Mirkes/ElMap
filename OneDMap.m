classdef OneDMap < MapGeometry
    %OneDMap is descendant for piecewise linear map
    %   Constructor contains one argument which is number of nodes
    
    methods
        function map = OneDMap(N)
            %Create map
            map@MapGeometry(1);
            %Calculate internal coordinates of nodes
            map.internal = (1:N)';
            %form array of links
            map.links = [map.internal(1:end-1), map.internal(2:end)];
            %form array of ribs
            map.ribs = [map.internal(1:end-2), map.internal(2:end-1), map.internal(3:end)];
            %Set mappedcoordinates to empty set
            map.mapped = [];
        end
    end
    
end

