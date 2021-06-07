classdef hemisphereMap < MapGeometry
    %hemisphereMap is descendant for piecewise linear map
    %   Constructor contains two argument which are number of parallel and
    %   number of meridians
    
    properties (SetAccess = protected)
        faces   %Array of faces to projections
    end
    
    methods
        function map = hemisphereMap(parallel, meridian)
	    % parallel will be circular component, while meridian will be
	    % linear (radii of circle) 
            if nargin < 2
                error('You MUST specify number of parallels and meridian for hemisphereMap');
            end
            % Create map
            map@MapGeometry(2);
            % Store size
            map.sizes = [parallel, meridian];
            % Calculate internal coordinates of nodes
            % Initial point is central point
            tmp = [0, 0];
            % Generate one circle
            ang = (1:meridian)' * 2 * pi / meridian;
            fragm = [cos(ang), sin(ang)];
            for m = 1:parallel - 1
                tmp = [tmp; m * fragm];
            end
            map.internal = tmp;
            % Form array of links
            % Radii
            tmp = [1, ((0:parallel - 2) * meridian) + 2];
            tmp = [tmp(1:end - 1)', tmp(2:end)'];
            temp = tmp;
            for m = 2:meridian
                tmp = tmp + 1;
                tmp(1, 1) = 1;
                temp = [temp; tmp];
            end
            
            patt = (1:parallel - 1) * meridian;
            N = size(map.internal, 1);
            tmp = [(2:N)', (3:N+1)'];
            tmp(patt, 2) = tmp(patt, 2) - meridian;
            map.links = [temp; tmp];
            % Form array of ribs
            % Radii
            temp = [];
            if parallel > 2
                tmp = [1, ((0:parallel - 2) * meridian) + 2]';
                tmp = [tmp(1:end - 2), tmp(2:end - 1), tmp(3:end)];
                temp = tmp;
                for m = 2:meridian
                    tmp = tmp + 1;
                    tmp(1, 1) = 1;
                    temp = [temp; tmp];
                end
            end
            % Circles
            tmp = (1:meridian)' + 1;
            tmp = [tmp; tmp(1:2)];
            tmp = [tmp(1:end - 2), tmp(2:end - 1), tmp(3:end)];
            for m = 2:parallel
                temp = [temp; tmp];
                tmp = tmp + meridian;
            end
            if mod(meridian, 2) == 0
                % In this case we have hemisphere
                tmp = [2, 1, meridian / 2 + 2];
                for m = 1:meridian / 2
                    temp = [temp; tmp];
                    tmp = tmp + [1, 0, 1];
                end
            end
            map.ribs = temp;
            % Form array of faces
            % Internal circle
            tmp = [(1:meridian) + 1, 2]';
            temp = [ones(meridian, 1), tmp(1:end -1), tmp(2:end)];
            % External circles
            for m = 2:parallel - 1
                tmp1 = tmp + meridian;
                temp = [temp; [tmp(1:end -1), tmp(2:end), tmp1(1:end -1)]];
                temp = [temp; [tmp1(1:end -1), tmp1(2:end), tmp(2:end)]];
                tmp = tmp1;
            end
    	    map.faces = temp;
            % Set mapped coordinates to empty set
            map.mapped = [];
        end
        
        function face = getFaces(map)
            %Function to access to the faces of map.
            %face is k-by-3 matrix. Each row contains three numbers of
            %nodes which form one face.
            face = map.faces;
        end
        
        function newMap = extendPrim(map)
            % Extension simply add one more parallel
            % Get size of existing map
            par = map.sizes(1);
            mer = map.sizes(2);
            % Create new map with greater size
            newMap = hemisphereMap(par + 1, mer);
            newMap.disp = map.disp;       % dispersion measure for PQSQ approach
            newMap.preproc = map.preproc; % true if data were preprocessed
            newMap.means = map.means;     % mean of data otherwise.
            newMap.PCs = map.PCs;         % set of PCs otherwise.
            % Form list of old nodes
            ind = 1:size(map.internal, 1);
            % Copy known positions of nodes into new map
            newMap.mapped(ind, :) = map.mapped;
            % Get indices of the last border
            bord = map.getBorder();
            newBord = newMap.getBorder();
            % Search previous layer to define new values of nodes
            if par < 3
                preBord = ones(length(bord), 1);
            else
                preBord = bord - mer;
            end
            % Calculate new positions of the added nodes 
            newMap.mapped(newBord, :) = 2 * map.mapped(bord, :) - map.mapped(preBord, :);
        end
        
        function res = getBorder(map)
            % Get map sizes
            n = size(map.internal, 1);
            res = n - map.sizes(2) + 1:n;
        end
    end
end

