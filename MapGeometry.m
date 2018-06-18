classdef MapGeometry  < handle
    %MapGeometry is class of map geometry for the elastic map and
    %self-organised map.
    %   MapGeometry contains description of underlied map.
    
    properties (SetAccess = protected)
        dimension   % map dimension
        internal    % internal coordinates of nodes
        mapped      % mapped coordinates of nodes
        links       % edges of map. Each edge contains numbers of two nodes 
                    % which are linked by this edge.
        ribs        % list of map's ribs. Each rib contains numbers of three
                    % nodes which form this rib.
        disp        % dispersion measure for PQSQ approach
        preproc     % true if data were preprocessed
        PCs         % empty if preproc is false and set of PCs otherwise.
    end
    
    methods
        function map = MapGeometry( dim )
            map.dimension = dim;
            map.preproc = false;
            map.PCs = [];
        end
        
        function dim = getDimension(map)
            %method to get map dimension.
            dim = map.dimension;
        end
        
        function coord = getInternalCoordinates(map)
            %method to access the internal coordinates of map
            coord = map.internal;
        end
        
        function coord = getMappedCoordinates(map)
            %method to access the mapped coordinates of map
            coord = map.mapped;
        end
        
        function links = getLinks(map)
            %method to access edges of map
            %links is k-by-2 matrix. Each row contains two numbers of nodes 
            %which form one edge.
            links = map.links;
        end
        
        function ribs = getRibs(map)
            %method to access ribs of map
            %ribs is k-by-3 matrix. Each row contains three numbers of 
            %nodes which form this rib.
            ribs = map.ribs;
        end
        
        function disp = getDisp(map)
            %method to access disp of map
            %disp is non-negative number which presents the maximal
            %distance from data point to nearest original node.
            disp = map.disp;
        end
        
        function init(map, data, type, reduce)
        %init is the method of map initialization. This method defines an 
        %initial mapped coordinates. In accordance of results of paper
        %[Akinduko, Ayodeji A., Evgeny M. Mirkes, and Alexander N. Gorban.
        %"SOM: Stochastic initialization versus principal components."
        %Information Sciences(2015),
        %http://dx.doi.org/10.1016/j.ins.2015.10.013] three methods have to
        %be implemented by each map: random initialization, random
        %selection and principal component initialization.
        %
        %Inputs:
        %   map is MapGeometry object to initialise map.
        %   data is n-by-m matrix with n data points and m coordinates for
        %       each point (each row is one data point)
        %   type is the type of initialization: 'random', 'randomSelection'
        %       or 'pci'. Default value 'pci'.
        %   reduce is nonnegative integer. If 'reduce' is positive and is
        %       less than n then specifien number of the first principal
        %       components is used. If 'reduce' is zero and m>n then the
        %       first n-1 principal components is used. If 'reduce' is
        %       positive and is greater or equal to n or 'reduce' is zero
        %       and n>m then dimensionality reduction nis not performed.
        %       
            if nargin < 2
                error('Input argument "data" MUST be specified');
            end
            if nargin < 3
                type = 'pci';
            end
            if nargin < 4
                reduce = 0;
            end
            
            data = map.preprocessDataInit(data, reduce);
            
            if strcmpi('random', type)
                %Random generation
                %Calculate intervals for each coordinates
                mini = min(data);
                maxi = max(data)-mini;
                %Generate random coordinates
                tmp = rand(size(map.internal,1),size(data,2));
                %Scale coordinates
                tmp = bsxfun(@times,tmp,maxi);
                %shift coordinates
                map.mapped = bsxfun(@plus,tmp,mini);
            elseif strcmpi('randomSelection',type)
                %Random selection
                %Generate vector of prototypes numbers
                tmp = randi(size(data,1),size(map.internal,1),1);
                %Get selected points and put as mapped coordinates
                map.mapped = data(tmp,:);
            elseif strcmp('pci',type)
                %Principal component initialization
                %Calculate requared number of PCs:
                [~, D, V] = svds(data, map.dimension);
                D = diag(D);
                [~, ind] = sort(D,'descend');
                V = V(:,ind);
                
                %Normalize PCs' direction
                for k=1:map.dimension
                    if V(k,k)<0
                        V(:,k)=-V(:,k);
                    end
                end
                %Calculate mean and dispersion along each PCs
                tmp = data*V;
                mini = min(tmp);
                maxi = max(tmp);
                means = sum(tmp)/size(data,1);
                disper = min([means-mini;maxi-means]);
                means = sum(data)/size(data,1);
                %Calculate mean and dispersion along internal coordinates
                minI = min(map.internal);
                maxI = max(map.internal);
                meanS = sum(map.internal)/size(map.internal,1);
                disP = min([meanS-minI;maxI-meanS]);
                %auxiliary calculations
                V = bsxfun(@times,V,disper./disP);
                %final values
                map.mapped=bsxfun(@plus,bsxfun(@minus, map.internal, meanS)*V',means);
            else
                error(['type "' type '" is not recognized as valid type of initialization']);
            end
            tmp = associate(map, map.mapped, data);
            map.disp = sqrt(max(tmp));
        end
        
        function data = preprocessDataInit(map, data, reduce)
        %Perform data preprocessing if necessary (see description of
        %'reduce')
        %
        %Inputs:
        %   map is MapGeometry object to initialise map.
        %   data is n-by-m matrix with n data points and m coordinates for
        %       each point (each row is one data point)
        %   reduce is nonnegative integer. If 'reduce' is positive and is
        %       less than n then specifien number of the first principal
        %       components is used. If 'reduce' is zero and m>n then the
        %       first n-1 principal components is used. If 'reduce' is
        %       positive and is greater or equal to n or 'reduce' is zero
        %       and n>m then dimensionality reduction nis not performed.
        %
            % Get sizes
            [n, m] = size(data);
            reduce = floor(reduce);
            if reduce >= m || n > m
                return;
            end
            % Define required number of coordinates
            k = n - 1;
            if reduce > 0 && reduce > k
                k = reduce;
            end
            
            % Search required number of PCs
            [~, D, V] = svds(data, k);
            D = diag(D);
            [~, ind] = sort(D,'descend');
            V = V(:,ind);
            % Standardise direction of PCs
            ind = diag(V) < 0;
            V(:, ind) = -V(:, ind);
            % Store results
            map.preproc = true;
            map.PCs = V;
            % Preprocess data
            data = map.preprocessData(data);
        end
        
        function data = preprocessData(map, data)
            if isempty(map.PCs)
                return;
            end
            data = data * map.PCs;
        end
        
        function coord = project(map, points, type, kind)
        %Project is the method to calculate projection of data point
        %(points) into map. There are d+1 types of projection for d
        %dimensional map: 0 means projection into nearest node of map, 1
        %means projection onto nearest edge of map, 2 means projection onto
        %nearest face of map. Projection can be calculated in the internal
        %or mapped coordinates. There are three input arguments for this
        %method: set of point to project, type of projection (integer
        %number) and coordinates space for projection: �internal� or
        %�mapped�.
        %points is n-by-m matrix where m is number of mapped coordinates
        %   and n is number of points to project.
        %type is type of projection: 0 or 1 for 1D maps, 0,1 or 2 for 2D
        %   maps.
        %kind is one of words 'internal' for internal coordinates and
        %   'mapped'for mapped coordinates.
            coord = projectPrime(map, map.mapped, points, type, kind);
        end
        
        function [coord, dist] = projectPrime(map, nodes, points, type, kind)
        %Project is the method to calculate projection of data point
        %(points) into map. There are d+1 types of projection for d
        %dimensional map: 0 means projection into nearest node of map, 1
        %means projection onto nearest edge of map, 2 means projection onto
        %nearest face of map. Projection can be calculated in the internal
        %or mapped coordinates. There are three input arguments for this
        %method: set of point to project, type of projection (integer
        %number) and coordinates space for projection: �internal� or
        %�mapped�.
        %coord is the set of requested projections.
        %dist is the vector of distances from points to map.
        %nodes is the current state of the mapped nodes. It can be diffed
        %   from map.mapped. It is useful for estimation of calculated
        %   mapped nodes without fixing it into map object
        %points is n-by-m matrix where m is number of mapped coordinates
        %   and n is number of points to project.
        %type is type of projection: 0 or 1 for 1D maps, 0,1 or 2 for 2D
        %   maps.
        %kind is one of words 'internal' for internal coordinates and
        %   'mapped'for mapped coordinates.
        
            %Check which type of coordinates is necessary to return
            cType = strcmpi('mapped',kind);
            N=size(points,1);
            
            switch type
                case 0
                    %Projection to the nearest node
                    %Search the nearest node
                    [dist, num] = map.associate(nodes, points);
                    if cType
                        coord = nodes(num,:);
                    else
                        coord = map.internal(num,:);
                    end
                case 1
                    %projection onto nearest edge
                    %Get array of edges end
                    V2 = (nodes(map.links(:,2),:))';
                    %Form matrix of edges directions
                    dir = (nodes(map.links(:,1),:))'-V2;
                    %Calculate squared length of edge directions
                    len = sum(dir.^2);
                    %Calculate projections length (l in documentation, matrix
                    %analogue of (2))
                    pr = bsxfun(@minus,points*dir,sum(V2.*dir));
                    %Copy projections to normalize (l* in documentation)
                    prn = bsxfun(@rdivide,pr,len);
                    %Non negativity
                    prn(prn<0) = 0;
                    %Cut too long projections (it is the same as step 3 of
                    %algorithm in documentation)
                    prn(prn>1) = 1;
                    %Calculate distances:
                    dist = bsxfun(@plus,sum(points.^2,2),sum(V2.^2))...
                        -2*points*V2+prn.*(bsxfun(@times,prn,len)-2*pr);
                    %Select the nearest edge
                    [dist, edge] = min(dist,[],2);
                    %form index to find length of projections
                    ind = sub2ind([N,size(map.links,1)],(1:N)',edge);
                    if cType
                        coord = bsxfun(@times,(1-prn(ind)),...
                            nodes(map.links(edge,2),:))...
                            +bsxfun(@times,prn(ind),...
                            nodes(map.links(edge,1),:));
                    else
                        coord = bsxfun(@times,(1-prn(ind)),...
                            map.internal(map.links(edge,2),:))...
                            +bsxfun(@times,prn(ind),...
                            map.internal(map.links(edge,1),:));
                    end
                case 2
                    %Projections onto face
                    if ~any(strcmp(methods(map), 'getFaces'))
                        error('request of the projection onto face for MapGeometry without methods "getFaces"');
                    end
                    %Get faces
                    face = map.getFaces;
                    %form auxiliary vectors
                    Y2 = (nodes(face(:,3),:))';
                    Y20 = Y2-(nodes(face(:,1),:))';
                    Y21 = Y2-(nodes(face(:,2),:))';
                    Y10 = (nodes(face(:,2),:)-nodes(face(:,1),:))';
                    Y20Y20 = sum(Y20.^2);
                    Y21Y20 = sum(Y20.*Y21);
                    Y21Y21 = sum(Y21.^2);
                    %Calculate projections
                    A20 = bsxfun(@minus,sum(Y2.*Y20),points*Y20);
                    A21 = bsxfun(@minus,sum(Y2.*Y21),points*Y21);
                    A10 = bsxfun(@rdivide,bsxfun(@minus,A20-A21,sum(Y21.*Y10)),sum(Y10.^2));
                    tmp = Y20Y20.*Y21Y21-Y21Y20.^2;
                    A0 = bsxfun(@rdivide,bsxfun(@times,A20,Y21Y21)...
                        -bsxfun(@times,A21,Y21Y20),tmp);
                    A1 = bsxfun(@rdivide,bsxfun(@times,A21,Y20Y20)...
                        -bsxfun(@times,A20,Y21Y20),tmp);
                    A20 = bsxfun(@rdivide,A20,Y20Y20);
                    A21 = bsxfun(@rdivide,A21,Y21Y21);
                    %Normalize projections
                    A20N = A20;
                    A20N(A20N<0) = 0;
                    A20N(A20N>1) = 1;
                    A21N = A21;
                    A21N(A21N<0) = 0;
                    A21N(A21N>1) = 1;
                    A10(A10<0) = 0;
                    A10(A10>1) = 1;
                    tmp = A0<0;
                    A0(tmp) = 0;
                    A1(tmp) = A21N(tmp);
                    tmp = A1<0;
                    A0(tmp) = A20N(tmp);
                    A1(tmp) = 0;
                    tmp = (1-(A0+A1))<0;
                    A0(tmp) = A10(tmp);
                    A1(tmp) = 1-A10(tmp);
                    %Calculate distances
                    dist = bsxfun(@plus,sum(points.^2,2),sum(Y2.^2))-2*points*Y2...
                        +bsxfun(@times,A0,Y20Y20).*(A0-2*A20)...
                        +bsxfun(@times,A1,Y21Y21).*(A1-2*A21)...
                        +2*bsxfun(@times,A0.*A1,Y21Y20);
                    %Select the nearest face
                    [dist, tmp] = min(dist,[],2);
                    %form index to find length of projections
                    ind = sub2ind([N,size(face,1)],(1:N)',tmp);
                    if cType
                        coord = bsxfun(@times,A0(ind),...
                            nodes(face(tmp,1),:))...
                            +bsxfun(@times,A1(ind),...
                            nodes(face(tmp,2),:))...
                            +bsxfun(@times,1-A0(ind)-A1(ind),...
                            nodes(face(tmp,3),:));
                    else
                        coord = bsxfun(@times,A0(ind),...
                            map.internal(face(tmp,1),:))...
                            +bsxfun(@times,A1(ind),...
                            map.internal(face(tmp,2),:))...
                            +bsxfun(@times,1-A0(ind)-A1(ind),...
                            map.internal(face(tmp,3),:));
                    end
                otherwise
                    error('unacceptable type or projections');
            end
        end
        
        function [dist, klas] = associate(~, node, data)
        %associate identify the nearest node for each data point and
        %return the squared distance between selected node and data
        %point and number of nearest node.
        %Inputs:
        %   node is n-by-k matrix of mapped coordinates for tested state
        %       of map, where n is number of nodes and m is dimension of
        %       data space.
        %   data is m-by-k data points to test, where m is number of
        %       points and k is dimension of data space.
        %Outputs:
        %   dist is m-by-1 matrix of squared distances from data point to
        %       nearest node
        %   klass is m-by-1 vector which contains number of nearest node
        %       for each data point.
            
            dist=bsxfun(@plus,sum(data.^2,2),sum(node.^2,2)')-2*(data*node');
            [dist, klas] = min(dist,[],2);
        end
        
        function fvu = FVU(map, data, node, type)
        %Calculate fraction of variance unexplained for specified data and
        %nodes.
        %data is set of data points
        %node is the set of considered mapped nodes. If t is ommited or
        %   empty then the map.mapped is used.
        %type is the type of projection: 0 means projection into nearest
        %   node of map, 1 means projection onto nearest edge of map, 2
        %   means projection onto nearest face of map. If this argument is
        %   omitted then 1 is used.
        
            %Check the input attributes
            if nargin < 4
                type = 1;
            end
            if nargin<3 || isempty(node)
                node = map.mapped;
            end
            
            %Calculate base variance
            N = size(data,1);
            means = sum(data)/N;
            base = sum(data(:).^2)-N*sum(means.^2);
            
            %Get distances to map
            [~, dist] = map.projectPrime(node, data, type, 'mapped');
            
            %Calculate FVU
            fvu = sum(dist)/base;
        end
        
        function putMapped(map, newMapped)
            %This method is used for the putting the fitted mapped
            %coordinates of map.
            %   newMapped is new matrix of mapped coordinates. It must have
            %       the same size as previously defined matrix
            if ~all(size(map.mapped)==size(newMapped))
                error('Matrix newMapped must have the same size as matrix mapped');
            end
            map.mapped = newMapped;
        end
    end
        
    
%     methods (Abstract)
%         %Distance is method to calculate distances between two nodes in the
%         %internal coordinates. This method is useful for SOM fitting
%         %procedure.
%         dist = distance(map,nodes)
%     end
    
end

