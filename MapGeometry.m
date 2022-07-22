classdef MapGeometry  < handle
    %MapGeometry is class of map geometry for the elastic map and
    %self-organised map.
    %   MapGeometry contains description of underlying map.
    
    properties (SetAccess = protected)
        sizes       % size of map is map dependent
        dimension   % map dimension
        internal    % internal coordinates of nodes
        mapped      % mapped coordinates of nodes
        links       % edges of map. Each edge contains numbers of two nodes 
                    % which are linked by this edge.
        ribs        % list of map's ribs. Each rib contains numbers of three
                    % nodes which form this rib.
        disp        % dispersion measure for PQSQ approach
        preproc     % true if data were preprocessed
        means       % mean of data otherwise.
        PCs         % set of PCs otherwise. Number of PCs can be specified 
                    % in the call of Init function. 
    end
    
    methods
        function map = MapGeometry( dim )
            map.dimension = dim;
            map.preproc = false;
            map.PCs = [];
        end
        
        function dim = getDimension(map)
            % Function to get map dimension.
            dim = map.dimension;
        end
        
        function coord = getInternalCoordinates(map)
            % Function to access the internal coordinates of map
            coord = map.internal;
        end
        
        function coord = getMappedCoordinates(map)
            % Function to access the mapped coordinates of map
            coord = map.mapped;
        end
        
        function links = getLinks(map)
            % Function to access edges of map. 
            % links is k-by-2 matrix. Each row contains two numbers of
            % nodes which form one edge.
            links = map.links;
        end
        
        function ribs = getRibs(map)
            % Function to access ribs of map
            % ribs is k-by-3 matrix. Each row contains three numbers of 
            % nodes which form this rib.
            ribs = map.ribs;
        end
        
        function disp = getDisp(map)
            % Function to access disp of map
            % disp is non-negative number which presents the maximal
            % distance from data point to nearest original node.
            disp = map.disp;
        end
        
        function init(map, data, type, reduce)
        % init is the  function of map initialization. This  function
        % defines an initial mapped coordinates. In accordance of results
        % of paper [Akinduko, Ayodeji A., Evgeny M. Mirkes, and Alexander
        % N. Gorban. "SOM: Stochastic initialization versus principal
        % components." Information Sciences(2015), 
        % http://dx.doi.org/10.1016/j.ins.2015.10.013] three methods have
        % to be implemented by each map: random initialization, random
        % selection and principal component initialization.
        %
        %Inputs:
        %   map is MapGeometry object to initialise map.
        %   data is n-by-m matrix with n data points and m coordinates for
        %       each point (each row is one data point)
        %   type is the type of initialization:
        %           'random' is random generation
        %           'randomSelection' is random selection of data points as
        %               initial nodes location
        %           'pci' is initialisation along the first, first two or
        %               first three PCs.
        %           vector with 1, 2 or 3 nonzero elements. These positive
        %               elements are numbers of used coordinates. For
        %               example, vector [2, 1] means usage of the second
        %               coordinate as first dimension and the first
        %               coordinate as the second dimension. Negative
        %               elements are numbers of used PCs. For example,
        %               vector [-1, -2] means usage of the first and the
        %               second PCs and coincides with default 'pci'
        %               initialisation. The number of required elements 
        %               depends of map dimension: 1 for 1D map, 2 for 2D
        %               map. For circular 1D map it is necessary to specify
        %               2 vectors. It is possible to combine positive and
        %               negative components.
        %           matrix with 1, 2 or 3 columns and with the same number
        %               of rows as number of columns in data. Each column
        %               vector of the matrix is cinsidered as direction of
        %               corresponding axis for map initialisation. 
        %        Default value 'pci'.
        %   reduce is nonzero integer. If 'reduce' is positive and is
        %       less than n then specified number of the first principal
        %       components are used. If 'reduce' is zero and m>n then the
        %       first n-1 principal components is used. If 'reduce' is
        %       positive and is greater or equal to n or 'reduce' is zero
        %       and n>m then dimensionality reduction is not performed. If
        %       reduce is negative then -reduce PCs are calculated but
        %       dimensionality reduction is not performed.
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
            
            reserve = data;
            data = map.preprocessDataInit(data, reduce);
            
            if strcmpi('random', type)
                %Random generation
                %Calculate intervals for each coordinates
                mini = min(data);
                maxi = max(data)-mini;
                %Generate random coordinates
                data = rand(size(map.internal,1),size(data,2));
                %Scale coordinates
                data = bsxfun(@times,data,maxi);
                %shift coordinates
                map.mapped = bsxfun(@plus,data,mini);
            elseif strcmpi('randomSelection',type)
                %Random selection
                %Generate vector of prototypes numbers
                data_ind = randi(size(data,1),size(map.internal,1),1);
                %Get selected points and put as mapped coordinates
                map.mapped = data(data_ind,:);
            elseif strcmp('pci',type) || isvector(type) || ismatrix(type)
                % Customised initialisation 
                embedding_dim = size(map.internal,2);
                % Form processed data matrix and vector of shift
                if strcmp('pci',type)
                    %Principal component initialization
                    % Get mean and PCs
                    if map.preproc
                        % Data were preprocessed
                        V = eye(size(data, 2), embedding_dim);
                        tmp = data(:, 1:embedding_dim);
                        meanDat = zeros(1, size(data, 2));
                    else
                        % Get requared number of PCs:
                        meanDat = map.means;
                        V = map.PCs(:, 1:embedding_dim);
                        tmp = data * V;
                    end
                elseif isvector(type)
                    if sum(type == 0) > 0
                        error('All elements of vector type must be nonzero');
                    end
                    % Vector of coordinates or PCs
                    if length(type) < embedding_dim
                        error(['Argumrnt type is vector with ',...
                            num2str(length(type)), ' elements, but for ',...
                            'this map ',num2str(embedding_dim),...
                            ' is necessary']);
                    end
                    % Are coordinates necessary?
                    nCo = max(type);
                    if nCo > 0
                        if map.preproc
                            error(['It is impossible to use coordinate ',...
                                'vectors for preprocessed data.',...
                                'Data preprocessing is producing if ',...
                                'positive "reduce" is specified or ',...
                                'if "reduce" is zero and the number of ',...
                                'attributes is greater than number of ',...
                                'observations.']);
                        end
                        if nCo > size(reserve, 2)
                            error(['The maximal number of requested ',...
                                'coordinate ', num2str(nCo),... 
                                ' is greater then the number',...
                                ' of coordinatess ',...
                                num2str(size(reserve, 2)),...
                                '. Initialisation is not possible.']);
                        end
                        % Get auxiliary arrays
                        DV = eye(size(reserve, 2), nCo);
                        meanDat = mean(reserve);
                        Dtmp = reserve;
                    end
                    % Check the PCs necessity
                    nPCs = -min(type);
                    if nPCs > 0
                        if nPCs > size(map.PCs, 2)
                            error(['The maximal number of requested PC ',...
                                num2str(nPCs),... 
                                ' is greater then the number',...
                                ' of calculated PCs ',...
                                num2str(size(map.PCs, 2)),...
                                '. Initialisation is not possible.']);
                        end
                        % Get auxiliary arrays
                        if map.preproc
                            % Data were preprocessed
                            PV = eye(size(data, 2), nPCs);
                            Ptmp = data(:, 1:embedding_dim);
                            meanDat = zeros(1, size(data, 2));
                        else
                            % Get requared number of PCs:
                            meanDat = map.means;
                            PV = map.PCs(:, 1:nPCs);
                            Ptmp = data * PV;
                        end
                    end
                    % Create and complete required matrices
                    tmp = zeros(size(data, 1), embedding_dim);
                    V = zeros(size(data, 2), embedding_dim);
                    for k = 1:embedding_dim
                        if type(k) > 0
                            % Coordinate
                            tmp(:, k) = reserve(:, type(k));
                            V(:, k) = DV(:, type(k));
                        else
                            % PC
                            tmp(:, k) = Ptmp(:, -type(k));
                            V(:, k) = PV(:, -type(k));
                        end
                    end
                else
                    % Matrix. Each column vector is axis vector. 
                    if size(type, 2) < embedding_dim
                        error(['Current map required ',...
                            num2str(embedding_dim),...
                            ' vectors but only ',...
                            num2str(size(type, 2)),...
                            ' were presented. Initialisation ',...
                            'is not possible.']);
                    end
                    if size(type, 1) ~= size(reserve, 2)
                        error(['For ', num2str(size(reserve, 2)),...
                            '-dimensional data matrix type must have ',...
                            num2str(size(reserve, 2)),...
                            'elements instead of presented ',...
                            num2str(size(type, 1)),...
                            '. Initialisation is not possible.']);
                    end
                    % Calculate othonormal basis for presented vectors
                    V = type(:, 1:embedding_dim);
                    for k = 1:embedding_dim
                        % Normalise current vector
                        V(:, k) = V(:, k) / sqrt(sum(V(:, k) .^ 2));
                        % Subtract from the further vectors
                        for kk = k + 1:embedding_dim
                            V(:, kk) = V(:,kk) - (V(:, k)' * V(:,kk)) .* V(:, k);
                        end
                    end
                    % Calculate mean and projections
                    meanDat = mean(reserve);
                    tmp = reserve * V;
                end
                %Calculate mean and dispersion along each PCs
                mini = min(tmp);
                maxi = max(tmp);
                meant = mean(tmp);
                disper = min([meant - mini; maxi - meant]);
                %Calculate mean and dispersion along internal coordinates
                minI = min(map.internal);
                maxI = max(map.internal);
                meanI = sum(map.internal) / size(map.internal,1);
                disP = min([meanI - minI; maxI - meanI]);
                %auxiliary calculations
                V = bsxfun(@times, V, disper ./ disP);
                %final values
                map.mapped=bsxfun(@plus,...
                    bsxfun(@minus, map.internal, meanI) * V', meanDat);
            else
                error(['type "' type '" is not recognized as valid type',...
                    ' of initialization']);
            end
            data = associate(map, map.mapped, data);
            map.disp = sqrt(max(data));
        end
        
        function data = preprocessDataInit(map, data, reduce)
        %Perform data preprocessing if necessary (see description of
        %'reduce')
        %
        %Inputs:
        %   map is MapGeometry object to initialise map.
        %   data is n-by-m matrix with n data points and m coordinates for
        %       each point (each row is one data point)
        %   reduce is integer. If 'reduce' is positive and is
        %       less than n then specified number of the first principal
        %       components are used. If 'reduce' is zero and m>n then the
        %       first n-1 principal components is used. If 'reduce' is
        %       positive and is greater or equal to n or 'reduce' is zero
        %       and n>m then dimensionality reduction is not performed. If
        %       reduce is negative then -reduce PCs are calculated but
        %       dimensionality reduction is not performed.

            % Get sizes
            [n, m] = size(data);
            reduce = floor(reduce);
            if reduce >= m || (reduce == 0 && n > m)
                % Calculate 3 PCs and mean but do not apply preprocessing
                k = 3;
                if k > m
                    k = m;
                end
                map.preproc = false;
            elseif reduce < 0
                k = -reduce;
                if k > m
                    k = m;
                end
                map.preproc = false;
            else
                % Define required number of PCs
                k = n - 1;
                if reduce > 0 && reduce > k
                    k = reduce;
                end
                map.preproc = true;
            end
            
            % Search required number of PCs
            map.means = mean(data);
            [~, D, V] = svds(bsxfun(@minus, data, map.means), k);
            D = diag(D);
            [~, ind] = sort(D,'descend');
            V = V(:,ind);
            % Standardise direction of PCs
            ind = diag(V) < 0;
            V(:, ind) = -V(:, ind);
            % Store results
            map.PCs = V;
            
            % Preprocess data if it is required
            if map.preproc
                data = map.preprocessData(data);
            end
        end
        
        function data = preprocessData(map, data)
            if ~map.preproc
                return;
            end
            data = bsxfun(@minus, data, map.means) * map.PCs;
        end
        
        function data = deprocessData(map, data)
            if ~map.preproc
                return;
            end
            data = bsxfun(@plus, data  * map.PCs', map.means);
        end
        
        function coord = project(map, points, type, kind)
        %Project is the function to calculate projection of data point
        %(points) into map. There are d+1 types of projection for d
        %dimensional map: 0 means projection into nearest node of map, 1
        %means projection onto nearest edge of map, 2 means projection onto
        %nearest face of map. Projection can be calculated in the internal
        %or mapped coordinates. There are three input arguments for this
        %method: set of point to project, type of projection (integer
        %number) and coordinates space for projection: "internal" or
        %"mapped".
        %
        %Inputs:
        %   map is MapGeometry object to use
        %   points is n-by-m matrix where m is number of mapped coordinates
        %       and n is number of points to project.
        %   type is type of projection: 0 or 1 for 1D maps, 0,1 or 2 for 2D
        %       maps.
        %   kind is one of words 'internal' for internal coordinates and
        %       'mapped'for mapped coordinates.
            coord = projectPrime(map, map.mapped, points, type, kind);
        end
        
        function [coord, dist] = projectPrime(map, nodes, points, type, kind)
        %projectPrime is the function to calculate projection of data point
        %(points) into map. There are d+1 types of projection for d
        %dimensional map: 0 means projection into nearest node of map, 1
        %means projection onto nearest edge of map, 2 means projection onto
        %nearest face of map. Projection can be calculated in the internal
        %or mapped coordinates. There are three input arguments for this
        %method: set of point to project, type of projection (integer
        %number) and coordinates space for projection: "internal or
        %mapped.
        %
        %Inputs:
        %   map is MapGeometry object to use
        %   nodes is the current state of the mapped nodes. It can be
        %       diffed from map.mapped. It is useful for estimation of
        %       calculated mapped nodes without fixing it into map object
        %   points is n-by-m matrix where m is number of mapped coordinates
        %       and n is number of points to project.
        %   type is type of projection: 0 or 1 for 1D maps, 0,1 or 2 for 2D
        %       maps.
        %   kind is one of words 'internal' for internal coordinates and
        %       'mapped'for mapped coordinates.
        %
        %Outputs:
        %   coord is the set of requested projections.
        %   dist is the vector of distances from points to map.
        
            %Check which type of coordinates is necessary to return
            cType = strcmpi('mapped', kind);
            N = size(points, 1);
            
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
                    Y2 = (nodes(face(:, 3), :))';
                    Y20 = Y2-(nodes(face(:, 1), :))';
                    Y21 = Y2-(nodes(face(:, 2), :))';
                    Y10 = (nodes(face(:, 2), :) - nodes(face(:, 1), :))';
                    Y2_2 = sum(Y2 .^ 2);
                    Y2_20 = sum(Y2 .* Y20);
                    Y2_21 = sum(Y2 .* Y21);
                    Y2020 = sum(Y20 .^ 2);
                    Y2120 = sum(Y20 .* Y21);
                    Y2121 = sum(Y21 .^ 2);
                    Y2110 = sum(Y21 .* Y10);
                    Y1010 = sum(Y10 .^ 2);
                    denom = Y2020.*Y2121-Y2120.^2;
                    
                    %Create arrays for results
                    dist = zeros(N, 1);
                    if cType
                        coord = zeros(N, size(points, 2));
                    else
                        coord = zeros(N, size(map.internal, 2));
                    end
                    
                    % maxLeng is maximal number of elements in array to create
                    maxLeng = 100000000;
                    maxPoints = floor(maxLeng / size(face, 1));
                    dataLength = sum(points .^ 2, 2)';
                    nodeLength = sum(face .^2, 2);
                    for i = 1:maxPoints:N
                        % Define last element for calculation
                        last = i + maxPoints - 1;
                        if last > N
                            last = N;
                        end
                        % Prepare index
                        ind = i:last;
                    
                        % Calculate projections
                        A20 = bsxfun(@minus, Y2_20, points(ind, :) * Y20);
                        A21 = bsxfun(@minus, Y2_21, points(ind, :) * Y21);
                        A10 = bsxfun(@rdivide, bsxfun(@minus, A20 - A21,...
                            Y2110), Y1010);
                        A0 = bsxfun(@rdivide, bsxfun(@times, A20, Y2121)...
                            -bsxfun(@times, A21, Y2120), denom);
                        A1 = bsxfun(@rdivide, bsxfun(@times, A21, Y2020)...
                            -bsxfun(@times, A20, Y2120), denom);
                        A20 = bsxfun(@rdivide, A20, Y2020);
                        A21 = bsxfun(@rdivide, A21, Y2121);

                        % Normalize projections
                        A20N = A20;
                        A20N(A20N < 0) = 0;
                        A20N(A20N > 1) = 1;
                        A21N = A21;
                        A21N(A21N < 0) = 0;
                        A21N(A21N > 1) = 1;
                        A10(A10 < 0) = 0;
                        A10(A10 > 1) = 1;
                        tmp = A0 < 0;
                        A0(tmp) = 0;
                        A1(tmp) = A21N(tmp);
                        tmp = A1 < 0;
                        A0(tmp) = A20N(tmp);
                        A1(tmp) = 0;
                        tmp = (1 - (A0 + A1)) < 0;
                        A0(tmp) = A10(tmp);
                        A1(tmp) = 1 - A10(tmp);
                        
                        % Calculate distances
                        d = bsxfun(@plus, dataLength(ind)', Y2_2)...
                            - 2 * points(ind, :) * Y2...
                            + bsxfun(@times, A0, Y2020) .* (A0 - 2 * A20)...
                            + bsxfun(@times, A1, Y2121) .* (A1 - 2 * A21)...
                            + 2 * bsxfun(@times, A0 .* A1, Y2120);
                    
                        % Select the nearest face
                        [dist(ind), tmp] = min(d,[],2);
                        
                        % Form index to find coordinates of projections
                        ind2 = sub2ind(size(d),(1:last - i + 1)',tmp);
                        
                        % calculate coordinates
                        if cType
                            coord(ind, :) = bsxfun(@times, A0(ind2),...
                                nodes(face(tmp, 1), :))...
                                + bsxfun(@times, A1(ind2),...
                                nodes(face(tmp, 2), :))...
                                + bsxfun(@times, 1 - A0(ind2) - A1(ind2),...
                                nodes(face(tmp, 3), :));
                        else
                            coord(ind, :) = bsxfun(@times, A0(ind2),...
                                map.internal(face(tmp, 1), :))...
                                + bsxfun(@times, A1(ind2),...
                                map.internal(face(tmp, 2), :))...
                                + bsxfun(@times, 1 - A0(ind2) - A1(ind2),...
                                map.internal(face(tmp, 3), :));
                        end
                    end
                otherwise
                    error('unacceptable type or projections');
            end
        end
        
        function [dist, klas] = associate(~, node, data)
        %associate identify the nearest node for each data point and
        %return the squared distance between selected node and data
        %point and number of nearest node.
        %
        %Inputs:
        %   node is n-by-k matrix of mapped coordinates for tested state
        %       of map, where n is number of nodes and m is dimension of
        %       data space.
        %   data is m-by-k data points to test, where m is number of
        %       points and k is dimension of data space.
        %
        %Outputs:
        %   dist is m-by-1 matrix of squared distances from data point to
        %       nearest node
        %   klass is m-by-1 vector which contains number of nearest node
        %       for each data point.
            
            % maxLeng is maximal number of elements in array to create 
            maxLeng = 100000000; %800M
            maxPoints = floor(maxLeng / size(node, 1));
            dataLength = sum(data .^ 2, 2)';
            nodeLength = sum(node .^2, 2);
            n = size(data, 1);
            klas = zeros(1, n);
            dist = zeros(1, n);


            for i = 1:maxPoints:n
                % Define last element for calculation
                last = i + maxPoints - 1;
                if last > n
                    last = n;
                end
                % Prepare index
                ind = i:last;
                % Calculate distances
                d = bsxfun(@plus, dataLength(ind), nodeLength) - 2 * (node * data(ind,:)');
                [dist(ind), klas(ind)] = min(d);
            end
            dist = dist';
            klas = klas';
        end
        
        function newMap = extend(map, val, data)
        %extend create extended version of map to reduce/prevent border
        %effect.
        %
        %Inputs:
        %   map is MapGeometry object to extend
        %   val is optional parameter to customise process:
        %       is greater of equal 1 is used to add val ribbons to each
        %           side of map.
        %       is positive number between 0 and 1 means maximal acceptable
        %           fraction of points which are projected onto map border.
        %       default value is 1
        %   data is n-by-m data points to test, where n is number of
        %       points and m is dimension of data space.
        
            % Check the val value
            if nargin < 2
                val = 1;
            elseif val < 0
                error(['Value of val attribute must be positive value',...
                    ' between 0 and 1 for fraction restriction or',...
                    ' positive integer to add val ribbons to each',...
                    ' side of map']);
            end
            if val < 1
                % Restriction for fraction of border cases
                if nargin < 3
                    error('To use 0 < val < 1 data argument must be specified');
                end
                dat = map.preprocessData(data);
                newMap = map;
                while newMap.borderCases(dat, newMap.getBorder()) > val
                    newMap = newMap.extendPrim();
                end
            else
                val = floor(val);
                newMap = map.extendPrim();
                for k = 2:val
                    newMap = newMap.extendPrim();
                end
            end 
        end
        
        function fvu = FVU(map, data, node, type)
        %Calculate fraction of variance unexplained for specified data and
        %nodes.
        %
        %Inputs:
        %   map is MapGeometry object to use
        %   data is set of data points
        %   node is the set of considered mapped nodes. If t is omitted or
        %       empty then the map.mapped is used.
        %   type is the type of projection: 0 means projection into nearest
        %       node of map, 1 means projection onto nearest edge of map, 2
        %       means projection onto nearest face of map. If this argument
        %       is omitted then 1 is used.
        
            %Check the input attributes
            if nargin < 4
                type = 1;
            end
            if nargin<3 || isempty(node)
                node = map.mapped;
            end
            
            %Calculate base variance
            N = size(data,1);
            meanS = sum(data)/N;
            base = sum(data(:).^2)-N*sum(meanS.^2);
            
            %Get distances to map
            [~, dist] = map.projectPrime(node, data, type, 'mapped');
            
            %Calculate FVU
            fvu = sum(dist)/base;
        end
        
        function putMapped(map, newMapped)
            %This function is used for the putting the fitted mapped
            %coordinates of map.
            %
            %Inputs:
            %   newMapped is new matrix of mapped coordinates. It must have
            %       the same size as previously defined matrix
            if ~all(size(map.mapped)==size(newMapped))
                error('Matrix newMapped must have the same size as matrix mapped');
            end
            map.mapped = newMapped;
        end
        
        function frac = borderCases(map, data, list)
        %borderCases calculates fraction of border cases among all data
        %points.
        %
        %Inputs:
        %   map is MapGeometry object to use
        %   data is n-by-m data points to test, where n is number of
        %       points and m is dimension of data space.
        %   list is list of indices of border nodes.
        
            [~, ass] = map.associate(map.getMappedCoordinates, data);
            % Calculate number of points for each node
            N = size(map.getMappedCoordinates, 1);
            tmp = accumarray(ass + 1, 1, [N + 1, 1]);
            % Normalise and remove dummy element
            tmp = tmp(2:end);
            frac = sum(tmp(list)) / size(data, 1);
        end
    end

    methods (Abstract)
        % Primitive extension of map - addition of one ribbon of nodes to
        % map in each direction
        newMap = extendPrim(map)
        % Form list of border nodes
        borders = getBorder(map)
    end
end
