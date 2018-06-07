function EM(map, data, varargin)
%EM is function to fit map to data
%   Syntax
%   EM( map, data )
%   EM( __, Name, Value )
%   
%Inputs:
%   map is an object of MapGeometry class or subclass of this class.
%   data is n-by-m matrix which is the data set to fit. Each row contains
%       m coordinates of one data point.
%   Name can be one of the following:
%       'type' is one of the following strings:
%          'hard' is hard map with stretch = 1 and bend = 1
%           'medium' is more flexible map with stretch = 0.7 and bend = 0.7
%           'soft'  is soft map with stretch = 0.5 and bend = 0.5
%           If 'type', 'stretch' and 'bending' are omitted then 'medium' is
%           used. 
%       'stretch' is a positive numeric value which is represent the value
%           of stretching modulo or a function with syntax 
%               val = stretch( epoch )
%           where epoch is number of epoch (see epoch definition below) and
%           val is the nonnegative stretching modulo to use on specified
%           epoch. Epochs are numerated from 1.
%           Default value corresponds to type 'medium'
%       'bend' is a positive numeric value which is represent the value of
%           bending modulo or a function with syntax
%               val = bend( epoch )
%           where epoch is number of epoch (see epoch definition below) and
%           val is the bending modulo to use on specified epoch. Epochs are
%           numerated from 1. 
%           Default value corresponds to type 'medium'
%       'weights' is n-by-1 vector of weights for data points. Weights must
%           be nonnegative.
%
% One epoch is fitting of map with fixed values of stretching and bending
% modulo. THis process can include several iterations of two step
% algorothm:
%   1. associate each datapoint with nearest node.
%   2. recalculate node position.
% Process of map fitting is stopped if new values of stratching and bending
% modulo are the same as on previous epoch OR if both stratching and bending
% modulo are zero.


    % Check the number of input attributes and types of the two first
    % attributes.
    if nargin < 2
        error('At least map and data must be specified');
    end    
    if ~isa(map,'MapGeometry')
        error('Incorrect type of the first argument');
    end
    if ~ismatrix(data) || ~isnumeric(data)
        error('Incorrect type of the second argument');
    end
    
    % Get sizes of data
    [n, dim] = size(data);
    
    % Decode varargin
    strFun = @constStretch;
    constStretching = 0.7;
    bendFun = @constBend;
    constBending = 0.7;
    weights = [];
    
    for i=1:2:length(varargin)
        if strcmpi(varargin{i}, 'type')
            switch lower(varargin{i + 1})
                case 'hard'
                    strFun = @constStretch;
                    constStretching = 1;
                    bendFun = @constBend;
                    constBending = 1;
                case 'medium'
                    strFun = @constStretch;
                    constStretching = 0.7;
                    bendFun = @constBend;
                    constBending = 0.7;
                case 'soft'                    
                    strFun = @constStretch;
                    constStretching = 0.5;
                    bendFun = @constBend;
                    constBending = 0.5;
                otherwise
                    error('Incorrect value for type argument');
            end
        elseif strcmpi(varargin{i}, 'stretch')
            tmp = varargin{i + 1};
            if isa(tmp, 'function_handle')
                strFun = tmp;
            else
                strFun = @constStretch;
                constStretching = tmp;
            end;
        elseif strcmpi(varargin{i}, 'bend')
            tmp = varargin{i + 1};
            if isa(tmp, 'function_handle')
                bendFun = tmp;
            else
                bendFun = @constBend;
                constBending = tmp;
            end;
        elseif strcmpi(varargin{i}, 'weights')
            weights = varargin{i + 1};
        else
            if ischar(varargin{i})
                error(['Wrong name of argument "', varargin{i}, '"']);
            else
                error(['Wrong name of argument "', num2str(varargin{i}), '"']);
            end
        end
    end

    % Check type and length of weights
    if isempty(weights)
        weights = ones(n, 1);
    elseif ~ismatrix(data) || ~isnumeric(data)
        error('Incorrect type weights argument');
    elseif size(weights, 1) ~= n || size(weights, 2) ~= 1
        error('Incorrect size of weights argument');
    end
    
    % Define weights. Now it is trash but maybe...
    TotalWeight = sum(weights);
    
    %Get initial state of nodes
    nodes = map.getMappedCoordinates();
    if size(nodes, 2) ~= dim
        error('Dimensions of mapped nodes and data must be the same');
    end
    N = size(nodes, 1);

    %Form matrices B and C
    tmp = map.getLinks();
    B = diag(accumarray(tmp(:), 1, [N, 1]));
    tmp = accumarray(tmp, 1, [N, N]);
    B = B - tmp - tmp';
    
    tmp = map.getRibs();
    C = diag(accumarray([tmp(:, 1); tmp(:, 3)], 1, [N, 1])...
        +accumarray(tmp(:, 2), 4, [N, 1]));
    w = accumarray(tmp(:, [1, 3]), 1,[N, N]);
    tmp = accumarray([tmp(:, 1:2); tmp(:, 2:3)], 2, [N, N]);
    C = C + w + w' - tmp - tmp';

    % Start iterative process
    epoch = 1; % Number of iteration
    ass = zeros(n, 1); % Initial associations. It is impossible combination
    % Get initial modulo
    stretch = strFun(epoch);
    bend = bendFun(epoch);
    while true
        % Save old associations.
        oldAss = ass;
        % Find new associations
        [~, ass] = associate(map, nodes, data);
        % If nothing is changed then we have end of epoch
        if all(oldAss == ass)
            epoch = epoch + 1;
            tmp = strFun(epoch);
            tmp1 = bendFun(epoch);
            if tmp == 0 && tmp1 == 0
                break;
            end
            if abs(tmp - stretch) + abs(tmp1 - bend) == 0
                break;
            end
            stretch = tmp;
            bend = tmp1;
        end
        
        % Form matrix A
        % For further robust and so on we consider possibility of zeros in
        % ass and create dummy element
        ass = ass + 1;
        % Calculate number of points for each node
        tmp = accumarray(ass, weights, [N + 1, 1]);
        % Normalise and remove dummy element
        NodeClusterRelativeSize = tmp(2:end) / TotalWeight;
        % To prevent appearance of NaN
        tmp(tmp == 0) = 1;
        NodeClusterCenters = zeros(N + 1, dim);
        for k = 1:dim
            NodeClusterCenters(:, k) =...
                accumarray(ass, data(:, k), [N + 1, 1]) ./ tmp;
        end
        % Remove dummy element
        NodeClusterCenters = NodeClusterCenters(2:end,:);
        
        % form SLAE
        SLAUMatrix = diag(NodeClusterRelativeSize) + stretch * B + bend * C;
        nodes = SLAUMatrix...
            \bsxfun(@times, NodeClusterRelativeSize, NodeClusterCenters);
    
        % Restore ass
        ass = ass - 1;
    end
    
    % Put new nodes into map
    map.putMapped(nodes);
    
    function stretch = constStretch( ~ )
        stretch = constStretching;
    end

    function bend = constBend( ~ )
        bend = constBending;
    end
end


