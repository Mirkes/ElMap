function drawMapInt( map, data, projType, varargin )
%drawMapInt create figure and draw data and maps in the internal maps
%coordinates.
%   map is object of class MapGeometry (descendant of this class).
%   data is matrix of data points. each row is one point.
%   projType is the type of projection: 
%       0 is projection to the nearest node,
%       1 is projection to the nearest edge,
%       2 is projection to the nearest face.
%       default value is 0.
%   There are several additional parameters in form m'Name', value:
%       'nodeMarker' is one of the possible marker shape symbol ('o', '+',
%           '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h') or
%           'none'. Default value is 'o'.
%       'nodeMarkerSize' is positive number which is size of node marker.
%           This parameter is ignored for the projection of the first type.
%           Default value is 6;
%       'nodeColour' is one of the possible Matlab colours ('r', 'g', 'b',
%           'y', 'm', 'c', 'w', 'k'). Default value is 'r';
%       'lineWidth' is non-negative number for map line width. Default
%           value is 2.

    if nargin<3
        projType = 0;
    end
    
    % Parse optional parameters
    classes = [];
    markColour = [];
    markShape = [];
    markSize = [];
    nodeMarker = 'o';
    nodeMarkerSize = 6;
    nodeColour = 'r';
    lineWidth = 2;
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case 'classes'
                classes = varargin{i + 1};
            case 'markcolour'
                markColour = varargin{i + 1};
            case 'markshape'
                markShape = varargin{i + 1};
            case 'marksize'
                markSize = varargin{i + 1};
            case 'nodemarker'
                nodeMarker = varargin{i + 1};
            case 'nodemarkersize'
                nodeMarkerSize = varargin{i + 1};
            case 'nodecolour'
                nodeColour = varargin{i + 1};
            case 'linewidth'
                lineWidth = varargin{i + 1};
            otherwise
                error(['Unknown argument at position "', num2str(i + 2)]);
        end
    end

%     N = size(data, 1);

    % Sanity check of arguments.
    if isempty(classes)
        classes = ones(size(data, 1), 1);
        cls = 1;
    else
        cls = unique(classes);
    end
    
    nCls = length(cls);
    
    if isempty(markColour)
        markColour = repmat('b', nCls, 1);
    end

    if isempty(markShape)
        markShape = repmat('s', nCls, 1);
    end
    
    if isempty(markSize)
        markSize = repmat(6, nCls, 1);
    end
    
    % Create figure
    figure;
    % Get map coordinates
    maps = map.getInternalCoordinates;
    % Get limits of map
    lims = zeros(1, 4);
    lims([1, 3]) = min(maps) - 1;
    lims([2, 4]) = max(maps) + 1;
    % Get map links
    links = map.getLinks;
    % Project data onto map
    if ~isempty(data)
        data = map.project(data, projType, 'internal');
    end
    % Get data dimension
    dim = size(maps,2);
    if dim > 3
        
    elseif dim == 3
        %3d data
        
    elseif dim == 2
        %two dimensional data
        %Project data to map
        hold on
        %Draw edges
        %Prepare arrays
        X=[maps(links(:,1),1)';maps(links(:,2),1)'];
        Y=[maps(links(:,1),2)';maps(links(:,2),2)'];
        %Draw edges
        plot(X, Y, nodeColour, 'LineWidth', lineWidth);
        if projType == 0
            if ~isempty(data)
                %Search unique points
                [dat, ~, ic] = unique(data,'rows');
                count = accumarray(ic, 1);
                ma = max(count);
                %Draw map nodes
                if nCls == 1 % No classes
                    scatter(dat(:, 1), dat(:, 2), count/ma * 400,...
                        nodeColour, 'filled', nodeMarker);
                else
                    % We have classes We need to draw each marker by it's
                    % own colour or several colours
                    % First of all calculate number of points of each class
                    % in each node.
                    nDat = size(dat, 1);
                    props = zeros(nDat, nCls);
                    for k = 1:nCls
                        ind = classes == cls(k);
                        props(:, k) = accumarray(ic(ind), 1, [nDat, 1]);
                    end
                    % Now we put pie Chart in each necesssary node.
                    for k = 1:nDat
                        drawPieChart(dat(k, 1), dat(k, 2),...
                            0.5 * count(k) / ma, props(k, :), markColour);
                    end
                end
            end
            axis(lims);
        else
            %Draw maps nodes
            plot(maps(:,1), maps(:,2), 'Marker', nodeMarker,...
            'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
            'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');

            %Draw data points
            if ~isempty(data)
                for k = 1:nCls
                    ind = classes == cls(k);
                    plot(data(ind, 1), data(ind, 2),...
                        [markColour(k), markShape(k)],...
                        'MarkerFaceColor', markColour(k),...
                        'MarkerSize', markSize(k));
                    hold on
                end
            end
        end
    else
        %one dimensional data
        hold on
        %Draw edges
        %Prepare arrays
        X=[maps(links(:,1),1)';maps(links(:,2),1)'];
        Y=zeros(2,size(links,1));
        %Draw edges
        plot(X, Y, nodeColour, 'LineWidth', lineWidth);
        %Draw map nodes
        if projType == 0
            if ~isempty(data)
                %Search unique points
                [dat, ~, ic] = unique(data,'rows');
                count = accumarray(ic,1);
                ma=max(count);
                %Draw map nodes
                if nCls == 1 % No classes
                    scatter(dat(:,1), zeros(size(dat,1),1),...
                        count / ma * 400, nodeColour, 'filled', nodeMarker);
                else
                    % We have classes We need to draw each marker by it's
                    % own colour or several colours
                    % First of all calculate number of points of each class
                    % in each node.
                    nDat = size(dat, 1);
                    props = zeros(nDat, nCls);
                    for k = 1:nCls
                        ind = classes == cls(k);
                        props(:, k) = accumarray(ic(ind), 1, [nDat, 1]);
                    end
                    % Now we put pie Chart in each necesssary node.
                    for k = 1:nDat
                        drawPieChart(dat(k, 1), 0,...
                            0.5 * count(k) / ma, props(k, :), markColour);
                    end
                    axis([lims(1:2), lims(1:2) - sum(lims(1:2)) / 2]);
                end
            end
        else
            %Draw maps nodes
            plot(maps(:,1),0,'ro','MarkerFaceColor','r','MarkerSize',10);
            %Draw data points
            if ~isempty(data)
                for k = 1:nCls
                    ind = classes == cls(k);
                    plot(data(ind, 1), 0,...
                        [markColour(k), markShape(k)],...
                        'MarkerFaceColor', markColour(k),...
                        'MarkerSize', markSize(k));
                    hold on
                end
            end
        end
    end
end

