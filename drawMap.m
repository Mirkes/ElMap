function drawMap(map, data, varargin)
%drawMap create figure and draw data and maps. If data dimension (number of
%columns) is less than 3 then original data coordinate system is used. If
%data dimension is greater than 3 then projection into space of the first
%three PCs is used for drowing.
%
%Usage
%   drawMap(map, data);
%   drawMap(__, Name, Value);
%
%Inputs:
%   map is object of class MapGeometry (descendant of this class).
%   data is n-by-dim matrix of data points. Each row is one point.
%   There are several possible Name, value pairs:
%       'classes' is n-by-1 vector of class labels for points. Each label
%           is positive integer number which is the index of cell array
%           with marker descriptions. Marker descriptions can be specified
%           by user in the arguments markColour, markShape, markSize.
%           Otherwise standard marker is used for all points.
%       'markColour' is K-by-1 vector of standard Matlab colours ('r', 'g',
%           'b', 'y', 'm', 'c', 'w', 'k'). Value K is number of defined
%           markers. If markColour is omitted then 'b' (blue) is used for
%           all markers. 
%       'markShape' is K-by-1 vector of standard Matlab marker shapes ('o',
%           '+', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h').
%           Value K is number of defined markers. If markShape is omitted
%           then 's' (square) is used for all markers.
%       'markSize' is K-by-1 vector of positive numbers. Value K is number
%           of defined markers. If markSize is omitted then 6 is used for
%           all markers.
%       'nodeMarker' is one of the possible marker shape symbol ('o', '+',
%           '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h') or
%           'none'. Default value is 'o'.
%       'nodeMarkerSize' is positive number which is size of node marker.
%           Default value is 6;
%       'nodeColour' is one of the possible Matlab colours ('r', 'g', 'b',
%           'y', 'm', 'c', 'w', 'k'). Default value is 'r';
%       'lineWidth' is non-negative number for map line width. Zero means
%           absense of map at all (widht of line is zero and 'nodeMarker'
%           is 'none'. Default value is 2.
%       'newFigure' is logical argument. Value true (default) causes
%           creation of new figure. Value false causes usage of current
%           active figure. THis option can be used for subplots and for map
%           colouring.

    % Data preprocessing
    data = map.preprocessData(data);

    %Get data dimension
    [N, dim] = size(data);

    % Parse varargin
    classes = [];
    markColour = [];
    markShape = [];
    markSize = [];
    nodeMarker = 'o';
    nodeMarkerSize = 6;
    nodeColour = 'r';
    lineWidth = 2;
    newFigure = false;
    
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
            case 'newFigure'
                newFigure = varargin{i + 1};
            otherwise
                error(['Unknown argument at position "', num2str(i + 2)]);
        end
    end
    
    % Sanity check of arguments.
    if isempty(classes)
        classes = ones(N, 1);
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
    
    if lineWidth == 0
        mapDraw = 0;
    else
        mapDraw = 1;
    end
    
    %Create figure
    if newFigure
        figure;
    end
    %Get map coordinates
    maps = map.getMappedCoordinates;
    %get map links
    links = map.getLinks;
    if dim > 3 % Dimension is greater than 3 and we need to use PCs
        % Calculate three first PCs and project data onto PCs
        if ~map.preproc
            dat = bsxfun(@minus, data, map.means);
            V = map.PCs(:, 1:3);
            data = dat * V;
            maps = bsxfun(@minus, maps, map.means) * V;
        end
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot3(data(ind, 1), data(ind, 2), data(ind, 3),...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        if mapDraw
            %Draw edges
            %Prepare arrays
            X=[maps(links(:,1), 1)';maps(links(:,2), 1)'];
            Y=[maps(links(:,1), 2)';maps(links(:,2), 2)'];
            Z=[maps(links(:,1), 3)';maps(links(:,2), 3)'];
            %Draw edges
            plot3(X, Y, Z, nodeColour, 'LineWidth', lineWidth);
            %Draw map nodes
            plot3(maps(:, 1), maps(:, 2), maps(:, 3), 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');
        end
        axis equal;
    elseif dim == 3
        %3d data
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot3(data(ind, 1), data(ind, 2), data(ind, 3),...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        if mapDraw
            %Draw edges
            %Prepare arrays
            X=[maps(links(:,1), 1)';maps(links(:,2), 1)'];
            Y=[maps(links(:,1), 2)';maps(links(:,2), 2)'];
            Z=[maps(links(:,1), 3)';maps(links(:,2), 3)'];
            %Draw edges
            plot3(X, Y, Z, nodeColour, 'LineWidth', lineWidth);
            %Draw map nodes
            plot3(maps(:, 1), maps(:, 2), maps(:, 3), 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');
        end
        axis equal;
    elseif dim == 2
        %two dimensional data
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot(data(ind, 1), data(ind, 2),...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        if mapDraw
            %Draw edges
            %Prepare arrays
            X=[maps(links(:,1),1)';maps(links(:,2),1)'];
            Y=[maps(links(:,1),2)';maps(links(:,2),2)'];
            %Draw edges
            plot(X, Y, nodeColour, 'LineWidth', lineWidth);
            %Draw map nodes
            plot(maps(:,1),maps(:,2), 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');
        end
        axis equal;
    else
        %one dimensional data
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot(data(ind, 1), 0,...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        if mapDraw
            
            %Draw edges
            %Prepare arrays
            X=[maps(links(:,1),1)'; maps(links(:,2),1)'];
            Y=zeros(2,size(links,1));
            %Draw edges
            plot(X, Y, nodeColour, 'LineWidth', lineWidth);
            %Draw map nodes
            plot(maps(:,1),0, 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');
        end
    end
end

