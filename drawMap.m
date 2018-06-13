function drawMap( map, data, classes, markColour, markShape, markSize )
%drawMap create figure and draw data and maps. If data dimension (number of
%columns) is less than 3 then original data coordinate system is used. If
%data dimension is greater than 3 then projection into space of the first
%three PCs is used for drowing.
%   map is object of class MapGeometry (descendant of this class).
%   data is n-by-dim matrix of data points. Each row is one point.
%   classes is n-by-1 vector of class labels for points. Each label is
%       positive integer number which is the index of cell array with
%       marker descriptions. Marker descriptions can be specified by user
%       in the arguments markColour, markShape, markSize. Otherwise
%       standard markers are used.
%   markColour is K-by-1 vector of standard Matlab colours ('r', 'g', 'b',
%       'y', 'm', 'c', 'w', 'k'). Value K is number of defined markers. If
%       markColour is omitted then 'b' is used for all markers.
%   markShape is K-by-1 vector of standard Matlab marker shapes ('o', '+',
%       '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h'). Value K is
%       number of defined markers. If markShape is omitted then 's' is used
%       for all markers.
%   markSize is K-by-1 vector of positive numbers. Value K is number of
%       defined markers. If markSize is omitted then 6 is used 
%       for all markers.

    %Get data dimension
    [N, dim] = size(data);

    if nargin<3 || isempty(classes)
        classes = ones(N, 1);
        cls = 1;
    else
        cls = unique(classes);
    end
    
    nCls = length(cls);
    
    if nargin < 4 || isempty(markColour)
        markColour = repmat('b', nCls, 1);
    end

    if nargin < 5 || isempty(markShape)
        markShape = repmat('s', nCls, 1);
    end
    
    if nargin < 6 || isempty(markSize)
        markSize = repmat(6, nCls, 1);
    end
    
    %Create figure
    figure;
    %Get map coordinates
    maps = map.getMappedCoordinates;
    %get map links
    links = map.getLinks;
    if dim>3 % Dimension is greater than 3 and we need to use PCs
        % Calculate three first PCs and project data onto PCs
        [~, ~, V] = svds(data, 3);
        data = data * V;
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot3(data(ind, 1), data(ind, 2), data(ind, 3),...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        %Draw edges
        %Prepare arrays
        maps = maps * V;
        X=[maps(links(:,1), 1)';maps(links(:,2), 1)'];
        Y=[maps(links(:,1), 2)';maps(links(:,2), 2)'];
        Z=[maps(links(:,1), 3)';maps(links(:,2), 3)'];
        %Draw edges
        plot3(X, Y, Z, 'r', 'LineWidth', 2);
        %Draw map nodes
        plot3(maps(:, 1), maps(:, 2), maps(:, 3), 'ro', 'MarkerFaceColor', 'r');        
        axis equal;
    elseif dim==3
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
        %Draw edges
        %Prepare arrays
        X=[maps(links(:,1), 1)';maps(links(:,2), 1)'];
        Y=[maps(links(:,1), 2)';maps(links(:,2), 2)'];
        Z=[maps(links(:,1), 3)';maps(links(:,2), 3)'];
        %Draw edges
        plot3(X, Y, Z, 'r', 'LineWidth', 2);
        %Draw map nodes
        plot3(maps(:, 1), maps(:, 2), maps(:, 3), 'ro', 'MarkerFaceColor', 'r');        
        axis equal;
    elseif dim==2
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
        %Draw edges
        %Prepare arrays
        X=[maps(links(:,1),1)';maps(links(:,2),1)'];
        Y=[maps(links(:,1),2)';maps(links(:,2),2)'];
        %Draw edges
        plot(X,Y,'r','LineWidth',2);
        %Draw map nodes
        plot(maps(:,1),maps(:,2),'ro','MarkerFaceColor','r');        
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
        %Draw edges
        %Prepare arrays
        X=[maps(links(:,1),1)';maps(links(:,2),1)'];
        Y=zeros(2,size(links,1));
        %Draw edges
        plot(X,Y,'r','LineWidth',2);
        %Draw map nodes
        plot(maps(:,1),0,'ro','MarkerFaceColor','r');        
    end
end

