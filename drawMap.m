function drawMap( map, data )
%drawMap create figure and draw data and maps. If data dimension (number of
%columns) is less than 3 then original data coordinate system is used. If
%data dimension is greater than 3 then projection into space of the first
%three PCs is used for drowing.
%   map is object of class MapGeometry (descendant of this class).
%   data is matrix of data points. each row is one point.

    %Create figure
    figure;
    axis equal;
    %Get map coordinates
    maps = map.getMappedCoordinates;
    %get map links
    links = map.getLinks;
    %Get data dimension
    dim = size(data,2);
    if dim>3 % Dimension is greater than 3 and we need to use PCs
        % Calculate three first PCs and project data onto PCs
        [~, ~, V] = svds(data, 3);
        data = data * V;
        %Draw data
        plot3(data(:, 1), data(:, 2), data(:, 3), 'sb', 'MarkerFaceColor','b');
        hold on
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
    elseif dim==3
        %3d data
        %Draw data
        plot3(data(:, 1), data(:, 2), data(:, 3), 'sb', 'MarkerFaceColor','b');
        hold on
        %Draw edges
        %Prepare arrays
        X=[maps(links(:,1), 1)';maps(links(:,2), 1)'];
        Y=[maps(links(:,1), 2)';maps(links(:,2), 2)'];
        Z=[maps(links(:,1), 3)';maps(links(:,2), 3)'];
        %Draw edges
        plot3(X, Y, Z, 'r', 'LineWidth', 2);
        %Draw map nodes
        plot3(maps(:, 1), maps(:, 2), maps(:, 3), 'ro', 'MarkerFaceColor', 'r');        
    elseif dim==2
        %two dimensional data
        %Draw data
        plot(data(:,1),data(:,2),'sb','MarkerFaceColor','b');
        hold on
        %Draw edges
        %Prepare arrays
        X=[maps(links(:,1),1)';maps(links(:,2),1)'];
        Y=[maps(links(:,1),2)';maps(links(:,2),2)'];
        %Draw edges
        plot(X,Y,'r','LineWidth',2);
        %Draw map nodes
        plot(maps(:,1),maps(:,2),'ro','MarkerFaceColor','r');        
    else
        %one dimensional data
        %Draw data
        plot(data(:,1),0,'sb','MarkerFaceColor','b');
        hold on
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

