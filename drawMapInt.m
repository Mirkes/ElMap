function drawMapInt( map, data, projType )
%drawMapInt create figure and draw data and maps. If data dimension (number
%of columns) is less than 3 then original data coordinate system is used.
%If data dimension is greater than 3 then projection into space of the
%first 
%three PCs is used for drowing.
%   map is object of class MapGeometry (descendant of this class).
%   data is matrix of data points. each row is one point.
%   projType is the type of projection: 
%       0 is projection to the nearest node,
%       1 is projection to the nearest edge,
%       2 is projection to the nearest face.
%       default value is 0.

    if nargin<3
        projType = 0;
    end
    %Create figure
    figure;
    %Get map coordinates
    maps = map.getInternalCoordinates;
    %get map links
    links = map.getLinks;
    %Project data onto map
    if ~isempty(data)
        data = map.project(data,projType,'internal');
    end
    %Get data dimension
    dim = size(maps,2);
    if dim>3
        
    elseif dim==3
        %3d data
        
    elseif dim==2
        %two dimensional data
        %Project data to map
        hold on
        %Draw edges
        %Prepare arrays
        X=[maps(links(:,1),1)';maps(links(:,2),1)'];
        Y=[maps(links(:,1),2)';maps(links(:,2),2)'];
        %Draw edges
        plot(X,Y,'r');
        if projType==0
            if ~isempty(data)
                %Search unique points
                [dat, ~, ic] = unique(data,'rows');
                count = accumarray(ic,1);
                ma=max(count);
                %Draw map nodes
                scatter(dat(:,1),dat(:,2),count/ma*400,'r','filled','o');
            end
        else
            %Draw maps nodes
            plot(maps(:,1),maps(:,2),'ro','MarkerFaceColor','r','MarkerSize',10);
            %Draw data points
            if ~isempty(data)
                plot(data(:,1),data(:,2),'bo','MarkerFaceColor','b');
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
        plot(X,Y,'r','LineWidth',2);
        %Draw map nodes
        if projType==0
            if ~isempty(data)
                %Search unique points
                [dat, ~, ic] = unique(data,'rows');
                count = accumarray(ic,1);
                ma=max(count);
                %Draw map nodes
                scatter(dat(:,1),zeros(size(dat,1),1),count/ma*400,'r','filled','o');
            end
        else
            %Draw maps nodes
            plot(maps(:,1),0,'ro','MarkerFaceColor','r','MarkerSize',10);
            %Draw data points
            if ~isempty(data)
                plot(data(:,1),zeros(size(data,1),1),'bs','MarkerFaceColor','b');
            end
        end
    end
end

