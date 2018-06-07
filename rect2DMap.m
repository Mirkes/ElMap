classdef rect2DMap < MapGeometry
    %OneDMap is descendant for piecewise linear map
    %   Constructor contains two argument which are number of rows and
    %   number of columns
    
    properties (SetAccess = protected)
        faces   %Array of faces to projections
    end
    
    methods
        function map = rect2DMap(rows, cols)
            if nargin < 2
                error('You MUST specify number of rows and columns for rect2DMap');
            end
            %Create map
            map@MapGeometry(2);
            %Calculate internal coordinates of nodes
            N=rows*cols;
            a1=repmat(1:cols,rows,1);
            a2=repmat((1:rows)',1,cols);
            map.internal = [a1(:),a2(:)];
            %form array of links
            A=reshape(1:N,rows,cols);
            B=A(1:end-1,:);
            C=A(2:end,:);
            D=A(:,1:end-1);
            E=A(:,2:end);
            map.links = [B(:), C(:); D(:), E(:)];
            %form array of ribs
            B1=A(1:end-2,:);
            B2=A(2:end-1,:);
            B3=A(3:end,:);
            C1=A(:,1:end-2);
            C2=A(:,2:end-1);
            C3=A(:,3:end);
            map.ribs = [B1(:), B2(:), B3(:); C1(:), C2(:), C3(:)];
            %form array of faces
            B1 = A(1:end-1,1:end-1);
            B2 = A(2:end,1:end-1);
            B3 = A(1:end-1,2:end);
            C1 = A(1:end-1,2:end);
            C2 = A(2:end,1:end-1);
            C3 = A(2:end,2:end);
            map.faces=[B1(:), B2(:), B3(:); C1(:), C2(:), C3(:)];
            %Set mappedcoordinates to empty set
            map.mapped = [];
        end
        
        function face = getFaces(map)
            %method to access to the feces of map
            %face is k-by-3 matrix. Each row contains three numbers of
            %nodes which form one face.
            face=map.faces;
        end
    end
    
end

