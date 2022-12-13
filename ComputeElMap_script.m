x = load('tests/_input.txt');

map = rect2DMap(40,40);
init(map,x,'pci');
EM(map,x,'bend',0.2,'stretch',0.1); 
EM(map,x,'bend',0.2,'stretch',0.01);
EM(map,x,'bend',0.2,'stretch',0.001);
mape = extend(map,0.01,x);


% Maps drawing
drawMapInt(mape, x, 2, 'coloring', 'density', 'drawData', false);
print(gcf, '-dpng', '-noui', '-loose', 'MapInInternalCoordinates.png');
drawMap(mape, x, 'drawData', false);
print(gcf, '-dpng', '-noui', '-loose', 'MapEmbeddedToOriginalSpace.png');

To use map outside of current software we can use following code
Xp = project(mape,x,2,'internal');

% Write projections of data points to map
fid = fopen('tests/_elmap_proj.txt','w');
for i=1:size(Xp,1)
    for j=1:size(Xp,2)
        fprintf(fid,'%f\t',Xp(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% Write coordinates of nodes ebedded to original space
fid = fopen('tests/_elmap_coords.txt','w');
for i=1:size(mape.mapped,1)
    for j=1:size(mape.mapped,2)
        fprintf(fid,'%f\t',mape.mapped(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% write edges as pairs of connected nodes
fid = fopen('tests/_elmap_edges.txt','w');
for i=1:size(mape.links,1)
    for j=1:size(mape.links,2)
        fprintf(fid,'%i\t',mape.links(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

%exit();
