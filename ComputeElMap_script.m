x = load('tests/_input.txt');

map = rect2DMap(40,40);
init(map,x,'pci');
EM(map,x,'bend',0.2,'stretch',0.1); 
EM(map,x,'bend',0.2,'stretch',0.01);
EM(map,x,'bend',0.2,'stretch',0.001);
mape = extend(map,0.01,x);
Xp = project(mape,x,2,'internal');

fid = fopen('tests/_elmap_proj.txt','w');
for i=1:size(Xp,1)
    for j=1:size(Xp,2)
        fprintf(fid,'%f\t',Xp(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('tests/_elmap_coords.txt','w');
for i=1:size(mape.mapped,1)
    for j=1:size(mape.mapped,2)
        fprintf(fid,'%f\t',mape.mapped(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('tests/_elmap_edges.txt','w');
for i=1:size(mape.links,1)
    for j=1:size(mape.links,2)
        fprintf(fid,'%i\t',mape.links(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);


%exit();