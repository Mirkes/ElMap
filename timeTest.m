%test projection time
tic;
Xp2 = project(mape,x,2,'internal');
toc
tic;
Xp12 = project(mape,x,2,'mapped');
toc