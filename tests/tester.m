what = 2;

% What is parameter to select test
% 0 - generate data for test
% 1 - load data from file
% 2 - test oneDMap
% 3 - test rect2DMap
% 4 - test tri2DMap
% 5 - Generate randomly distorted arc of circle
% 6 - load arc of circle from arc
% 7 - Create two arcs with shift of z coordinate

switch what
    case 0
        data = createTestData(3, 100, [5, 3, 1], [30, 30]);
        save('data3D.mat', 'data', '-v7.3');
    case 1
        load('data3D.mat', 'data');
    case 2
        map = OneDMap(10);
        init(map, data, 'pci');
%         drawMap(map, data);
%         drawMapInt(map, data);
%         EM(map, data, 'stretch', 0.01, 'bend', 0.1 );
%         EM(map, data, 'stretch', 0.01, 'bend', 0.1, 'potential', @LLog,...
%             'Number_of_intervals', 5);
        EM(map, data, 'stretch', 0.01, 'bend', 0.1, 'potential', @L2,...
            'Number_of_intervals', 2, 'intshrinkage', 0.7);
        drawMap(map, data);
        drawMapInt(map, data);
    case 3
        map = rect2DMap(4, 4);
        init(map, data, 'pci');
        drawMap(map, data);
        drawMapInt(map, data);
        EM(map, data, 'stretch', 0.001, 'bend', 0.01);
        drawMap(map, data);
        drawMapInt(map, data);
    case 4
        map = tri2DMap(4, 4);
        init(map, data, 'pci');
        drawMap(map, data);
        drawMapInt(map, data);
        EM(map, data, 'stretch', 0.001, 'bend', 0.01);
        drawMap(map, data);
        drawMapInt(map, data);
    case 5 
        data = createTestCurve(3, 100, pi / 4, 20, 1);
        save('data3DArc.mat', 'data', '-v7.3');
    case 6
        load('data3DArc.mat', 'data');
    case 7
        data = createTestCurve(3, 500, pi / 4, 20, 1);
        data(1:5:500, 3) = data(1:5:500, 3) + 5;
        save('data3DArc2.mat', 'data', '-v7.3');
        
end
