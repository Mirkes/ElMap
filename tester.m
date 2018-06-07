what = 4;

% What is parameter to select test
% 0 - generate data for test
% 1 - load data from file
% 2 - test oneDMap
% 3 - test rect2DMap
% 4 - test tri2DMap

switch what
    case 0
        data = createTestData(3, 100, [5, 3, 1], [30, 30]);
        save('data3D.mat', 'data', '-v7.3');
    case 1
        load('data3D.mat', 'data', '-v7.3');
    case 2
        map = OneDMap(4);
        init(map, data, 'pci');
        drawMap(map, data);
        drawMapInt(map, data);
        EM(map, data, 'stretch', 0.01, 'bend', 0.1);
        drawMap(map, data);
        drawMapInt(map, data);
    case 3
        map = rect2DMap(4, 4);
        init(map, data, 'pci');
        drawMap(map, data);
        drawMapInt(map, data);
        EM(map, data, 'stretch', 0.01, 'bend', 0.1);
        drawMap(map, data);
        drawMapInt(map, data);
    case 4
        map = tri2DMap(4, 4);
        init(map, data, 'pci');
        drawMap(map, data);
        drawMapInt(map, data);
        EM(map, data, 'stretch', 0.01, 'bend', 0.1);
        drawMap(map, data);
        drawMapInt(map, data);
end
