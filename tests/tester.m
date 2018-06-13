what = 12;

% What is parameter to select test
% 0 - generate data for test
% 1 - Generate randomly distorted arc of circle
% 2 - Create two arcs with shift of z coordinate
% 3 - Create two arcs with shift of x coordinate
% 4 - Create fragment of sphere
% 5 - Create two fragments of sphere with shift

% 11 - test oneDMap
% 12 - test rect2DMap
% 13 - test tri2DMap

switch what
    case 0
        data = createTestData(3, 100, [5, 3, 1], [30, 30]);
        save('data3D.mat', 'data', '-v7.3');
    case 1
        data = createTestCurve(3, 100, pi / 4, 20, 1);
        save('data3DArc.mat', 'data', '-v7.3');
    case 2
        data = createTestCurve(3, 500, pi / 4, 20, 1);
        data(1:5:500, 3) = data(1:5:500, 3) + 5;
        save('data3DArc2.mat', 'data', '-v7.3');
    case 3
        data = createTestCurve(3, 500, pi / 4, 20, 1);
        data(1:5:500, 1) = data(1:5:500, 1) + 5;
        save('data3DArc3.mat', 'data', '-v7.3');
    case 4
        data = createTestSurf(3, 200, pi / 4, 20, 1);
        save('data3DSurf.mat', 'data', '-v7.3');
    case 5
        data = createTestSurf(3, 500, pi / 4, 20, 1);
        data(1:5:500, 3) = data(1:5:500, 3) + 5;
        classes = ones(500, 1);
        classes(1:5:500) = classes(1:5:500) + 1;
        save('data3DSurf2.mat', 'data', 'classes', '-v7.3');
    case 11
        map = OneDMap(10);
        init(map, data, 'pci');
%         drawMap(map, data);
%         drawMapInt(map, data);
%         EM(map, data, 'stretch', 0.01, 'bend', 0.1 );
        EM(map, data, 'stretch', 0.001, 'bend', 0.1, 'potential', @L1,... @LLog,...
            'Number_of_intervals', 5);
%         EM(map, data, 'stretch', 0.01, 'bend', 0.1, 'potential', @L2,...
%             'Number_of_intervals', 2, 'intshrinkage', 0.5);
        drawMap(map, data);
        drawMapInt(map, data);
    case 12
        map = rect2DMap(10, 10);
        init(map, data, 'pci');
        if exist('classes', 'var') == 1
            drawMap(map, data, classes, ['b','g']);
        else
            drawMap(map, data);
        end
        drawMapInt(map, data);
%         EM(map, data, 'stretch', 0.01, 'bend', 0.1);
%         EM(map, data, 'stretch', 0.0001, 'bend', 0.1, 'potential', @LLog,...
%             'Number_of_intervals', 5);
        EM(map, data, 'stretch', 0.01, 'bend', 0.1, 'potential', @L2,...
            'Number_of_intervals', 2, 'intshrinkage', 0.7);
        if exist('classes', 'var') == 1
            drawMap(map, data, classes, ['b','g']);
        else
            drawMap(map, data);
        end
        drawMapInt(map, data);
    case 13
        map = tri2DMap(4, 4);
        init(map, data, 'pci');
        drawMap(map, data);
        drawMapInt(map, data);
        EM(map, data, 'stretch', 0.001, 'bend', 0.01);
        drawMap(map, data);
        drawMapInt(map, data);
end
