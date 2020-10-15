function [dist, klas] = associate4(maxLeng, node, data)
%associate identify the nearest node for each data point and
%return the squared distance between selected node and data
%point and number of nearest node.
%
%Inputs:
%   node is n-by-k matrix of mapped coordinates for tested state
%       of map, where n is number of nodes and m is dimension of
%       data space.
%   data is m-by-k data points to test, where m is number of
%       points and k is dimension of data space.
%
%Outputs:
%   dist is m-by-1 matrix of squared distances from data point to
%       nearest node
%   klass is m-by-1 vector which contains number of nearest node
%       for each data point.

    % maxLeng is maximal number of elements in array to create
    % maxLeng = 10000000;
    maxPoints = floor(maxLeng / size(node, 1));
    dataLength = sum(data .^ 2, 2);
    node = node';
    nodeLength = sum(node .^2);
    n = size(data, 1);
    klas = zeros(n, 1);
    dist = zeros(n, 1);
        
    for i = 1:maxPoints:n
        % Define last element for calculation
        last = i + maxPoints - 1;
        if last > n
            last = n;
        end
        % Prepare index
        ind = i:last;
        % Calculate distances
        d = bsxfun(@plus, dataLength(ind), nodeLength) - 2 * (data(ind,:) * node);
        [dist(ind), klas(ind)] = min(d,[],2);
    end
end
