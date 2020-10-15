function [dist, klas] = associate2(~, node, data)
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

    dist = bsxfun(@plus, sum(data .^ 2, 2)',...
        sum(node .^2, 2)) - 2 * (node * data');
    [dist, klas] = min(dist);
    dist = dist';
    klas = klas';
end
