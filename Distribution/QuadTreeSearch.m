function [index,tree] = QuadTreeSearch(X,xq)
% Quadtree binary search algorithm for 2D data
% INPUTS:
%   X = reference data (rows = observations, columns = coordinates)
%   xq = query data
% OUTPUTS:
%   index = index of closest point in X to xq
%   dist = distance between X(index) and xq

indX = [1:size(X,1)]';

minX = min(X(:,1)); maxX = max(X(:,1));
minY = min(X(:,2)); maxY = max(X(:,2));

tree = [minX minY; maxX minY; maxX maxY; minX maxY];
branch = 1;
% Assume XN = [x1,x2,x3,x4,x12,x14,x23,x24,xC]
caseXN = [9 7 3 8; 5 2 7 9; 6 9 8 4; 1 5 9 6];
while true
    % Pull out 4 corner points
    x1 = tree(branch,:);
    x2 = tree(branch+1,:);
    x3 = tree(branch+2,:);
    x4 = tree(branch+3,:);
    % Compute halfway points
    x12 = 0.5*(x1+x2);
    x14 = 0.5*(x1+x4);
    x23 = 0.5*(x2+x3);
    x34 = 0.5*(x3+x4);
    xC = 0.5*(x1+x3);
    % Logical left/down testing
    Left = xq(1)<x12(1);
    Down = xq(2)<x14(2);
    % Determine what case we are dealing with
    XYN = [x1;x2;x3;x4;x12;x14;x23;x34;xC];
    Case = 2*Left + Down + 1;
    CaseXY = caseXN(Case,:);
    % Set new 4 corner points
    branch = branch + 4;
    tree(branch,:) = XYN(CaseXY(1),:);
    tree(branch+1,:) = XYN(CaseXY(2),:);
    tree(branch+2,:) = XYN(CaseXY(3),:);
    tree(branch+3,:) = XYN(CaseXY(4),:);
    % Find out how many data points are in the new domain
    logicX = (X(:,1)>tree(branch,1)) & (X(:,1)<tree(branch+2,1));
    logicY = (X(:,2)>tree(branch,2)) & (X(:,2)<tree(branch+2,2));
    logicXY = logicX & logicY;
    TOTAL = sum(logicXY);
    if TOTAL < 2
        branch = branch-4;
        logicX = (X(:,1)>=tree(branch,1)) & (X(:,1)<=tree(branch+2,1));
        logicY = (X(:,2)>=tree(branch,2)) & (X(:,2)<=tree(branch+2,2));
        logicXY = logicX & logicY;
        index = indX(logicXY);
        break;
    end
end
% Find minimum distance
IND = []; DIST = []; tic;
for i=1:512
    [index,tree] = QuadTreeSearch([xs(:) log(ysN(:)+1e-5)],[xs2(i,1) log(ys2N(i,1)+1e-5)]);
    dist = sqrt((xs2(i,1)-xs(index)).^2 + (log(ys2N(i,1)+1e-5)-log(ysN(index)+1e-5)).^2);
    [themin,minind] = min(dist);
    IND = [IND; index(minind)];
    DIST = [DIST; themin];
end
toc;

end

