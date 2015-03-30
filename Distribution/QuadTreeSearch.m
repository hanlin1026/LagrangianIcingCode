function [index,tree] = QuadTreeSearch(X,xq,BB)
% Quadtree binary search algorithm for 2D data
% INPUTS:
%   X = reference data (rows = observations, columns = coordinates)
%   xq = query data
% OUTPUTS:
%   index = index of closest point in X to xq
%   dist = distance between X(index) and xq

indX = [1:size(X,1)]';

minX = BB(1); maxX = BB(2);
minY = BB(3); maxY = BB(4);

tree = [minX minY; maxX minY; maxX maxY; minX maxY];
branch = 1;
% Assume XN = [x1,x2,x3,x4,x12,x14,x23,x24,xC]
caseXN = [9 7 3 8; 5 2 7 9; 6 9 8 4; 1 5 9 6];
XDOM = X;
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
    % Save old logic
    XDOMprev = XDOM;
    indXprev = indX;
    % Find out how many data points are in the new domain
    %logicX = (XDOM(:,1)>=tree(branch,1)) & (XDOM(:,1)<=tree(branch+2,1));
    %logicY = (XDOM(:,2)>=tree(branch,2)) & (XDOM(:,2)<=tree(branch+2,2));
    logicXY = (XDOM(:,1)>=tree(branch,1)) & (XDOM(:,1)<=tree(branch+2,1)) & (XDOM(:,2)>=tree(branch,2)) & (XDOM(:,2)<=tree(branch+2,2));
    XDOM = XDOM(logicXY,:);
    indX = indX(logicXY);
    TOTAL = sum(logicXY);
    % If number of data points inside bucket is less than threshold, roll 
    % back to the previous bucket
    if TOTAL < 50
        ind1 = indXprev;
        XDOM = XDOMprev;
        break;
    end
end
% Transform coordinates by rotation corresponding to local tang/normal
%{
cent = X(min(ind1),:);
cent_IP1 = X(min(ind1)+1,:);
t = (cent_IP1-cent)/norm(cent_IP1-cent);
n = [-t(2); t(1)];
%}
%a = t(1)^2 + n(1)^2;
%b = t(1)*t(2) + n(1)*n(2);
%c = t(2)^2 + n(2)^2;
% Find minimum distance (1-norm) in transformed coordinates
%dist = abs(xq(1)*t(1) + xq(2)*t(2)) + abs(xq(1)*n(1) + xq(2)*n(2));
Xbin = X(ind1,:);
dist = sqrt((xq(1)-Xbin(:,1)).^2 + (xq(2)-Xbin(:,2)).^2);
[themin,theind] = min(dist);
index = ind1(theind);
%[~,ind2] = min(dist);
%index = ind1(ind2);
end

