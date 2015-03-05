function [pBL, pBR, pTR, pTL, surfaceFlag] = findCell(x,y,xq,yq,ind)
% Function to find the cell occupied by the query point
% 'ind' is a vector of indices of the nearest grid point neighbors of the
% query points

% NEW CODE (BARYCENTRIC SEARCH)
%{
MEANX = 0.25*(x(1:end-1,1:end-1)+x(2:end,1:end-1)+x(1:end-1,2:end)+x(2:end,2:end));
MEANY = 0.25*(y(1:end-1,1:end-1)+y(2:end,1:end-1)+y(1:end-1,2:end)+y(2:end,2:end));
NS = createns([MEANX(:),MEANY(:)]);
ind = knnsearch(NS,[xq,yq]);
%}

n = size(x,1); m = size(x,2);
total = length(x(:));

pCENT = [x(ind); y(ind)];
xLOC = [xq;yq] - pCENT;

surfaceFlag=0;
if ind-2*n <= 0 % Are we on the airfoil surface?
% YES...
    surfaceFlag=1;
    % Is the query point above or below the nearest grid point?
    UPvec = [x(ind+1); y(ind+1)] - pCENT;
    uptest = dot(xLOC,UPvec);
    if uptest<0
        % Below
        pBR = pCENT;
        pTR = [x(ind+1); y(ind+1)];
        pBL = [x(ind+n); y(ind+n)];
        pTL = [x(ind+n+1); y(ind+n+1)];
    else
        % Above
        pTR = pCENT;
        pBR = [x(ind-1); y(ind-1)];
        pTL = [x(ind+n); y(ind+n)];
        pBL = [x(ind+n-1); y(ind+n-1)];
    end
else
% NO...
if ind+n <= total % Points to the left?
    % YES...
    if mod(ind,n)==0 % No points above?
        % NO (pts above)...
        % Points to left, no points above
        LEFTvec = [x(ind+n); y(ind+n)] - pCENT;
        lefttest = dot(xLOC,LEFTvec);
        if lefttest > 0
            % Points to left, no points above, lefttest>0 --> pCENT is top right
            pTR = pCENT;
            pTL = [x(ind+n); y(ind+n)];
            pBR = [x(ind-1); y(ind-1)];
            pBL = [x(ind+n-1); y(ind+n-1)];
        else
            % Points to left, no points above, lefttest<0 --> pCENT is top left
            pTL = pCENT;
            pTR = [x(ind-n); y(ind-n)];
            pBR = [x(ind-n-1); y(ind-n-1)];
            pBL = [x(ind-1); y(ind-1)];
        end
    else
        % YES (pts above)...
        % Points to the left, points above
        LEFTvec = [x(ind+n); y(ind+n)] - pCENT;
        UPvec = [x(ind+1); y(ind+1)] - pCENT;
        lefttest = dot(xLOC,LEFTvec);
        uptest = dot(xLOC,UPvec);
        if (lefttest>0 && uptest>0)
            % Point left, point above, lefttest>0 & uptest>0 --> pCENT is bottom right
            pBR = pCENT;
            pTR = [x(ind+1); y(ind+1)];
            pBL = [x(ind+n); y(ind+n)];
            pTL = [x(ind+n+1); y(ind+n+1)];
        elseif (lefttest>0 && uptest<0)
            % Point left, point above, lefttest>0 & uptest<0 --> pCENT is top right
            pTR = pCENT;
            pTL = [x(ind+n); y(ind+n)];
            pBL = [x(ind+n-1); y(ind+n-1)];
            pBR = [x(ind-1); y(ind-1)];
        elseif (lefttest<0 && uptest<0)
            % Point left, point above, lefttest<0 & uptest<0 --> pCENT is top left
            pTL = pCENT;
            pBL = [x(ind-1); y(ind-1)];
            pTR = [x(ind-n); y(ind-n)];
            pBR = [x(ind-n-1); y(ind-n-1)];
        else
            % Point left, point above, lefttest<0 & uptest>0 --> pCENT is bottom left
            pBL = pCENT;
            pTL = [x(ind+1); y(ind+1)];
            pBR = [x(ind-n); y(ind-n)];
            pTR = [x(ind-n+1); y(ind-n+1)];
        end        
        
    end
else
    % NO...
    if mod(ind,n)==0 % No points above?
        % NO (pts above)...
        % No points to the left, no points above --> pCENT is top left
        pTL = [x(ind); y(ind)];
        pTR = [x(ind-n); y(ind-n)];
        pBL = [x(ind-1); y(ind-1)];
        pBR = [x(ind-n-1); y(ind-n-1)];
    else
        % YES (pts above)...
        % No points to the left, point above
        UPvec = [x(ind+1); y(ind+1)] - pCENT;
        uptest = dot(xLOC,UPvec);
        if uptest > 0
            % No points left, point above, uptest>0 --> pCENT is bottom left
            pBL = pCENT;
            pBR = [x(ind-n); y(ind-n)];
            pTL = [x(ind+1); y(ind+1)];
            pTR = [x(ind-n+1); y(ind-n+1)];
        else
            % No points left, point above, uptest<0 --> pCENT is top left
            pTL = pCENT;
            pTR = [x(ind-n); y(ind-n)];
            pBL = [x(ind-1); y(ind-1)];
            pBR = [x(ind-n-1); y(ind-n-1)];
        end
    end
end

end
