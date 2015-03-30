classdef Airfoil < hgsetget
    % Defines an airfoil and associated operations
    
    properties
        X; Y; % Spline points of the airfoil (provided at initialization)
        PANELx; PANELy; % (x,y) locations of panel center points (for normal vectors)
        TANG; NORM; % Tangent and normal vectors
        s; % s-coordinate locations of panel center points
        FILMsplash = []; % Height of surface water film accumulation in s-coordinates
        FILMspread = [];
        FILM = [];
        LIMup; LIMdown; % Upper/lower limits of impingement
        interpXYtoS; % Scattered interpolant for x,y to s-coordinates
        originalImpingeScoord = []; % S coordinates of parent particles impinged
        originalImpingeScoordSplash = [];
        originalImpingeScoordSpread = [];
        stagPt; % Stagnation point of the airfoil
        % Tree-searcher object for panel centroids
        NS;
    end
    
    methods
        function airfoil = Airfoil(vars)
            % Constructs an airfoil by feeding in vars = [x, y]
            airfoil.X = vars(:,1);
            airfoil.Y = vars(:,2);
            % Calculate tangent/normal vectors at panel center points
            airfoil.setTangNorm();
            % Initialize s-coord mapping
            airfoil.calculateSCoords();
            % Initialize tree-searcher for panel center points
            airfoil.NS = createns([airfoil.PANELx, airfoil.PANELy]);
        end
        
        function setTangNorm(airfoil)
            % Function to compute tangent and normal vectors of airfoil
            
            x = airfoil.X; y = airfoil.Y;
            
            ds_x = diff(x); ds_y = diff(y);
            panelx = x(1:end-1) + 0.5*ds_x;
            panely = y(1:end-1) + 0.5*ds_y;
            ds = [ds_x ds_y];
            tangvec = ds./repmat(sqrt(ds_x.^2 + ds_y.^2),1,2);
            normvec = -[ds_y -ds_x]./repmat(sqrt(ds_x.^2 + ds_y.^2),1,2);
            
            airfoil.PANELx = panelx;
            airfoil.PANELy = panely;
            airfoil.TANG = tangvec;
            airfoil.NORM = normvec;
            
        end
        
        function [px,py,nx,ny,tx,ty] = findPanel(airfoil,xq,yq)
            % Function to determine the closest panel to a query point
            % Returns closest panel point and normal vector
            
            normvec = airfoil.NORM; tangvec = airfoil.TANG;
            % Tree-search for nearest panel centroid
            ind = knnsearch(airfoil.NS,[xq,yq]);
            px = airfoil.PANELx(ind);
            py = airfoil.PANELy(ind);
            % Tangent and normal vectors
            tx = tangvec(ind,1);
            ty = tangvec(ind,2);
            nx = normvec(ind,1);
            ny = normvec(ind,2);
            
            %{
            x = airfoil.PANELx; y = airfoil.PANELy;
            normvec = airfoil.NORM; tangvec = airfoil.TANG;
            
            for i=1:size(xq,1)
                repXY = repmat([xq(i),yq(i)],length(x),1);
                L2norm = sqrt((repXY(:,1)-x).^2 + (repXY(:,2)-y).^2);
                [val,index] = min(L2norm);
                px(i,1) = x(index(1)); py(i,1) = y(index(1));
                nx(i,1) = normvec(index(1),1); ny(i,1) = normvec(index(1),2);
                tx(i,1) = tangvec(index(1),1); ty(i,1) = tangvec(index(1),2);
            end
            %}
        end
        
        function calculateSCoords(airfoil)
            % Calculate s-coordinates of panel center points
            
            x = airfoil.PANELx; y = airfoil.PANELy;
            s = zeros(length(x),1);
            DX = diff(x); DY = diff(y);
            DS = sqrt(DX.^2 + DY.^2);
            s(2:end) = cumsum(DS);
            
            airfoil.s = s;
            airfoil.interpXYtoS = scatteredInterpolant([x,y],s);
            
        end
        
        function TH = findTH(airfoil,xq,yq,uq,vq)
            % Function to determine the angle of incidence of a droplet
            % impinging on the airfoil surface
            
            % Find closest panel
            [px,py,nx,ny,~,~] = airfoil.findPanel(xq,yq);
            % Find angle between normal vectors and query velocities
            VelNorm = sqrt(uq.^2 + vq.^2);
            Vel = [uq vq]./repmat(VelNorm,1,2);
            proj = sum(Vel.*[nx,ny],2);
            TH = acos(proj);
            
        end
        
        function airfoil = set.FILMsplash(airfoil,vars)
            if ~isempty(vars)
                % Set particular elements of film
                ind1 = size(airfoil.FILMsplash,1)+1;
                ind2 = ind1+size(vars,1)-1;
                airfoil.FILMsplash(ind1:ind2,1) = vars(:,1);
                airfoil.FILMsplash(ind1:ind2,2) = vars(:,2);
            else
                airfoil.FILMsplash = [];
            end
        end
        
        function airfoil = set.FILMspread(airfoil,vars)
            % Set particular elements of film
            airfoil.FILMspread = vars;
        end
        
        function airfoil = set.FILM(airfoil,vars)
            % Set particular elements of film
            airfoil.FILM = vars;
        end
        
        function airfoil = set.originalImpingeScoordSplash(airfoil,vars)
            if ~isempty(vars)
                % Set particular elements of film
                ind1 = size(airfoil.originalImpingeScoordSplash,1)+1;
                ind2 = ind1+size(vars,1)-1;
                airfoil.originalImpingeScoordSplash(ind1:ind2,1) = vars;
            else
                airfoil.originalImpingeScoordSplash = [];
            end
        end
        
        function airfoil = set.originalImpingeScoordSpread(airfoil,vars)
            % Set particular elements of film
            airfoil.originalImpingeScoordSpread = vars;
        end
        
        function airfoil = set.originalImpingeScoord(airfoil,vars)
            % Set particular elements of film
            airfoil.originalImpingeScoord = vars;
        end
        
        function ind = findPanelIndex(airfoil,px,py)
            % Function to find index of a panel point
            
            x = airfoil.PANELx; y = airfoil.PANELy;
            ind = find(x==px & y==py);
        end
        
        function setLim(airfoil,xq,yq,str)
            % Function to set either upper or lower limit of impingement
            
            if strcmp(str,'UP')
                airfoil.LIMup = airfoil.XYtoScoords(xq,yq);
            elseif strcmp(str,'DOWN')
                airfoil.LIMdown = airfoil.XYtoScoords(xq,yq);
            end
        end
        
        function s = XYtoScoords(airfoil,xq,yq)
            % Function to convert (x,y) surface coordinates to s-coordinates
            
            s = airfoil.interpXYtoS([xq,yq]);
            
        end
        
        function [xq,yq] = interpStoXY(airfoil,sq)
            % Function to convert s-coordinates to x,y
            
            x = airfoil.PANELx; y = airfoil.PANELy;
            s = airfoil.s;
            
            xq = spline(s,x,sq);
            yq = spline(s,y,sq);
        end
        
        function calcStagPt(airfoil,fluid)
            % Calculate the stagnation point of the airfoil
            
            % Find stagnation pt
            I = size(fluid.x,1);
            indWrap1 = I+1; indWrap2 = I*2;
            RHOU = fluid.RHOU(indWrap1:indWrap2); RHOV = fluid.RHOV(indWrap1:indWrap2);
            [themin,theind] = min(RHOU(1:floor(I/2)).^2 + RHOV(1:floor(I/2)).^2);
            airfoil.stagPt = airfoil.XYtoScoords(fluid.x(I+theind),fluid.y(I+theind));
            
        end
        
        function clearFilm(airfoil)
            % Function to clear airfoil film
            
            set(airfoil,'FILM',[]);
            set(airfoil,'FILMsplash',[]);
            set(airfoil,'FILMspread',[]);
            set(airfoil,'originalImpingeScoord',[]);
            set(airfoil,'originalImpingeScoordSplash',[]);
            set(airfoil,'originalImpingeScoordSpread',[]);
        end
        
    end
    
end

