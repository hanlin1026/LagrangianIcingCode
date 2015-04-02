classdef Fluid < hgsetget
    % Defines a fluid class for holding all properties and interpolating
    % functions of a fluid
    
    properties
        % Scalar flow quantities
        pinf;
        R;
        Tinf;
        rhoinf;
        Ubar;
        Uinf;
        rhol;
        % Scattered data interpolating functions
        RHO;
        RHOU;
        RHOV;
        E;
        % Mesh and solution filenames
        meshfile;
        solnfile;
        % Mesh coordinates
        x; y;
        % Tree-search object for nearest neighbor search of the grid
        NS;
        refinedWrapInd;
        wrapRefinement;
        indexNN;
        % Cell stuff
        cellarea;
        MEANx; MEANy; % Centers
        xx; xy; yx; yy; % Grid metrics
        rhoC; rhouC; rhovC; eC; % Flow variables at cell centers
        Lmin; % Minimum cell lengths
        % SQRT plane mapping of grid corners
        xSQRT; ySQRT;
        % Other solution quantities
        mach;
        alpha;
        Re;
        
    end
    
    methods
        function fluid = Fluid(scalars,mesh,soln)
            % Constructor for fluid
            % scalars = [pinf,R,Tinf,rhoinf,Uinf,rhol,Rd,Tinf]
            
            fluid.pinf =   scalars(1);
            fluid.R =      scalars(2);
            fluid.Tinf =   scalars(3);
            fluid.rhoinf = scalars(4);
            fluid.Ubar =   scalars(5);
            fluid.rhol =   scalars(6);
            
            fluid.meshfile = mesh;
            fluid.solnfile = soln;
            [x,y,rho,rhou,rhov,e,mach,alpha,Re,~] = readp3d(mesh,soln);
            fluid.x = x;
            fluid.y = y;
            fluid.RHO = rho(:);
            fluid.RHOU = rhou(:);
            fluid.RHOV = rhov(:);
            fluid.E = e(:);
            fluid.mach = mach;
            fluid.alpha = alpha;
            fluid.Re = Re;
            fluid.Uinf = mach*340;
            % Calculate cell properties
            fluid.computeCellAreas();      
            fluid.computeCellCenters();
            fluid.computeGridMetrics();
            % Create tree-search object
            fluid.createTreeSearcher();
        end
        
        function [pg,ug,vg] = getNNFluidProps(fluid,xq,yq)
            % Function to return fluid values at nearest neighbor cell
            % center pts
            
            rhoinf = fluid.rhoinf;
            Ubar = fluid.Ubar;
            RHO = fluid.rhoC;
            RHOU = fluid.rhouC;
            RHOV = fluid.rhovC;
            
            % Use values at closest grid points
            ind = fluid.searchTree([xq,yq]);
            pg = rhoinf*RHO(ind);
            ug = RHOU(ind)*Ubar*rhoinf./pg;
            vg = RHOV(ind)*Ubar*rhoinf./pg;
        end
        
        function computeCellAreas(fluid)
            % Function to compute cell areas
            
            x = fluid.x; y = fluid.y;
            DI_x = x(2:end,1:end-1)-x(1:end-1,1:end-1);
            DI_y = y(2:end,1:end-1)-y(1:end-1,1:end-1);
            DI = sqrt(DI_x(:).^2 + DI_y(:).^2);
            DJ_x = x(1:end-1,2:end)-x(1:end-1,1:end-1);
            DJ_y = y(1:end-1,2:end)-y(1:end-1,1:end-1);
            DJ = sqrt(DJ_x(:).^2 + DJ_y(:).^2);
            fluid.cellarea = DI.*DJ;
        end

        function sqrtMap(fluid)
            % Function to compute sqrt mapping of grid coordinates
            
            CENT = [4.885655747050746E-003  6.112252758408689E-003];
            x = fluid.x - CENT(1); y = fluid.y - CENT(2);
            MEANX = 0.25*(x(1:end-1,1:end-1)+x(2:end,1:end-1)+x(1:end-1,2:end)+x(2:end,2:end));
            MEANY = 0.25*(y(1:end-1,1:end-1)+y(2:end,1:end-1)+y(1:end-1,2:end)+y(2:end,2:end));
            x2 = MEANX; y2 = MEANY;
            scal=1;
            bx=0;
            for j=1:129
                u=1; v=0; angl= pi + pi;
             	for k = 1:513
                    angl      = angl + atan2((u*y(k,j)-v*x(k,j)),(u*x(k,j)+v*y(k,j)));
                    r         = scal*sqrt(x(k,j)^2  +y(k,j)^2);
                    u         = x(k,j); 
                    v         = y(k,j);
                    r         = bx*r + sqrt((bx*r)^2  +2.*r);
                    xs(k,j)     = r*cos(.5*angl);
                    ys(k,j)     = r*sin(.5*angl);
                end
            end
            
            for j=1:128
                angl2 = pi + pi;
                u2=1; v2=0;
                for k=1:512
                    angl2     = angl2  +atan2((u2*y2(k,j)  -v2*x2(k,j)),(u2*x2(k,j)  +v2*y2(k,j)));
                    r2         = scal*sqrt(x2(k,j)^2 + y2(k,j)^2);
                    u2         = x2(k,j);
                    v2         = y2(k,j);
                    r2         = bx*r2  +sqrt((bx*r2)^2  +2.*r2);
                    xs2(k,j)     = r2*cos(.5*angl2);
                    ys2(k,j)     = r2*sin(.5*angl2);
                end
            end
            % Subtract off the airfoil
            ys = ys - repmat(ys(:,1),1,129);
            ys2 = ys2 - repmat(ys2(:,1),1,128);
            % Take logarithm of y-scale
            %ys = atan(ys + 1e-5);
            %ys2 = atan(ys2 + 1e-5);
            
        end
        
        function createTreeSearcher(fluid)
            % Function to create tree-searcher object on refined grid
            
            MEANx = fluid.MEANx; MEANy = fluid.MEANy;
            fluid.NS = createns([MEANx(:),MEANy(:)]);
            %{
            x = fluid.x; y = fluid.y;
            I = size(x,1); J = size(x,2);
            % Calculate aspect ratio of each wrap
            AR = [];
            for i=1:J-1
                ind = (i-1)*I + floor(I/2);
                lJ = norm([x(ind),y(ind)]-[x(ind+I),y(ind+I)]);
                lI = norm([x(ind),y(ind)]-[x(ind+1),y(ind+1)]);
                AR(i) = lI/lJ;
            end
            RES = floor(AR/1);
            RES(RES <= 2) = 2;
            % Refine grid
            XNEW = [];
            tmp = find(RES>2);
            lastWrap = tmp(end);
            for j=1:lastWrap
                for i=1:I-1
                    XNEW = [XNEW; [linspace(x(i,j),x(i+1,j),RES(j))' linspace(y(i,j),y(i+1,j),RES(j))']];
                end
            end
            XNEW = [XNEW; [x(1+lastWrap*size(x,1):end)',y(1+lastWrap*size(x,1):end)']];
            % Create tree-searcher
            fluid.NS = createns(XNEW);
            % Create look-up table to convert refined indices to original
            % grid indices
            convertRefOrigInd = [];
            for i=1:lastWrap
                if i==1
                    convertRefOrigInd(i,1) = 1;
                else
                    convertRefOrigInd(i,1) = convertRefOrigInd(i-1,2)+1;
                end
                convertRefOrigInd(i,2) = convertRefOrigInd(i,1) + (RES(i)*(I-1)-1);
            end
            for i=lastWrap+1:J
                convertRefOrigInd(i,1) = convertRefOrigInd(i-1,2)+1;
                convertRefOrigInd(i,2) = convertRefOrigInd(i,1) + I-1;
            end
            fluid.refinedWrapInd = convertRefOrigInd;
            RES(RES==2) = 1;
            fluid.wrapRefinement = RES;
            %}
        end
        
        function index = searchTree(fluid,xq)
            % Function to search the tree object
            
            % Search the refined grid
            index = knnsearch(fluid.NS,xq);
            % Convert refined grid index to original grid index
            %{
            I = size(fluid.x,1);
            RES = fluid.wrapRefinement;
            refinedWrapInd = fluid.refinedWrapInd;
            index = zeros(size(xq,1),1);
            for i=1:size(xq,1)
                ind1 = ind(i) >= refinedWrapInd(:,1);
                ind2 = ind(i) <= refinedWrapInd(:,2);
                wrap(i,1) = find(ind1&ind2 == 1);
                cell(i,1) = ceil((ind(i)-refinedWrapInd(wrap(i),1))/RES(wrap(i)));
                index(i) = (wrap(i)-1)*I + cell(i);
            end
            %}
        end
        
        function [pg,ug,vg] = interpFluid(fluid,xq,yq)
            % Function to interpolate fluid at (xq,yq)
            
            % Implement by simply finding nearest neighbor mesh pt
            [pg,ug,vg] = fluid.getNNFluidProps(xq,yq);
            
        end
        
        function area = getCellArea(fluid,indGridPt)
            % Function to return cell area corresponding to closest grid pt
            
            x = fluid.x; y = fluid.y;
            I = size(x,1); J = size(x,2);
            % Assume that closest grid pt is the bottom left of cell in
            % which particle is located
            indWrap = ceil(indGridPt./I);
            indCell = indGridPt - (indWrap-1).*I;
            indCenter = (indWrap-1).*(I-1) + indCell;
            area = fluid.cellarea(indCenter);
        end
        
        function computeCellCenters(fluid)
            % Function to compute cell centers
            
            x = fluid.x; y = fluid.y;
            nI = size(x,1); nJ = size(x,2);
            fluid.MEANx = 0.25*(x(1:end-1,1:end-1)+x(2:end,1:end-1)+x(1:end-1,2:end)+x(2:end,2:end));
            fluid.MEANy = 0.25*(y(1:end-1,1:end-1)+y(2:end,1:end-1)+y(1:end-1,2:end)+y(2:end,2:end));
            % Compute flow variables at cell centers
            R = reshape(fluid.RHO,nI,nJ); RU = reshape(fluid.RHOU,nI,nJ); 
            RV = reshape(fluid.RHOV,nI,nJ); E = reshape(fluid.E,nI,nJ);
            I = 1:(nI-1); J = 1:(nJ-1);
            fluid.rhoC = 0.25*(R(I,J)+R(I+1,J)+R(I,J+1)+R(I+1,J+1));
            fluid.rhouC = 0.25*(RU(I,J)+RU(I+1,J)+RU(I,J+1)+RU(I+1,J+1));
            fluid.rhovC = 0.25*(RV(I,J)+RV(I+1,J)+RV(I,J+1)+RV(I+1,J+1));
            fluid.eC = 0.25*(E(I,J)+E(I+1,J)+E(I,J+1)+E(I+1,J+1));
        end
        
        function computeGridMetrics(fluid)
            % Function to calculate grid metrics
            
            x = fluid.x; y = fluid.y;
            xC = fluid.MEANx; yC = fluid.MEANy;
            % Calculate Jacobian metric tensor components at cell centers
            I = 1:(size(x,1)-1); J = 1:(size(x,2)-1);
            fluid.xx = 0.5*((x(I+1,J+1)-x(I,J+1)) + (x(I+1,J)-x(I,J)));
            fluid.xy = 0.5*((x(I,J+1)-x(I,J)) + (x(I+1,J+1)-x(I+1,J)));
            fluid.yx = 0.5*((y(I+1,J+1)-y(I,J+1)) + (y(I+1,J)-y(I,J)));
            fluid.yy = 0.5*((y(I,J+1)-y(I,J)) + (y(I+1,J+1)-y(I+1,J)));
            % Calculate minimum cell lengths (for use in CFL condition)
            LI = sqrt(fluid.xx.^2 + fluid.yx.^2);
            LJ = sqrt(fluid.xy.^2 + fluid.yy.^2);
            fluid.Lmin = min(LI,LJ);
            
        end
        
        function [I,J] = transformXYtoIJ(fluid,IND,xq)
            % Function to transform global (x,y) coords to local (I,J)
            % computational domain coords (centered at cell center IND)
            % NOTE: this routine is vectorized with the assumption that
            % element xq(k) is being transformed using Jacobian(IND(k))
            
            x = xq(:,1); y = xq(:,2);
            XC = fluid.MEANx(IND); YC = fluid.MEANy(IND);
            A = fluid.cellarea(IND);
            
            xx = fluid.xx(IND); xy = fluid.xy(IND);
            yx = fluid.yx(IND); yy = fluid.yy(IND);
            
            % Inverse transformation
            X = x-XC; Y = y-YC;
            I = (X.*yy - Y.*xy)./A;
            J = (-X.*yx + Y.*xx)./A;
            
        end
    end
    
end

