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
        % Mesh x,y coordinates
        x; y;
        % Tree-search object for nearest neighbor search of the grid
        NS;
        % Cell areas
        cellarea;
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
            fluid.RHO = scatteredInterpolant(x(:),y(:),rho(:));
            fluid.RHOU = scatteredInterpolant(x(:),y(:),rhou(:));
            fluid.RHOV = scatteredInterpolant(x(:),y(:),rhov(:));
            fluid.E = scatteredInterpolant(x(:),y(:),e(:));
            fluid.mach = mach;
            fluid.alpha = alpha;
            fluid.Re = Re;
            fluid.Uinf = mach*340;
            % Create tree-search object of cell centroids
            MEANX = 0.25*(x(1:end-1,1:end-1)+x(2:end,1:end-1)+x(1:end-1,2:end)+x(2:end,2:end));
            MEANY = 0.25*(y(1:end-1,1:end-1)+y(2:end,1:end-1)+y(1:end-1,2:end)+y(2:end,2:end));
            fluid.NS = createns([MEANX(:),MEANY(:)]);
            % Calculate cell areas
            fluid.computeCellAreas();            
        end
        
        function [pg,ug,vg] = interpFluid(fluid,xq,yq)
            % Function to interpolate fluid at query points
            
            rhoinf = fluid.rhoinf;
            Ubar = fluid.Ubar;
            RHO = fluid.RHO;
            RHOU = fluid.RHOU;
            RHOV = fluid.RHOV;
            
            pg = rhoinf*RHO(xq,yq);
            ug = RHOU(xq,yq)*Ubar*rhoinf./pg;
            vg = RHOV(xq,yq)*Ubar*rhoinf./pg;
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
            %CENT = [0 0];
            x = fluid.x - CENT(1); y = fluid.y - CENT(2);
            MEANX = 0.25*(x(1:end-1,1:end-1)+x(2:end,1:end-1)+x(1:end-1,2:end)+x(2:end,2:end));
            MEANY = 0.25*(y(1:end-1,1:end-1)+y(2:end,1:end-1)+y(1:end-1,2:end)+y(2:end,2:end));
            x2 = MEANX; y2 = MEANY;
            COMP = x + i*y;
            SQRT = sqrt(COMP);
            xSQRT = real(SQRT); ySQRT = imag(SQRT);
            scal=1;
        
                  
                  bx=0;
                  for j=1:129
                      u=1; v=0; angl= pi + pi;
                  for k = 1:513
                        angl      = angl  +atan2((u*y(k,j)  -v*x(k,j)),(u*x(k,j)  +v*y(k,j)));
                        r         = scal*sqrt(x(k,j)^2  +y(k,j)^2);
                        u         = x(k,j); 
                        v         = y(k,j);
                        r         = bx*r  +sqrt((bx*r)^2  +2.*r);
                        xs(k,j)     = r*cos(.5*angl);
                        ys(k,j)     = r*sin(.5*angl);
                      end
                  end
       
                                    
                  for j=1:128
                  angl2 = pi + pi;
                  u2=1; v2=0;
                  for k=1:512
                      angl2     = angl2  +atan2((u2*y2(k,j)  -v2*x2(k,j)),(u2*x2(k,j)  +v2*y2(k,j)));
                      r2         = scal*sqrt(x2(k,j)^2  +y2(k,j)^2);
                      u2         = x2(k,j);
                      v2         = y2(k,j);
                      r2         = bx*r2  +sqrt((bx*r2)^2  +2.*r2);
                      xs2(k,j)     = r2*cos(.5*angl2);
                      ys2(k,j)     = r2*sin(.5*angl2);
                  end
                  end


            
        end
    end
    
end

