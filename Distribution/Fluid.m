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
            fluid.NS = createns([x(:),y(:)]);
            fluid.RHO = scatteredInterpolant(x(:),y(:),rho(:));
            fluid.RHOU = scatteredInterpolant(x(:),y(:),rhou(:));
            fluid.RHOV = scatteredInterpolant(x(:),y(:),rhov(:));
            fluid.E = scatteredInterpolant(x(:),y(:),e(:));
            fluid.mach = mach;
            fluid.alpha = alpha;
            fluid.Re = Re;
            fluid.Uinf = mach*340;
            
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
    end
    
end

