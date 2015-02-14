classdef InjectionDomain < hgsetget
    % Class defining the injection domain
    
    properties
        % Injection domain bounds and area
        XY_bounds;
        area;
        rhol;
        % Parameters for the number density PDF defined over the domain f(x,u,R,e)
        % Separate PDFs for xy,u,v,R,e (assuming these are all independent)
        fxy; % Space
        fu; fv; % Velocity
        fR; % Radius
        fe; % Energy
        % Total mass of domain
        mass;
        % Parameters for sampling f(x,u,R,e) over the domain
        samples;
        x = []; y = [];
        u = []; v = [];
        R = []; e = [];
        n = [];
    end
    
    methods
        
        function domain = InjectionDomain(dx,strPDFTypes,PDFparams,fluid,airfoil)
            % Constructor
            % INPUTS:
            %   dx: length scale for 'calcInjectionDomain'
            %   strPDFTypes: Vector of strings for 'setPDF'
            %   PDFparams: 2-column array of parameters for 'setPDF'
            
            % Set PDFs for u,v,R,e
            domain.fu = domain.setPDF(strPDFTypes(1),PDFparams(1,:));
            domain.fv = domain.setPDF(strPDFTypes(2),PDFparams(2,:));
            domain.fR = domain.setPDF(strPDFTypes(3),PDFparams(3,:));
            domain.fe = domain.setPDF(strPDFTypes(4),PDFparams(4,:));
            
            % Calculate bounds of injection domain
            if strcmp(strPDFTypes(3),'Gaussian')==1
                Ravg = PDFparams(3,1);
            elseif strcmp(strPDFTypes(3),'Uniform')==1
                Ravg = 0.5*(PDFparams(3,1)+PDFparams(3,2));
            end
            domain.calcInjectionDomain(dx,Ravg,fluid,airfoil);
            
            xy = domain.XY_bounds;
            xv = [xy(1:2,1); xy(4,1); xy(3,1)];
            yv = [xy(1:2,2); xy(4,2); xy(3,2)];
            % Calculate area of trapezoidal domain
            xREC = linspace(min(xy(:,1)),max(xy(:,1)),500)';
            yREC = linspace(min(xy(:,2)),max(xy(:,2)),500)';
            dx = xREC(2)-xREC(1);
            dy = yREC(2)-yREC(1);
            [XX,YY] = meshgrid(xREC,yREC);
            X = XX(:); Y = YY(:);
            area = sum(inpolygon(X,Y,xv,yv))*(dx*dy);
            domain.area = area;
            
            % Set spatial PDF as uniform in space over the trapezoidal domain
            domain.fxy = @(xq,yq) inpolygon(xq,yq,xv,yv)./area;
            % Set fluid density of droplets
            domain.rhol = fluid.rhol;
            
        end
        
        function calcInjectionDomain(domain,dx,Ravg,fluid,airfoil)
            % Function to calculate a box which contains the droplet clusters
            % INPUTS: 
            %   dx: length scale (desired domain length in the x-direction)
            %   Ravg: average droplet radius
            %   fluid,airfoil: fluid and airfoil objects

            % Determine impingement limits in y-direction at x0
            xL = -0.5; Yhit = -0.1;
            Ymiss = 0.3;
            yLimUP = impingementLimitsSLD(Ravg,fluid,airfoil,xL,Ymiss,Yhit,'UP');
            Ymiss = -0.3;
            yLimDOWN = impingementLimitsSLD(Ravg,fluid,airfoil,xL,Ymiss,Yhit,'DOWN');

            % Initialize cloud object of 2 particles placed at impingement limits
            x0 = [xL; xL];
            y0 = [yLimUP; yLimDOWN];
            [pg,ug,vg] = interpFluid(fluid,x0,y0);
            u0 = ug + 0.01*ug;
            v0 = vg + 0.01*vg;
            rd0 = [Ravg; Ravg];
            time0 = [0; 0];
            particles = 2;
            cloud = SLDcloud([x0 y0 u0 v0 rd0 time0 [1:particles]'],fluid.rhol,particles);

            % Advect particles backwards for time dx/U
            Uinf = fluid.Uinf;
            tStar = -dx/Uinf;
            while max(cloud.time) > tStar
                calcDtandImpinge(cloud,airfoil,fluid);
                set(cloud,'dt',-cloud.dt);
                transportSLD(cloud,fluid);
            end

            % Save new boundaries of particles
            domain.XY_bounds = [x0 y0; cloud.x cloud.y];
        end
        
        function pdf = setPDF(domain,strType,params)
            % Function to return PDF of a specified type
            % INPUTS:
            %   'strType' = string specifying PDF type
            %   'params' = PDF dependent shape parameters
            
            if strcmp(strType,'Gaussian')
                % 'params' = [mu sigma]
                mu = params(1);
                sigma = params(2);
                xsamp = linspace(mu-3*sigma,mu+3*sigma,1000)';
                func = 1/sigma/sqrt(2*pi)*exp(-0.5*(xsamp-mu).^2./(sigma^2));
            elseif strcmp(strType,'Uniform')
                % 'params' = [minx maxx]
                minx = params(1);
                maxx = params(2);
                normalize = maxx-minx;
                xsamp = linspace(minx,maxx,1000)';
                func = (heaviside(xsamp-minx) - heaviside(xsamp-maxx)).*(1/normalize);
            end
            pdf = [xsamp, func];
        end
        
        function integral = integratePDF(domain,strI,fluid)
            % Function to compute an integral over the domain, weighted by
            % the number density function
            % INPUTS:
            %   'strI': string defining type of integrand
            %   'fluid': object defining properties of fluid
                
            rhol = domain.rhol;
            if strcmp(strI,'Mass')
                R = domain.fR(:,1); fR = domain.fR(:,2);
                Integrand = 4/3*pi*rhol*(R.^3);
                I_mass = trapz(R,Integrand.*fR);
                u = domain.fu(:,1); fu = domain.fu(:,2); I_u = trapz(u,fu);
                v = domain.fv(:,1); fv = domain.fv(:,2); I_v = trapz(v,fv);
                e = domain.fe(:,1); fe = domain.fe(:,2); I_e = trapz(e,fe);
                integral = domain.NDOM.*I_mass.*I_u.*I_v.*I_e;
            end
            
        end
        
        
        
    end
    
end

