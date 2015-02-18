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
        ft; % Time
        % Ensemble average quantities of the domain
        LWC_avg;
        mass_avg;
        R_avg;
        numDroplets_avg;
        dtTraverse_avg;
        % Parameters for sampling f(x,u,R,e) over the domain
        samples;
        simTime;
        x = []; y = [];
        u = []; v = [];
        R = []; e = [];
        nDroplet = [];
        t = [];
    end
    
    methods
        
        function domain = InjectionDomain(strPDFTypes,PDFparams,fluid,airfoil,LWC_avg,simTime)
            % Constructor
            % INPUTS:
            %   strPDFTypes: Vector of strings for 'setPDF'
            %   PDFparams: 2-column array of parameters for 'setPDF'
            %   LWC_avg: ensemble average for LWC
            
            % Set PDFs for u,v,R,e
            domain.fu = domain.setPDF(strPDFTypes(1),PDFparams(1,:));
            domain.fv = domain.setPDF(strPDFTypes(2),PDFparams(2,:));
            domain.fR = domain.setPDF(strPDFTypes(3),PDFparams(3,:));
            domain.fe = domain.setPDF(strPDFTypes(4),PDFparams(4,:));
            domain.ft = domain.setPDF(strPDFTypes(5),PDFparams(5,:));
            domain.LWC_avg = LWC_avg;
            domain.simTime = simTime;
            % Set fluid density of droplets
            domain.rhol = fluid.rhol;
            % Calculate bounds of injection domain and avg time to traverse
            if strcmp(strPDFTypes(3),'Gaussian')==1
                Ravg = PDFparams(3,1);
            elseif strcmp(strPDFTypes(3),'Uniform')==1
                Ravg = 0.5*(PDFparams(3,1)+PDFparams(3,2));
            end
            domain.R_avg = Ravg;
            domain.calcInjectionDomain(fluid,airfoil);
            % Calculate area of trapezoidal injection domain
            domain.calcDomainArea();
            % Calculate the average droplet mass from the distribution f(R)
            R = domain.fR(:,1); fR = domain.fR(:,2);
            Integrand = 4/3*pi*domain.rhol*(R.^3);
            md_avg = trapz(R,Integrand.*fR);
            % Calculate average number of total droplets in the simulation
            % (assumes differential thickness of domain = one MVD)
            ThicknessSlice = 2*domain.R_avg;
            VolSlice = ThicknessSlice*domain.area;
            I_t = simTime/domain.dtTraverse_avg;
            domain.numDroplets_avg = I_t*LWC_avg*VolSlice/md_avg;
            % Set spatial PDF as uniform in space over the trapezoidal
            % domain, scaled by the average number of particles
            domain.fxy = @(xq,yq) domain.num_avg*inpolygon(xq,yq,xv,yv)./domain.area;
        end
        
        function sampleRealization(numClumps)
            % Function to create a realization based on desired number of
            % particle clumps and simulation time
            % INPUTS:
            %   numClumps: desired number of particle clumps
            %   simTime: desired time of simulation
            
            % Code to draw random samples from (x,u,R,e;t)
            % ...
            % Code to evaluate f(x,u,R,e;t) at sample points
            % ...
            
            dV = domain.numDroplets_avg/sum(fEval.*mEval);
            domain.nDroplet = fEval.*dV;
            
        end
        
        function calcInjectionDomain(domain,fluid,airfoil)
            % Function to calculate a box which contains the droplet clusters
            % INPUTS: 
            %   fluid,airfoil: fluid and airfoil objects
            
            Ravg = domain.R_avg;
            % Determine impingement limits in y-direction at x0
            xL = -0.5; Yhit = -0.1;
            Ymiss = 0.3;
            yLimUP = impingementLimitsSLD(Ravg,fluid,airfoil,xL,Ymiss,Yhit,'UP');
            Ymiss = -0.3;
            yLimDOWN = impingementLimitsSLD(Ravg,fluid,airfoil,xL,Ymiss,Yhit,'DOWN');
            x1 = [xL; yLimUP];
            x2 = [xL; yLimDOWN];
            % Determine impingement limits in y-direction at x = x0-1
            xL = xL-1; Yhit = -0.25;
            Ymiss = 0.3;
            yLimUP = impingementLimitsSLD(Ravg,fluid,airfoil,xL,Ymiss,Yhit,'UP');
            Ymiss = -0.35;
            yLimDOWN = impingementLimitsSLD(Ravg,fluid,airfoil,xL,Ymiss,Yhit,'DOWN');
            x3 = [xL; yLimUP];
            x4 = [xL; yLimDOWN];
            % Save new boundaries of particles
            domain.XY_bounds = [x1'; x2'; x3'; x4'];
            % Estimate time to traverse domain
            domain.dtTraverse_avg = 1/fluid.Uinf;
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
                xsamp = linspace(mu-5*sigma,mu+5*sigma,3000)';
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
        
        function calcDomainArea(domain)
            % Function to compute the area of the injection domain
            
            xy = domain.XY_bounds;
            xv = [xy(1:2,1); xy(4,1); xy(3,1)];
            yv = [xy(1:2,2); xy(4,2); xy(3,2)];
            xREC = linspace(min(xy(:,1)),max(xy(:,1)),500)';
            yREC = linspace(min(xy(:,2)),max(xy(:,2)),500)';
            dx = xREC(2)-xREC(1);
            dy = yREC(2)-yREC(1);
            [XX,YY] = meshgrid(xREC,yREC);
            X = XX(:); Y = YY(:);
            area = sum(inpolygon(X,Y,xv,yv))*(dx*dy);
            domain.area = area;
            
        end
        
        function calcAvgMass(domain)
            % Function to compute the average mass of the injection domain
            % INPUTS:
            %   'strI': string defining type of integrand
            
            rhol = domain.rhol;
            I_xy = domain.num_avg;
            I_u = 1; I_v = 1; I_e = 1;
            R = domain.fR(:,1); fR = domain.fR(:,2);
            Integrand = 4/3*pi*rhol*(R.^3);
            I_mass = trapz(R,Integrand.*fR);
            domain.mass_avg = I_mass.*I_xy.*I_u.*I_v.*I_e;
            
        end
        
        
        
    end
    
end

