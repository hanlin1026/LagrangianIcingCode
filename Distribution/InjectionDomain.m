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
        nDroplet = [];
        numParcels;
    end
    
    methods
        
        function domain = InjectionDomain(strPDFTypes,PDFparams,fluid,airfoil,LWC_avg,simTime)
            % Constructor
            % INPUTS:
            %   strPDFTypes: Vector of strings for 'setPDF'
            %   PDFparams: cell of parameters for 'setPDF'
            %   LWC_avg: ensemble average for LWC
            
            % Set PDFs for u,v,R,e
            if (strcmp(strPDFTypes(1),'Implicit')==0 || strcmp(strPDFTypes(2),'Implicit')==0)
                domain.fu = domain.setPDF(strPDFTypes(1),PDFparams{1});
                domain.fv = domain.setPDF(strPDFTypes(2),PDFparams{2});
            end
            domain.fR = domain.setPDF(strPDFTypes(3),PDFparams{3});
            domain.fe = domain.setPDF(strPDFTypes(4),PDFparams{4});
            domain.ft = domain.setPDF(strPDFTypes(5),PDFparams{5});
            domain.LWC_avg = LWC_avg;
            domain.simTime = simTime;
            % Set fluid density of droplets
            domain.rhol = fluid.rhol;
            % Calculate bounds of injection domain and avg time to traverse
            if strcmp(strPDFTypes(3),'Gaussian')==1
                Ravg = PDFparams{3}(1);
            elseif strcmp(strPDFTypes(3),'Uniform')==1
                Ravg = 0.5*(PDFparams{3}(1)+PDFparams{3}(2));
            elseif strcmp(strPDFTypes(3),'Custom')==1
                Ravg = trapz(domain.fR(:,1),domain.fR(:,1).*domain.fR(:,2));
            end
            domain.R_avg = Ravg;
            %
            domain.calcInjectionDomain(fluid,airfoil);
            % Calculate area of trapezoidal injection domain
            domain.calcDomainArea();
            % Calculate the average droplet mass from the distribution f(R)
            R = domain.fR(:,1); fR = domain.fR(:,2);
            Integrand = 4/3*pi*domain.rhol*(R.^3);
            md_avg = trapz(R,Integrand.*fR);
            % Calculate average number of droplets in the injection domain
            % (assumes differential thickness of domain = one MVD)
            ThicknessSlice = 2*domain.R_avg;
            VolSlice = ThicknessSlice*domain.area;
            domain.numDroplets_avg = LWC_avg*VolSlice/md_avg;
            % Set spatial PDF as uniform in space over the trapezoidal
            % domain, scaled by the average number of particles
            xy = domain.XY_bounds;
            xv = [xy(1:2,1); xy(4,1); xy(3,1)];
            yv = [xy(1:2,2); xy(4,2); xy(3,2)];
            domain.fxy = @(xq,yq) domain.numDroplets_avg*inpolygon(xq,yq,xv,yv)./domain.area;
            % Scale temporal PDF s.t. integral is tSim/dtTraverse
            domain.ft(:,2) = domain.ft(:,2)*simTime/domain.dtTraverse_avg;
            %}
        end
        
        function sampleRealization(domain,numClumps,fluid)
            % Function to create a realization based on desired number of
            % particle clumps and simulation time
            % INPUTS:
            %   numClumps: desired number of particle clumps
            %   simTime: desired time of simulation
            
            % Code to draw random samples from (x,u,R,e;t)
            minx = min(domain.XY_bounds(:,1)); maxx = max(domain.XY_bounds(:,1));
            miny = min(domain.XY_bounds(:,2)); maxy = max(domain.XY_bounds(:,2));
            samplesXY = []; fxyEval = [];
            % Sample xy domain
            while size(samplesXY,1) < numClumps
                xy = unifrnd([minx miny],[maxx maxy]);
                xyEval = domain.fxy(xy(1),xy(2));
                if xyEval ~= 0
                    fxyEval = [fxyEval; xyEval];
                    samplesXY = [samplesXY; xy];
                end
            end
            if (isempty(domain.fu) || isempty(domain.fv))
                % (u,v) PDFs are some perturbation of the gas velocities
                [~,ug,vg] = fluid.interpFluid(samplesXY(:,1),samplesXY(:,2));
                meanU = mean(ug); sigU = max(1e-10,sqrt(var(ug)));
                meanV = mean(vg); sigV = max(1e-10,sqrt(var(vg)));
                domain.fu = domain.setPDF('Gaussian',[meanU sigU]);
                domain.fv = domain.setPDF('Gaussian',[meanV sigV]);
            end
            fu1 = domain.fu(1,1); fu2 = domain.fu(end,1);
            fv1 = domain.fv(1,1); fv2 = domain.fv(end,1);
            fR1 = domain.fR(1,1); fR2 = domain.fR(end,1);
            fe1 = domain.fe(1,1); fe2 = domain.fe(end,1);
            ft1 = domain.ft(1,1); ft2 = domain.ft(end,1);
            s0 = repmat([fu1 fv1 fR1 fe1 ft1],numClumps,1);
            sf = repmat([fu2 fv2 fR2 fe2 ft2],numClumps,1);
            samples = [samplesXY, unifrnd(s0,sf,numClumps,5)];
            % Code to evaluate f(x,u,R,e;t) at sample points
            fuEval = interp1(domain.fu(:,1),domain.fu(:,2),samples(:,3));
            fvEval = interp1(domain.fv(:,1),domain.fv(:,2),samples(:,4));
            fREval = interp1(domain.fR(:,1),domain.fR(:,2),samples(:,5));
            feEval = interp1(domain.fe(:,1),domain.fe(:,2),samples(:,6));
            ftEval = interp1(domain.ft(:,1),domain.ft(:,2),samples(:,7));
            fEval = fxyEval.*fuEval.*fvEval.*fREval.*feEval.*ftEval;
            mEval = 4/3*pi*domain.rhol.*samples(:,5).^3;
            % Determine number of droplets per clump by balancing mass
            domain.calcAvgMass();
            dV = domain.mass_avg/sum(fEval.*mEval);
            nDroplet = fEval.*dV;
            samples = samples;
            % Round the number of particles per clump to nearest integer
            ind = find(nDroplet<1);
            nDroplet(ind) = 1;
            % TEMPORARY OVERRIDE *************
            nDroplet(:) = 1;
            % ********************************
            domain.nDroplet = round(nDroplet);
            domain.samples = samples;
            domain.numParcels = numClumps;
        end
        
        function calcInjectionDomain(domain,fluid,airfoil)
            % Function to calculate a box which contains the droplet clusters
            % INPUTS: 
            %   fluid,airfoil: fluid and airfoil objects
            
            % CALCULATE DOMAIN BOUNDS *************************************
            x0 = -5; dx = 0.1;
            yLOW = -0.5; yHIGH = 0;
            [yLIM1,yLIM2] = domain.calcImpingementLimits(fluid,airfoil,x0,yLOW,yHIGH);
            [yLIM3,yLIM4] = domain.calcImpingementLimits(fluid,airfoil,x0-dx,yLOW,yHIGH);
            domain.XY_bounds = [x0 yLIM2; x0 yLIM1; x0-dx yLIM4; x0-dx yLIM3];
            % *************************************************************
            % OR
            % LOAD FROM PRECOMPUTATIONS ***********************************
            %{
            % Save new boundaries of particles
            %domain.XY_bounds = [x1'; x2'; x3'; x4'];
            load('XY_bounds.mat');
            domain.XY_bounds = XY_bounds;
            % Estimate time to traverse domain
            dx = 0.1;
            %}
            domain.dtTraverse_avg = dx/fluid.Uinf;
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
            elseif strcmp(strType,'Custom')
                % 'params' = [X f(X)]
                minx = params(1,1);
                maxx = params(end,1);
                xsamp = linspace(minx,maxx,1000)';
                ysamp = interp1(params(:,1),params(:,2),xsamp);
                normalize = trapz(xsamp,ysamp);
                func = ysamp./normalize;
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
            % Function to compute the average total mass of the simulation
            % INPUTS:
            %   'strI': string defining type of integrand
            
            rhol = domain.rhol;
            I_xy = domain.numDroplets_avg;
            I_u = 1; I_v = 1; I_e = 1;
            I_t = domain.simTime/domain.dtTraverse_avg;
            R = domain.fR(:,1); fR = domain.fR(:,2);
            Integrand = 4/3*pi*rhol*(R.^3);
            I_mass = trapz(R,Integrand.*fR);
            domain.mass_avg = I_mass.*I_xy.*I_u.*I_v.*I_e.*I_t;
            
        end
        
        function dispSampleStatistics(domain)
            % Function to display sample statistics
            
            samples = domain.samples;
            % xy
            maxC = max(samples(:,5)); minC = min(samples(:,5));
            C = 9*(samples(:,5)-minC)./(maxC-minC) + 1;
            figure(10); subplot(2,4,1); plot(domain.XY_bounds(:,1),domain.XY_bounds(:,2),'o');
            hold on; scatter(samples(:,1),samples(:,2),[],C,'filled');
            % u,v,R,e,t
            fX = {domain.fu(:,1),domain.fv(:,1),domain.fR(:,1),domain.fe(:,1),domain.ft(:,1)};
            fY = {domain.fu(:,2),domain.fv(:,2),domain.fR(:,2),domain.fe(:,2),domain.ft(:,2)};
            for i=3:7
                figure(10); subplot(2,4,i-1); bar(samples(:,i),domain.nDroplet);
                hold on; plot(fX{i-2},fY{i-2});
            end
            % mR
            figure(10); subplot(2,4,7); bar(samples(:,5).^3*4/3*pi*domain.rhol,domain.nDroplet);
            hold on; plot(fX{3}.^3*4/3*pi*domain.rhol,fY{3})
            % Cumulative distribution of nDroplets
            [N,bins] = hist(domain.nDroplet,1:max(domain.nDroplet));
            figure(10); subplot(2,4,8); plot(bins,cumsum(N)./sum(N))
            
            
        end
        
        function [yLIM1,yLIM2] = calcImpingementLimits(domain,fluid,airfoil,x0,yLOW,yHIGH)
            % Function to calculate the impingement limits of the airfoil
            % [yLIM1, yLIM2] = lower and upper impingement limits at x0
            
            disp('Begin calculation of impingement limits');
            % Set up screen of particles at domain outlet
            ySamp = 400;
            ySCREEN = linspace(yLOW,yHIGH,ySamp)';
            yTOL = 1e-5;
            % Initialize cloud using the screen
            [pg,ug,vg] = fluid.interpFluid(repmat(x0,ySamp,1),ySCREEN);
            x = fluid.x; y = fluid.y; rhol = fluid.rhol;
            Rd = domain.R_avg*ones(ySamp,1);
            Unew = ug + 0.001*ug*unifrnd(-1,1);
            Vnew = vg + 0.001*vg*unifrnd(-1,1);
            cloud = SLDcloud([repmat(x0,ySamp,1) ySCREEN Unew Vnew Rd zeros(ySamp,3)],rhol,ySamp,fluid,'NoTResolve');
            STATE = [cloud.x cloud.y];
            % Find initial bounds for upper and lower impingement limits
            iter = 1; maxiter = 3000;
            while((iter<maxiter) && (length(cloud.impingeTotal)<ySamp))
                % Call subroutine to calculate local timesteps and impinging particles
                calcDtandImpinge(cloud,airfoil,fluid);
                % Advect particles one time step
                transportSLD(cloud,fluid);
                if (~isempty(cloud.impinge))
                    cloud.computeImpingementParams(airfoil);
                end
                % Save new positions
                Xnew = cloud.x; Ynew = cloud.y;
                STATE = [STATE; [Xnew Ynew]];
                iter = iter+1;
            end
            figure(12); clf; plot(STATE(:,1),STATE(:,2),'b.');
            hold on; plot(airfoil.PANELx,airfoil.PANELy,'k'); axis([-.5,.3,-.3,.3]); drawnow;
            impinge = cloud.impingeTotal;
            indLOW = min(impinge);
            indHIGH = max(impinge);
            yLOW1 = ySCREEN(indLOW-1);
            yLOW2 = ySCREEN(indLOW);
            yHIGH1 = ySCREEN(indHIGH-1);
            yHIGH2 = ySCREEN(indHIGH);
            % Iteratively advect and refine lower screen limits
            fprintf('Calculating lower impingement limit...');
            ySCREEN = linspace(yLOW1,yLOW2,ySamp)';
            diffLOW = 10*yTOL;
            [pg,ug,vg] = fluid.interpFluid(repmat(x0,ySamp,1),ySCREEN);
            uSCREEN = ug + 0.001*ug*unifrnd(-1,1);
            vSCREEN = vg + 0.001*vg*unifrnd(-1,1);
            while(diffLOW > yTOL)
                STATE = [];
                iter = 1;
                clear cloud;
                cloud = SLDcloud([repmat(x0,ySamp,1) ySCREEN uSCREEN vSCREEN Rd zeros(ySamp,3)],rhol,ySamp,fluid,'NoTResolve');
                % Advect the sheet of particles
                while((iter<maxiter) && (length(cloud.impingeTotal)<ySamp))
                    % Call subroutine to calculate local timesteps and impinging particles
                    calcDtandImpinge(cloud,airfoil,fluid);
                    % Advect particles one time step
                    transportSLD(cloud,fluid);
                    if (~isempty(cloud.impinge))
                        cloud.computeImpingementParams(airfoil);
                    end
                    % Update field quantities at new particle positions
                    Xnew = cloud.x; Ynew = cloud.y;
                    STATE = [STATE; [Xnew,Ynew]];
                    iter = iter+1;
                end
                figure(12); clf; plot(STATE(:,1),STATE(:,2),'b.');
                hold on; plot(airfoil.PANELx,airfoil.PANELy,'k'); axis([-.5,.3,-.3,.3]); drawnow;
                % Find the topmost and bottommost impingements
                impinge = cloud.impingeTotal;
                indHIT = min(impinge);
                % Update screen properties
                yHIT = ySCREEN(indHIT);
                if (indHIT == 1)
                    yMISS = yHIT - diffLOW;
                else
                    indMISS = indHIT-1;
                    yMISS = ySCREEN(indMISS);
                end
                ySCREEN = linspace(yMISS,yHIT,ySamp)';
                [pg,ug,vg] = fluid.interpFluid(repmat(x0,ySamp,1),ySCREEN);
                uSCREEN = ug + 0.001*ug*unifrnd(-1,1);
                vSCREEN = vg + 0.001*vg*unifrnd(-1,1);
                % Calculate diffLOW and diffHIGH
                diffLOW = abs(yHIT-yMISS)/(ySamp-1);
            end
            yLIM1 = yHIT;
            fprintf('Done.\n');
            % Iteratively advect and refine upper screen limits
            fprintf('Calculating upper impingement limit...');
            ySCREEN = linspace(yHIGH1,yHIGH2,ySamp)';
            diffHIGH = 10*yTOL;
            [pg,ug,vg] = fluid.interpFluid(repmat(x0,ySamp,1),ySCREEN);
            uSCREEN = ug + 0.001*ug*unifrnd(-1,1);
            vSCREEN = vg + 0.001*vg*unifrnd(-1,1);
            while(diffHIGH > yTOL)
                STATE = [];
                iter = 1;
                clear cloud;
                cloud = SLDcloud([repmat(x0,ySamp,1) ySCREEN uSCREEN vSCREEN Rd zeros(ySamp,3)],rhol,ySamp,fluid,'NoTResolve');
                % Advect the sheet of particles
                while((iter<maxiter) && (length(cloud.impingeTotal)<ySamp))
                    % Call subroutine to calculate local timesteps and impinging particles
                    calcDtandImpinge(cloud,airfoil,fluid);
                    % Advect particles one time step
                    transportSLD(cloud,fluid);
                    if (~isempty(cloud.impinge))
                        cloud.computeImpingementParams(airfoil);
                    end
                    % Update field quantities at new particle positions
                    Xnew = cloud.x; Ynew = cloud.y;
                    STATE = [STATE; [Xnew,Ynew]];
                    iter = iter+1;
                end
                figure(12); clf; plot(STATE(:,1),STATE(:,2),'b.');
                hold on; plot(airfoil.PANELx,airfoil.PANELy,'k'); axis([-.5,.3,-.3,.3]); drawnow;
                % Find the topmost and bottommost impingements
                impinge = cloud.impingeTotal;
                indHIT = max(impinge);
                yHIT = ySCREEN(indHIT);
                if (indHIT == ySamp)
                    yMISS = yHIT + diffHIGH;
                else
                    indMISS = indHIT+1;
                    yMISS = ySCREEN(indMISS);
                end
                ySCREEN = linspace(yMISS,yHIT,ySamp)';
                [pg,ug,vg] = fluid.interpFluid(repmat(x0,ySamp,1),ySCREEN);
                uSCREEN = ug + 0.001*ug*unifrnd(-1,1);
                vSCREEN = vg + 0.001*vg*unifrnd(-1,1);
                % Calculate diffLOW and diffHIGH
                diffHIGH = abs(yHIT-yMISS)/(ySamp-1);
            end
            yLIM2 = yHIT;
            fprintf('Done.\n');
            
        end

        
    end
    
end

