classdef SLDcloud < hgsetget
    % Defines a cloud of SLD droplets
    
    properties
        % Vectors (NOT arrays) of current (NOT also past) state information
        x=[]; y=[]; % Position
        u=[]; v=[]; % Velocity
        rd=[]; % Radius
        Temp=[]; % Temperature
        time=[]; % Times at which parcels enter the injection domain
        numDroplets=[]; % Number of droplets per clump
        % Other parameters (time step, fracture, impingement, etc.)
        FLAGtimeResolve; % '0' if no time resolution, '1' if time resolved
        tGLOB; % Current global time of the simulation
        dt; % Timestep (single global value for time-resolved simulation)
        indT; % Indices of particles currently in the simulation
        indAdv; % Indices of particles currently in simulation being advected
        rhol=[]; % Density
        fracture=[]; % '0' if not, '1' if yes
        impinge=[]; % Indices of currently impinging particles
        impingeTotal=[]; % Indices of all impinged particles, past and present
        bounce=[]; spread = []; splash = []; % Impingement regime indices for impinged particles
        parentind=[]; % Numeric indices of parents to splash particles
        childind={}; % Numeric indices of child particles of a splash
        particles; % Number of total current particles
        originalNumParticles; % Number of original particles
        % Impingement Parameters
        sigma;
        T;
        K = [];
        fs = [];
        fb = [];
        Ks0;
        Kb0;
        normvelsq=[];
        tangvel = [];
        % Current cells occupied by particles
        indCell;

    end 
    
    methods
        function cloud = SLDcloud(state0,rhol,particles,fluid,FLAGtimeResolve)
            % Constructor: state = (x0,y0,u0,v0,r0,temp0,t0,nDrop0)
            
            cloud.x(:,1) = state0(:,1);
            cloud.y(:,1) = state0(:,2);
            cloud.u(:,1) = state0(:,3);
            cloud.v(:,1) = state0(:,4);
            cloud.rd(:,1) = state0(:,5);
            cloud.Temp(:,1) = state0(:,6);
            cloud.time(:,1) = state0(:,7);
            cloud.numDroplets(:,1) = state0(:,8);
            
            if strcmp(FLAGtimeResolve,'TimeResolved')
                cloud.FLAGtimeResolve = 1;
            else
                cloud.FLAGtimeResolve = 0;
            end
            cloud.tGLOB = min(cloud.time);
            cloud.rhol = rhol;
            cloud.particles = particles;
            cloud.originalNumParticles = particles;
            
            cloud.sigma = 75.64e-3; % Surface tension of water at 0 deg C against air
            cloud.T = 270; % Temperature (K) for SLD's (estimate, slightly below freezing)
            
            % Initialize particle positions
            indInit = fluid.searchTree([cloud.x,cloud.y]);
            set(cloud,'indCell',indInit);
            
        end
        
        function cloud = addParticle(cloud,state,indCellP)
            % Function to add a particle to the SLD cloud
            
            i = size(cloud.x,1)+1;
            numnew = size(state,1);
            cloud.x(i:i+numnew-1,1) = state(:,1);
            cloud.y(i:i+numnew-1,1) = state(:,2);
            cloud.u(i:i+numnew-1,1) = state(:,3);
            cloud.v(i:i+numnew-1,1) = state(:,4);
            cloud.rd(i:i+numnew-1,1) = state(:,5);
            cloud.Temp(i:i+numnew-1,1) = state(:,6);
            cloud.time(i:i+numnew-1,1) = state(:,7);
            cloud.numDroplets(i:i+numnew-1,1) = state(:,8);
            cloud.indCell(i:i+numnew-1,1) = indCellP;
            % Update total number of particles
            cloud.particles = cloud.particles+numnew;
        end
        
        function cloud = computeImpingementParams(cloud,airfoil)
            % Function to compute and store the various parameters needed
            % in the impingement module calculations
            
            % Pull out state variables of all particles which have impinged
            indimp = cloud.impinge;
            x = cloud.x(indimp); y = cloud.y(indimp); u = cloud.u(indimp); v = cloud.v(indimp);
            rd = cloud.rd(indimp); t = cloud.time(indimp); rhol = cloud.rhol;
            T = cloud.Temp(indimp); % Temperature (K) for SLD's (estimate, slightly below freezing)
            
            sigma = cloud.sigma; % Surface tension of water at 0 deg C against air
            mul = (2.414e-5)*10.^(247.8./(T-140)); % Estimate of mu_water (N*s/m^2)

            % Find local points of impingement, normal vectors
            [~,~,nx,ny,tx,ty] = airfoil.findPanel(x,y);
            % Compute (vectorized) dot products of droplet velocities with normal vectors
            vnormsq = ((u.*nx) + (v.*ny)).^2;
            vtang = sum(([u,v].*[tx,ty]),2);
            We = 2*rhol.*rd.*(vnormsq)./sigma;
            % TEMPORARY PLOTTING ***********************************
            s = airfoil.XYtoScoords(x,y);
            figure(15); hold on; plot(s-airfoil.stagPt,vnormsq,'b.');
            % ******************************************************
            % Compute Ohnesorge number
            Oh = sqrt(mul./2./rhol./sigma./rd);
            % Compute Cossali number
            K = We.*(Oh.^(-0.4));
            % Parameters used in impingement regime calculation
            Ks0 = 3000;
            Kb0 = 600;
            wsr = 20/3;
            wsf = 5/6;
            wbr = 32;
            wbf = 1;
            hr = 20e-6; % Mean surface roughness height (m)
            hf = 0; % Mean surface liquid film height (m)
            % Calculated parameters that depend on elemental ones above
            R = hr./2./rd; % Dimensionless wall roughness height
            delta = hf./2./rd; % Dimensionless film thickness
            R_tilda = (R.^2)./(R+delta); % Modified wall roughness due to presence of film
            fs = (1+R_tilda.^2).*(1+delta.^2)./(1+wsr*R_tilda.^2)./(1+wsf*delta.^2); % Correction factor
            fb = (1+R_tilda.^2).*(1+delta.^2)./(1+wsr*R_tilda.^2)./(1+wbr*R_tilda.^4)./(1+wbf*delta.^2); % Correction factor
            % Set impingement parameter properties in the cloud
            cloud.K = K;
            cloud.Ks0 = Ks0;
            cloud.Kb0 = Kb0;
            cloud.fs = fs;
            cloud.fb = fb;
            % Calculate impingement regime
            ind1 = find(K < Kb0*fb);
            ind2 = find((K > Kb0*fb) & (K < Ks0*fs));
            ind3 = find(K > Ks0*fs);
            bounce = ind1;
            %spread = ind2;
            spread = [];
            splash = [ind2; ind3];
            % TEMPORARY PLOTTING ***********************************
            s = airfoil.XYtoScoords(x,y);
            figure(16); hold on; plot(s(ind3)-airfoil.stagPt,K(ind3)./(Ks0*fs(ind3)),'b.');
            % ******************************************************
            % Set index trackers in the cloud
            set(cloud,'bounce',[]); set(cloud,'bounce',bounce);
            set(cloud,'spread',[]); set(cloud,'spread',spread);
            set(cloud,'splash',[]); set(cloud,'splash',splash);
            set(cloud,'normvelsq',[]); set(cloud,'normvelsq',vnormsq);
            set(cloud,'tangvel',[]); set(cloud,'tangvel',vtang);
            % Update impingeTotal
            splashSpread = [indimp(spread); indimp(splash)];
            impingeNew = setdiff(splashSpread,cloud.impingeTotal);
            set(cloud,'impingeTotal',impingeNew);
        end
        
        function computeNewCellLocations(cloud,fluid)
            % Function to compute new cells occupied by particles after
            % transporting them one timestep
            
            indAdv = cloud.indAdv;
            indCell = cloud.indCell(indAdv);
            % Particle positions, cell centers
            xP = cloud.x(indAdv); yP = cloud.y(indAdv);
            xC = fluid.MEANx; yC = fluid.MEANy;
            ni = size(fluid.MEANx,1); nj = size(fluid.MEANx,2);
            C = indCell;
            N = C+ni; S = C-ni; E = C+1; W = C-1;
            SW = S-1; SE = S+1; NW = N-1; NE = N+1;
            % Transformed particle position
            [PI,PJ] = fluid.transformXYtoIJ(C,[xP,yP]);
            % Transformed neighbor positions
            [CI,CJ] = fluid.transformXYtoIJ(C,[xC(C),yC(C)]);
            [NI,NJ] = fluid.transformXYtoIJ(C,[xC(N),yC(N)]);
            [EI,EJ] = fluid.transformXYtoIJ(C,[xC(E),yC(E)]);
            [WI,WJ] = fluid.transformXYtoIJ(C,[xC(W),yC(W)]);
            [NWI,NWJ] = fluid.transformXYtoIJ(C,[xC(NW),yC(NW)]);
            [NEI,NEJ] = fluid.transformXYtoIJ(C,[xC(NE),yC(NE)]);
            % Catch any particles trying to "glitch" through the airfoil surface
            noGlitch = S>0; glitch = S<=0;
            SI = zeros(length(indAdv),1); SJ = zeros(length(indAdv),1);
            SWI = zeros(length(indAdv),1); SWJ = zeros(length(indAdv),1);
            SEI = zeros(length(indAdv),1); SEJ = zeros(length(indAdv),1);
            [SI(noGlitch),SJ(noGlitch)] = fluid.transformXYtoIJ(C(noGlitch),[xC(S(noGlitch)),yC(S(noGlitch))]);
            [SWI(noGlitch),SWJ(noGlitch)] = fluid.transformXYtoIJ(C(noGlitch),[xC(SW(noGlitch)),yC(SW(noGlitch))]);
            [SEI(noGlitch),SEJ(noGlitch)] = fluid.transformXYtoIJ(C(noGlitch),[xC(SE(noGlitch)),yC(SE(noGlitch))]);
            % 9NN search in the transformed plane
            dC = (CI-PI).^2 + (CJ-PJ).^2;
            dN = (NI-PI).^2 + (NJ-PJ).^2;
            dE = (EI-PI).^2 + (EJ-PJ).^2;
            dW = (WI-PI).^2 + (WJ-PJ).^2;
            dNW = (NWI-PI).^2 + (NWJ-PJ).^2;
            dNE = (NEI-PI).^2 + (NEJ-PJ).^2;
            % Catch any particles trying to "glitch" through airfoil surface
            dS = zeros(length(indAdv),1); dSW = zeros(length(indAdv),1); dSE = zeros(length(indAdv),1);
            dS(noGlitch) = (SI(noGlitch)-PI(noGlitch)).^2 + (SJ(noGlitch)-PJ(noGlitch)).^2;
            dSW(noGlitch) = (SWI(noGlitch)-PI(noGlitch)).^2 + (SWJ(noGlitch)-PJ(noGlitch)).^2;
            dSE(noGlitch) = (SEI(noGlitch)-PI(noGlitch)).^2 + (SEJ(noGlitch)-PJ(noGlitch)).^2;
            dS(glitch) = inf;
            dSW(glitch) = inf;
            dSE(glitch) = inf;
            % Find minimum
            [~,indMin] = min([dC,dN,dS,dE,dW,dSW,dSE,dNW,dNE]');
            indNN = [C,N,S,E,W,SW,SE,NW,NE];
            indCellNew = zeros(length(indAdv),1);
            for i=1:length(indAdv)
                indCellNew(i) = indNN(i,indMin(i));
            end
            set(cloud,'indCell',[indCellNew,indAdv]);
            % Plot
            %{
            for i=1:length(indAdv)
                figure(1); hold on; scatter(fluid.MEANx(indCell(i)),fluid.MEANy(indCell(i)),'ro');
            end
            %}
            %{
            for i=1:length(indAdv)
                figure(1+i); hold on; scatter([CI(i);NI(i);SI(i);EI(i);WI(i);SWI(i);SEI(i);NWI(i);NEI(i)],...
                    [CJ(i);NJ(i);SJ(i);EJ(i);WJ(i);SWJ(i);SEJ(i);NWJ(i);NEJ(i)],'b');
                hold on; scatter(PI(i),PJ(i),'r');
            end
            %}
        end
        
        function cloud = deleteParticles(cloud,ind)
            % Function to delete a particle from the SLD cloud
            
            % Delete appropriate elements from current state vectors
            cloud.x(ind) = [];
            cloud.y(ind) = [];
            cloud.u(ind) = [];
            cloud.v(ind) = [];
            cloud.rd(ind) = [];
            cloud.Temp(ind) = [];
            cloud.numDroplets(ind) = [];
            cloud.time(ind) = [];
            cloud.indT(ind) = [];
            set(cloud,'indCell',[ind;-1]);
            % Delete elements from current particles being advected
            [commonInd,indIndAdv,~] = intersect(cloud.indAdv,ind);
            set(cloud,'dt',[indIndAdv;-1]);
            cloud.indAdv(indIndAdv) = [];
            % Delete elements from currently impinging particles
            [commonImpinge,indImp,~] = intersect(cloud.impinge,ind);
            set(cloud,'impinge',[indImp;-1]);
            set(cloud,'K',[indImp;-1]);
            set(cloud,'fs',[indImp;-1]);
            set(cloud,'fb',[indImp;-1]);
            set(cloud,'normvelsq',[indImp;-1]);
            set(cloud,'tangvel',[indImp;-1]);
            [~,indBounce,~] = intersect(cloud.bounce,indImp);
            [~,indSpread,~] = intersect(cloud.spread,indImp);
            [~,indSplash,~] = intersect(cloud.splash,indImp);
            set(cloud,'bounce',[indBounce;-1]);
            set(cloud,'spread',[indSpread;-1]);
            set(cloud,'splash',[indSplash;-1]);
            % Update total number of particles and indices
            numdel = size(ind,1);
            cloud.particles = cloud.particles - numdel;
        end
        
        function cloud = removeIrrelevantParticles(cloud,fluid)
            % Function to remove particles that won't impinge
            
            nWrap = size(fluid.MEANx,1);
            ind1 = find(cloud.indCell > 60*nWrap);
            ind2 = find(cloud.x > 0);
            ind = intersect(ind1,ind2);
            if (~isempty(cloud.impingeTotal))
                cloud.deleteParticles(ind);
            end
        end

        function cloud = set.x(cloud,x)
            cloud.x = x;
        end
        
        function cloud = set.y(cloud,y)
            cloud.y = y;
        end
        
        function cloud = set.u(cloud,u)
            if size(u,2) == 2
                % Set only those indices provided
                ind = u(:,2);
                cloud.u(ind) = u(:,1);
            else
                % Set all indices of u
                cloud.u = u;
            end
        end
        
        function cloud = set.v(cloud,v)
            if size(v,2) == 2
                % Set only those indices provided
                ind = v(:,2);
                cloud.v(ind) = v(:,1);
            else
                % Set all indices of v
                cloud.v = v;
            end
        end
        
        function cloud = set.rd(cloud,vars)
            if isempty(vars)
                % Clear rd
                cloud.rd = [];
            elseif size(vars,2)==2
                % Set specific indexed elements of rd
                cloud.rd(vars(:,1)) = vars(:,2);
            else
                % Set all of rd
                cloud.rd = vars;
            end
        end
        
        function cloud = set.time(cloud,time)
            cloud.time = time;
        end

        function cloud = set.fracture(cloud,fracture)
            if isempty(fracture)
                % Clear fracture
                cloud.fracture = [];
            elseif length(fracture)==1
                % Append fracture
                ind=length(cloud.fracture)+1;
                cloud.fracture(ind) = fracture;
            else
                % Set fracture
                cloud.fracture=fracture;
            end
        end
        
        function cloud = set.impinge(cloud,vars)
            if isempty(vars)
                % Clear impinge
                cloud.impinge = [];
            elseif size(vars,1)==1
                % Append a single element to impinge
                ind=length(cloud.impinge)+1;
                cloud.impinge(ind,1) = vars;
            elseif vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.impinge(ind) = [];
            else
                % Append multiple elements
                cloud.impinge = [cloud.impinge; vars];
            end
        end
        
        function cloud = set.impingeTotal(cloud,vars)
            if isempty(vars)
                cloud.impingeTotal = cloud.impingeTotal;
            else
                cloud.impingeTotal = [cloud.impingeTotal; vars];
            end
        end
        
        function cloud = set.bounce(cloud,vars)
            if isempty(vars)
                % Clear bounce
                cloud.bounce = [];
            elseif vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.bounce(ind) = [];
            else
                % Case to append new index trackers to cloud.bounce
                cloud.bounce = [cloud.bounce; vars];
            end
        end
        
        function cloud = set.spread(cloud,vars)
            if isempty(vars)
                cloud.spread = [];
            elseif vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.spread(ind) = [];
            else
                cloud.spread = [cloud.spread; vars];
            end
        end
        
        function cloud = set.splash(cloud,vars)
            if isempty(vars)
                cloud.splash = [];
            elseif vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.splash(ind) = [];
            else
                cloud.splash = [cloud.splash; vars];
            end
        end
        
        function cloud = set.normvelsq(cloud,vars)
            if isempty(vars)
                cloud.normvelsq = [];
            elseif vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.normvelsq(ind) = [];
            else
                cloud.normvelsq = [cloud.normvelsq; vars];
            end
        end
        
        function cloud = set.tangvel(cloud,vars)
            if isempty(vars)
                cloud.tangvel = [];
            elseif vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.tangvel(ind) = [];
            else
                cloud.tangvel = [cloud.tangvel; vars];
            end
        end
        
        function cloud = set.dt(cloud,vars)
            if isempty(vars)
                % Clear dt
                cloud.dt = [];
            elseif size(vars,2)==2
                % Set specific indexed elements of dt
                cloud.dt(vars(:,1)) = vars(:,2);
            elseif vars(end) == -1
                % Delete
                ind = vars(1:end-1);
                cloud.dt(ind) = [];
            else
                % Set all of dt
                cloud.dt = vars;
            end
        end
        
        function cloud = set.parentind(cloud,ind)
            % Add a parent index for splash tracking
            
            sizeP = size(cloud.parentind,1);
            numnew = size(ind,1);
            ind1 = sizeP+1;
            ind2 = ind1+numnew-1;
            cloud.parentind(ind1:ind2,1) = ind;
        end
        
        function cloud = set.childind(cloud,ind)
            % Add indices of child particles for splash tracking
            
            sizeC = size(cloud.childind,1);
            cloud.childind{sizeC+1,1} = ind;
            
        end
        
        function cloud = set.tGLOB(cloud,val)
            cloud.tGLOB = val;
        end
        
        function cloud = set.indCell(cloud,val)
            % Set indCell
            
            if size(val,2)==1
                if val(end) == -1
                    % Delete
                    ind = val(1:end-1);
                    cloud.indCell(ind) = [];
                else
                    % Set all elements
                    cloud.indCell = val;
                end
            elseif size(val,2)==2
                % Set specified elements
                cloud.indCell(val(:,2)) = val(:,1);
            end
        end
        
        function cloud = set.K(cloud,vars)
            if vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.K(ind) = [];
            else
                cloud.K = vars;
            end
            
        end
        
        function cloud = set.fs(cloud,vars)
            if vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.fs(ind) = [];
            else
                cloud.fs = vars;
            end
            
        end
        
        function cloud = set.fb(cloud,vars)
            if vars(end) == -1
                % Delete specified members
                ind = vars(1:end-1);
                cloud.fb(ind) = [];
            else
                cloud.fb = vars;
            end
        end
        
        function cloud = set.indAdv(cloud,vars)
            cloud.indAdv = vars;
        end
        
        function cloud = set.indT(cloud,vars)
            cloud.indT = vars;
        end
        
        function state = getState(cloud)
            % Function to return all current state variables
            indT = cloud.indT;
            state = [cloud.x(indT) cloud.y(indT) cloud.u(indT) cloud.v(indT) cloud.rd(indT) cloud.Temp(indT) cloud.time(indT) cloud.numDroplets(indT)];
        end
        
        function newVector = resetVector(cloud,oldVector,deleteInd)
            % Function that returns new indices after some indices have
            % been deleted. For use with 'deleteParticle' routine and index
            % tracking. Assumes that oldVector and deleteInd have no
            % elements in common.
            
            if isempty(oldVector)
                newVector = [];
            else
                newVector = zeros(length(oldVector),1);
                for i=1:length(oldVector)
                    % Find total number of deleted indices less in cardinality
                    % than the tracker index
                    subtract = sum(deleteInd<oldVector(i));
                    newVector(i) = oldVector(i) - subtract;
                end
            end
            
        end
        
        
    end
    
end

