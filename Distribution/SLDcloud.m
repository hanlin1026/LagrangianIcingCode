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

    end 
    
    methods
        function cloud = SLDcloud(state0,rhol,particles)
            % Constructor: state = (x0,y0,u0,v0,r0,temp0,t0,nDrop0)
            
            cloud.x(:,1) = state0(:,1);
            cloud.y(:,1) = state0(:,2);
            cloud.u(:,1) = state0(:,3);
            cloud.v(:,1) = state0(:,4);
            cloud.rd(:,1) = state0(:,5);
            cloud.Temp(:,1) = state0(:,6);
            cloud.time(:,1) = state0(:,7);
            cloud.numDroplets(:,1) = state0(:,8);
            
            cloud.tGLOB = min(cloud.time);
            cloud.rhol = rhol;
            cloud.particles = particles;
            cloud.originalNumParticles = particles;
            
            cloud.sigma = 75.64e-3; % Surface tension of water at 0 deg C against air
            cloud.T = 270; % Temperature (K) for SLD's (estimate, slightly below freezing)
            
        end
        
        function cloud = addParticle(cloud,state)
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
            vnormsq = sum(([u,v].*[nx,ny]).^2,2);
            vtang = sum(([u,v].*[tx,ty]),2);
            We = 2*rhol.*rd.*(vnormsq)./sigma;
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
            spread = ind2;
            splash = ind3;
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
        
        function cloud = deleteParticle(cloud,ind)
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
            % Update total number of particles and indices
            numdel = length(ind);
            cloud.particles = cloud.particles - numdel;
            % Fix index trackers due to particles being deleted
            newBounce = cloud.resetVector(cloud.bounce,ind);
            set(cloud,'bounce',[]); set(cloud,'bounce',newBounce);
            newSpread = cloud.resetVector(cloud.spread,ind);
            set(cloud,'spread',[]); set(cloud,'spread',newSpread);
            newImpinge = cloud.resetVector(cloud.impinge,ind);
            set(cloud,'impinge',[]); set(cloud,'impinge',newImpinge);
            
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
            elseif length(vars)==1
                % Set a particular element of impinge
                ind=length(cloud.impinge)+1;
                cloud.impinge(ind,1) = vars;
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
            else
                % Case to append new index trackers to cloud.bounce
                cloud.bounce = [cloud.bounce; vars];
            end
        end
        
        function cloud = set.spread(cloud,vars)
            if isempty(vars)
                cloud.spread = [];
            else
                cloud.spread = [cloud.spread; vars];
            end
        end
        
        function cloud = set.splash(cloud,vars)
            if isempty(vars)
                cloud.splash = [];
            else
                cloud.splash = [cloud.splash; vars];
            end
        end
        
        function cloud = set.normvelsq(cloud,vars)
            if isempty(vars)
                cloud.normvelsq = [];
            else
                cloud.normvelsq = [cloud.normvelsq; vars];
            end
        end
        
        function cloud = set.tangvel(cloud,vars)
            if isempty(vars)
                cloud.tangvel = [];
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

