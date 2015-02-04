%% Initialization

meshfile = 'MESH.P3D';
solnfile = 'q103.0.50E+01.bin';

[x,y,rho,rhou,rhov,e,mach,alpha,reynolds,time] = readp3d(meshfile,solnfile);

% Dimensional reference quantities
pinf = 1.01325e5; % N/m^2
R = 287.058; % J/kg/K
Tinf = 300; % K
rhoinf = pinf/R/Tinf;
Uinf = sqrt(pinf/rhoinf);
rhol = 1000; % kg/m^3
Rd = 250e-6; % m
Tinf = 300;  % K

clear xpart ypart;
tsteps = 1000; particles = 100;
xpart = cell(tsteps,1);
ypart = cell(tsteps,1);
vx = cell(tsteps,1);
vy = cell(tsteps,1);
rd = cell(tsteps,1);
time = cell(tsteps,1);

xpart{1} = unifrnd(-1,-0.5,[particles 1]);
ypart{1} = unifrnd(-0.25,0,[particles 1]);
rd{1} = Rd*ones(particles,1);
time{1} = zeros(particles,1);
%figure(1); hold on; plot(x,y,'k');

RHO = scatteredInterpolant(x(:),y(:),rho(:));
RHOU = scatteredInterpolant(x(:),y(:),rhou(:));
RHOV = scatteredInterpolant(x(:),y(:),rhov(:));
E = scatteredInterpolant(x(:),y(:),e(:));

% Initialize velocities to local cell velocities
xq = xpart{1}; yq = ypart{1};
pg = rhoinf*RHO(xq,yq);
ug = RHOU(xq,yq)*Uinf*rhoinf./pg;
vg = RHOV(xq,yq)*Uinf*rhoinf./pg;
vx{1} = ug + 0.1*ug.*unifrnd(-1,1,[particles,1]);
vy{1} = vg + 0.1*vg.*unifrnd(-1,1,[particles,1]);

% Declare vector to keep track of particles which have hit the surface
IMPINGE = cell(tsteps,1);
IMPINGE{1} = zeros(particles,1);

% Initialize the airfoil surface
ind = find(x(:,1)<=1);
ax = x(ind,1); ay = y(ind,1);
airfoil = Airfoil([ax,ay]);

%% Simulation
for t = 1:tsteps
    i = 1;
    % Pull out vector of current state variables
    Xt = xpart{t}; Yt = ypart{t}; VXt = vx{t}; VYt = vy{t}; RDt = rd{t}; TIMEt = time{t}; IMPINGEt = IMPINGE{t};
    DTt = [];
    % Declare new state variable vectors
    Xnew = []; Ynew = []; VXnew = []; VYnew = []; RDnew = []; TIMEnew = []; IMPINGEnew = [];
    while i<=particles
        % Pull out individual particle
        xq = Xt(i); yq = Yt(i); uq = VXt(i); vq = VYt(i); rq = RDt(i); impinge = IMPINGEt(i);
        [pBL, pBR, pTR, pTL] = findCell(x,y,xq,yq);
        if (uq==0 && vq==0) || (size(pBL,1)==1) % Are we on the airfoil surface?
            % YES
            Xnew(i,1) = xq;
            Ynew(i,1) = yq;
            VXnew(i,1) = 0;
            VYnew(i,1) = 0;
            RDnew(i,1) = rq;
            IMPINGEnew(i,1) = 1;
        else
            % NO
            L1 = norm(pBL-pBR);
            L2 = norm(pTR-pBR);
            area = L1*L2;
            vel = [uq; vq];
            dt = 0.2*sqrt(area)/norm(vel); DTt(i,1) = dt;
            pg = rhoinf*RHO(xq,yq);
            ug = RHOU(xq,yq)*Uinf*rhoinf/pg;
            vg = RHOV(xq,yq)*Uinf*rhoinf/pg;
            % Advect particle one time step
            [Xnew(i,1), Ynew(i,1), VXnew(i,1), VYnew(i,1)] = transportSLD(xq,yq,uq,vq,ug,vg,pg,dt,rhol,rq,Tinf);
            IMPINGEnew(i,1) = 0;
            TIMEnew(i,1) = TIMEt(i) + dt;
        end
        i=i+1;
    end
    RDnew = RDt;
    
    % Check for fragmentation
    % Interpolate flow field over all particles
    RHOfield = rhoinf*RHO(Xnew,Ynew);
    Ufield = RHOU(Xnew,Ynew)*Uinf*rhoinf./RHOfield;
    Vfield = RHOV(Xnew,Ynew)*Uinf*rhoinf./RHOfield;
    [rFRAG, xFRAG, yFRAG, uFRAG, vFRAG, timeFRAG, impingeFRAG, indFRAG, numFRAG] = fragmentSLD(Xnew,Ynew,VXnew,VYnew,TIMEnew,Ufield,Vfield,RHOfield,RDt,DTt,rhol,IMPINGEnew);
    if ~isempty(indFRAG)
        % Replace original vectors with new ones
        Xnew = replaceVector(Xnew,xFRAG,indFRAG);
        Ynew = replaceVector(Ynew,yFRAG,indFRAG);
        VXnew = replaceVector(VXnew,uFRAG,indFRAG);
        VYnew = replaceVector(VYnew,vFRAG,indFRAG);
        RDnew = replaceVector(RDnew,rFRAG,indFRAG);
        TIMEnew = replaceVector(TIMEnew,timeFRAG,indFRAG);
        IMPINGEnew = replaceVector(IMPINGEnew,impingeFRAG,indFRAG);
        % Update number of particles
        particles = particles + numFRAG - length(indFRAG);
    end
    % Set new variables
    xpart{t+1} = Xnew; ypart{t+1} = Ynew; vx{t+1} = VXnew; vy{t+1} = VYnew; rd{t+1} = RDnew; time{t+1} = TIMEnew; IMPINGE{t+1} = IMPINGEnew;
    % Runtime output, plots, etc.
    figure(10); hold on; plot(t,particles,'o','Color','b');
    t
end

%% Plotting

%Find max time
tmax=0;
for i=1:tsteps
    [val,ind] = max(time{i});
    if val>tmax
        tmax=val;
    end
end
% Plot streamlines
cmap = jet(1000);
tRES = linspace(0,tmax,1000)';
X = []; Y = []; IND = [];
for i=1:tsteps+1
    part = size(xpart{i},1);
    for j=1:part
        xtmp = xpart{i}(j); ytmp = ypart{i}(j);
        X = [X; xtmp]; Y = [Y; ytmp];
        tsamp = time{i}(j);
        [val,ind] = min(abs(tRES-tsamp)); ind = ind(1);
        IND = [IND; ind];
    end
end
figure(1); hold on; scatter(X,Y,10,cmap(IND,:),'filled');
hold on; plot(x(:,1),y(:,1),'k'); axis equal