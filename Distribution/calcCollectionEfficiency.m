
% Function to calculate the collection efficiency of an airfoil
rng('shuffle');
% Addpath basedir
strImpMod = 'Impingement';
workdir = pwd;
basedir = '/home/adegenna/LagrangianIcingCode/Distribution/';
batchdir = '/home/adegenna/LagrangianIcingCode/Distribution/BatchJobs/';
addpath(basedir);
% Load stuff
cd(batchdir);
tmpF = load('fluid.mat'); tmpA = load('airfoil.mat');
fluid = tmpF.fluid; airfoil = tmpA.airfoil;

% Dimensional reference quantities
mach = 0.2287; 
alpha = 2.5;
pinf = 1.01325e5; % N/m^2
R = 287.058; % J/kg/K
Tinf = 300; % K
rhoinf = pinf/R/Tinf;
Ubar = sqrt(1.4*pinf/rhoinf);
rhol = 1000; % kg/m^3
Rd = 0.5*(111e-6); % m
LWC = 0.73e-3; % kg/m^3
Uinf = mach*340; % m/s
% Initialize domain and associated distribution functions
strPDFTypes = {'Implicit','Implicit','Uniform','Gaussian','Uniform'};
simTime = 60*30;
uDIST = 0*Uinf*cos(alpha/180)*[1 0.00001];
vDIST = 0*Uinf*sin(alpha/180)*[1 0.00001];
tmpB = load('MVD52.mat'); 
DISTr = tmpB.MVD52;
%BIN27(:,1) = BIN27(:,1);
%rDIST = BIN27;
% Hacky crap to use a monodistributed PDF for R
dirnum = [];
for i=1:length(workdir)
    tmp = str2num(workdir(i));
    if (~isempty(tmp)) && (isreal(tmp))
        dirnum = [dirnum workdir(i)];
    end
end
dirNUM = str2num(dirnum);
meanR = DISTr(dirNUM,1);
rDIST = [0.99*meanR 1.01*meanR];

PDFparams = {};
PDFparams{1} = uDIST; 
PDFparams{2} = vDIST; 
PDFparams{3} = rDIST;
PDFparams{4} = [0 2]+273.15;
PDFparams{5} = [0 simTime];
LWC = 0.73e-3;
domain = InjectionDomain(strPDFTypes,PDFparams,fluid,airfoil,LWC,simTime);
nClumps = 1000;
domain.sampleRealization(nClumps,fluid);
dlmwrite([workdir '/domainXY.dat'],domain.samples(:,1:2),'\t');
%}

cd(basedir);
airfoil.FILM = [];
% Initialize a cloud using domain realization
particles = domain.numParcels;
rhol = fluid.rhol;
cloud = SLDcloud([domain.samples, domain.nDroplet],rhol,particles,fluid,'NoTResolve');

% Advect particles (with/without fracture or impingement submodules)
disp('Advecting test particles between limits...');
totalImpinge = 0;
impinge = [];
STATE = {};
t = cloud.tGLOB; simTime = domain.simTime;
tSAMP = 1; tREFRESH = 0.03; tRATE = 20; tMARK = 1;
iter = 1; maxiter = 5000;
[xtmp,ytmp] = airfoil.interpStoXY(airfoil.LIMup);
% Calculate total mass flux through the injection domain screen
R = cloud.rd; nDrop = cloud.numDroplets;
mTOT = sum((4/3*pi*rhol)*(R.^3).*nDrop);
dy = domain.XY_bounds(1,2) - domain.XY_bounds(2,2);
if strcmp(strImpMod,'NoImpingement')
    % Calculate collection efficiency without impingement module
    while totalImpinge<cloud.particles && iter<maxiter
        % Fragmentation module
        %fragmentSLD(cloud,fluid);
        % Call subroutine to calculate local timesteps and impinging particles
        calcDtandImpinge(cloud,airfoil,fluid);
        % Advect particles one time step
        transportSLD(cloud,fluid);
        % Compute bookkeeping of which parcels have impinged
        if ~isempty(cloud.impinge)
            cloud.computeImpingementParams(airfoil);
            bounceDynamics(cloud,airfoil);
        end
        totalImpinge = size(cloud.impingeTotal,1);
        % Save state variables every 0.1 sec
        %
        maxC = max(cloud.rd); minC = min(cloud.rd);
        C = 9*(cloud.rd-minC)./(maxC-minC) + 1;
        if mod(iter,100)==0
            figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.1 .5 -.3 .3]); %axis([-.02 .15 -.1 .1])
        end
        %
        if floor(tSAMP/tRATE)>=1
            state = cloud.getState();
            part = [state(:,1), state(:,2)];
            %figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.5 1 -.3 .3]);
            figure(1); hold on; scatter(part(:,1),part(:,2),10,C,'filled');
            hold on; plot(xtmp,ytmp,'ko');
            drawnow;
            tMARK = iter;
        end
        tSAMP = t-tMARK;
        %}
        iter = iter+1;
        t = iter;
    end
    % Record mass deposited on the airfoil surface
    indStick = cloud.impingeTotal;
    xStick = cloud.x(indStick); yStick = cloud.y(indStick); 
    RStick = cloud.rd(indStick); nStick = cloud.numDroplets(indStick);
    sStick = airfoil.XYtoScoords(xStick,yStick);
    mStick = (4/3*pi*rhol)*(RStick.^3).*nStick;
    set(airfoil,'FILM',[sStick, mStick]);
elseif strcmp(strImpMod,'Impingement')
    countMovie = 1;
    % Calculate collection efficiency with impingement module
    while totalImpinge<cloud.particles && t<maxiter
        % Check for fracture
        %fragmentSLD(cloud,fluid);
        % Before starting, remove any particles which are extraneous
        %cloud.removeIrrelevantParticles(fluid);
        % Call subroutine to calculate local timesteps and impinging particles
        calcDtandImpinge(cloud,airfoil,fluid);
        % Call subroutine to calculate impingement regimes for impinging particles
        if ~isempty(cloud.impinge)
            impingementRegimeSLD(cloud,airfoil);
        end
        totalImpinge = size(cloud.impingeTotal,1);
        % Advect particles one time step
        transportSLD(cloud,fluid);
        % Save state variables every 0.1 sec
        %
        maxC = max(cloud.rd); minC = min(cloud.rd);
        C = 9*(cloud.rd-minC)./(maxC-minC) + 1;
        if mod(iter,tRATE)==tRATE-1
            figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.1 .5 -.3 .3]); %axis([-.02 .15 -.1 .1])
        end
        %
        if floor(tSAMP/tRATE)>=1
            state = cloud.getState();
            part = [state(:,1), state(:,2)];
            %figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.5 1 -.3 .3]);
            figure(1); hold on; scatter(part(:,1),part(:,2),10,C,'filled');
            drawnow;
            hold on; plot(xtmp,ytmp,'ko');
            % Movie
            %{
            if iter>700
                figure(1); F(countMovie) = getframe(gca);
                countMovie = countMovie+1;
            end
            %}
            tMARK = iter;
        end
        tSAMP = t-tMARK;
        %}
        iter = iter+1; t = iter;
        [t, cloud.particles]
    end
    % Record mass deposited on the airfoil surface
    sStick = airfoil.FILM(:,1);
    mStick = airfoil.FILM(:,2);
    % Save movie
    %save('MOVIE.mat','F');
end

% Calculate collection efficiency
disp('Calculating collection efficiency...');
STATE = cloud.getState();
% Plot airfoil and impinged parcels
%{
impinge = sort(cloud.impingeTotal);
maxC = max(cloud.rd); minC = min(cloud.rd);
C = 9*(cloud.rd-minC)./(maxC-minC) + 1;
figure(11); plot(airfoil.X,airfoil.Y,'k');
figure(11); hold on; scatter(STATE(impinge,1),STATE(impinge,2),[],C(impinge),'filled');
drawnow;
%}
sDrop = airfoil.interpXYtoS(STATE(impinge,1),STATE(impinge,2));
% Calculate total impinged mass in each of the surface bins
bins = 35;
sBin = linspace(min(sStick),max(sStick),bins+1)';
ds = sBin(2)-sBin(1);
for i=1:bins
    mBin(i,1) = sum(mStick(sStick>sBin(i) & sStick<sBin(i+1)));
end
% Calculate collection efficiency
beta = (mBin/ds).*(dy/mTOT);
s = sBin;
%figure(12); plot(sDrop,cloud.rd(impinge),'bo'); drawnow;
%figure(14); hist(cloud.rd(impinge)); drawnow;

% Output beta
airfoil.calcStagPt(fluid);
sCENT = 0.5*(s(2:end)+s(1:end-1)) - airfoil.stagPt;
BETA = [sCENT,beta]
dlmwrite([workdir '/BETA.dat'],BETA,'\t');
%}