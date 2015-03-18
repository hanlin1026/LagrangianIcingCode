%% Initialization

% Read in flow solution
meshfile = 'NACA23012/MESH.P3D';
solnfile = 'NACA23012/q103.0.25E+01.bin';
[~,~,~,~,~,~,mach,alpha,~,~] = readp3d(meshfile,solnfile);
% Dimensional reference quantities
pinf = 1.01325e5; % N/m^2
R = 287.058; % J/kg/K
Tinf = 300; % K
rhoinf = pinf/R/Tinf;
Ubar = sqrt(1.4*pinf/rhoinf);
rhol = 1000; % kg/m^3
Rd = 0.5*(111e-6); % m
LWC = 0.73e-3; % kg/m^3
Uinf = mach*340; % m/s
% Initialize fluid object
scalars = [pinf;R;Tinf;rhoinf;Ubar;rhol];
fluid = Fluid(scalars,meshfile,solnfile);
% Initialize the airfoil surface
x = fluid.x; y = fluid.y;
ind = find(x(:,1)<=1);
ax = x(ind,1); ay = y(ind,1);
airfoil = Airfoil([ax,ay]);
airfoil.calcStagPt(fluid);
%%
% Initialize domain and associated distribution functions
strPDFTypes = {'Implicit','Implicit','Custom','Gaussian','Uniform'};
simTime = 60*30;
uDIST = 0*Uinf*cos(alpha/180)*[1 0.00001];
vDIST = 0*Uinf*sin(alpha/180)*[1 0.00001];
load('BIN27.mat');
BIN27(:,1) = BIN27(:,1);
rDIST = BIN27;
PDFparams = {};
PDFparams{1} = uDIST; 
PDFparams{2} = vDIST; 
PDFparams{3} = rDIST;
PDFparams{4} = [0 2];
PDFparams{5} = [0 simTime];
domain = InjectionDomain(strPDFTypes,PDFparams,fluid,airfoil,LWC,simTime);
nClumps = 1000;
domain.sampleRealization(nClumps,fluid);
domain.dispSampleStatistics();

%% Collection efficiency

STATE = {};
strImpMod = 'NoImpingement';
for i=1:100
    i
    domain.nDroplet = [];
    domain.samples = [];
    domain.numParcels = [];
    domain.sampleRealization(nClumps);
    
    [STATE,totalImpinge,impinge,s,beta,sStick,mStick] = calcCollectionEfficiency(airfoil,fluid,domain,strImpMod);
    airfoil.calcStagPt(fluid);
    sCENT = 0.5*(s(2:end)+s(1:end-1)) - airfoil.stagPt;
    figure(2); hold on; plot(sCENT,beta);
    [xSMOOTH,ySMOOTH] = smoothGaussianKernel1D(sCENT,beta,0.0075);
    figure(3); hold on; plot(xSMOOTH,ySMOOTH,'k');
    BETA{i} = [sCENT,beta];
end
save('BETA_23012_BIN27_SPLASH.mat','BETA');

%% Simulation

tic;
STATE = {};
% Create tree-search object for nearest neighbor searches
IMP = 1;
for t=1:tsteps
    % Subroutine to calculate new particles entering the domain (used if a
    % continuous stream of particles is desired)
    cloud.calcNewParticles(fluid);
    % Call subroutine to calculate local timesteps and impinging particles
    calcDtandImpinge(cloud,airfoil,fluid);
    % Advect particles one time step
    transportSLD(cloud,fluid);
    % Check for fracture
    fragmentSLD(cloud,fluid);
    % Call subroutine to calculate impingement regimes for impinging particles
    if ~isempty(cloud.impinge)
        %break;
        impingementRegimeSLD(cloud,airfoil);
    end
    % Save state variables
    STATE{t} = cloud.getState();
    t
    
end
toc

%% Visualization (Video)

indend = 1500;
indstop = floor(indend/10);
impimp = cloud.impingeTotal(cloud.impingeTotal<=particles);
impsplash = cloud.parentind;

% Colormap scheme
cmap = jet(100);
maxt = max(STATE{indend}(:,6));
GAIN = 130/maxt;
F(indstop) = struct('cdata',[],'colormap',[]);
for i=1:indstop
    ind=10*i;
    figure(10); hold on; plot(x(:,1),y(:,1),'k');
    scatter(STATE{ind}(:,1),STATE{ind}(:,2),3,'b','filled');
    scatter(STATE{ind}(impimp,1),STATE{ind}(impimp,2),3,'go');
    scatter(STATE{ind}(impsplash,1),STATE{ind}(impsplash,2),3,'ro');
    xlim([-.1 .4]); ylim([-.25 .25]); 
    F(i) = getframe;
    pause(.01); clf;
end

%% Visualization (Plot)

X = []; Y = []; CMAP = [];
for i=1:indstop
    ind = 10*i;
    X = [X; STATE{ind}(:,1)]; Y = [Y; STATE{ind}(:,2)];
    %colorind = round(GAIN*STATE{ind}(:,6));
    %padcmap = size(STATE{ind},1);
    %CMAP = [CMAP; repmat(cmap(i,:),padcmap,1)];
end
tic;
figure(1); hold on; scatter(X,Y,10,'k','filled');
hold on; plot(x(:,1),y(:,1),'b','LineWidth',3); axis equal;
xlim([-.25 1.1]); ylim([-.6 .6]);
toc
%{
indimp = cloud.impinge;
hold on; plot(cloud.x(indimp),cloud.y(indimp),'o','Color','r');

hold on; scatter(cloud.x(cloud.bounce),cloud.y(cloud.bounce),'filled','MarkerFaceColor','g');
hold on; scatter(cloud.x(cloud.spread),cloud.y(cloud.spread),'filled','MarkerFaceColor','r')
hold on; scatter(cloud.x(cloud.splash),cloud.y(cloud.splash),'filled','MarkerFaceColor','k')

% Highlight current state
X = STATE{indend}(:,1); Y = STATE{indend}(:,2);
figure(1); hold on; scatter(X,Y,10,'o','MarkerFaceColor','g','MarkerEdgeColor','g');
%}







