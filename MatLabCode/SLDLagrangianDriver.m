%% Initialization

% Read in flow solution
meshfile = 'MESH.P3D';
solnfile = 'q103.0.50E+01.bin';

% Dimensional reference quantities
pinf = 1.01325e5; % N/m^2
R = 287.058; % J/kg/K
Tinf = 300; % K
rhoinf = pinf/R/Tinf;
Uinf = sqrt(pinf/rhoinf);
rhol = 1000; % kg/m^3
Rd = 100e-6; % m

% Initialize fluid object
scalars = [pinf;R;Tinf;rhoinf;Uinf;rhol;Rd];
fluid = Fluid(scalars,meshfile,solnfile);

% Random initial positions
tsteps = 1500; particles = 50;
x0 = unifrnd(-1,-0.5,[particles 1]);
y0 = unifrnd(-0.25,0,[particles 1]);
% Initialize radius, time
rd0 = Rd*ones(particles,1);
time0 = zeros(particles,1);
% Initialize velocities to local cell velocities
xq = x0; yq = y0;
[pg,ug,vg] = interpFluid(fluid,xq,yq);
u0 = ug + 0.01*ug.*unifrnd(-1,1,[particles,1]);
v0 = vg + 0.01*vg.*unifrnd(-1,1,[particles,1]);

% Initialize SLD cloud
cloud = SLDcloud([x0 y0 u0 v0 rd0 time0 [1:particles]'],rhol,particles);

% Initialize the airfoil surface
x = fluid.x; y = fluid.y;
ind = find(x(:,1)<=1);
ax = x(ind,1); ay = y(ind,1);
airfoil = Airfoil([ax,ay]);

%% Collection efficiency

STATE = {};
strImpMod = 'Impingement';
[STATE,totalImpinge,impinge,s,beta] = calcCollectionEfficiency(airfoil,fluid,strImpMod);
airfoil.calcStagPt(fluid);
figure; plot(s-airfoil.stagPt*ones(length(s),1),beta);

%% Simulation

tic;
STATE = {};
% Create tree-search object for nearest neighbor searches
IMP = 1;
for t=1:tsteps
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







