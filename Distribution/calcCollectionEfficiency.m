function [STATE,totalImpinge,impinge,s,beta] = calcCollectionEfficiency(rdAvg,airfoil,fluid,strImpMod)
% Function to calculate the collection efficiency of an airfoil

% Calculate impingement limits of the airfoil
disp('Calculating impingement limits...');
xL = -0.5; Yhit = -0.1;
Ymiss = 0.3;
yLimUP = impingementLimitsSLD(rdAvg,fluid,airfoil,xL,Ymiss,Yhit,'UP');
Ymiss = -0.3;
yLimDOWN = impingementLimitsSLD(rdAvg,fluid,airfoil,xL,Ymiss,Yhit,'DOWN');

% Create cloud of particles between limits
particles = 100;
dy = (yLimUP-yLimDOWN)/(particles-1);
y0 = linspace(yLimDOWN,yLimUP,particles)';
x0 = xL*ones(particles,1);
[pg,ug,vg] = fluid.interpFluid(x0,y0);
u0 = ug + 0*ug.*unifrnd(-1,1,[particles,1]);
v0 = vg + 0*vg.*unifrnd(-1,1,[particles,1]);
rd0 = fluid.Rd*ones(particles,1);
time0 = zeros(particles,1);
rhol = fluid.rhol;
cloud = SLDcloud([x0 y0 u0 v0 rd0 time0 [1:particles]'],rhol,particles);

% Advect particles (with/without fracture or impingement submodules)
disp('Advecting test particles between limits...');
totalImpinge = 0;
impinge = [];
STATE = {};
t = 1;
maxiter = 2000;
if strcmp(strImpMod,'NoImpingement')
    % Calculate collection efficiency without impingement module
    while totalImpinge<particles && t<maxiter
        % Call subroutine to calculate local timesteps and impinging particles
        calcDtandImpinge(cloud,airfoil,fluid);
        % Advect particles one time step
        transportSLD(cloud,fluid);
        % Save state variables
        totalImpinge = size(cloud.impinge,1);
        STATE{t} = cloud.getState();
        t = t+1
    end
elseif strcmp(strImpMod,'Impingement')
    % Calculate collection efficiency with impingement module
    while totalImpinge<particles && t<maxiter
        % Call subroutine to calculate local timesteps and impinging particles
        calcDtandImpinge(cloud,airfoil,fluid);
        % Advect particles one time step
        transportSLD(cloud,fluid);
        % Check for fracture
        fragmentSLD(cloud,fluid);
        % Call subroutine to calculate impingement regimes for impinging particles
        if ~isempty(cloud.impinge)
            impingementRegimeSLD(cloud,airfoil);
        end
        % Save state variables
        STATE{t} = cloud.getState();
        t = t+1
    end
end

% Calculate collection efficiency
disp('Calculating collection efficiency...');
if strcmp(strImpMod,'NoImpingement')
    impinge = sort(cloud.impinge);
    xq = cloud.x(impinge);
    yq = cloud.y(impinge);
    s = airfoil.XYtoScoords(xq,yq);
    figure(1); hold on; scatter(STATE{1}(impinge,1),STATE{1}(impinge,2),'MarkerFaceColor','k');
    figure(1); hold on; scatter(cloud.x(impinge),cloud.y(impinge),'MarkerFaceColor','k');
    beta = dy./(s(2:end)-s(1:end-1));
    beta(end+1) = beta(end);
elseif strcmp(strImpMod,'Impingement')
    % Original mass that was in each of the advected tubes
    m_dy = 2*4/3*pi*fluid.Rd^3*rhol;
    % Obtain impingement locations of original droplets
    sImpOrig = sort(airfoil.originalImpingeScoord);
    % Calculate total impinged mass in each of the surface bins
    s_imp = airfoil.FILM(:,1);
    m_imp = airfoil.FILM(:,2);
    m_ds = zeros(length(sImpOrig)-1,1);
    for i=1:length(sImpOrig)-1
        ind = find(s_imp <= sImpOrig(i+1) & s_imp >= sImpOrig(i));
        m_ds(i) = sum(m_imp(ind));
    end
    % Calculate m_ds/m_dy ratio
    mds_mdy = m_ds./m_dy;
    % Calculate dy/ds ratio
    dy_ds = dy./(sImpOrig(2:end)-sImpOrig(1:end-1));
    % Collection efficiency
    s = sImpOrig;
    beta = mds_mdy.*dy_ds;
    beta(end+1) = beta(end);
    
end

% Plot in s-coordinates and x,y coordinates
%
figure(10); hold on; plot(s,beta,'r');
cmap = jet(totalImpinge);
figure(11); plot(airfoil.PANELx,airfoil.PANELy,'k');
[val,ind] = sort(beta);
[xq,yq] = airfoil.interpStoXY(s(ind));
for i=1:totalImpinge
    hold on; scatter(xq(i),yq(i),'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k');
end
%}











end