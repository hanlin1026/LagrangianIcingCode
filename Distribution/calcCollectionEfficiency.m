function [STATE,totalImpinge,impinge,s,beta] = calcCollectionEfficiency(airfoil,fluid,domain,strImpMod)
% Function to calculate the collection efficiency of an airfoil

% Initialize a cloud using domain realization
particles = domain.numParcels;
rhol = fluid.rhol;
cloud = SLDcloud([domain.samples, domain.nDroplet],rhol,particles,'NoTResolve');

% Advect particles (with/without fracture or impingement submodules)
disp('Advecting test particles between limits...');
totalImpinge = 0;
impinge = [];
STATE = {};
t = cloud.tGLOB; simTime = domain.simTime;
tSAMP = 1; tREFRESH = 0.03; tRATE = 20; tMARK = 1;
iter = 1; maxiter = 1500;
figure(1); plot(airfoil.X,airfoil.Y); axis([-.5 1 -.3 .3]);
[xtmp,ytmp] = airfoil.interpStoXY(airfoil.LIMup);
hold on; plot(xtmp,ytmp,'ko');
if strcmp(strImpMod,'NoImpingement')
    % Calculate collection efficiency without impingement module
    while totalImpinge<particles && iter<maxiter
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
        %{
        if mod(iter,100)==0
            figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.5 1 -.3 .3]);
        end
        %
        if floor(tSAMP/tRATE)>=1
            state = cloud.getState();
            part = [state(:,1), state(:,2)];
            %figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.5 1 -.3 .3]);
            figure(1); hold on; scatter(part(:,1),part(:,2),[],'b','filled');
            hold on; plot(xtmp,ytmp,'ko');
            tMARK = iter;
        end
        tSAMP = t-tMARK;
        %}
        iter = iter+1;
        t = iter
    end
    % Record mass deposited on the airfoil surface
    indStick = cloud.impingeTotal;
    xStick = cloud.x(indStick); yStick = cloud.y(indStick); 
    RStick = cloud.rd(indStick); nStick = cloud.numDroplets(indStick);
    sStick = airfoil.XYtoScoords(xStick,yStick);
    mStick = (4/3*pi*rhol)*(RStick.^3).*nStick;
    set(airfoil,'FILM',[sStick, mStick]);
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
    STATE = cloud.getState();
    % Plot airfoil and impinged parcels
    impinge = sort(cloud.impingeTotal);
    xq = cloud.x(impinge);
    yq = cloud.y(impinge);
    s = airfoil.XYtoScoords(xq,yq);
    maxC = max(domain.nDroplet(impinge)); minC = min(domain.nDroplet(impinge));
    C = 9*(domain.nDroplet(impinge)-minC)./(maxC-minC) + 1;
    figure(1); hold on; scatter(STATE(impinge,1),STATE(impinge,2),[],C,'filled');
    % Calculate total impinged mass in each of the surface bins
    bins = 50;
    sBin = linspace(min(sStick),max(sStick),bins+1)';
    ds = sBin(2)-sBin(1);
    for i=1:bins
        mBin(i,1) = sum(mStick(sStick>sBin(i) & sStick<sBin(i+1)));
    end
    % Calculate total mass flux through the injection domain screen
    R = cloud.rd; nDrop = cloud.numDroplets;
    mTOT = sum((4/3*pi*rhol)*(R.^3).*nDrop);
    dy = domain.XY_bounds(1,2) - domain.XY_bounds(2,2);
    % Calculate collection efficiency
    beta = (mBin/ds).*(dy/mTOT);
    s = sBin;
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
%{
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