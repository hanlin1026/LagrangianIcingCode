function [STATE,totalImpinge,impinge,s,beta,sStick,mStick] = calcCollectionEfficiency(airfoil,fluid,domain,strImpMod)
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
iter = 1; maxiter = 1700;
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
        %
        maxC = max(cloud.rd); minC = min(cloud.rd);
        C = 9*(cloud.rd-minC)./(maxC-minC) + 1;
        if mod(iter,100)==0
            figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.02 .04 -.04 .04]);
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
    while totalImpinge<cloud.particles && t<maxiter
        % Call subroutine to calculate local timesteps and impinging particles
        calcDtandImpinge(cloud,airfoil,fluid);
        % Advect particles one time step
        transportSLD(cloud,fluid);
        % Check for fracture
        %fragmentSLD(cloud,fluid);
        % Call subroutine to calculate impingement regimes for impinging particles
        if ~isempty(cloud.impinge)
            impingementRegimeSLD(cloud,airfoil);
        end
        totalImpinge = size(cloud.impingeTotal,1);
        % Save state variables every 0.1 sec
        %
        maxC = max(cloud.rd); minC = min(cloud.rd);
        C = 9*(cloud.rd-minC)./(maxC-minC) + 1;
        if mod(iter,tRATE)==tRATE-1
            figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.02 .075 -.06 .06]);
        end
        %
        if floor(tSAMP/tRATE)>=1
            state = cloud.getState();
            part = [state(:,1), state(:,2)];
            %figure(1); clf; plot(airfoil.X,airfoil.Y); axis([-.5 1 -.3 .3]);
            figure(1); hold on; scatter(part(:,1),part(:,2),10,C,'filled');
            drawnow;
            hold on; plot(xtmp,ytmp,'ko');
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
end

% Calculate collection efficiency
disp('Calculating collection efficiency...');
STATE = cloud.getState();
% Plot airfoil and impinged parcels
impinge = sort(cloud.impingeTotal);
maxC = max(cloud.rd); minC = min(cloud.rd);
C = 9*(cloud.rd-minC)./(maxC-minC) + 1;
figure(11); plot(airfoil.X,airfoil.Y,'k');
figure(11); hold on; scatter(STATE(impinge,1),STATE(impinge,2),[],C(impinge),'filled');
drawnow;
sDrop = airfoil.interpXYtoS(STATE(impinge,1),STATE(impinge,2));
% Calculate total impinged mass in each of the surface bins
bins = 35;
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
figure(12); plot(sDrop,cloud.rd(impinge),'bo'); drawnow;
figure(14); hist(cloud.rd(impinge)); drawnow;

end