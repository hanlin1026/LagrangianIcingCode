function impingementRegimeSLD(cloud,airfoil)
% Function to compute the impingement regime for those drops which have
% struck the airfoil surface. Returns indices for each of three cases.
% INPUTS: SLDcloud object, gas velocity and density at particle positions, Airfoil object


cloud.computeImpingementParams(airfoil);
%{
% CALCULATE IMPINGEMENT REGIMES *******************************************
% Pull out state variables of all particles which have impinged
indimp = cloud.impinge;
x = cloud.x(indimp); y = cloud.y(indimp); u = cloud.u(indimp); v = cloud.v(indimp);
rd = cloud.rd(indimp); t = cloud.time(indimp); dt = cloud.dt(indimp); rhol = cloud.rhol;

sigma = 75.64e-3; % Surface tension of water at 0 deg C against air
T = 270; % Temperature (K) for SLD's (estimate, slightly below freezing)
mul = (2.414e-5)*10^(247.8/(T-140)); % Estimate of mu_water (N*s/m^2)

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
% Calculate impingement regime
ind1 = find(K < Kb0*fb);
ind2 = find((K > Kb0*fb) & (K < Ks0*fs));
ind3 = find(K > Ks0*fs);
bounce = indimp(ind1);
spread = indimp(ind2);
splash = indimp(ind3);
% Set index trackers in the cloud
set(cloud,'bounce',[]); set(cloud,'bounce',cloud.index(bounce));
set(cloud,'spread',[]); set(cloud,'spread',cloud.index(spread));
set(cloud,'splash',[]); set(cloud,'splash',cloud.index(splash));
set(cloud,'normvel',[]); set(cloud,'normvel',vnormsq);
% END IMPINGEMENT REGIME **************************************************
%}


% CALCULATE SPLASHING MODE DYNAMICS ***************************************
indSplash = cloud.splash; % Indices of cloud.impinge which have splashed
indStateSplash = cloud.impinge(indSplash); % Splash indices for state variables
indexTrackSplash = cloud.index(indStateSplash); % Index trackers

origSplash = [];
if ~isempty(indSplash)
    % Pull out quantities for the splashing mode
    xSplash = cloud.x(indStateSplash); ySplash = cloud.y(indStateSplash); uSplash = cloud.u(indStateSplash); vSplash = cloud.v(indStateSplash);
    rdSplash = cloud.rd(indStateSplash); tSplash = cloud.time(indStateSplash); rhol = cloud.rhol;
    
    KSplash = cloud.K(indSplash); KsSplash = cloud.fs(indSplash)*cloud.Ks0;
    vnormsqSplash = cloud.normvelsq(indSplash); vtangSplash = cloud.tangvel(indSplash);
    [~,~,~,~,txSplash,tySplash] = airfoil.findPanel(xSplash,ySplash);
    sCoordSplash = airfoil.XYtoScoords(xSplash,ySplash);
    % Calculate impingement mass loss parameters for splashing mode
    TH = abs(pi/2 - airfoil.findTH(xSplash,ySplash,uSplash,vSplash));
    a = 1-0.3*sin(TH);
    b = (1/8)*(1+3*cos(TH));
    % Calculate splashing ejection mass (ms) and sticking mass (mStick)
    ms_m0 = max(a - (KsSplash./KSplash).^b,0);
    m0 = (4/3)*pi*rdSplash.^3;
    ms = ms_m0.*m0;
    mStick = rhol*(m0 - ms);
    % Calculate splashed droplet size
    rnew = {};
    state = [];
    NUMNEWTOTAL = 0;
    for i=1:length(indSplash)
        if ms(i) ~= 0
            % Interpolate analytical expression for the CDF to get droplet size
            var = 0.2;
            A0 = 0.09; A1 = 0.51; delK = 1500;
            rm_rd = A0 + A1*exp(-KSplash(i)/delK);
            rm = rm_rd*rdSplash(i);
            mu = log(rm);
            dropsize = linspace(0.05*rm,rdSplash(i),1000)';
            CUMDIST = 0.5 + 0.5*erf((1/sqrt(2*var))*(log(dropsize)-mu));
            mCHILD = 0; 
            mPARENT = ms(i)*rhol;
            numnew = 0;
            % Draw child particles until mass is conserved
            while mCHILD < mPARENT
                numnew = numnew+1;
                dsamp = unifrnd(0.05,0.95,[1,1]);
                rnew{i}(numnew,1) = interp1(CUMDIST,dropsize,dsamp);
                mCHILD = mCHILD + 4/3*pi*(rnew{i}(numnew,1))^3*rhol;
            end
            % Calculate post splashing droplet velocities
            % Interpolate analytical expression for the CDF to get magnitude of v2 for
            % splash droplets, where v_new = v1 + v2
            vratio = linspace(0,1,1000)';
            CUMDIST = 1 - exp(-13.7984.*vratio.^2.5);
            xsamp = unifrnd(0,1,[numnew,1]);
            vrat = interp1(CUMDIST,vratio,xsamp);
            v2mag = vrat*sqrt(vnormsqSplash(i));
            % Calculate elevation (relative to surface tangent) of splash droplet rebounds
            e1 = unifrnd(0,25*pi/180,[numnew,1]);
            e2 = unifrnd(pi-25*pi/180,pi,[numnew,1]);
            elevation = randsample([e1;e2],numnew);
            foilAngle = atan2(tySplash(i),txSplash(i));
            v2 = [v2mag,v2mag].*[cos(foilAngle+elevation), sin(foilAngle+elevation)];
            % Tangential splashing
            v1panelframe = 0.8*vtangSplash(i);
            v1 = v1panelframe*[cos(foilAngle), sin(foilAngle)];
            % Total splashed velocity = v1 + v2
            vSplashDrop = repmat(v1,numnew,1) + v2;
            % Save state of new splashed child particle
            sx = repmat(xSplash(i),numnew,1);
            sy = repmat(ySplash(i),numnew,1);
            sr = rnew{i}(:);
            st = repmat(tSplash(i),numnew,1);
            state = [state; [sx sy vSplashDrop(:,1) vSplashDrop(:,2) sr st]];
            NUMNEWTOTAL = NUMNEWTOTAL + numnew;
        else
            % No splashing actually occurs, per model limits
            
        end
        % Add mass which has "stuck" to airfoil
        set(airfoil,'FILMsplash',[sCoordSplash(i), mStick(i)]);
    end
    if ~isempty(state)
        % Add new splashed child particles
        indS1 = max(cloud.index)+1;
        indS2 = indS1+NUMNEWTOTAL-1;
        indSnew = [indS1:1:indS2]';
        cloud.addParticle(state);
        set(cloud,'parentind',indexTrackSplash);
        set(cloud,'childind',indSnew);
        % Record splashing of original droplets
        indtmp = find(indexTrackSplash <= cloud.originalNumParticles);
        origSplash = sCoordSplash(indtmp);
        set(airfoil,'originalImpingeScoordSplash',origSplash);
    end
end
% END SPLASHING MODE DYNAMICS *********************************************

% CALCULATE BOUNCE MODE DYNAMICS ******************************************
indBounce = cloud.bounce; % Indices of cloud.impinge which have bounced
indStateBounce = cloud.impinge(indBounce); % Bounce indices for state variables
indexTrackBounce = cloud.index(indStateBounce); % Index trackers

if ~isempty(indBounce)
    % Pull out quantities for the bounce mode
    xBounce = cloud.x(indStateBounce); yBounce = cloud.y(indStateBounce); uBounce = cloud.u(indStateBounce); vBounce = cloud.v(indStateBounce);
    rdBounce = cloud.rd(indStateBounce); rhol = cloud.rhol;
    
    KBounce = cloud.K(indBounce); KsBounce = cloud.fs(indBounce)*cloud.Ks0; vnormsqBounce = cloud.normvelsq(indBounce);
    [~,~,nxBounce,nyBounce,txBounce,tyBounce] = airfoil.findPanel(xBounce,yBounce);
    KbBounce = cloud.fb(indBounce)*cloud.Kb0;
    vnormBounce = sqrt(vnormsqBounce); vtangBounce = cloud.tangvel(indBounce);
    % Calculate post impact velocity for bouncing droplets
    vn = 4*vnormBounce.*(sqrt(KBounce./KbBounce)-KBounce./KbBounce);
    vt = 0.8*vtangBounce;
    vnew = [vn,vn].*[nxBounce,nyBounce] + [vt,vt].*[txBounce,tyBounce];
    % Reset velocities of bouncing droplets
    set(cloud,'u',[vnew(:,1), indStateBounce]);
    set(cloud,'v',[vnew(:,2), indStateBounce]);
end
% END BOUNCE MODE DYNAMICS ************************************************

% CALCULATE SPREAD MODE DYNAMICS ******************************************
% Treat as pure "sticking" and simply add mass to airfoil film
indSpread = cloud.spread;
indStateSpread = cloud.impinge(indSpread);
indexTrackSpread = cloud.index(indStateSpread);
origSpread = [];
if ~isempty(indSpread)
    pxSpread = cloud.x(indStateSpread); pySpread = cloud.y(indStateSpread);
    rdSpread = cloud.rd(indStateSpread);
    rhol = cloud.rhol;
    mSpread = (4/3)*pi*rdSpread.^3*rhol;
    sCoordSpread = airfoil.XYtoScoords(pxSpread,pySpread);
    set(airfoil,'FILMspread',[sCoordSpread, mSpread]);
    % Record splashing of original droplets
    indtmp = find(indexTrackSpread <= cloud.originalNumParticles);
    origSpread = sCoordSpread(indtmp);
    set(airfoil,'originalImpingeScoordSpread',origSpread);
end
% END SPREAD MODE DYNAMICS ************************************************


% Delete old parent splashed and spread droplets
set(airfoil,'FILM',[airfoil.FILMsplash; airfoil.FILMspread]);
set(airfoil,'originalImpingeScoord',[airfoil.originalImpingeScoordSplash; origSpread]);
if ~isempty(indStateSplash)
    cloud.deleteParticle(indStateSplash);
end

end