function splashDynamics(cloud,airfoil)

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



end