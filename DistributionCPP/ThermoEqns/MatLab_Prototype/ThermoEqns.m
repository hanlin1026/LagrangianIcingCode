%% Example solver (nonlinear iteration)

% INPUT ****************************
SURF = 'LOWER';
% Import skin friction coefficient data
CF = importdata('heatflux');
stagPt = 1.008358;
NPts = 1000;
if (strcmp(SURF,'LOWER')==1)
    s_min = stagPt - 0.4; s_max = stagPt + 0.0;
    [val,indFirst] = min(abs(CF(:,1)-s_min));
    [val,indLast] = min(abs(CF(:,1)-s_max));
    % Find stagPt w.r.t. skin friction
    [val,ind] = min(CF(indLast-5:indLast+5,3));
    indLast = indLast-5+(ind-1);
    stagPt = CF(indLast,1);
elseif (strcmp(SURF,'UPPER')==1)
    s_min = stagPt; s_max = stagPt + 0.4;
    [val,indFirst] = min(abs(CF(:,1)-s_min));
    [val,indLast] = min(abs(CF(:,1)-s_max));
    % Find stagPt w.r.t. skin friction
    [val,ind] = min(CF(indFirst-5:indFirst+5,3));
    indFirst = indFirst-5+(ind-1);
    stagPt = CF(indFirst,1);
end
% Interpolation
s = linspace(CF(indFirst,1),CF(indLast,1),NPts)';
ch = interp1(CF(indFirst:indLast,1),CF(indFirst:indLast,2),s);
tau_wall = interp1(CF(indFirst:indLast,1),CF(indFirst:indLast,3),s);
s = s - stagPt;
% Piecewise interpolation *************
% indDisc = find(abs(diff(CF(:,3)))/mean(abs(diff(CF(:,3)))) > .5e2);
% stagPt = 1.008358;
% s_min = stagPt; s_max = stagPt + 0.4;
% [val,indFirst] = min(abs(CF(:,1)-s_min));
% [val,indLast] = min(abs(CF(:,1)-s_max));
% indDisc = [indFirst-1;indDisc;indLast];
% NPts = 1000;
% s = []; tau_wall = []; ch = [];
% 
% for i=1:(length(indDisc)-1)
%     startpt = indDisc(i)+1;
%     breakpt = indDisc(i+1);
%     if startpt == breakpt
%         s = [s; CF(startpt,1)];
%         ch = [ch; CF(startpt,2)];
%         tau_wall = [tau_wall; CF(startpt,3)];
%     else
%         npts = NPts*(breakpt-startpt)/(indLast-indFirst);
%         snew = linspace(CF(startpt,1),CF(breakpt,1),npts)';
%         s = [s; snew];
%         ch = [ch; interp1(CF(startpt:breakpt,1),CF(startpt:breakpt,2),snew)];
%         tau_wall = [tau_wall; interp1(CF(startpt:breakpt,1),CF(startpt:breakpt,3),snew)];
%     end
% end
% s = s - stagPt;
% **************************************
ch = -1*1.412*(1.01e5/1.412)^1.5/(273.15-250)*ch;
ds = zeros(length(s),1);
ds(1:end-1) = diff(s); ds(end) = s(end)-s(end-1);
pw = 1000;
uw = 1.787e-3;
% Input incoming liquid mass k(s)
BETA = importdata('BetaXY.dat',',');
Uinf = 100;
LWC = 0.55e-3;
mimp = Uinf*LWC*interp1(BETA(:,1),BETA(:,2),s);
mimp(isnan(mimp)) = 0;
% Scaling for skin friction
tau_wall = tau_wall.*ds*0.5*1.412*Uinf^2;
% Guess ice profile z(s)
Z = 0*mimp;
if (strcmp(SURF,'LOWER')==1)
    s = -flipud(s);
    ch = flipud(ch);
    tau_wall = flipud(tau_wall);
    mimp = flipud(mimp);
end
% **********************************
% Set up structure for parameters
scalars.s_ = s;
scalars.ds_ = ds;
scalars.pw_ = pw;
scalars.uw_ = uw;
scalars.cw_ = 4217.6; % J/(kg C) at T = 0 C and P = 100 kPa
scalars.Td_ = -12;
scalars.ud_ = 80;
scalars.cice_ = 2093; % J/(kg C) at T = 0
scalars.Lfus_ = 334774; % J/kg
scalars.ch_ = ch; % W/(m^2 C)
scalars.mimp_ = mimp;
scalars.tau_wall_ = tau_wall;
scalars.Z_ = Z;
% Set convergence tolerances for water and ice constraints
epsWATER = -1e-8;
epsICE = 1e-8;
% Iterate on mass/energy eqns until physical solution attained
C_filmPos = true; C_icePos = true; C_waterWarm = false; C_iceCold = false;
iter = 1;
%%
%while (((C_filmPos && C_icePos && C_waterWarm && C_iceCold) == false) && (iter < 11) )
con = 1;
figure(13); plot(s,mimp,'k--');
while ((iter<3) )
    iter
    con = 0;
    % MASS (solve for X)
    %x0 = linspace(0,10e-3,length(s))'; % Initial guess
    %eps = 1e-4;
    %X = NewtonKrylovIteration(@MassBalance,scalars,x0,eps);
    X = sqrt((2*uw/pw./tau_wall).*cumtrapz(s,mimp-Z));
    %{
    if (iter==1)
        X = zeros(length(s),1);
    end
    [X,err,con] = explicitSolver(@MassBalance,X,scalars,1e2,1e-4,'MASS');
    if (con == 1)
        disp('CONSTRAINT: Mass balance violated (negative film height)');
    end
    %}
    scalars.X_ = X;
    % Constraint: check that conservation of mass is not violated
    scalars = correctFilmHeight(scalars);
    X = scalars.X_; Z = scalars.Z_;
    % ENERGY (solve for Y)
    eps = 1e-4;
    if (iter == 1)
        %Y = scalars.Td_*ones(length(s),1);
        Y = 0*ones(length(s),1);
    end
    %Ynew = NewtonKrylovIteration(@EnergyBalance,scalars,Y,eps); Y = Ynew;
    [Y,~,~] = explicitSolver(@EnergyBalance,Y,scalars,1e1,1e-4,'ENERGY');
    scalars.Y_ = Y;
    % CONSTRAINTS
    %
    YZ = Y.*Z;
    % Ice cannot be warm
    indICE = find(YZ>100*epsICE);
    if (isempty(indICE))
        C_iceCold = true;
    else
        disp('CONSTRAINT: Ice above freezing detected');
        con = 1;
        % If we have warm ice, cool it down using epsICE
        C_iceCold = false;
        Z(indICE) = epsICE./Y(indICE);
        %scalars.Z_ = Z;
        % Resolve for ice profile
        Z(Z<0) = 0;
        scalars.Z_ = Z;
    end
    %}
    
    %
    XY = X.*Y;
    indWATER = find(XY<100*epsWATER);
    % Water cannot be cold
    if (isempty(indWATER))
        C_waterWarm = true;
    else
        disp('CONSTRAINT: Water below freezing detected');
        con = 1;
        % If we have freezing water, warm it up using epsWATER
        C_waterWarm = false;
        Ytmp = Y;
        Ytmp(indWATER) = 0;
        Z = SolveThermoForIceRate(X,Ytmp,scalars);
        scalars.Z_ = Z;
    end
    scalars.Z_ = Z;
    %}
    
    %}
    % Plot
    figure(1); hold on; plot(s,X); drawnow;
    figure(2); hold on; plot(s,Y); drawnow;
    figure(3); hold on; plot(s,Z); drawnow;
    iter = iter + 1;
    
end
%%
if (con==0)
    disp('All compatibility relations satisfied');
else
    % Check to see if need refinement of ice profile
    % Need to work out exactly where glaze accretion stops and rime begins
    % Yupper = glaze ice everywhere; Ylower = rime profile everywhere
    massSurplus = trapz(s,mimp-Z);
    Yzero = zeros(length(s),1);
    Zupper = SolveThermoForIceRate(X,Yzero,scalars);
    Zlower = mimp;
    if (massSurplus<0)
        indLOWER = find(Z==Zupper);
        % Water mass deficit (too much ice)
        % Iteratively convert glaze to rime accretion until mass balance is fixed
        % Start at current glaze/rime interface and iteratively march forward
        for i=length(indLOWER):-1:1
            indSTMP = indLOWER(i);
            Ztmp = Z; Ytmp = Y;
            Ztmp(indSTMP:end) = Zlower(indSTMP:end);
            %Ytmp(indSTMP:end) = 0;
            massSurplus = trapz(s,mimp-Ztmp);
            figure(3); hold on; plot(s,Ztmp,'k');
            if (massSurplus>=0)
                X = real(sqrt((2*uw/pw./tau_wall).*cumtrapz(s,mimp-Z)));
                X(indSTMP:end) = 0;
                scalars.X_ = X;
                scalars.Z_ = Ztmp;
                [Ytmp,~,~] = explicitSolver(@EnergyBalance,Y,scalars,1e1,1e-4,'ENERGY');
                scalars.Y_ = Ytmp;
                break;
            end
        end
    else
        indRAISE = [find(Z==Zlower); find(Z<1e-4)];
        if (indRAISE(1)==1)
            indRAISE = indRAISE(2:end);
        end
        % Water mass surplus (too much water)
        % Iteratively convert rime to glaze accretion until mass balance is fixed
        % Start at current glaze/rime intergace and iteratively march aft
        for i=1:length(indRAISE)
            indSTMP = indRAISE(i);
            Ztmp = Z; Ytmp = Y;
            Ztmp(1:indSTMP) = Zupper(1:indSTMP);
            Ytmp(1:indSTMP) = 0;
            massSurplus = trapz(s,mimp-Ztmp);
            figure(3); hold on; plot(s,Ztmp,'k');
            if (massSurplus<=0)
                X = real(sqrt((2*uw/pw./tau_wall).*cumtrapz(s,mimp-Z)));
                break;
            end
        end
    end
    % Accept final solution
    %
    Z = Ztmp; Y = Ytmp;
    %}
end

if (strcmp(SURF,'LOWER')==1)
    s = -flipud(s);
    ch = flipud(ch);
    tau_wall = flipud(tau_wall);
    mimp = flipud(mimp);
    X = flipud(X);
    Y = flipud(Y);
    Z = flipud(Z);
end






