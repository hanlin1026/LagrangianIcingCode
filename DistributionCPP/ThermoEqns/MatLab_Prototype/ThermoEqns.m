%% Example solver (nonlinear iteration)

% INPUT ****************************
% Import skin friction coefficient data
CF = importdata('SkinFrictionXY.dat',',');
CF(:,2) = 1;
% Piecewise interpolation *************
indDisc = find(abs(diff(CF(:,2)))/mean(abs(diff(CF(:,2)))) > .5e2);
s_min = 0; s_max = 0.4;
[val,indFirst] = min(abs(CF(:,1)-s_min));
[val,indLast] = min(abs(CF(:,1)-s_max));
indDisc = [indFirst-1;indDisc;indLast];
NPts = 1000;
s = []; tau_wall = [];

for i=1:(length(indDisc)-1)
    startpt = indDisc(i)+1;
    breakpt = indDisc(i+1);
    if startpt == breakpt
        s = [s; CF(startpt,1)];
        tau_wall = [tau_wall; CF(startpt,2)];
    else
        npts = NPts*(breakpt-startpt)/(indLast-indFirst);
        snew = linspace(CF(startpt,1),CF(breakpt,1),npts)';
        s = [s; snew];
        tau_wall = [tau_wall; interp1(CF(startpt:breakpt,1),CF(startpt:breakpt,2),snew)];
    end
end
% **************************************
tau_wall(abs(tau_wall) < 0.05*max(abs(tau_wall))) = 0;
ds = zeros(length(s),1);
ds(1:end-1) = diff(s); ds(end) = s(end)-s(end-1);
pw = 1000;
uw = 1.787e-3;
% Input incoming liquid mass k(s)
BETA = importdata('BetaXY.dat',',');
Uinf = 100;
LWC = 1;
mimp = Uinf*LWC*interp1(BETA(:,1),BETA(:,2),s);
mimp(isnan(mimp)) = 0;
% Guess ice profile z(s)
Z = 0*mimp;
% **********************************
% Set up structure for parameters
scalars.s_ = s;
scalars.ds_ = ds;
scalars.pw_ = pw;
scalars.uw_ = uw;
scalars.cw_ = 4217.6; % J/(kg C) at T = 0 C and P = 100 kPa
scalars.Td_ = -20;
scalars.ud_ = 80;
scalars.cice_ = 2093; % J/(kg C) at T = 0
scalars.Lfus_ = 334774; % J/kg
scalars.ch_ = 100; % W/(m^2 C)
scalars.mimp_ = mimp;
scalars.tau_wall_ = tau_wall;
scalars.Z_ = Z;
% Set convergence tolerances for water and ice constraints
epsWATER = -1e-4;
epsICE = 1e-4;
% Iterate on mass/energy eqns until physical solution attained
C_filmPos = true; C_icePos = true; C_waterWarm = false; C_iceCold = false;
iter = 1;
%%
%while (((C_filmPos && C_icePos && C_waterWarm && C_iceCold) == false) && (iter < 11) )
while (iter<5)
    iter
    % MASS (solve for X)
    %x0 = linspace(0,10e-3,length(s))'; % Initial guess
    %eps = 1e-4;
    %X = NewtonKrylovIteration(@MassBalance,scalars,x0,eps);
    %{
    x0 = zeros(length(s),1);
    X = x0; epsT = 1e-8; iterMASS = 1; ERR = 1;
    DX = MassBalance(X,scalars);
    X = X - epsT*DX;
    X(X<0) = 0;
    ERR0 = max(abs(DX));
    while ((ERR > 1e-2*ERR0) && (iterMASS < 20000))
        iterMASS = iterMASS+1;
        DX = MassBalance(X,scalars);
        X = X - epsT*DX;
        X(X<0) = 0;
        ERR = max(abs(DX));
    end
    iterMASS
    %}
    X = sqrt((2*uw/pw./tau_wall).*cumtrapz(s,mimp-Z));
    scalars.X_ = X;
    % Constraint: check that conservation of mass is not violated
    scalars = correctFilmHeight(scalars); 
    X = scalars.X_;
    % ENERGY (solve for Y)
    eps = 1e-4;
    if (iter == 1)
        Y = scalars.Td_*ones(length(s),1);
    end
    Ynew = NewtonKrylovIteration(@EnergyBalance,scalars,Y,eps);
    scalars.Y_ = Y;
    % Check constraints
    XY = X.*Y;
    indWATER = find(XY<epsWATER);
    % CONSTRAINTS
    % Water cannot be cold
    if (isempty(indWATER))
        C_waterWarm = true;
    else
        disp('Water below freezing detected');
        % If we have freezing water, warm it up using epsWATER
        C_waterWarm = false;
        Y(indWATER) = epsWATER./X(indWATER);
        scalars.Y_ = Y;
        % Resolve for ice profile (Z)
        Z = SolveThermoForIceRate(X,Y,scalars);
        Z(Z<0) = 0;
        scalars.Z_ = Z;
    end
    YZ = Y.*Z;
    % Ice cannot be warm
    indICE = find(YZ>epsICE);
    if (isempty(indICE))
        C_iceCold = true;
    else
        disp('Ice above freezing detected');
        % If we have warm ice, cool it down using epsICE
        C_iceCold = false;
        Y(indWATER) = epsICE./Z(indWATER);
        scalars.Y_ = Y;
        % Resolve for ice profile
        Z = SolveThermoForIceRate(X,Y,scalars);
        Z(Z<0) = 0;
        scalars.Z_ = Z;
    end
    
    % Plot
    figure(11); hold on; plot(s,X); drawnow;
    figure(12); hold on; plot(s,Y); drawnow;
    figure(13); hold on; plot(s,Z); drawnow;
    iter = iter + 1;
    
end
disp('All compatibility relations satisfied');










