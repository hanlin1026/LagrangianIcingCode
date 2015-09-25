%% Example solver (nonlinear iteration)

% INPUT ****************************
% Import skin friction coefficient data
CF = importdata('../heatflux');
% Piecewise interpolation *************
indDisc = find(abs(diff(CF(:,3)))/mean(abs(diff(CF(:,3)))) > .5e2);
stagPt = 1.008358;
s_min = stagPt; s_max = stagPt + 0.4;
[val,indFirst] = min(abs(CF(:,1)-s_min));
[val,indLast] = min(abs(CF(:,1)-s_max));
indDisc = [indFirst-1;indDisc;indLast];
NPts = 1000;
s = []; tau_wall = []; ch = [];

for i=1:(length(indDisc)-1)
    startpt = indDisc(i)+1;
    breakpt = indDisc(i+1);
    if startpt == breakpt
        s = [s; CF(startpt,1)];
        ch = [ch; CF(startpt,2)];
        tau_wall = [tau_wall; CF(startpt,3)];
    else
        npts = NPts*(breakpt-startpt)/(indLast-indFirst);
        snew = linspace(CF(startpt,1),CF(breakpt,1),npts)';
        s = [s; snew];
        ch = [ch; interp1(CF(startpt:breakpt,1),CF(startpt:breakpt,2),snew)];
        tau_wall = [tau_wall; interp1(CF(startpt:breakpt,1),CF(startpt:breakpt,3),snew)];
    end
end
s = s - stagPt;
% **************************************
%tau_wall(abs(tau_wall) < 0.05*max(abs(tau_wall))) = 0;
ch = -1*1.412*(1.01e5/1.412)^1.5/(273.15-250)*ch;
ds = zeros(length(s),1);
ds(1:end-1) = diff(s); ds(end) = s(end)-s(end-1);
pw = 1000;
uw = 1.787e-3;
% Input incoming liquid mass k(s)
BETA = importdata('../BetaXY.dat',',');
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
scalars.Td_ = -1;
scalars.ud_ = 80;
scalars.cice_ = 2093; % J/(kg C) at T = 0
scalars.Lfus_ = 334774; % J/kg
scalars.ch_ = ch; % W/(m^2 C)
scalars.mimp_ = mimp;
scalars.tau_wall_ = tau_wall;
scalars.Z_ = Z;
% Set convergence tolerances for water and ice constraints
epsWATER = -1e-4;
epsICE = 1e-4;
% Iterate on mass/energy eqns until physical solution attained
C_filmPos = true; C_icePos = true; C_waterWarm = false; C_iceCold = false;
iter = 1;
%% Solve

%u0 = linspace(0,10e-3,1000)';
%u0 = 1e-3*ones(1000,1);
%X = NewtonKrylovIteration(@massBalance,@JX,u0,scalars);
X = sqrt((2*uw/pw./tau_wall).*cumtrapz(s,mimp-Z));
scalars.X_ = X;
scalars = correctFilmHeight(scalars); 
X = scalars.X_;
% ENERGY (solve for Y)
eps = 1e-4;
if (iter == 1)
    Y = scalars.Td_*ones(length(s),1);
end
Ynew = NewtonKrylovIteration(@EnergyBalance,@JX,Y,scalars);















