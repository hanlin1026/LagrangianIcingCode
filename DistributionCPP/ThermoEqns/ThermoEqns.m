%% Example solver (nonlinear iteration)

% INPUT ****************************
ds = 1;
s = [0:ds:999]';
pw = 1000;
uw = 1.787e-3;
% Input incoming liquid mass k(s)
mimp = exp(-0.5*(s-mean(s)).^2/25^2);
% Guess ice profile z(s)
Z = 0.1*mimp;
% Input droplet impingement energy terms
Eimp = 0.5*100^2;
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
scalars.Z_ = Z;
% Define exact solution
xEXACT = sqrt((2*uw/pw)*cumsum(mimp-Z)*ds);
% MASS
x0 = linspace(0,10e-3,length(s))'; % Initial guess
eps = 1e-4;
xn = NewtonKrylovIteration(@MassBalance,scalars,x0,eps);
scalars.X_ = xn;
% ENERGY
%
y0 = scalars.Td_*ones(length(s),1);
yn = NewtonKrylovIteration(@EnergyBalance,scalars,y0,eps);
%}
X = xn; Y = yn;
figure(1); plot(s,X); hold on; plot(s,xEXACT,'g');
figure(2); plot(s,Y);
figure(3); plot(s,Z);

XY = X.*Y;
YZ = Y.*Z;
if (isempty(X(X<-1e-10)) && isempty(Z(Z<-1e-10)) && isempty(XY(XY<-1e-10)) && isempty(YZ(YZ>1e-10)))
    disp('All compatibility relations satisfied');
end












