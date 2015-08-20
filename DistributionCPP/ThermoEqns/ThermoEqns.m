%% Example solver (nonlinear iteration)

% INPUT ****************************
ds = 1;
s = [0:ds:999]';
pw = 1000;
uw = 1.787e-3;
% Input incoming liquid mass k(s)
mimp = exp(-0.5*(s-mean(s)).^2/25^2);
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
scalars.Z_ = Z;
% Set convergence tolerances for water and ice constraints
epsWATER = -2.5e-5;
epsICE = 2.5e-5;
% Iterate on mass/energy eqns until physical solution attained
C_filmPos = true; C_icePos = true; C_waterWarm = false; C_iceCold = false;
iter = 1;
while (((C_filmPos && C_icePos && C_waterWarm && C_iceCold) == false) && (iter < 6) )
    iter
    % MASS (solve for X)
    %x0 = linspace(0,10e-3,length(s))'; % Initial guess
    %eps = 1e-4;
    %xn = NewtonKrylovIteration(@MassBalance,scalars,x0,eps);
    X = sqrt((2*uw/pw)*cumsum(mimp-Z)*ds);
    scalars.X_ = X;
    % Constraint: check that conservation of mass is not violated
    tolIMAG = 1e-6;
    indX = find(abs(imag(X)) > tolIMAG);
    if (~isempty(indX))
        % Lower ice accretion rate in affected areas
        indFix = find(Z(indX) > mimp(indX));
        Z(indX(indFix)) = mimp(indX(indFix));
        scalars.Z_ = Z;
        X = real(X); scalars.X_ = X;
    end
    % ENERGY (solve for Y)
    eps = 1e-4;
    if (iter == 1)
        Y = scalars.Td_*ones(length(s),1);
    end
    Ynew = NewtonKrylovIteration(@EnergyBalance,scalars,Y,eps);
    Y = Ynew; scalars.Y_ = Y;
    % Check constraints
    XY = X.*Y;
    indWATER = find(XY<epsWATER);
    % CONSTRAINTS
    % Water cannot be cold
    if (isempty(indWATER))
        C_waterWarm = true;
    else
        % If we have freezing water, warm it up using epsWATER
        C_waterWarm = false;
        Y(indWATER) = epsWATER./X(indWATER);
        scalars.Y_ = Y;
        % Resolve for ice profile (Z)
        Z = SolveThermoForIceRate(X,Y,scalars);
        scalars.Z_ = Z;
    end
    YZ = Y.*Z;
    % Ice cannot be warm
    indICE = find(YZ>epsICE);
    if (isempty(indICE))
        C_iceCold = true;
    else
        % If we have warm ice, cool it down using epsICE
        C_iceCold = false;
        Y(indWATER) = epsICE./Z(indWATER);
        scalars.Y_ = Y;
        % Resolve for ice profile
        Z = SolveThermoForIceRate(X,Y,scalars);
        scalars.Z_ = Z;
    end
    
    % Plot
    figure(1); hold on; plot(s,X); drawnow;
    figure(2); hold on; plot(s,Y); drawnow;
    figure(3); hold on; plot(s,Z); drawnow;
    iter = iter + 1;
    
end
disp('All compatibility relations satisfied');










