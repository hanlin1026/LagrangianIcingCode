%% Example solver for thermo eqns (quasilinearization method)

% INPUT ****************************
ds = 1;
s = [0:ds:999]';
pw = 1000;
% Input incoming liquid mass k(s)
mimp = exp(-0.5*(s-mean(s)).^2/25^2);
% Guess ice profile z(s)
Z = 0.5*mimp;
% Input droplet impingement energy terms
Eimp = 0.5*100^2;
% **********************************

dX = 0.1*ones(1000,1); dY = 0.1*ones(1000,1);
for iter = 1:100
    iter
    figure(3); hold on; plot(s,Z); drawnow;
    % MASS EQN *************************
    % Solve mass eqn explicitly (RK4)
    f = mimp - Z;
    steps = 999;
    ya = 0;
    if (iter>1)
        Xprev = X;
        Yprev = Y;
        Zprev = Z;
    end
    [T,X] = rk4([s,f],ya,steps);
    T = T'; X = X';
    % Plot
    figure(1); hold on; plot(T,X); drawnow;
    % **********************************

    % THERMO EQN ***********************
    % Set up mass matrix
    k0 = 2; kf = length(s)-1;
    C = (-1/2/ds)*X(k0:kf);
    D = (1/2/ds)*(X(k0+1:kf+1)-X(k0-1:kf-1)) + 1 - Z(k0:kf);
    E = (1/2/ds)*X(k0:kf);
    A1 = (1/ds)*(X(2)-X(1)) - (1/ds)*X(1) + 1 - Z(1);
    B1 = (1/ds)*X(1);
    Af = (1/ds)*(X(end)-X(end-1)) - (1/ds)*X(end) + 1 - Z(end);
    Bf = (1/ds)*X(end);
    M = diag([A1;D;Bf]) + diag([B1;E],1) + diag([C;Af],-1);
    % Set up RHS
    RHS = Eimp + Z;
    % Solve implicit linear system
    Y = M\RHS;
    % Plot
    figure(2); hold on; plot(T,Y); drawnow;
    % **********************************

    % VARIATION ************************
    DY = zeros(length(Y),1);
    DY(1) = (Y(2)-Y(1))./(T(2)-T(1));
    DY(2:end-1) = (Y(3:end)-Y(1:end-2))./(T(3:end)-T(1:end-2));
    DY(end) = (Y(end)-Y(end-1))./(T(end)-T(end-1));

    DX = zeros(length(X),1);
    DX(1) = (X(2)-X(1))./(T(2)-T(1));
    DX(2:end-1) = (X(3:end)-X(1:end-2))./(T(3:end)-T(1:end-2));
    DX(end) = (X(end)-X(end-1))./(T(end)-T(end-1));
    
%     if (iter > 1)
%         dX = X - Xprev;
%         dY = Y - Yprev;
%     end
    eps = 0.01;
    dZ = ((-DY).*dX + (-DX-1+Z).*dY)./(eps-Y);
    dZ(dZ>2) = 2;
    dZ(dZ<-2) = -2;
    figure(4); hold on; scatter(iter,trapz(s,abs(dZ./dX) + abs(dZ./dY)),'filled'); set(gca,'yscale','log')
    Z = Z - dZ;
    Z(Z<0) = 0;
    Z = smooth(T,Z);
    % **********************************
    
end







%% Example solver (nonlinear iteration)

% INPUT ****************************
ds = 1;
s = [0:ds:999]';
pw = 1000;
uw = 1.787e-3;
% Input incoming liquid mass k(s)
mimp = exp(-0.5*(s-mean(s)).^2/25^2);
% Guess ice profile z(s)
Z = 0.5*mimp;
% Input droplet impingement energy terms
Eimp = 0.5*100^2;
% **********************************
% Set up structure for parameters
scalars.s_ = s;
scalars.ds_ = ds;
scalars.pw_ = pw;
scalars.uw_ = uw;
scalars.mimp_ = mimp;
scalars.Z_ = Z;
scalars.Eimp_ = Eimp;
% Define exact solution
xEXACT = sqrt((2*uw/pw)*cumsum(mimp-Z)*ds);
% Initial guess
x0 = linspace(0,10e-3,length(s))';
% Newton-Krylov iteration
eps = 1e-4;
xn = NewtonKrylovIteration(@MassBalance,scalars,x0,eps);
scalars.X_ = xn;

figure; plot(s,xn);
hold on; plot(s,mimp,'r');
hold on; plot(s,xEXACT,'g');
















