function [XY_bounds] = calcInjectionDomain(dx,Ravg,fluid,airfoil)
% Function to calculate a box which contains the droplet clusters
% INPUTS: 
%   dx: length scale (desired domain length in the x-direction)
%   Ravg: average droplet radius
%   fluid,airfoil: fluid and airfoil objects
% OUTPUTS:
%   XY_bounds: xy coordinates of domain boundaries

% Determine impingement limits in y-direction at x0
xL = -0.5; Yhit = -0.1;
Ymiss = 0.3;
yLimUP = impingementLimitsSLD(Ravg,fluid,airfoil,xL,Ymiss,Yhit,'UP');
Ymiss = -0.3;
yLimDOWN = impingementLimitsSLD(Ravg,fluid,airfoil,xL,Ymiss,Yhit,'DOWN');

% Initialize cloud object of 2 particles placed at impingement limits
x0 = [xL; xL];
y0 = [yLimUP; yLimDOWN];
[pg,ug,vg] = interpFluid(fluid,x0,y0);
u0 = ug + 0.01*ug;
v0 = vg + 0.01*vg;
rd0 = [Ravg; Ravg];
time0 = [0; 0];
particles = 2;
cloud = SLDcloud([x0 y0 u0 v0 rd0 time0 [1:particles]'],fluid.rhol,particles);

% Advect particles backwards for time dx/U
Uinf = fluid.Uinf;
tStar = -dx/Uinf;
while max(cloud.time) > tStar
    calcDtandImpinge(cloud,airfoil,fluid);
    set(cloud,'dt',-cloud.dt);
    transportSLD(cloud,fluid);
end

% Save new boundaries of particles
XY_bounds = [x0 y0; cloud.x cloud.y];

end

