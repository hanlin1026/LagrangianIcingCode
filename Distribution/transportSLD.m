function transportSLD(cloud,fluid)
% Lagrangian advection difference eqns (Villedieu, Guffond and Bobo, AIAA 2012)
% Particles have positions (x,y) and velocities (u,v)
% Gas has velocity (ug,vg) (at the cell locations occupied by the particles)

Tinf = fluid.Tinf;
% Pull out state variables of all particles
x = cloud.x; y = cloud.y; u = cloud.u; v = cloud.v;
rd = cloud.rd; t = cloud.time; dt = cloud.dt; rhol = cloud.rhol;

% Interpolate fluid properties at particle positions
[pg,ug,vg] = fluid.interpFluid(x,y);

% Sutherland's law
C1 = 1.458e-6; % kg/(ms*sqrt(K))
S = 110.4; % K
mug = C1*(Tinf^1.5)/(Tinf+S);

g = [0; -9.81];
Red = 2.*pg.*rd./mug.*sqrt((u-ug).^2 + (v-vg).^2);
CD = 24./Red.*(1+0.15.*(Red.^0.687));
tau = 24./Red./CD.*(2*rhol.*rd.^2./9./mug);
md = 4/3*rhol.*pi.*rd.^3;
if Red==0
    tau = 0;
end

% Advection equations
xnp1 = x + ug.*dt + (u-ug).*(1-exp(-dt./tau)).*tau + (dt-(1-exp(-dt./tau)).*tau).*tau.*g(1);
ynp1 = y + vg.*dt + (v-vg).*(1-exp(-dt./tau)).*tau + (dt-(1-exp(-dt./tau)).*tau).*tau.*g(2);
unp1 = ug + exp(-dt./tau).*(u-ug) + (1-exp(-dt./tau)).*tau.*g(1);
vnp1 = vg + exp(-dt./tau).*(v-vg) + (1-exp(-dt./tau)).*tau.*g(2);

% Update particle states
set(cloud,'x',xnp1); set(cloud,'y',ynp1); set(cloud,'u',unp1); set(cloud,'v',vnp1);
set(cloud,'time',t+dt);

end