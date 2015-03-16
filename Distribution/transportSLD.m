function transportSLD(cloud,fluid)
% Lagrangian advection difference eqns (Villedieu, Guffond and Bobo, AIAA 2012)
% Particles have positions (x,y) and velocities (u,v)
% Gas has velocity (ug,vg) (at the cell locations occupied by the particles)

% Check that there are actually particles to advect
if ~isempty(cloud.indAdv)
    Tinf = fluid.Tinf;
    % Find parcels which are have already entered the injection domain
    indAdv = cloud.indAdv;

    % Pull out state variables of those particles currently in the simulation
    x = cloud.x(indAdv); y = cloud.y(indAdv); u = cloud.u(indAdv); v = cloud.v(indAdv); 
    rd = cloud.rd(indAdv); dt = cloud.dt; rhol = cloud.rhol;

    % Interpolate fluid properties at particle positions
    [pg,ug,vg] = fluid.interpFluid(x,y);

    % Sutherland's law
    C1 = 1.458e-6; % kg/(ms*sqrt(K))
    S = 110.4; % K
    mug = C1*(Tinf^1.5)/(Tinf+S);

    % Force parameter calculations
    g = [0; -9.81];
    Red = 2.*pg.*rd./mug.*sqrt((u-ug).^2 + (v-vg).^2);
    CD = 24./Red.*(1+0.15.*(Red.^0.687));
    tau = 24./Red./CD.*(2*rhol.*rd.^2./9./mug);
    md = 4/3*rhol.*pi.*rd.^3;
    tau(Red==0) = 0;

    % Advection equations
    xnp1 = x + ug.*dt + (u-ug).*(1-exp(-dt./tau)).*tau + (dt-(1-exp(-dt./tau)).*tau).*tau.*g(1);
    ynp1 = y + vg.*dt + (v-vg).*(1-exp(-dt./tau)).*tau + (dt-(1-exp(-dt./tau)).*tau).*tau.*g(2);
    unp1 = ug + exp(-dt./tau).*(u-ug) + (1-exp(-dt./tau)).*tau.*g(1);
    vnp1 = vg + exp(-dt./tau).*(v-vg) + (1-exp(-dt./tau)).*tau.*g(2);

    % Update particle states
    xNEW = cloud.x; xNEW(indAdv) = xnp1; yNEW = cloud.y; yNEW(indAdv) = ynp1;
    uNEW = cloud.u; uNEW(indAdv) = unp1; vNEW = cloud.v; vNEW(indAdv) = vnp1;
    set(cloud,'x',xNEW); set(cloud,'y',yNEW); set(cloud,'u',uNEW); set(cloud,'v',vNEW);
    if cloud.FLAGtimeResolve==1
        dtGLOB = max(dt);
        set(cloud,'tGLOB',cloud.tGLOB+dtGLOB);
    end
end

end