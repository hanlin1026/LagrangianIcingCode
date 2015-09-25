function zn = SolveThermoForIceRate(x,y,scalars)
% Solve thermo eqn for temperature

s = scalars.s_;
ds = scalars.ds_;
pw = scalars.pw_;
uw = scalars.uw_;
cw = scalars.cw_;
Td = scalars.Td_;
ud = scalars.ud_;
cice = scalars.cice_;
Lfus = scalars.Lfus_;
ch = scalars.ch_;
mimp = scalars.mimp_;
tau_wall = scalars.tau_wall_;

% Implementation using finite volume with Roe scheme calculation of fluxes
% Calculate body centered fluxes
F = zeros(length(x),1);
F = (0.5*cw/uw)*x.^2.*y.*tau_wall;
% Calculate fluxes at cell faces (Roe scheme upwinding)
k = [1:length(x)-1];
f = zeros(length(x)-1,1);
xFACE = 0.5*(x(k)+x(k+1));
sFACE = 0.5*(s(k)+s(k+1));
tauFACE = 0.5*(tau_wall(k)+tau_wall(k+1));
DF = (cw/2/uw).*tauFACE.*xFACE.^2;
f = 0.5*(F(k)+F(k+1)) - 0.5*abs(DF).*(y(k+1)-y(k));
% Solve discretization for ice accretion rate
k = [2:length(x)-1];
D_flux = f(k)-f(k-1);
dsFACE = sFACE(k)-sFACE(k-1);
RHS = mimp(k)*(cw*Td + 0.5*ud^2) + ch(k).*(Td - y(k));

zn = zeros(length(x),1);
zn(2:end-1) = ((pw./dsFACE).*D_flux - RHS)./(Lfus - cice*y(k));
zn(1) = zn(2);
zn(end) = zn(end-1);

% Implementation using finite differences
%{
XY = x.^2.*y;
DXY = zeros(length(x),1);
DXY(2:end-1) = (XY(3:end)-XY(1:end-2))./(s(3:end)-s(1:end-2));
DXY(1) = (XY(2)-XY(1))./(s(2)-s(1));
DXY(end) = (XY(end)-XY(end-1))./(s(end)-s(end-1));

LHS = (pw*cw/2/uw)*DXY;
RHS = mimp.*(cw*Td + 0.5*ud^2) + ch.*(Td - y);

zn = (LHS-RHS)./(Lfus - cice*y);
%}
end