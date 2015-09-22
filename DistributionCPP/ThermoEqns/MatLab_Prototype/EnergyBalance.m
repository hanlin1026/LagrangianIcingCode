function err = EnergyBalance(y,scalars)
% Energy balance in the form E(x) = 0

% Get relevant parameters
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
z = scalars.Z_;
x = scalars.X_;

% Calculate body centered fluxes
F = zeros(length(x),1);
F = (0.5*cw/uw)*x.^2.*y.*tau_wall;
% Calculate fluxes at cell faces (Roe scheme upwinding)
k = [1:length(x)-1];
f = zeros(length(x)-1,1);
xFACE = 0.5*(x(k)+x(k+1));
tauFACE = 0.5*(tau_wall(k)+tau_wall(k+1));
DF = (cw/2/uw).*tauFACE.*xFACE.^2;
f = 0.5*(F(k)+F(k+1)) - 0.5*abs(DF).*(y(k+1)-y(k));
% Calculate error for internal cells
err = zeros(length(x),1);
k = [2:length(x)-1];
D_flux = f(k)-f(k-1);
RHS = (1/pw)*(mimp(k)*(cw*Td + 0.5*ud^2) + z(k).*(Lfus - cice.*y(k)) + ch(k).*(Td - y(k)));
I_sources = ds(k).*RHS;
err(k) = D_flux-I_sources;

% Boundary conditions
if (z(1) < 0.01*max(z))
    err(1) = y(1)-Td;
else
    err(1) = y(1) - 0;
end
err(end) = err(end-1);


end