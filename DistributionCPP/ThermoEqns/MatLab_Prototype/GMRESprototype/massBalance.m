function err = massBalance(x,scalars)
% Mass balance in form F(x) = err

% Get relevant parameters
s = scalars.s_;
ds = scalars.ds_;
pw = scalars.pw_;
uw = scalars.uw_;
mimp = scalars.mimp_;
Z = scalars.Z_;
tau_wall = scalars.tau_wall_;

err = zeros(length(x),1);
% Calculate body centered fluxes
F = zeros(length(x),1);
F = (0.5/uw)*x.^2.*tau_wall;
% Calculate fluxes at cell faces (Roe scheme upwinding)
k = [1:length(x)-1];
f = zeros(length(x)-1,1);
xFACE = 0.5*(x(k)+x(k+1));
tauFACE = 0.5*(tau_wall(k)+tau_wall(k+1));
DF = (1/uw)*(xFACE.*tauFACE);
f = 0.5*(F(k)+F(k+1)) - 0.5*abs(DF).*(x(k+1)-x(k));
% Calculate error for internal cells
k = [2:length(x)-1];
D_flux = f(k)-f(k-1);
I_sources = (1/pw).*ds(k).*(mimp(k)-Z(k));
err(k) = (D_flux-I_sources);
% Boundary conditions
err(1) = x(1)-0;
err(end) = err(end-1);


end