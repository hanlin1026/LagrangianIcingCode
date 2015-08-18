function err = MassBalance(x,scalars)
% Mass balance in form F(x) = err

% Get relevant parameters
s = scalars.s_;
ds = scalars.ds_;
pw = scalars.pw_;
uw = scalars.uw_;
mimp = scalars.mimp_;
Z = scalars.Z_;

err = zeros(length(x),1);
% Calculate body centered fluxes
F = zeros(length(x),1);
F = (0.5*pw/uw)*x.^2;
% Calculate fluxes at cell faces (Roe scheme upwinding)
k = [1:length(x)-1];
f = zeros(length(x)-1,1);
xFACE = 0.5*(x(k)+x(k+1));
f = 0.5*(F(k)+F(k+1)) - 0.5*(pw/uw)*abs(xFACE).*(x(k+1)-x(k));
% Calculate error for internal cells
k = [2:length(x)-1];
err(k) = f(k)-f(k-1) - ds*(mimp(k)-Z(k)); 
% Boundary conditions
err(1) = x(1)-0;
err(end) = err(end-1);

end