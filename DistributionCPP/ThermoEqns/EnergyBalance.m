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
z = scalars.Z_;
x = scalars.X_;

% Calculate body centered fluxes
F = zeros(length(x),1);
F = (0.5*pw*cw/uw)*x.^2.*y;
% Calculate fluxes at cell faces (Roe scheme upwinding)
k = [1:length(x)-1];
f = zeros(length(x)-1,1);
xFACE = 0.5*(x(k)+x(k+1));
f = 0.5*(F(k)+F(k+1)) - (0.5*pw*cw/uw)*xFACE.^2.*(y(k+1)-y(k));
% Calculate error for internal cells
err = zeros(length(x),1);
k = [2:length(x)-1];
df = f(k)-f(k-1);
RHS = mimp(k)*(cw*Td + 0.5*ud^2) + z(k).*(Lfus - cice.*y(k)) + ch*(Td - y(k));
err(k) = df - ds*RHS;

%LHS = (pw*cw/16/uw)*((x(k+1)+x(k)).^2.*(y(k+1)+y(k)) - (x(k)+x(k-1)).^2.*(y(k)+y(k-1)));
%RHS = mimp(k)*(cw*Td + 0.5*ud^2) + z(k).*(Lfus - cice.*y(k)) + ch*(Td - y(k));
%err(2:end-1) = LHS - ds*RHS;
% Boundary conditions
if (z(1) < 1e-4)
    err(1) = y(1)-Td;
else
    err(1) = y(1) - 0;
end
err(end) = err(end-1);


end