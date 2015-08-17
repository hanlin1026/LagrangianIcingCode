function err = EnergyBalance(x,scalars)
% Energy balance in the form E(x) = 0

% Get relevant parameters
s = scalars.s_;
ds = scalars.ds_;
pw = scalars.pw_;
uw = scalars.uw_;
cw = scalars.cw_;
mimp = scalars.mimp_;
z = scalars.Z_;
x = scalars.X_;

err = zeros(length(x),1);
k = [2:length(x)-1];
% Calculate error for internal cells
err(2:end-1) = (pw*cw/16/uw)*((x(k+1)+x(k)).^2.*(y(k+1)+y(k)) - (x(k)+x(k-1)).^2.*(y(k)+y(k-1)));
% Boundary conditions
err(1) = x(1)-0;
err(end) = err(end-1);


end