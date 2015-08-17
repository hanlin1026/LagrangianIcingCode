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
k = [2:length(x)-1];
% Calculate error for internal cells
err(2:end-1) = (0.5*pw/uw)*0.25*((x(k+1)+x(k)).^2 - (x(k)+x(k-1)).^2) - ds*(mimp(k)-Z(k));
% Boundary conditions
err(1) = x(1)-0;
err(end) = err(end-1);

end