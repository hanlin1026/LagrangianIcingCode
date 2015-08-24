function Y = IntegrateEnergyEqn(scalars)
% Function to explicitly integrate energy eqn for temperature Y(s)

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
tau_wall = 1;

% Define parameters
K1 = pw*cw/2/uw;
K2 = tau_wall.*(x.^2);
K3 = -cice*z - ch;
K4 = mimp.*(cw*Td + (ud.^2)/2) + Lfus*z + ch*Td;

DK2 = zeros(length(s),1);
DK2(2:end) = (1/ds)*diff(K2);
%DK2(2:end-1) = (1/2/ds)*(K2(3:end)-K2(1:end-2));
%DK2(1) = (1/ds)*(K2(2)-K2(1));
%DK2(end) = (1/ds)*(K2(end)-K2(end-1));
C1 = zeros(length(s),1);
C2 = zeros(length(s),1);
% Correction for locations where K2 = 0
tolK2 = 0.01;
indK2 = find(abs(K2)./max(abs(K2)) < tolK2);
K2(indK2) = tolK2*max(abs(K2))*sign(K2(indK2));
ind2 = find(K2 == 0);
K2(ind2) = tolK2*max(abs(K2));
C1 = K3./K1./K2 - DK2./K2;
C2 = K4./K1./K2;

% Compute integrals
integralC1 = zeros(length(s),1);
integralC2 = zeros(length(s),1);
yINT = zeros(length(s),1);
for i=2:length(s)
    integralC1(i) = trapz(s(1:i),C1(1:i));
    yINT(i) = C2(i)*exp(-integralC1(i));
end
for i=2:length(s)
    integralC2(i) = trapz(s(1:i),yINT(1:i));
end
% Boundary Condition
if (z(1) < 1e-4)
    C3 = Td;
else
    C3 = 0;
end
% Solve for Y(s)
Y = exp(integralC1).*integralC2 + C3.*exp(integralC1);

Y(indK2) = K4(indK2)./(K1*DK2(indK2)-K3(indK2));
figure(10); hold on; plot(s,DK2);

end