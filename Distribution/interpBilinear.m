function [pg,ug,vg] = interpBilinear(x,y,rho,rhou,rhov,Ubar,rhoinf,IND,xq,yq)
% Function to perform fast, vectorized bilinear interpolation

NUM = length(IND);
IND_ip1_j = IND+1;
IND_i_jp1 = IND+513;
IND_ip1_jp1 = IND_ip1_j+513;

X = zeros(NUM,4);
Y = zeros(NUM,4);
V = zeros(3*NUM,4);

X(:,1) = x(IND);
X(:,2) = x(IND_ip1_j);
X(:,3) = x(IND_ip1_jp1);
X(:,4) = x(IND_i_jp1);

Y(:,1) = y(IND);
Y(:,2) = y(IND_ip1_j);
Y(:,3) = y(IND_ip1_jp1);
Y(:,4) = y(IND_i_jp1);

V(1:NUM,1) = rho(IND);         V(NUM+1:2*NUM,1) = rhou(IND);         V(2*NUM+1:end,1) = rhov(IND);
V(1:NUM,2) = rho(IND_ip1_j);   V(NUM+1:2*NUM,2) = rhou(IND_ip1_j);   V(2*NUM+1:end,2) = rhov(IND_ip1_j);
V(1:NUM,3) = rho(IND_ip1_jp1); V(NUM+1:2*NUM,3) = rhou(IND_ip1_jp1); V(2*NUM+1:end,3) = rhov(IND_ip1_jp1);
V(1:NUM,4) = rho(IND_i_jp1);   V(NUM+1:2*NUM,4) = rhou(IND_i_jp1);   V(2*NUM+1:end,4) = rhov(IND_i_jp1);

%X = [x(i,j); x(i+1,j); x(i+1,j+1); x(i,j+1)];
%Y = [y(i,j); y(i+1,j); y(i+1,j+1); y(i,j+1)];
%V = [  rho(i,j)     rhou(i,j)     rhov(i,j); ...
%       rho(i+1,j)   rhou(i+1,j)   rhov(i+1,j); ...
%       rho(i+1,j+1) rhou(i+1,j+1) rhov(i+1,j+1); ...
%       rho(i,j+1)   rhou(i,j+1)   rhov(i,j+1)];

% Algorithm to do biliear interpolation
%AI = [  1     0     0     0; ...
%       -1     1     0     0; ...
%       -1     0     0     1; ...
%        1    -1     1    -1];
    
a = zeros(NUM,4);
b = zeros(NUM,4);
a(:,1) = X(:,1);
a(:,2) = -X(:,1)+X(:,2);
a(:,3) = -X(:,1)+X(:,4);
a(:,4) = X(:,1)-X(:,2)+X(:,3)-X(:,4);

b(:,1) = Y(:,1);
b(:,2) = -Y(:,1)+Y(:,2);
b(:,3) = -Y(:,1)+Y(:,4);
b(:,4) = Y(:,1)-Y(:,2)+Y(:,3)-Y(:,4);

%a = AI*X;
%b = AI*Y;
% Quadratic equation coeffs, aa*mm^2+bb*m+cc=0
aa = a(:,4).*b(:,3) - a(:,3).*b(:,4);
bb = a(:,4).*b(:,1) - a(:,1).*b(:,4) + a(:,2).*b(:,3) - a(:,3).*b(:,2) + xq.*b(:,4) - yq.*a(:,4);
cc = a(:,2).*b(:,1) - a(:,1).*b(:,2) + xq.*b(:,2) - yq.*a(:,2);
% Compute mapped coordinates (l,m)
det = sqrt(bb.*bb - 4.*aa.*cc);
m = (-bb+det)./(2*aa);
l = (xq-a(:,1)-a(:,3).*m)./(a(:,2)+a(:,4).*m);
% Convert (l,m) indices to physical interpolated coordinates
Xint = a(:,1) + a(:,2).*l + a(:,3).*m + a(:,4).*l.*m;
Yint = b(:,1) + b(:,2).*l + b(:,3).*m + b(:,4).*l.*m;
% Interpolate function values
dl = repmat(l,3,1);
dm = repmat(m,3,1);
%Vq = (1-dl).*(1-dm).*V(:,1) + ...
%		   dl.*(1-dm).*V(:,2) + ...
%		   dl.*dm.*V(:,3) +...
%		   (1-dl).*dm.*V(:,4);	
Vq = 0.25*(V(:,1)+V(:,2)+V(:,3)+V(:,4));
       
pg = rhoinf*Vq(1:NUM);
ug = Vq(NUM+1:2*NUM)*Ubar*rhoinf./pg;
vg = Vq(2*NUM+1:end)*Ubar*rhoinf./pg;

end