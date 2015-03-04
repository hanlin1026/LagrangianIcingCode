function [xCONV,yCONV] = smoothGaussianKernel1D(x,y,S)
% Function to smooth a 1D signal using a Gaussian Kernel
% INPUTS:
%   x,y: coordinates of the signal
%   S: standard deviation of the Gaussian smoothing kernel

% Set smoothing kernel
G = @(xq,mu) 1/sqrt(2*pi)/S*exp(-0.5*(xq-mu).^2./(S^2));
% 1D interpolation over fine grid
N = length(x);
xINT = linspace(min(x),max(x),500)';
yINT = interp1(x,y,xINT);
% Convolution
xCONV = xINT;
yCONV = zeros(length(xINT),1);
for i=1:length(xINT)
    yCONV(i) = trapz(xINT,yINT.*G(xINT(i)-xINT,0));
end

end

