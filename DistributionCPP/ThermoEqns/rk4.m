function [T,Y] = rk4(f,y0,m)
%---------------------------------------------------------------------------
%RK4   Runge-Kutta solution for y' = f(t,y) with y(a) = ya.
% Sample call
%   [T,Y] = rk4('f',a,b,ya,m)
% Inputs
%   f    [x,y] points of function
%   y0   initial value
%   m    number of steps
% Return
%   T    solution: vector of abscissas
%   Y    solution: vector of ordinates
%
% NUMERICAL METHODS: MATLAB Programs, (c) John H. Mathews 1995
% To accompany the text:
% NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
% Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
% Prentice Hall, Inc.; USA, Canada, Mexico ISBN 0-13-624990-6
% Prentice Hall, International Editions:   ISBN 0-13-625047-5
% This free software is compliments of the author.
% E-mail address:      in%"mathews@fullerton.edu"
%
% Algorithm 9.4 (Runge-Kutta Method of Order 4).
% Section	9.5, Runge-Kutta Methods, Page 460
%---------------------------------------------------------------------------

a = f(1,1); b = f(end,1);
h = (b - a)/m;
T = zeros(1,m+1);
Y = zeros(1,m+1);
T(1) = a;
Y(1) = y0;
for j=1:m,
  tj = T(j);
  yj = Y(j);
  k1 = h*interp1(f(:,1),f(:,2),tj);
  k2 = h*interp1(f(:,1),f(:,2),tj+h/2);
  k3 = h*interp1(f(:,1),f(:,2),tj+h/2);
  k4 = h*interp1(f(:,1),f(:,2),tj+h);
  Y(j+1) = yj + (k1 + 2*k2 + 2*k3 + k4)/6;
  T(j+1) = a + h*j;
end