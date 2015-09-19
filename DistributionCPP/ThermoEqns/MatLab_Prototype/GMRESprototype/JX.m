function jx = JX(x,u0,scalars)

eps = 1e-4;
x2 = u0 + eps*x;
f2 = massBalance(x2,scalars);
f1 = massBalance(u0,scalars);
jx = (1/eps)*(f2-f1);

end