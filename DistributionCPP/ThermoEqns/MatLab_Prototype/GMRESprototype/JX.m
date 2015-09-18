function jx = JX(x,u0)

eps = 1e-4;
x2 = u0 + eps*x;
f2 = testBalance(x2);
f1 = testBalance(u0);
jx = (1/eps)*(f2-f1);

end