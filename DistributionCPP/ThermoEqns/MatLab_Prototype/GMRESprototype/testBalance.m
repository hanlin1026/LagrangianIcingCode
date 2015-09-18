function err = testBalance(x)

A = [8 1 6; 3 5 7; 4 9 2];
RHS = [1;-1;2];

err = A*x-RHS;

end