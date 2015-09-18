function un = NewtonKrylovIteration(F,J,u0)

nitermax = 30;
restart = 10;
maxit = 100;
tol = 1e-3;
globaltol = 1e-5;

%r = -F(u0,scalars);
un = u0;

disp('******* Start JFNK *******');
for i = 1:nitermax
    u0 = un;
    x = u0;
    F_k = F(u0);
    du0 = gmresCustom(J, -F_k, restart, tol, maxit, x, u0 );
    %[du0,flag] = bicgstab( @(du) J(du), -F_k, tol, maxit );
    %[du0,flag] = bicgstab( @(du) J(du), -F_k, tol, maxit,[],[],du0 );
    un = u0 + du0;
    globalerr = F(un);
    r = -globalerr - J(un,u0);
    disp(  sprintf('JFNK error: %0.5g \t global error: %0.5g',norm(r),norm(globalerr)) );
    if ( norm(globalerr) < globaltol )
        return;
    end
   
end

end