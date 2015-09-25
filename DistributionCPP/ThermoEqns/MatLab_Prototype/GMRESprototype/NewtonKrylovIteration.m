function un = NewtonKrylovIteration(F,J,u0,scalars)

nitermax = 20;
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
    F_k = F(u0,scalars);
    [du0,flag] = gmres(@(du) J(du,u0,scalars),-F_k,restart,tol,maxit);
    %du0 = gmres(J, -F_k, restart, tol, maxit, x, u0, scalars);
    %[du0,flag] = bicgstab( @(du) J(du), -F_k, tol, maxit );
    %[du0,flag] = bicgstab( @(du) J(du), -F_k, tol, maxit,[],[],du0 );
    un = u0 + du0;
    globalerr = F(un,scalars);
    r = -globalerr - J(un,u0,scalars);
    disp(  sprintf('JFNK error: %0.5g \t global error: %0.5g',norm(r),norm(globalerr)) );
    if ( norm(globalerr) < globaltol )
        return;
    end
   
end

end