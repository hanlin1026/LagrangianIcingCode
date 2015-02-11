meshfile = 'MESH.P3D';
solnfile = 'q103.0.50E+01.bin';

[x,y,rho,rhou,rhov,E,mach,alpha,reynolds,time] = readp3d(meshfile,solnfile);
disp(mach)
disp(alpha)
disp(reynolds)

figure(1)
plot(x,y,'b-',x',y','b-'); 
axis equal
xlim([-2 2])
ylim([-2 2])

figure(2);

subplot(2,2,1)
contourf(x,y,rho)
axis equal
xlim([-2 2])
ylim([-2 2])

subplot(2,2,2)
contourf(x,y,rhou)
axis equal
xlim([-2 2])
ylim([-2 2])

subplot(2,2,3)
contourf(x,y,rhov)
axis equal
xlim([-2 2])
ylim([-2 2])

subplot(2,2,4)
contourf(x,y,E)
axis equal
xlim([-2 2])
ylim([-2 2])

u = rhou./rho;
v = rhov./rho;
figure(3);
quiver(x,y,u,v,0.4)
axis equal
xlim([-2 2])
ylim([-2 2])

