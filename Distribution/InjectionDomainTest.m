% Test script for InjectionDomain

% Read in flow solution
meshfile = 'MESH.P3D';
solnfile = 'q103.0.50E+01.bin';
[~,~,~,~,~,~,mach,alpha,~,~] = readp3d(meshfile,solnfile);
% Dimensional reference quantities
pinf = 1.01325e5; % N/m^2
R = 287.058; % J/kg/K
Tinf = 300; % K
rhoinf = pinf/R/Tinf;
Ubar = sqrt(pinf/rhoinf);
rhol = 1000; % kg/m^3
Rd = 50e-6; % m
LWC = 0.4e-3; % kg/m^3
Uinf = mach*340; % m/s
% Initialize fluid object
scalars = [pinf;R;Tinf;rhoinf;Ubar;rhol];
fluid = Fluid(scalars,meshfile,solnfile);
% Initialize the airfoil surface
x = fluid.x; y = fluid.y;
ind = find(x(:,1)<=1);
ax = x(ind,1); ay = y(ind,1);
airfoil = Airfoil([ax,ay]);
% Initialize domain and associated distribution functions
strPDFTypes = {'Gaussian','Gaussian','Gaussian','Gaussian','Uniform'};
simTime = 60*10;
PDFparams = [20 2; 10 1; 50e-6 10e-6; 0 2; 0 simTime];
domain = InjectionDomain(strPDFTypes,PDFparams,fluid,airfoil,LWC,simTime);
nClumps = 100;
domain.sampleRealization(nClumps);
domain.dispSampleStatistics();