%{
STATE = {};
strImpMod = 'Impingement';
realizations = 10;
xtmp = linspace(-.3,.15,1000)';
BETA = {};
for i=1:realizations
    i
    domain.nDroplet = [];
    domain.samples = [];
    domain.numParcels = [];
    domain.sampleRealization(nClumps);
    
    [STATE,totalImpinge,impinge,s,beta,sStick,mStick] = calcCollectionEfficiency(airfoil,fluid,domain,strImpMod);
    airfoil.calcStagPt(fluid);
    sCENT = 0.5*(s(2:end)+s(1:end-1)) - airfoil.stagPt;
    figure(2); hold on; plot(sCENT,beta); drawnow;
    [xSMOOTH,ySMOOTH] = smoothGaussianKernel1D(sCENT,beta,0.0075);
    figure(3); hold on; plot(xSMOOTH,ySMOOTH,'k'); drawnow;
    BETA{i} = [sCENT,beta];
    
    % Plot moving average
    ytmp = zeros(1000,i);
    for j=1:i
        ytmp(:,j) = interp1(BETA{j}(:,1),BETA{j}(:,2),xtmp);
    end
    MEAN = sum(ytmp,2)/i;
    figure(22); plot(-xtmp,MEAN); drawnow;
end
save('BETA_23012_BIN27_SPLASH.mat','BETA');
%}
%{
% Initialize domain and associated distribution functions
strPDFTypes = {'Implicit','Implicit','Custom','Gaussian','Uniform'};
simTime = 60*30;
uDIST = 0*Uinf*cos(alpha/180)*[1 0.00001];
vDIST = 0*Uinf*sin(alpha/180)*[1 0.00001];
load('BIN27.mat');
BIN27(:,1) = BIN27(:,1);
rDIST = BIN27;
PDFparams = {};
PDFparams{1} = uDIST; 
PDFparams{2} = vDIST; 
PDFparams{3} = rDIST;
PDFparams{4} = [0 2];
PDFparams{5} = [0 simTime];
domain = InjectionDomain(strPDFTypes,PDFparams,fluid,airfoil,LWC,simTime);
nClumps = 1000;
domain.sampleRealization(nClumps,fluid);
domain.dispSampleStatistics();
%}

x = fluid.x; y = fluid.y;
% Sample point
I = 50; J = 200;
IND = 512*(I-1) + J;
[Ig,Jg] = transformXYtoIJ(fluid,IND,[x(:),y(:)]);
[Ic,Jc] = transformXYtoIJ(fluid,IND,[fluid.MEANx(:),fluid.MEANy(:)]);
II = reshape(Ig,513,129); JJ = reshape(Jg,513,129);
figure(2); plot(II,JJ,'b',II',JJ','b');
hold on; scatter(Ic,Jc,'ro');
xlim([-5,5]); ylim([-5,5]);
% Draw circle around point
xC = fluid.MEANx(IND); yC = fluid.MEANy(IND);
th = linspace(0,2*pi,100)'; R = 1e-5;
xt = xC + R*sin(th);
yt = yC + R*cos(th);
figure(1); hold on; plot(x,y,'k',x',y','k');
hold on; plot(xC,yC,'ro',xt,yt,'r');
% Draw transformed circle
[It,Jt] = transformXYtoIJ(fluid,IND,[xt,yt]);
figure(2); hold on; plot(It,Jt,'r');





