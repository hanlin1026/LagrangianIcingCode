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