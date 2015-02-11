function [X,Y,rho,rhou,rhov,E,mach,alpha,reynolds,time] = readp3d(meshfile,solnfile)

fid = fopen(meshfile,'r');

xy = importdata(meshfile,' ',1);
fclose(fid);

n = sscanf(char(xy.textdata),'%i %i');
ni = n(1); nj = n(2);
n = ni*nj;
xy = transpose(xy.data);
xy = xy(:);

% clear out error from reading end value of non-square array
while true
    if isnan(xy(end))
        xy = xy(1:end-1);
    else
        break
    end
end

x = xy(1:n);
y = xy(n+1:end);


solnfid = fopen(solnfile,'r');
fread(solnfid,[2],'int');  % mesh size
params = fread(solnfid,[4],'float'); % mach, alpha, reynolds, time
mach = params(1);
alpha = params(2);
reynolds = params(3);
time = params(4);

rho  = fread(solnfid,[n],'float'); rho = reshape(rho,ni,nj);
rhou = fread(solnfid,[n],'float'); rhou = reshape(rhou,ni,nj);
rhov = fread(solnfid,[n],'float'); rhov = reshape(rhov,ni,nj);
E    = fread(solnfid,[n],'float'); E = reshape(E,ni,nj);

X = reshape(x,ni,nj);
Y = reshape(y,ni,nj);

end
