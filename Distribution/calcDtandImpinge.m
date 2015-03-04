function calcDtandImpinge(cloud,airfoil,fluid)
% Function which does the following:
% (1) calculate which cell is occupied by all droplets
% (2) set global time step based on cell volumes
% (3) determine and track which droplets have impinged on the airfoil

x = fluid.x; y = fluid.y; NS = fluid.NS;
dt = [];t = cloud.time;
if cloud.FLAGtimeResolve==1
    % Find parcels which are have already entered the injection domain
    tGLOB = cloud.tGLOB;
    indT = find(t<=tGLOB);
    cloud.indT = indT;
else
    % Consider all particles simultaneously
    indT = [1:cloud.particles]';
    cloud.indT = indT;
end
% Find parcels which have already impinged/stuck to airfoil and except them
indStick = cloud.impingeTotal;
indAdv = setdiff(indT,indStick);
cloud.indAdv = indAdv;
% Nearest neighbor grid point search over those particles currently in the
% simulation which have not impinged/stuck to airfoil
xq = cloud.x(indAdv); yq = cloud.y(indAdv); 
uq = cloud.u(indAdv); vq = cloud.v(indAdv);
particles = size(xq,1);
ind = knnsearch(NS,[xq,yq]);
% Clear impinge index tracker
set(cloud,'impinge',[]);
for i=1:particles
    [pBL, pBR, pTR, pTL, surfaceFlag] = findCell(x,y,xq(i),yq(i),ind(i));
    % Determine local time step
    L1 = norm(pBL-pBR);
    L2 = norm(pTR-pBR);
    area = L1*L2;
    vel = [uq(i); vq(i)];
    if surfaceFlag==1
        % Calculate normal velocity to make sure impingement is occuring
        [~,~,nx,ny,~,~] = airfoil.findPanel(xq(i),yq(i));
        normvel = vel(1)*nx+vel(2)*ny;
        if normvel<0.01
            % NOTE: dt will be reset to dt=0 for splash/spread impingement
            % modes separately in their respective modules!
            dt(i,1) = 0.2*sqrt(area)/norm(vel);
            % Only retain new impingement indices
            impinge = setdiff(indAdv(i),cloud.impingeTotal);
            if ~isempty(impinge)
                set(cloud,'impinge',impinge);
            end
        else
            % If not, set local time step
            dt(i,1) = 0.2*sqrt(area)/norm(vel);
        end
    elseif xq(i)>1
        % If we are past the airfoil, stop advecting
        dt(i,1) = 0;
    else
        % If not, set local time step
        dt(i,1) = 0.2*sqrt(area)/norm(vel);
    end
end

if cloud.FLAGtimeResolve==1
    % Set global timestep as minimum of dt
    dtNZ = dt(dt~=0);
    if isempty(dtNZ)
    % If we are in 'dead region' between next impinging particle, fast forward
    % in time at a specified rate while advecting no particles
        dt = zeros(length(indAdv),1);
        dtINTER = 0.01;
        set(cloud,'tGLOB',cloud.tGLOB+dtINTER);
    else
        dt(:) = min(dtNZ);
    end
end
set(cloud,'dt',dt);

end