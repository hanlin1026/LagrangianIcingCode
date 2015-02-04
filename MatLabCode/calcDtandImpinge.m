function calcDtandImpinge(cloud,airfoil,fluid)
% Function which does the following:
% (1) calculate which cell is occupied by all droplets
% (2) set local time steps based on the cell volumes
% (3) determine and track which droplets have impinged on the airfoil

x = fluid.x; y = fluid.y; NS = fluid.NS;
dt = [];
% Get total number of current particles
particles = cloud.particles;
% Nearest neighbor grid point search over all particles
xq = cloud.x; yq = cloud.y; uq = cloud.u; vq = cloud.v;
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
        if normvel<0
            % If we have impinged on the airfoil surface, stop advecting
            dt(i,1) = 0;
            set(cloud,'impinge',i);
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
% Set local timesteps
set(cloud,'dt',dt);

end