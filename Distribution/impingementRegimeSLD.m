function impingementRegimeSLD(cloud,airfoil)
% Function to compute the impingement regime for those drops which have
% struck the airfoil surface.
% INPUTS: SLDcloud object, gas velocity and density at particle positions, Airfoil object
% OUTPUT: sets airfoil.FILM

% Compute impingement parameters
cloud.computeImpingementParams(airfoil);
% Compute impingement regime dynamics
splashDynamics(cloud,airfoil);
bounceDynamics(cloud,airfoil);
spreadDynamics(cloud,airfoil);

% Set airfoil "film"
set(airfoil,'FILM',[airfoil.FILMsplash; airfoil.FILMspread]);

% Reset number of particles in simulation
t = cloud.time;
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

end