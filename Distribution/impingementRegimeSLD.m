function impingementRegimeSLD(cloud,airfoil)
% Function to compute the impingement regime for those drops which are
% currently striking the airfoil surface.
% INPUTS: SLDcloud object, gas velocity and density at particle positions, Airfoil object
% OUTPUT: sets airfoil.FILM

% Compute impingement parameters
cloud.computeImpingementParams(airfoil);

% Clear airfoil.FILM
airfoil.clearFilm();

% Compute impingement regime dynamics
splashDynamics(cloud,airfoil);
bounceDynamics(cloud,airfoil);
spreadDynamics(cloud,airfoil);

% Set airfoil "film"
%[length(airfoil.FILMsplash), length(airfoil.FILMspread), length(airfoil.FILM)]
if ~isempty([airfoil.FILMsplash; airfoil.FILMspread])
    set(airfoil,'FILM',[airfoil.FILMsplash; airfoil.FILMspread]);
end

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