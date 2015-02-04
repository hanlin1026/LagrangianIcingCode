function impingementRegimeSLD(cloud,airfoil)
% Function to compute the impingement regime for those drops which have
% struck the airfoil surface. Returns indices for each of three cases.
% INPUTS: SLDcloud object, gas velocity and density at particle positions, Airfoil object

% Compute impingement parameters
cloud.computeImpingementParams(airfoil);
% Compute impingement regime dynamics
splashDynamics(cloud,airfoil);
bounceDynamics(cloud,airfoil);
spreadDynamics(cloud,airfoil);

% Delete old parent splashed and spread droplets
set(airfoil,'FILM',[airfoil.FILMsplash; airfoil.FILMspread]);
indSplash = cloud.splash; % Indices of cloud.impinge which have splashed
indStateSplash = cloud.impinge(indSplash); % Splash indices for state variables
if ~isempty(indStateSplash)
    cloud.deleteParticle(indStateSplash);
end

end