function formInjectionBox(cloud,airfoil,fluid)
% Calculate length of box needed for continuous injection of particles

% Find largest initial time step
cloud.calcDtandImpinge(airfoil,fluid);
dt0_dtSpace = ceil(max(cloud.dt./cloud.dtSpacingAvg));
% Advect all initial particles backwards by the largest local timestep
thecloud = cloud;
set(thecloud,'dt',-thecloud.dt);
transportSLD(thecloud,fluid);
% Place enough particles between initial and backwards integration
% locations to satisfy constraint on LWC (mean separation distance)



end