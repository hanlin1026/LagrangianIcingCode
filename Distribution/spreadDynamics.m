function spreadDynamics(cloud,airfoil)
% Functionality:
% (1) Add all of the spreading mass to airfoil.FILMspread

% CALCULATE SPREAD MODE DYNAMICS ******************************************
% Treat as pure "sticking" and simply add mass to airfoil film
indSpread = cloud.spread;
indStateSpread = cloud.impinge(indSpread);

origSpread = [];
if ~isempty(indSpread)
    % Pull out relevant properties
    pxSpread = cloud.x(indStateSpread); pySpread = cloud.y(indStateSpread);
    rdSpread = cloud.rd(indStateSpread);
    nDropSpread = cloud.numDroplets(indStateSpread);
    rhol = cloud.rhol;
    mSpread = (4/3)*pi*rdSpread.^3*rhol;
    % Deposit mass at s-coordinate location
    sCoordSpread = airfoil.XYtoScoords(pxSpread,pySpread);
    set(airfoil,'FILMspread',[]);
    set(airfoil,'FILMspread',[sCoordSpread, nDropSpread.*mSpread]);
    % Record splashing of original droplets
    indtmp = find(indStateSpread <= cloud.originalNumParticles);
    origSpread = sCoordSpread(indtmp);
    set(airfoil,'originalImpingeScoordSpread',origSpread);
    set(airfoil,'originalImpingeScoord',[airfoil.originalImpingeScoordSplash; origSpread]);
end
% END SPREAD MODE DYNAMICS ************************************************

end