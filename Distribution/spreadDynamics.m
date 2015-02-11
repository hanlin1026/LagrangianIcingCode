function spreadDynamics(cloud,airfoil)

% CALCULATE SPREAD MODE DYNAMICS ******************************************
% Treat as pure "sticking" and simply add mass to airfoil film
indSpread = cloud.spread;
indStateSpread = cloud.impinge(indSpread);
indexTrackSpread = cloud.index(indStateSpread);

origSpread = [];
if ~isempty(indSpread)
    % Pull out relevant properties
    pxSpread = cloud.x(indStateSpread); pySpread = cloud.y(indStateSpread);
    rdSpread = cloud.rd(indStateSpread);
    rhol = cloud.rhol;
    mSpread = (4/3)*pi*rdSpread.^3*rhol;
    % Update spread particle properties: set dt=0
    set(cloud,'dt',[indexTrackSpread, zeros(length(indexTrackSpread),1)]);
    % Deposit mass at s-coordinate location
    sCoordSpread = airfoil.XYtoScoords(pxSpread,pySpread);
    set(airfoil,'FILMspread',[sCoordSpread, mSpread]);
    % Record splashing of original droplets
    indtmp = find(indexTrackSpread <= cloud.originalNumParticles);
    origSpread = sCoordSpread(indtmp);
    set(airfoil,'originalImpingeScoordSpread',origSpread);
    set(airfoil,'originalImpingeScoord',[airfoil.originalImpingeScoordSplash; origSpread]);
end
% END SPREAD MODE DYNAMICS ************************************************

end