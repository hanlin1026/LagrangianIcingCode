function bounceDynamics(cloud,airfoil)

% CALCULATE BOUNCE MODE DYNAMICS ******************************************
indBounce = cloud.bounce; % Indices of cloud.impinge which have bounced
indStateBounce = cloud.impinge(indBounce); % Bounce indices for state variables
indexTrackBounce = cloud.index(indStateBounce); % Index trackers

if ~isempty(indBounce)
    % Pull out quantities for the bounce mode
    xBounce = cloud.x(indStateBounce); yBounce = cloud.y(indStateBounce); uBounce = cloud.u(indStateBounce); vBounce = cloud.v(indStateBounce);
    rdBounce = cloud.rd(indStateBounce); rhol = cloud.rhol;
    
    KBounce = cloud.K(indBounce); KsBounce = cloud.fs(indBounce)*cloud.Ks0; vnormsqBounce = cloud.normvelsq(indBounce);
    [~,~,nxBounce,nyBounce,txBounce,tyBounce] = airfoil.findPanel(xBounce,yBounce);
    KbBounce = cloud.fb(indBounce)*cloud.Kb0;
    vnormBounce = sqrt(vnormsqBounce); vtangBounce = cloud.tangvel(indBounce);
    % Calculate post impact velocity for bouncing droplets
    vn = 4*vnormBounce.*(sqrt(KBounce./KbBounce)-KBounce./KbBounce);
    vt = 0.8*vtangBounce;
    vnew = [vn,vn].*[nxBounce,nyBounce] + [vt,vt].*[txBounce,tyBounce];
    % Reset velocities of bouncing droplets
    set(cloud,'u',[vnew(:,1), indStateBounce]);
    set(cloud,'v',[vnew(:,2), indStateBounce]);
end
% END BOUNCE MODE DYNAMICS ************************************************

end