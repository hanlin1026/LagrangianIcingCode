function bounceDynamics(cloud,airfoil)

% CALCULATE BOUNCE MODE DYNAMICS ******************************************
indBounce = cloud.impinge(cloud.bounce); % Indices of cloud.impinge which have bounced
bounce = cloud.bounce;

if ~isempty(cloud.bounce)
    % Pull out quantities for the bounce mode
    xBounce = cloud.x(indBounce); yBounce = cloud.y(indBounce); uBounce = cloud.u(indBounce); vBounce = cloud.v(indBounce);
    rdBounce = cloud.rd(indBounce); rhol = cloud.rhol;
    
    KBounce = cloud.K(bounce); KsBounce = cloud.fs(bounce)*cloud.Ks0; vnormsqBounce = cloud.normvelsq(bounce);
    [~,~,nxBounce,nyBounce,txBounce,tyBounce] = airfoil.findPanel(xBounce,yBounce);
    KbBounce = cloud.fb(bounce)*cloud.Kb0;
    vnormBounce = sqrt(vnormsqBounce); vtangBounce = cloud.tangvel(bounce);
    % Calculate post impact velocity for bouncing droplets
    vn = 4*vnormBounce.*(sqrt(KBounce./KbBounce)-KBounce./KbBounce);
    vt = 0.8*vtangBounce;
    vnew = [vn,vn].*[nxBounce,nyBounce] + [vt,vt].*[txBounce,tyBounce];
    % Reset velocities of bouncing droplets
    set(cloud,'u',[vnew(:,1), indBounce]);
    set(cloud,'v',[vnew(:,2), indBounce]);
end
% END BOUNCE MODE DYNAMICS ************************************************

end