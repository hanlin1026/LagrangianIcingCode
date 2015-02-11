function fragmentSLD(cloud,fluid)
% Function to compute fragmentation of SLD's (before surface impact)

% Pull out state variables of all particles
x = cloud.x; y = cloud.y; u = cloud.u; v = cloud.v;
rd = cloud.rd; t = cloud.time; dt = cloud.dt; rhol = cloud.rhol;
[rhog,ug,vg] = fluid.interpFluid(x,y);

sigma = 75.64e-3; % Surface tension of water at 0 deg C against air
% Compute Weber number
We = 2*rhog.*((u-ug).^2 + (v-vg).^2).*rd./sigma;
% Compute mean break-up duration (from empirical correlation)
tau = 5*sqrt(rhol./rhog)*2.*rd./sqrt((u-ug).^2 + (v-vg).^2);

% If greater than threshold, compute break-up probability
prob = zeros(size(u,1),1);
ind = find(We>12);
prob(ind) = 1 - exp(-dt(ind)./tau(ind));

% Choose uniform random variable on [0,1]. Fragmentation occurs if
% alpha<prob
alpha = unifrnd(0,1,[size(u,1),1]);
ind_prob = find(alpha<prob);
indFRAG = ind_prob;

% Discount particles which have already impinged on the airfoil surface
% from fracturing
indIMPINGE = cloud.impingeTotal;
tmpind = setdiff(ind_prob,indIMPINGE);
ind_prob = tmpind;

if ~isempty(ind_prob)
    % Select only those particles which have fractured for processing
    we = We(ind_prob);
    % Choose radius of particles according to number density function of
    % resulting droplets
    ind1 = find(we>12 & we<=18);
    ind2 = find(we>18 & we<=45);
    ind3 = find(we>45);
    r32_r0 = zeros(size(ind_prob,1),1);
    r32_r0(ind1) = 0.32*we(ind1).^(-1/3).*(4.1./(we(ind1)-12).^(1/4)).^(2/3);
    r32_r0(ind2) = 0.32*we(ind2).^(-1/3).*((2.45*(we(ind2)-12).^(1/2)-1.9)./(we(ind2)-12).^(1/4)).^(2/3);
    r32_r0(ind3) = 0.32*we(ind3).^(-1/3).*(12.2./(we(ind3)-12).^(1/4)).^(2/3);
    
    r32 = r32_r0.*rd(ind_prob);
    rstar = r32./3;
    
    numFRAG = 0;
    % For each fracturing droplet, calculate children's sizes
    for i=1:size(ind_prob,1)
        xFRAG = []; yFRAG = []; uFRAG = []; vFRAG = [];
        rFRAG = []; timeFRAG = [];
        num_children = randi([1,5],1,1); % Number of break-up particle children
        F_children = unifrnd(0,1,[num_children,1]); % Random CDF sampling
        % Interpolate analytical expression for the CDF to find children's size
        rsamp = linspace(0,rd(ind_prob(i)),1000)';
        CUMDISTbup = 1 - 1/6./rstar(i).^3.*(6*rstar(i).^3 + 6*rstar(i).^2.*rsamp + 3*rstar(i).*rsamp.^2 + rsamp.^3).*exp(-rsamp./rstar(i));
        rFRAG = interp1(CUMDISTbup,rsamp,F_children);
        numFRAG = numFRAG + length(rFRAG);
        % Set relevant state properties
        uFRAG = ug(ind_prob(i)) - (ug(ind_prob(i))-u(ind_prob(i)))./(1 + 2.7*sqrt(rhog(ind_prob(i))./rhol).*rFRAG./rd(ind_prob(i)));
        vFRAG = vg(ind_prob(i)) - (vg(ind_prob(i))-v(ind_prob(i)))./(1 + 2.7*sqrt(rhog(ind_prob(i))./rhol).*rFRAG./rd(ind_prob(i)));
        xFRAG = x(ind_prob(i))*ones(length(rFRAG),1);
        yFRAG = y(ind_prob(i))*ones(length(rFRAG),1);
        timeFRAG = cloud.time(ind_prob(i))*ones(length(rFRAG),1);
        % Add new SLD particles to the cloud with relevant state properties
        for j=1:num_children-1
            state = [xFRAG(j) yFRAG(j) uFRAG(j) vFRAG(j) rFRAG(j) timeFRAG(j)];
            cloud.addParticle(state);
        end
        % Record last child particle in the place of the original one
        set(cloud,'u',[uFRAG(num_children), ind_prob(i)]);
        set(cloud,'v',[vFRAG(num_children), ind_prob(i)]);
        set(cloud,'rd',[ind_prob(i), rFRAG(num_children)]);
    end
    % Record indices of old fractured parent droplets
    set(cloud,'fracture',[]); set(cloud,'fracture',ind_prob);
    % Delete old fractured parent droplets
    %cloud.deleteParticle(ind_prob);

end


end