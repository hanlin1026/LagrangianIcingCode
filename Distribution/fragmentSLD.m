function fragmentSLD(cloud,fluid)
% Function to compute fragmentation of SLD's (before surface impact)

indAdv = cloud.indAdv;
if ~isempty(cloud.indAdv)
    % Pull out state variables of all advecting particles
    x = cloud.x(indAdv); y = cloud.y(indAdv); u = cloud.u(indAdv); v = cloud.v(indAdv);
    rd = cloud.rd(indAdv); temp = cloud.Temp(indAdv);
    t = cloud.time(indAdv); numDrop = cloud.numDroplets(indAdv);
    indCell = cloud.indCell(indAdv);
    dt = cloud.dt; rhol = cloud.rhol;
    % Get fluid properties at nearest neighbor grid pts
    [pg,ug,vg] = fluid.getNNFluidProps(x,y);

    sigma = 75.64e-3; % Surface tension of water at 0 deg C against air
    % Compute Weber number
    We = 2*pg.*((u-ug).^2 + (v-vg).^2).*rd./sigma;
    % Compute mean break-up duration (from empirical correlation)
    tau = 5*sqrt(rhol./pg)*2.*rd./sqrt((u-ug).^2 + (v-vg).^2);

    % If greater than threshold, compute break-up probability
    prob = zeros(size(u,1),1);
    ind = find(We>12);
    prob(ind) = 1 - exp(-dt(ind)./tau(ind));

    % Choose uniform random variable on [0,1]. Fragmentation occurs if
    % alpha<prob
    alpha = unifrnd(0,1,[size(u,1),1]);
    ind_prob = find(alpha<prob);
    indFRAG = ind_prob;

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

        % For each fracturing droplet, calculate children's sizes
        state = []; indCellNew = []; NewTotal = 0;
        for i=1:size(ind_prob,1)
            xFRAG = []; yFRAG = []; uFRAG = []; vFRAG = [];
            rFRAG = []; timeFRAG = []; tempFRAG = [];
            num_children = randi([1,5],1,1); % Number of break-up particle children
            F_children = unifrnd(0,1,[num_children,1]); % Random CDF sampling
            % Interpolate analytical expression for the CDF to find children's size
            rsamp = linspace(0,rd(ind_prob(i)),1000)';
            CUMDISTbup = 1 - 1/6./rstar(i).^3.*(6*rstar(i).^3 + 6*rstar(i).^2.*rsamp + 3*rstar(i).*rsamp.^2 + rsamp.^3).*exp(-rsamp./rstar(i));
            rFRAG = interp1(CUMDISTbup,rsamp,F_children);
            numFRAG = repmat(numDrop(ind_prob(i)),num_children,1);
            % Set relevant state properties
            uFRAG = ug(ind_prob(i)) - (ug(ind_prob(i))-u(ind_prob(i)))./(1 + 2.7*sqrt(pg(ind_prob(i))./rhol).*rFRAG./rd(ind_prob(i)));
            vFRAG = vg(ind_prob(i)) - (vg(ind_prob(i))-v(ind_prob(i)))./(1 + 2.7*sqrt(pg(ind_prob(i))./rhol).*rFRAG./rd(ind_prob(i)));
            xFRAG = x(ind_prob(i))*ones(num_children,1);
            yFRAG = y(ind_prob(i))*ones(num_children,1);
            tempFRAG = temp(ind_prob(i))*ones(num_children,1);
            timeFRAG = cloud.time(ind_prob(i))*ones(num_children,1);
            % Add new SLD particles to the cloud with relevant state properties
            for j=1:num_children
                state = [state; [xFRAG(j) yFRAG(j) uFRAG(j) vFRAG(j) rFRAG(j) tempFRAG(j) timeFRAG(j) numFRAG(j)]];
            end
            iCell = repmat(indCell(ind_prob(i)),num_children,1);
            indCellNew = [indCellNew; iCell];
            % Reset mass of parent particle (by conservation of mass)
            massCHILD = sum((4/3*pi*rhol)*rFRAG.^3);
            massPARENT = (4/3*pi*rhol)*rd(ind_prob(i))^3;
            rNEW = ((massPARENT-massCHILD)/(4/3*pi*rhol))^(1/3);
            if ~isreal(rNEW) || (rNEW < 2.5e-6) || isnan(rNEW)
                rNEW = 2.5e-6;
            end
            NewTotal = NewTotal + num_children;
            set(cloud,'rd',[indAdv(ind_prob(i)), rNEW]);
            %}
        end
        % Record indices of old fractured parent droplets
        set(cloud,'fracture',[]); set(cloud,'fracture',indAdv(ind_prob));
        % Add new chidl particles to the cloud
        cloud.addParticle(state,indCellNew);

    end
    if ~isempty(ind_prob)
        [indAdv(ind_prob) cloud.rd(indAdv(ind_prob))]
    end

end


end