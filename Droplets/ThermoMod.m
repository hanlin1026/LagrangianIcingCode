classdef ThermoMod < hgsetget
    % Defines the airfoil grid used by the thermodynamic module
    
    properties
        % Grid properties
        cellCenterPts; % S-coords of airfoil grid cell center points
        
        % Mass balance quantities (quasi-steady assumption)
        mImp;
        mEvapSub;
        mIce;
        mRunIn;
        mRunOut;
               
        % Energy balance properties (quasi-steady assumption)
        
        % Flow properties
        Uinf;
        
        % Airfoil
        airfoil;
        
        % Mass/energy balance parameters
        ch; % Convective heat transfer coefficient
        cpAir; % Specific heat capacity of air
        Ts; % Surface temperature
        Pvp; % Saturation vapor pressure at the surface
        Pvinf; % Saturation vapor pressure of water in ambient air
        Ps; % Absolute pressure above the control volume outside the boundary layer
        Hrinf; % Relative humidity
        Tinf; % Free-stream temperature
        Ts; % Equilibrium temperature at the air/film/ice/wall interface
        
    end
    
    methods
        function thermoMod = ThermoMod(CFDModAirfoil,fluid,scalars)
            % Constructor for thermoFoil
            
            thermoMod.airfoil = CFDModAirfoil;
            % Calculate impinging water mass
            thermoFoil.cellCenterPts = CFDModAirfoil.s;
            ds = diff(s); ds(end+1) = ds(1);
            UFS = fluid.UFS;
            LWC = fluid.LWC;
            beta = CFDModAirfoil.beta;
            thermoFoil.mIMP = beta.*UFS*LWC.*ds;
            % Mass balance parameters
            thermoMod.cpAir = 1.006; % kJ/kg/K
            thermoMod.Tinf = fluid.Tinf; % K
            thermoMod.Ts = scalars(1); % K
            thermoMod.Pvp = thermoMod.calcSatVapPres();
            
            
        end
        
        function Pv = calcSatVapPres(thermoMod)
            % Calculate saturation vapor pressure using saturated steam tables
            
            Ts = thermoMod.Ts;
            That = 72 + 1.8*(Ts-273.15);
            Pv = 3386*(0.0039 + (6.80961e-6)*That.^2 + (3.55791e-7)*That.^3);
        end
        
    end
    
end

