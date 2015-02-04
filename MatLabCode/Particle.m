classdef Particle < hgsetget
    % Defines a droplet (SLD) and associated operations
    
    properties
        % These properties are time histories of all state quantities
        x=[]; y=[]; % Position
        u=[]; v=[]; % Velocity
        time=[]; % Time
        rd; % Radius
        dt; % Timestep
        rhol; % Density
        index; % Numeric index
        fracture; % '0' if not, '1' if yes
        impinge; % '0' if not, '1' if yes
        childind=[]; % Numeric indices of child particles
                       
    end
    
    methods
        function particle = Particle(state0)
            % Constructor: state = (x0,y0,u0,v0,r0,t0,rho,ind0)
            
            particle.x(1,1) = state0(1);
            particle.y(1,1) = state0(2);
            particle.u(1,1) = state0(3);
            particle.v(1,1) = state0(4);
            particle.rd = state0(5);
            particle.time(1,1) = state0(6);
            particle.rhol = state0(7);
            particle.index = state0(8);
            
            particle.fracture = 0;
            particle.impinge = 0;
        end
        
        function particle = set.x(particle,x)
            % Function to update state
            count = length(particle.x)+1;
            particle.x(count,1) = x;
            
        end
        
        function particle = set.y(particle,y)
            % Function to update state
            count = length(particle.y)+1;
            particle.y(count,1) = y;
            
        end
        
        function particle = set.u(particle,u)
            % Function to update state
            count = length(particle.u)+1;
            particle.u(count,1) = u;
            
        end
        
        function particle = set.v(particle,v)
            % Function to update state
            count = length(particle.v)+1;
            particle.v(count,1) = v;
            
        end
        
        function particle = set.time(particle,time)
            % Function to update state
            count = length(particle.time)+1;
            particle.time(count,1) = time;
            
        end
        
        function particle = set.dt(particle,dt)
            % Function to update state
            particle.dt = dt;
            
        end
        
        function particle = set.fracture(particle,fracture)
            % Function to set fracture indicator
            particle.fracture = fracture;
        end
        
        function particle = set.impinge(particle,impinge)
            % Function to set impingement indicator
            particle.impinge = impinge;
        end
        
        function particle = set.childind(particle,childind)
            % Function to set indices for child fracture particles
            particle.childind = childind;
        end
        
        
    end
    
end

