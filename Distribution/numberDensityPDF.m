classdef numberDensityPDF < hgsetget
    % Defines properties/methods for a PDF describing the number density
    % function's dependence on x,u,R,e
    
    properties
        % Separate PDFs for x,u,R,e (assuming these are all independent)
        fx; fy; % Space
        fu; fv; % Velocity
        fR; % Radius
        fe; % Energy
        
    end
    
    methods
        function fn = numberDensityPDF(strTypes,params)
            % Constructor
            
            fn.fx = fn.computePDF(strTypes(1),params(1,:));
            fn.fy = fn.computePDF(strTypes(2),params(2,:));
            fn.fu = fn.computePDF(strTypes(3),params(3,:));
            fn.fv = fn.computePDF(strTypes(4),params(4,:));
            fn.fR = fn.computePDF(strTypes(5),params(5,:));
            fn.fe = fn.computePDF(strTypes(6),params(6,:));
        
        end
        
        function func = computePDF(fn,strType,params)
            % Function to return PDF of a specified type at locations
            % specified by user
            % 'x' = evaluation pts
            % 'params' = PDF dependent shape parameters
            
            if strcmp(strType,'Gaussian')
                % 'params' = [mu sigma]
                mu = params(1);
                sigma = params(2);
                func = @(x) 1/sigma/sqrt(2*pi)*exp(-0.5*(x-mu).^2./(sigma^2));
            elseif strcmp(strType,'Uniform')
                % 'params' = [minx maxx]
                minx = params(1);
                maxx = params(2);
                normalize = maxx-minx;
                func = @(x) (heaviside(x-minx) - heaviside(x-maxx)).*(1/normalize);
            end
        end
        
        
        
    end
    
end

