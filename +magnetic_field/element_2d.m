%{
This is the abstract parent class for all 2d field generator elements.
Useful only for axisymmetric or planar fields with no Btheta (or By)
component.

This class defines the basic interfaces and utilities

[psi,Bz,Br] = field_2d(h,z,r) queries the field components in local axes
[dBz_dz,dBz_dr,dBr_dz,dBr_dr] = derivatives_2d(h,z,r) queries the
derivatives 
h_line = plot_generator_2d(h,varargin) plots the position of the
generator

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121011
%----------------------------------------------------------------------
%}
classdef element_2d < magnetic_field.element 
%----------------------------------------------------------------------    
    properties (Abstract = true)
        is_axi % whether axisimetric (1) or planar (0) geometry is used
    end
%----------------------------------------------------------------------    
    methods (Abstract = true)
        [psi,Bz,Br] = field_2d(h,z,r)
        [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = derivatives_2d(h,z,r)
        h_line = plot_generator_2d(h,varargin)
    end
%----------------------------------------------------------------------    
    methods % Computation methods
        function psi = psi_2d(h,z,r)
        % Shortcut function for psi
            [psi,~,~] = h.field_2d(z,r);
        end
        function Bz = Bz_2d(h,z,r)
        % Shortcut function for Bz
            [~,Bz,~] = h.field_2d(z,r);
        end
        function Br = Br_2d(h,z,r)
        % Shortcut function for Br
            [~,~,Br] = h.field_2d(z,r);
        end
        function B = B_2d(h,z,r)
        % Shortcut function for B
            [~,Bz,Br] = h.field_2d(z,r);
            B = sqrt(Bz.^2 + Br.^2);
        end
        function dBz_dz = dBz_dz_2d(h,z,r)
        % Shortcut function for one of the derivatives
            [dBz_dz,~,~,~] = h.derivatives_2d(z,r);
        end
        function dBz_dr = dBz_dr_2d(h,z,r)
        % Shortcut function for one of the derivatives
            [~,dBz_dr,~,~] = h.derivatives_2d(z,r);
        end
        function dBr_dz = dBr_dz_2d(h,z,r)
        % Shortcut function for one of the derivatives
            [~,~,dBr_dz,~] = h.derivatives_2d(z,r);
        end
        function dBr_dr = dBr_dr_2d(h,z,r)
        % Shortcut function for one of the derivatives
            [~,~,~,dBr_dr] = h.derivatives_2d(z,r);
        end
        function kB = curvature_2d(h,ZZ,RR)
        % Computes curvature of magnetic lines. The sign is according to
        % the direction of n_vector
            [~,Bz,Br] = h.field_2d(ZZ,RR);
            [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = h.derivatives_2d(ZZ,RR);            
            dz = Br.*dBz_dr + Bz.*dBz_dz;
            dr = Br.*dBr_dr + Bz.*dBr_dz;            
            kB = (-dz.*Br+dr.*Bz)./ sqrt(Bz.^2+Br.^2).^3;
        end
        function B_vector = B_vector_2d(h,z,r)
        % Returns the magnetic field vector at z,r. B_vector is laid
        % along a new last dimension  
            input_size = size(z); % find non-singleton size of input
            for i = length(input_size):-1:1
                if input_size(i) == 1
                    input_size = input_size(1:i-1);
                else
                    break
                end
            end
            B_vector(cat(2,input_size,2)) = 0; % allocate
            [~,bz,br] = h.field_2d(z,r);             
            B_vector(1:numel(z)) = bz;
            B_vector(1+numel(z):end) = br; 
        end
        function b_vector = b_vector_2d(h,z,r)
        % Returns the normalized (i.e. unitary) magnetic field vector at
        % z,r. B_vector is laid along a new last dimension  
            [b_vector] = h.B_vector_2d(z,r);
            % normalize each vector
            B = sqrt(b_vector(1:numel(z)).^2+b_vector(1+numel(z):end).^2);
            b_vector(1:numel(z)) = b_vector(1:numel(z))./B;
            b_vector(1+numel(z):end) = b_vector(1+numel(z):end)./B; 
        end
        function n_vector = n_vector_2d(h,z,r)
        % Returns the normalized (i.e. unitary) normal vector rotated 90
        % deg from b_vector
            [n_vector] = h.b_vector_2d(z,r);
            % rotate 90 deg            
            temp = n_vector(1+numel(z):end);
            n_vector(1+numel(z):end) = n_vector(1:numel(z));
            n_vector(1:numel(z)) = -temp;
        end
        function [Z,R] = streamline_2d(h,z,r,ds,n_steps)
        % Calculates the streamline(s) through z,r, which can be a single
        % point or an array of points. The streamline is calculated on steps
        % ds and for a total n_steps
        % The output are the arrays of points of the streamline(s), added in
        % a new last dimension wrt z,r arrays.
            % Allocate
            input_size = size(z); % find non-singleton size of input
            for i = length(input_size):-1:1
                if input_size(i) == 1
                    input_size = input_size(1:i-1);
                else
                    break
                end
            end            
            Z(cat(2,input_size,n_steps)) = 0; 
            R = Z; 
            % Initial points
            for i = 1:numel(z)
                Z(i) = z(i);
                R(i) = r(i);
            end
            % Rest of points
            for i = 1:numel(Z)-numel(z) % This way it always uses points that have been already computed                
                [Z(i+numel(z)),R(i+numel(z))] = h.next_point_2d(Z(i),R(i),ds);
            end            
        end        
        function [z,r] = next_point_2d(h,z,r,ds)
        % High-accuracy computation of next point on streamline(s)
        % through z,r, based on streamfunction psi. Output is the new
        % point(s), a straight distance ds downstream
            for i = 1:numel(z);
                if r(i) == 0 % separate axis case
                    z(i) = z(i) + ds;
                    continue;
                end
                % Initial approximation to new point
                [psi0,Bz,Br] = h.field_2d(z(i),r(i)); 
                % Compute theta (angle to new point)
                theta = fzero(@(theta)psi0-h.field_2d(z(i)+ ds*cos(theta),r(i)+ds*sin(theta)),atan2(Br,Bz));
                % Assign
                z(i) = z(i)+ ds*cos(theta);
                r(i) = r(i)+ ds*sin(theta);  
            end
        end        
    end
%----------------------------------------------------------------------    
    methods % Plotting methods 
        function [h_surface] = plot_2d(h,varargin)
        % Plots 'var' on a plane in current axes. Can pass extra
        % name:value pairs with surface options
            % Validate input and defaults
            p = inputParser;               
            p.addParameter('var','sqrt(Bz.^2 + Br.^2)',@ischar);            
            p.addParameter('ZLim',[0,10],@isnumeric);
            p.addParameter('RLim',[0,8],@isnumeric);
            p.addParameter('zpoints',300,@isnumeric);
            p.addParameter('rpoints',200,@isnumeric); 
            % Parse and assign input
            p.KeepUnmatched = true;
            p.parse(varargin{:});            
            var = p.Results.var;
            ZLim = p.Results.ZLim;
            RLim = p.Results.RLim;
            zpoints = p.Results.zpoints;
            rpoints = p.Results.rpoints; 
            surfaceopts = p.Unmatched; % anything else are surface options   
            % Meshgrid and field arrays 
            [z,r] = meshgrid(linspace(ZLim(1),ZLim(2),zpoints),linspace(RLim(1),RLim(2),rpoints));                
            [psi,Bz,Br] = h.field_2d(z,r);
            % Compute variable to plot and call surface
            quantity = eval(var); % !!! Potential security risk
            h_surface = surface(z,r,quantity*0,quantity,'Displayname',[var,' of ',inputname(1)],'linestyle','none',surfaceopts);     
        end
    end
%----------------------------------------------------------------------
end

