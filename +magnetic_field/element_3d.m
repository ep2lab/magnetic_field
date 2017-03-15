%{
This is the abstract parent class for all 3d field generator elements.

It defines the basic interfaces and utilities
 
%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121011
%----------------------------------------------------------------------
%}
classdef element_3d < magnetic_field.element 
%----------------------------------------------------------------------    
    methods (Abstract = true)
        [Bx,By,Bz] = field_3d(h,x,y,z)
        [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = derivatives_3d(h,x,y,z)
        h_line = plot_generator_3d(h,varargin)
    end
%----------------------------------------------------------------------
    methods % Computation methods 
        function Bx = Bx_3d(h,x,y,z)
        % Shortcut function for Bx
            [Bx,~,~] = h.field_3d(x,y,z);
        end 
        function By = By_3d(h,x,y,z)
        % Shortcut function for Bx
            [~,By,~] = h.field_3d(x,y,z);
        end 
        function Bz = Bz_3d(h,x,y,z)
        % Shortcut function for Bx
            [~,~,Bz] = h.field_3d(x,y,z);
        end
        function B = B_3d(h,x,y,z)
        % Shortcut function for B
            [Bx,By,Bz] = h.field_3d(x,y,z);
            B = sqrt(Bx.^2 + By.^2 + Bz.^2);
        end        
        function dBx_dx = dBx_dx_3d(h,z,r)
            % Shortcut function for one of the derivatives
            [dBx_dx,~,~,~,~,~,~,~,~] = h.derivatives_3d(z,r);
        end
        function dBx_dy = dBx_dy_3d(h,z,r)
            % Shortcut function for one of the derivatives
            [~,dBx_dy,~,~,~,~,~,~,~] = h.derivatives_3d(z,r);
        end
        function dBx_dz = dBx_dz_3d(h,z,r)
            % Shortcut function for one of the derivatives
            [~,~,dBx_dz,~,~,~,~,~,~] = h.derivatives_3d(z,r);
        end
        function dBy_dx = dBy_dx_3d(h,z,r)
            % Shortcut function for one of the derivatives
            [~,~,~,dBy_dx,~,~,~,~,~] = h.derivatives_3d(z,r);
        end
        function dBy_dy = dBy_dy_3d(h,z,r)
            % Shortcut function for one of the derivatives
            [~,~,~,~,dBy_dy,~,~,~,~] = h.derivatives_3d(z,r);
        end
        function dBy_dz = dBy_dz_3d(h,z,r)
            % Shortcut function for one of the derivatives
            [~,~,~,~,~,dBy_dz,~,~,~] = h.derivatives_3d(z,r);
        end
        function dBz_dx = dBz_dx_3d(h,z,r)
            % Shortcut function for one of the derivatives
            [~,~,~,~,~,~,dBz_dx,~,~] = h.derivatives_3d(z,r);
        end
        function dBz_dy = dBz_dy(h,z,r)
            % Shortcut function for one of the derivatives
            [~,~,~,~,~,~,~,dBz_dy,~] = h.derivatives_3d(z,r);
        end
        function dBz_dz = dBz_dz_3d(h,z,r)
        % Shortcut function for one of the derivatives
            [~,~,~,~,~,~,~,~,dBz_dz] = h.derivatives_3d(z,r);
        end     
        function kB = curvature_3d(h,XX,YY,ZZ)
        % Computes curvature of magnetic lines. The sign is according to
        % the direction of n_vector
            h;
            kB = XX+YY+ZZ*NaN;
            error('Not yet implemented!');
        end
        function B_vector = B_vector_3d(h,x,y,z)
        % Returns the magnetic field vector at x,y,z. B_vector is laid
        % along a new last dimension  
            input_size = size(x); % find non-singleton size of input
            for i = length(input_size):-1:1
                if input_size(i) == 1
                    input_size = input_size(1:i-1);
                else
                    break
                end
            end
            B_vector(cat(2,input_size,3)) = 0; % allocate
            [bx,by,bz] = h.field_3d(x,y,z);             
            B_vector(1:numel(x)) = bx;
            B_vector(1+numel(x):2*numel(x)) = by; 
            B_vector(1+2*numel(x):end) = bz; 
        end
        function b_vector = b_vector_3d(h,x,y,z)
        % Returns the normalized (i.e. unitary) magnetic field vector at
        % x,y,z. B_vector is laid along a new last dimension  
            [b_vector] = h.B_vector_3d(x,y,z);
            % normalize each vector
            B = sqrt(b_vector(1:numel(x)).^2+b_vector(1+numel(x):2*numel(x)).^2+b_vector(1+2*numel(x):end).^2);
            b_vector(1:numel(x)) = b_vector(1:numel(x))./B;
            b_vector(1+numel(x):end) = b_vector(1+numel(x):end)./B; 
        end
        function n_vector = n_vector_3d(h,XX,YY,ZZ)
        % Computes curvature of magnetic lines
            h;
            n_vector = XX+YY+ZZ*NaN;
            error('Not yet implemented!');
        end
        function [X,Y,Z] = streamline_3d(h,x,y,z,ds,n_steps,odeoptions)
        % Calculates the streamline(s) through x,y,z, which can be a single
        % point or an array of points. The streamline is calculated on steps
        % ds and for a total n_steps with ode45.
        % The output are the arrays of points of the streamline(s), added in
        % a new last dimension wrt x,y,z arrays.
            % Input parsing
            if ~exist('odeoptions','var')
                odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'InitialStep',0.1);
            end
            % Allocate output
            input_size = size(x); % find non-singleton size of input
            for i = length(input_size):-1:1
                if input_size(i) == 1
                    input_size = input_size(1:i-1);
                else
                    break
                end
            end
            X(cat(2,input_size,n_steps)) = 0; 
            Y = X; 
            Z = X;
            % Integration
            b = @(t,p) h.b_vector_3d(p(1),p(2),p(3)).';
            if n_steps>2
                s_span = linspace(0,ds*(n_steps-1),n_steps).';
            else
                s_span = [0,ds,2*ds]; % add an extra, ignored point because otherwise ode45 does not work
            end            
            for i = 1:numel(x)
                [~,p] = ode45(b,s_span,[x(i);y(i);z(i)],odeoptions); 
                for j = 1:n_steps
                    X(i+(j-1)*numel(x)) = p(j,1);
                    Y(i+(j-1)*numel(x)) = p(j,2);
                    Z(i+(j-1)*numel(x)) = p(j,3);                
                end
            end 
        end       
        function [X,Y,Z] = next_point_3d(h,x,y,z,ds,odeoptions)
        % Computes next point along streamline(s) that pass by x,y,z
        % using streamline_3d. Output is the new point(s), a curved
        % distance ds downstream 
            % Input parsing
            if ~exist('odeoptions','var')
                odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'InitialStep',0.1);
            end
            % Call streamline with just one new step
            [X,Y,Z] = h.streamline_3d(x,y,z,ds,2,odeoptions);
            % Keep only new step
            X = X(1+numel(x):end);
            Y = Y(1+numel(x):end);
            Z = Z(1+numel(x):end);
        end       
    end
%----------------------------------------------------------------------    
    methods % Plotting methods           
        function [h_surface] = plot_3d(h,varargin)
        % Plots 'var' on a plane defined by origin vector o, and tangent
        % vectors v1,v2, in current axes. Can pass extra name:value
        % pairs with surface options 
            % Validate input and defaults
            p = inputParser;               
            p.addParameter('var','sqrt(Bx.^2 + By.^2 + Bz.^2)',@ischar);            
            p.addParameter('o',[0,0,0],@isnumeric);
            p.addParameter('v1',[1,0,0],@isnumeric);
            p.addParameter('v2',[0,0,1],@isnumeric);
            p.addParameter('t1Lim',[-5,5],@isnumeric);
            p.addParameter('t2Lim',[-5,5],@isnumeric);
            p.addParameter('t1points',300,@isnumeric);
            p.addParameter('t2points',320,@isnumeric); 
            % Parse and assign input
            p.KeepUnmatched = true;
            p.parse(varargin{:});            
            var = p.Results.var;
            o = p.Results.o;
            v1= p.Results.v1;
            v2= p.Results.v2;            
            t1Lim = p.Results.t1Lim;
            t2Lim = p.Results.t2Lim;
            t1points = p.Results.t1points;
            t2points = p.Results.t2points; 
            surfaceopts = p.Unmatched; % anything else are surface options
            % Meshgrid and field arrays             
            [x,y,z] = utilities.meshgrid_inclined_plane(o,v1,v2,t1Lim,t2Lim,t1points,t2points); % !!! requires matlabtools                        
            [Bx,By,Bz] = h.field_3d(x,y,z); 
            % Compute variable to plot and call surface
            quantity = eval(var); % !!! Potential security risk
            h_surface = surface(x,y,z,quantity,'Displayname',[var,' of ',inputname(1)],'linestyle','none',surfaceopts);     
        end
    end
%----------------------------------------------------------------------
end

