%{
This class describes a conductor with a electric current I, given by a series of
node points (x,y,z). 

The field at any (array) position can be calculated with method 'field_3d', with
the Biot-Savart integral.

SI units have to be used throughout.

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121015
%----------------------------------------------------------------------
%} 
classdef path_3d < magnetic_field.element_3d 
%----------------------------------------------------------------------    
    properties
        I; % intensity
        x; % xyz position of nodes
        y;
        z;
    end
    properties (Dependent = true)
        nodes; % n x 3 matrix with all node coordinates
        n_nodes; % total number of nodes
    end 
%----------------------------------------------------------------------
    methods % Constructor, set/get
        function h = path_3d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('I',1,@isnumeric);            
            p.addParameter('x',1:10,@isnumeric);
            p.addParameter('y',1:10,@isnumeric);
            p.addParameter('z',1:10,@isnumeric); 
            p.addParameter('n',[1:10;1:10;1:10].',@isnumeric); 
            % Parse and assign input
            p.parse(varargin{:});            
            h.I = p.Results.I;
            h.x = p.Results.x;
            h.y = p.Results.y;
            h.z= p.Results.z;
            if ~ismember('n',p.UsingDefaults)
                % Only if n has been given as an actual input assign it,
                % else use x,y,z
                h.n= p.Results.n;
            end 
        end
        function v = get.nodes(h)
            v = [h.x,h.y,h.z];
        end
        function set.nodes(h,v)
        % assumes that x,y,z are each one column
            h.x = v(:,1);
            h.y = v(:,1);
            h.z = v(:,1);
        end
        function v = get.n_nodes(h)
            v = length(h.x);
        end        
        function set.x(h,v) % force column vectors:
            h.x = v(:);
        end
        function set.y(h,v)
            h.y = v(:);
        end
        function set.z(h,v)
            h.z = v(:);
        end
    end
%----------------------------------------------------------------------        
    methods % Calculation methods
        function set_B0(h,B0)
        % Scales the whole field to make the field at the origin equal
        % to B0
            if ~exist('B0','var')
                B0 = 1; % if nothing provided, normalize to 1
            end
            h.I = h.I * B0/h.B_3d(0,0,0);        
        end
        function [Bx,By,Bz] = field_3d(h,x,y,z) 
        % calculates field at position with Biot-Savart            
            Bx = x.*0; % allocate
            By = Bx;
            Bz = Bx;            
            if h.I == 0    % Skip null-current paths
                return;
            end            
            C = h.I*h.const.mu0/(4*pi); % the multiplying constant            
            for i = 1:h.n_nodes-1 % number of intermediate points = h.n_nodes-1                
                dx = x-(h.x(i+1)+h.x(i))/2; % distance from xyz to intermediate point between nodes i+1 and i
                dy = y-(h.y(i+1)+h.y(i))/2;
                dz = z-(h.z(i+1)+h.z(i))/2;                
                dlx = h.x(i+1)-h.x(i); % length vector between nodes i+1 and i
                dly = h.y(i+1)-h.y(i);
                dlz = h.z(i+1)-h.z(i);                                
                dr_3_2 = (dx.^2 + dy.^2 + dz.^2).^(3/2);                
                good = (dr_3_2 ~= 0);                
                Bx(good) = Bx(good) + (dly*dz(good) - dlz*dy(good))./dr_3_2(good);
                By(good) = By(good) + (dlz*dx(good) - dlx*dz(good))./dr_3_2(good);
                Bz(good) = Bz(good) + (dlx*dy(good) - dly*dx(good))./dr_3_2(good);                
            end            
            Bx = C*Bx;
            By = C*By;
            Bz = C*Bz;            
        end
        function  [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = derivatives_3d(h,x,y,z)
            % Allocate
            dBx_dx = x.*0;
            dBx_dy = dBx_dx;
            dBx_dz = dBx_dx;
            dBy_dx = dBx_dx;
            dBy_dy = dBx_dx;
            dBy_dz = dBx_dx;
            dBz_dx = dBx_dx;
            dBz_dy = dBx_dx;
            dBz_dz = dBx_dx;            
            if h.I == 0    % Skip null-current paths
                return;
            end                        
            for i = 1:h.n_nodes-1 % number of intermediate points = h.n_nodes-1                
                dx = x-(h.x(i+1)+h.x(i))/2; % distance from xyz to intermediate point between nodes i+1 and i
                dy = y-(h.y(i+1)+h.y(i))/2;
                dz = z-(h.z(i+1)+h.z(i))/2;
                dd = sqrt(dx.^2+dy.^2+dz.^2);                
                dlx = h.x(i+1)-h.x(i); % length vector between nodes i+1 and i
                dly = h.y(i+1)-h.y(i);
                dlz = h.z(i+1)-h.z(i);
                dl = sqrt(dlx^2+dly^2+dlz^2);                
                % auxiliar magnitudes
                zeta = (dx.*dlx + dy.*dly + dz.*dlz)/dl; % zeta coordinate (axial) (see notebook II, p. 71)
                xi = sqrt(dd.^2 - zeta.^2); % xi coord (radial)
                etax1 = (dly*dz-dlz*dy)./(dl.*xi); % escalar product 1eta*1x
                etay1 = (dlz*dx-dlx*dz)./(dl.*xi); % escalar product 1eta*1y
                etaz1 = (dlx*dy-dly*dx)./(dl.*xi); % escalar product 1eta*1z
                xix1 = (etay1*dlz-etaz1*dly)/dl; % escalar product 1xi*1x
                xiy1 = (etaz1*dlx-etax1*dlz)/dl; % escalar product 1xi*1y
                xiz1 = (etax1*dly-etay1*dlx)/dl; % escalar product 1xi*1z
                % Auxiliar terms in cylindrical coords
                dBeta_dxi = dl*(1./dd.^3 - 3*xi.^2./dd.^5);
                dBeta_dzeta = -3*dl.*xi.*zeta./dd.^5;                
                % Projection of the grad tensor in xyz coords
                dBx_dx_ = (dBeta_dxi.*xix1 + dBeta_dzeta.*dlx./dl).*etax1 - ...
                         dl./dd.^3.*xix1.*etax1;
                dBx_dy_ = (dBeta_dxi.*xiy1 + dBeta_dzeta.*dly./dl).*etax1 - ...
                         dl./dd.^3.*xix1.*etay1;
                dBx_dz_ = (dBeta_dxi.*xiz1 + dBeta_dzeta.*dlz./dl).*etax1 - ...
                         dl./dd.^3.*xix1.*etaz1;
                dBy_dx_ = dBx_dy_;
                dBy_dy_ = (dBeta_dxi.*xiy1 + dBeta_dzeta.*dly./dl).*etay1 - ...
                         dl./dd.^3.*xiy1.*etay1;
                dBy_dz_ = (dBeta_dxi.*xiz1 + dBeta_dzeta.*dlz./dl).*etay1 - ...
                         dl./dd.^3.*xiy1.*etaz1;                     
                dBz_dx_ = dBx_dz_; 
                dBz_dy_ = dBy_dz_;
                dBz_dz_ = (dBeta_dxi.*xix1 + dBeta_dzeta.*dlx./dl).*etax1 - ...
                         dl./dd.^3.*xix1.*etax1;                     
                good = (dd ~= 0);     
                dBx_dx(good) = dBx_dx(good) + dBx_dx_(good);
                dBx_dy(good) = dBx_dy(good) + dBx_dy_(good);
                dBx_dz(good) = dBx_dz(good) + dBx_dz_(good);
                dBy_dx(good) = dBy_dx(good) + dBy_dx_(good);
                dBy_dy(good) = dBy_dy(good) + dBy_dy_(good);
                dBy_dz(good) = dBy_dz(good) + dBy_dz_(good);
                dBz_dx(good) = dBz_dx(good) + dBz_dx_(good);
                dBz_dy(good) = dBz_dy(good) + dBz_dy_(good);
                dBz_dz(good) = dBz_dz(good) + dBz_dz_(good);
            end                        
            C = h.I*h.const.mu0/(4*pi); % the multiplying constant            
            dBx_dx = C*dBx_dx;
            dBx_dy = C*dBx_dy;
            dBx_dz = C*dBx_dz;
            dBy_dx = C*dBy_dx;
            dBy_dy = C*dBy_dy;
            dBy_dz = C*dBy_dz;
            dBz_dx = C*dBz_dx;
            dBz_dy = C*dBz_dy;
            dBz_dz = C*dBz_dz;            
        end
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_3d(h,varargin)
            h_line = line('XData',h.x,'YData',h.y,'ZData',h.z,'DisplayName',inputname(1),varargin{:});
        end
    end
%----------------------------------------------------------------------        
end

