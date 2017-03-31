%{
This class describes the magnetic field of a circular loop analytically
in 3d. It inherits not only from element_3d but also from loop_2d, to
use (internally) the same methods for field computation as in there.

Look local axees are referred to the global axes by the local origin
'origin' and the loop axis unit vector 'axis'. You are adviced to keep
ZL = 0 to avoid reference frame transformation confusions

SI units have to be used throughout.

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121015
%----------------------------------------------------------------------
%}
classdef loop_3d < magnetic_field.element_3d & magnetic_field.loop_2d
%----------------------------------------------------------------------    
    properties         
        axis % axis unitary vector of the loop
        origin % origin of local axes wrt global.
    end
    properties (Hidden = true)
        plot_generator_points
    end 
%----------------------------------------------------------------------
    methods % Constructor, set/get           
        function h = loop_3d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('I',1,@isnumeric);            
            p.addParameter('ZL',0,@isnumeric);
            p.addParameter('RL',3.5,@isnumeric);
            p.addParameter('axis',[0;0;1],@isnumeric);
            p.addParameter('origin',[0;0;0],@isnumeric);
            p.addParameter('plot_generator_points',100,@isnumeric);
            % Parse and assign input
            p.parse(varargin{:});            
            h.I = p.Results.I;
            h.ZL = p.Results.ZL;
            h.RL = p.Results.RL;
            h.axis = p.Results.axis;
            h.origin = p.Results.origin;
            h.plot_generator_points = p.Results.plot_generator_points;
        end 
        function set.axis(h,v)
            h.axis = v(:)/norm(v,2); % just normalize it, so it is unitary. Force column
        end
        function set.origin(h,v)
            h.origin= v(:); % force column.
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
        % Compute the magnetic field of the loop analytically with
        % field_2d, then rotate and translate it to 3d axes
            % Allocate
            Bx = x.*0;
            By = Bx;
            Bz = Bx;                                    
            if h.I == 0    % Skip null-current loop
                return;
            end            
            % transform to local axes
            x = x(:)-h.origin(1);
            y = y(:)-h.origin(2);
            z = z(:)-h.origin(3);            
            points = [x, y, z]; % construct array
            Z = points * h.axis; % Z, distance along axis, is just the projection on the unitary vector h.axis
            points = points - Z * h.axis.'; % here points becomes the perpe vectors from the axis to the points
            R = sqrt(sum(points.*points,2));            
            % Calculate field in local axes
            [~,Bz_,Br_] = h.field_2d(Z,R);            
            % Express Bz_, Br_ in global axes. Br_ only when R .ne. 0:
            good = logical(R ~= 0);
            Bx(good) = Br_(good,1) .* points(good,1) ./ R(good,1); % the extra indices are necessary for matlab
            By(good) = Br_(good,1) .* points(good,2) ./ R(good,1); 
            Bz(good) = Br_(good,1) .* points(good,3) ./ R(good,1);            
            % now add the Bz_ contribution
            Bx(:) = Bx(:) + Bz_(:) .* h.axis(1);
            By(:) = By(:) + Bz_(:) .* h.axis(2);
            Bz(:) = Bz(:) + Bz_(:) .* h.axis(3);                      
        end
        function [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = derivatives_3d(h,x,y,z)
        % Calculates dervatives of the field using the 2d code, then
        % rotating and translating
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
            if h.I == 0    % Skip null-current loop
                return;
            end            
            % transform to local axes
            points = [x(:)-h.origin(1), y(:)-h.origin(2), z(:)-h.origin(3)];
            Z = points * h.axis; % Z, distance along axis, is just the projection on the unitary vector h.axis
            points = points - Z * h.axis.'; % here points becomes the perpe vectors from the axis to the points
            R = sqrt(sum(points.*points,2));            
            good = logical(R ~= 0);
            points(good,1) = points(good,1) ./ R(good); % this gives the unit vectors in the R direction. At the axis, it is 0
            points(good,2) = points(good,2) ./ R(good); 
            points(good,3) = points(good,3) ./ R(good); 
            points(~good,:) = 0;            
            % calculate 2d derivatives and BR field (in local axes)
            [dBZ_dZ,dBZ_dR,dBR_dZ,dBR_dR] = h.derivatives_2d(Z,R);
            [~,~,BR] = h.field_2d(Z,R);
            BR_R(good) = BR(good)./R(good);
            BR_R(~good) = dBR_dR(~good);            
            % Transform into global axes derivatives of the field. 
            Zx1 = h.axis(1); % escalar product 1Z*1x. These are the projections of Z and R versors along X,Y,Z.
            Zy1 = h.axis(2); % escalar product 1Z*1y
            Zz1 = h.axis(3); % escalar product 1Z*1z            
            Rx1 = points(:,1); % escalar product 1R*1x
            Ry1 = points(:,2); % escalar product 1R*1y
            Rz1 = points(:,3); % escalar product 1R*1z            
            Tx1 = h.axis(2)*points(:,3)-h.axis(3)*points(:,2); % escalar product 1THETA*1x
            Ty1 = h.axis(3)*points(:,1)-h.axis(1)*points(:,3); % escalar product 1THETA*1y
            Tz1 = h.axis(1)*points(:,2)-h.axis(2)*points(:,1); % escalar product 1THETA*1z
            % Compute 3d derivatives
            dBx_dx(:) = (dBZ_dZ.*Zx1 + dBZ_dR.*Rx1).*Zx1 +...
                        (dBR_dZ.*Zx1 + dBR_dR.*Rx1).*Rx1 +...
                        BR_R.*Tx1.*Tx1;
            dBx_dy(:) = (dBZ_dZ.*Zy1 + dBZ_dR.*Ry1).*Zx1 +...
                        (dBR_dZ.*Zy1 + dBR_dR.*Ry1).*Rx1 +...
                        BR_R.*Ty1.*Tx1;
            dBx_dz(:) = (dBZ_dZ.*Zz1 + dBZ_dR.*Rz1).*Zx1 +...
                        (dBR_dZ.*Zz1 + dBR_dR.*Rz1).*Rx1 +...
                        BR_R.*Tz1.*Tx1;                    
            dBy_dx = dBx_dy; % use symmetries (rot(B) = 0)
            dBy_dy(:) = (dBZ_dZ.*Zy1 + dBZ_dR.*Ry1).*Zy1 +...
                        (dBR_dZ.*Zy1 + dBR_dR.*Ry1).*Ry1 +...
                        BR_R.*Ty1.*Ty1;
            dBy_dz(:) = (dBZ_dZ.*Zz1 + dBZ_dR.*Rz1).*Zy1 +...
                        (dBR_dZ.*Zz1 + dBR_dR.*Rz1).*Ry1 +...
                        BR_R.*Tz1.*Ty1;                    
            dBz_dx = dBx_dz;
            dBz_dy = dBy_dz;
            dBz_dz(:) = (dBZ_dZ.*Zz1 + dBZ_dR.*Rz1).*Zz1 +...
                        (dBR_dZ.*Zz1 + dBR_dR.*Rz1).*Rz1 +...
                        BR_R.*Tz1.*Tz1;            
        end 
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_3d(h,varargin)
            % First, construct a vector which is perpendicular to h.axis
            [~,i] = max((abs(h.axis))); % take the maximum abs value, better
            p(1:3) = h.axis(i); 
            p(i) = -h.axis(mod(i,3)+1) -h.axis(mod(i+1,3)+1); % from dotproduct = 0
            p = p(:)/norm(p,2); % normalize
            % Second perpe vector
            q = cross(h.axis,p);
            % Find the points of the loop
            theta = linspace(0,2*pi,h.plot_generator_points);
            points(3,h.plot_generator_points) = 0; % allocate
            for i = 1:h.plot_generator_points
            	points(:,i) = h.RL*(p*cos(theta(i)) + q*sin(theta(i))) + h.origin;
            end
            % Plot                       
            h_line = line('XData',points(1,:),'YData',points(2,:),'ZData',points(3,:),'DisplayName',inputname(1));
            set(h_line,varargin{:});
        end
    end   
%----------------------------------------------------------------------
end
