%{
This class describes the magnetic field of an infinite straight
wire analytically in the 3d case. It inherits from element_3d and
wire_2d

The wire can be positioned and oriented with a point vector and a
direction vector. ZW and RW should be kept equal to zero.

SI units have to be used throughout.

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170314
%----------------------------------------------------------------------
%}
classdef wire_3d < magnetic_field.element_3d & magnetic_field.wire_2d
%----------------------------------------------------------------------    
    properties         
        direction % direction unitary vector of the wire
        point % a point of the wire
    end
    properties (Hidden = true)
        plot_generator_tLim % how much distance, from point, to plot
    end 
%----------------------------------------------------------------------
    methods % Constructor, set/get           
        function h = wire_3d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('I',1,@isnumeric);        
            p.addParameter('direction',[0;0;1],@isnumeric);
            p.addParameter('point',[0;3.5;0],@isnumeric);
            p.addParameter('plot_generator_tLim',[-10,10],@isnumeric);
            % Parse and assign input
            p.parse(varargin{:});            
            h.I = p.Results.I; 
            h.direction = p.Results.direction;
            h.point = p.Results.point;
            h.plot_generator_tLim = p.Results.plot_generator_tLim;
            h.RW = 0; % meaningless otherwise in the 3d case
            h.ZW = 0;
        end 
        function set.direction(h,v)
            h.direction = v(:)/norm(v,2); % just normalize it, so it is unitary. Force column
        end
        function set.point(h,v)
            h.point= v(:); % force column.
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
        % Compute the magnetic field of the wire analytically with
        % field_2d, then rotate and translate it to 3d axes
            % Allocate
            Bx = x.*0;
            By = Bx;
            Bz = Bx;                                    
            if h.I == 0    % Skip null-current wire
                return;
            end         
            if h.ZW ~= 0 || h.RW ~= 0
                error('ZW and RW should be 0 in wire_3d')
            end
            % transform to local axes
            x = x(:)-h.point(1);
            y = y(:)-h.point(2);
            z = z(:)-h.point(3);
            points = [x, y, z]; % construct array
            THETA = points * h.direction; % Z, distance in the direction parallel to the wire
            points = points - THETA * h.direction.'; % here points becomes the perpe vectors from the wire to the points
            RHO = sqrt(sum(points.*points,2)); 
            DIRECTION = points.*0;
            DIRECTION(:,1) = h.direction(1);
            DIRECTION(:,2) = h.direction(2);
            DIRECTION(:,3) = h.direction(3);
            CROSS = cross(points,DIRECTION,2); % direction of B field at each point
            [~,Bz_,~] = h.field_2d(RHO.*0,RHO);                 
            % Translate
            Bx(:) = Bz_(:) .* CROSS(:,1)./sqrt(sum(CROSS.*CROSS,2));
            By(:) = Bz_(:) .* CROSS(:,2)./sqrt(sum(CROSS.*CROSS,2));
            Bz(:) = Bz_(:) .* CROSS(:,3)./sqrt(sum(CROSS.*CROSS,2));
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
            if h.I == 0    % Skip null-current wire
                return;
            end               
            if h.ZW ~= 0 || h.RW ~= 0
                error('ZW and RW should be 0 in wire_3d')
            end    
            % transform to local axes
            x = x(:)-h.point(1);
            y = y(:)-h.point(2);
            z = z(:)-h.point(3);
            points = [x, y, z]; % construct array
            THETA = points * h.direction; % Z, distance in the direction parallel to the wire
            points = points - THETA * h.direction.'; % here points becomes the perpe vectors from the wire to the points
            RHO = sqrt(sum(points.*points,2));            
            DIRECTION = points.*0;
            DIRECTION(:,1) = h.direction(1);
            DIRECTION(:,2) = h.direction(2);
            DIRECTION(:,3) = h.direction(3);
            CROSS = cross(points,DIRECTION); % direction of B field at each point          
            % calculate 2d derivatives and BR field (in local axes)
            [~,dBZ_dR,~,~] = h.derivatives_2d(RHO.*0,RHO);    
            % Transform into global axes derivatives of the field. 
            Zx1 = CROSS(:,1); % escalar product 1Z*1x. These are the projections of Z and R versors along X,Y,Z.
            Zy1 = CROSS(:,2); % escalar product 1Z*1y
            Zz1 = CROSS(:,3); % escalar product 1Z*1z            
            Rx1 = points(:,1); % escalar product 1R*1x
            Ry1 = points(:,2); % escalar product 1R*1y
            Rz1 = points(:,3); % escalar product 1R*1z            
            % Compute 3d derivatives
            dBx_dx(:) = 2*dBZ_dR.*Rx1.*Zx1;
            dBx_dy(:) = 2*dBZ_dR.*Ry1.*Zx1;
            dBx_dz(:) = 2*dBZ_dR.*Rz1.*Zx1; 
            dBy_dx = dBx_dy; % use symmetries (rot(B) = 0)
            dBy_dy(:) = dBZ_dR.*Ry1.*Zy1;
            dBy_dz(:) = dBZ_dR.*Rz1.*Zy1; 
            dBz_dx = dBx_dz;
            dBz_dy = dBy_dz;
            dBz_dz = -dBx_dx -dBy_dy; % Div(B) = 0
        end 
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_3d(h,varargin) 
            h_line = line('XData',h.point(1)+h.direction(1)*h.plot_generator_tLim,...
                          'YData',h.point(2)+h.direction(2)*h.plot_generator_tLim,...
                          'ZData',h.point(3)+h.direction(3)*h.plot_generator_tLim,...
                          'DisplayName',inputname(1));
            set(h_line,varargin{:});
        end
    end   
%----------------------------------------------------------------------
end
