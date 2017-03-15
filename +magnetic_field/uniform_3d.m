%{ 
This class describes a uniform magnetic field oriented according to a 3d
vector 'direction'

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121015
%----------------------------------------------------------------------
%}
classdef uniform_3d < magnetic_field.element_3d
%----------------------------------------------------------------------    
    properties   
        B0 % magnetic field intensity
        direction % direction of the field
    end
%----------------------------------------------------------------------
    methods % Constructor, set/get
        function h = uniform_3d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('B0',1,@isnumeric);            
            p.addParameter('direction',[0,0,1],@isnumeric);
            % Parse and assign input
            p.parse(varargin{:});            
            h.B0 = p.Results.B0;
            h.direction = p.Results.direction;
        end
        function set.direction(h,v)
            h.direction = v(:)/norm(v,2); % normalize it, so it is unitary. Force column
        end 
    end
%----------------------------------------------------------------------    
    methods % Calculation methods
        function set_B0(h,B0)
        % Scales the whole field to make the field equal to B0
            if ~exist('B0','var')
                B0 = 1; % if nothing provided, normalize to 1
            end
            h.B0 = B0;
        end
        function [Bx,By,Bz] = field_3d(h,~,~,z)            
            Bx = h.B0*h.direction(1) + z.*0;
            By = h.B0*h.direction(2) + z.*0;
            Bz = h.B0*h.direction(3) + z.*0;                                   
        end
        function [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = derivatives_3d(~,x,~,~)
            dBx_dx = x.*0;
            dBx_dy = dBx_dx;
            dBx_dz = dBx_dx;
            dBy_dx = dBx_dx;
            dBy_dy = dBx_dx;
            dBy_dz = dBx_dx;
            dBz_dx = dBx_dx;
            dBz_dy = dBx_dx;
            dBz_dz = dBx_dx;               
        end 
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_3d(~,varargin) 
            % Plot nothing
            h_line = [];
        end
    end           
%----------------------------------------------------------------------
end
