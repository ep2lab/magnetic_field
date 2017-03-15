%{ 
This class describes a uniform magnetic field. 

MMM20121015
%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121015
%----------------------------------------------------------------------
%}
classdef uniform_2d < magnetic_field.element_2d
%----------------------------------------------------------------------
    properties 
        B0 % magnetic field intensity
        is_axi
    end
%----------------------------------------------------------------------
    methods % Constructor methods
        function h = uniform_2d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('B0',1,@isnumeric);            
            p.addParameter('is_axi',1,@isnumeric);
            % Parse and assign input
            p.parse(varargin{:});            
            h.B0 = p.Results.B0;
            h.is_axi = p.Results.is_axi;
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
        function [psi,Bz,Br] = field_2d(h,z,r) 
        % calculates field in local axes            
            Br = z.*0;
            Bz = h.B0 + Br;            
            if h.is_axi % axisymmetric
                psi = 0.5 * r.^2 * h.B0;          
            else % planar
                psi = r * h.B0;
            end
        end 
        function [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = derivatives_2d(~,z,~) 
        % Calculates dervatives of Bz, Br on z,r
            dBz_dz = z.*0;
            dBz_dr = dBz_dz;
            dBr_dz = dBz_dz;
            dBr_dr = dBz_dz;                        
        end 
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_2d(~,varargin) 
            % Plot nothing
            h_line = [];
        end
    end   
%----------------------------------------------------------------------
end
