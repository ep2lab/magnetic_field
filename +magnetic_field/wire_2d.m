%{
This class describes the magnetic field of an infinite straight
wire analytically. It inherits from element_2d.

The wire has intensity I, and coordinates in the plane ZW, RW

SI units have to be used throughout.

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170314
%----------------------------------------------------------------------
%}
classdef wire_2d < magnetic_field.element_2d
%----------------------------------------------------------------------    
    properties
        I  % intensity
        ZW % axial position of loop
        RW % radius of loop         
        is_axi = 0 % whether axisimetric (1) or planar (0) geometry is used
    end 
%----------------------------------------------------------------------    
    methods % Constructor
        function h = wire_2d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('I',1,@isnumeric);            
            p.addParameter('ZW',0,@isnumeric);
            p.addParameter('RW',3.5,@isnumeric); 
            % Parse and assign input
            p.parse(varargin{:});            
            h.I = p.Results.I;
            h.ZW = p.Results.ZW;
            h.RW = p.Results.RW; 
        end
        function set.is_axi(h,v)
            if v ~= 0
                error('Wires cannot be set to axisymmetric mode, they are always planar');
            end
            h.is_axi = v;               
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
            h.I = h.I * B0/h.B_2d(0,0);
        end
        function [psi,Bz,Br] = field_2d(h,z,r) 
        % calculates field of single wire analytically. The code
        % has been optimized for speed
            % allocate
            psi = z.*0; 
            Bz = psi;
            Br = psi;            
            if h.I == 0    % Skip null-current wire
                return;
            end
            % Multiplying constant
            C = h.I*h.const.mu0/(4*pi);     
            % Computation
            z=z-h.ZW;
            r=r-h.RW;
            rho2 = z.^2+r.^2;
            psi = -C/2*log(rho2);
            Bz = -C*r./rho2;
            Br = C*z./rho2;           
            % Return Inf, zeros if (z,r) coincides with the wire
            psi(z==0 & r==0) = Inf;
            Bz(z==0 & r==0)  = 0;
            Br(z==0 & r==0)  = 0;  
        end 
        function [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = derivatives_2d(h,z,r) 
        % Calculates dervatives of Bz, Br on z,r            
            % allocate
            dBz_dz = z.*0;
            dBz_dr = dBz_dz;
            dBr_dz = dBz_dz;
            dBr_dr = dBz_dz;                        
            if h.I == 0    % Skip null-current wire
                return;
            end            
            % Multiplying constant
            C = h.I*h.const.mu0/(2*pi);    
            % Computation
            z=z-h.ZW;
            r=r-h.RW;
            rho2 = z.^2+r.^2;
            dBz_dz = 2*C*z.*r./rho2;
            dBz_dr = 2*C*r.^2./rho2 - C./sqrt(rho2);  
            % Return zeros if (z,r) coincides with the wire
            dBz_dz(z==0 & r==0) = 0;
            dBz_dr(z==0 & r==0)  = 0;
            % Compute the othres from conditions 
            dBr_dr = - dBz_dz; % Divergence-free condition
            dBr_dz = dBz_dr; % Curl-free condition     
        end 
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_2d(h,varargin)  
            h_line = line('XData',h.ZW,'YData',h.RW,'DisplayName',inputname(1));
            set(h_line,'Marker','s','linestyle','none',varargin{:});
        end
    end    
%----------------------------------------------------------------------    
end
