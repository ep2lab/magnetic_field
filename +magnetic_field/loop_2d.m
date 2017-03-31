%{
This class describes the magnetic field of a circular loop analytically. 
It inherits from element_2d.

The loop has intensity I, radius RL, centered at ZL,
and is oriented positively along the +Z axis.  
You can define ellipKE = 1 for quick (interpolated) computation of the
elliptic integrals K,E

SI units have to be used throughout.

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121015
%----------------------------------------------------------------------
%}
classdef loop_2d < magnetic_field.element_2d
%----------------------------------------------------------------------    
    properties
        I  % intensity
        ZL % axial position of loop
        RL % radius of loop         
        is_axi = 1 % whether axisimetric (1) or planar (0) geometry is used
    end
%----------------------------------------------------------------------    
    properties (Hidden = true)
        interpKE % whether to use interpolated ellipke for speed (1/0)
    end 
%----------------------------------------------------------------------    
    methods % Constructor
        function h = loop_2d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('I',1,@isnumeric);            
            p.addParameter('ZL',0,@isnumeric);
            p.addParameter('RL',3.5,@isnumeric);
            p.addParameter('interpKE',0,@isnumeric);
            % Parse and assign input
            p.parse(varargin{:});            
            h.I = p.Results.I;
            h.ZL = p.Results.ZL;
            h.RL = p.Results.RL;
            h.interpKE = p.Results.interpKE; 
        end
        function set.is_axi(h,v)
            if v == 0
                error('Loops cannot be set to planar mode, they are always axisymmetric');
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
        % calculates field of single current loop analytically. The code
        % has been optimized for speed
            % allocate
            psi = z.*0; 
            Bz = psi;
            Br = psi;            
            if h.I == 0    % Skip null-current loop
                return;
            end
            % Multiplying constant
            C = h.I*h.const.mu0/(4*pi);  
            % Intermediate values
            z = z - h.ZL; % subtract coil position
            sgn = sign(r);
            r = abs(r);
            R2 = r.^2;
            Z2 = z.^2;
            SQRT = (h.RL+r).^2+Z2;
            k2 = 4.*h.RL.*r./SQRT;
            SQRT = sqrt(SQRT);
            if h.interpKE % test for faster computation of KE
                [K,E] = h.ellipke(k2);
            else
                [K,E] = ellipke(k2);
            end
            E_FACTOR = E./((r-h.RL).^2+Z2);          
            % Calculate psi, Br, Bz
            psi =  sgn .* C .* SQRT .* ((2-k2).*K-2.*E);
            Bz  =  C.*2 ./ SQRT .* (K + (h.RL^2-R2-Z2).*E_FACTOR);
            Br  = - sgn .* C.*2 .* z./r ./ SQRT .* (K-(h.RL^2+R2+Z2).*E_FACTOR);           
            % L'hopital for Br near the axis, where expression becomes 0/0
            bad_points = (k2 < 1e-4 & r < 1e-4*h.RL); % 1e-4 is arbitrary, but works fine for the match of solutions
            Br(bad_points) = C.*z(bad_points)./SQRT(bad_points) ./((r(bad_points)-h.RL).^2+Z2(bad_points)) .* ...
                              3*pi.*r(bad_points).*h.RL^2./(h.RL^2 + Z2(bad_points));           
            % Return Inf, zeros if (z,r) == (ZL,RL) (the 0 is because z
            % is z-ZL)
            psi(z==h.ZL & r==h.RL) = Inf;
            Bz(z==h.ZL & r==h.RL)  = 0;
            Br(z==h.ZL & r==h.RL)  = 0;            
        end 
        function [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = derivatives_2d(h,z,r) 
        % Calculates dervatives of Bz, Br on z,r            
            % allocate
            dBz_dz = z.*0;
            dBz_dr = dBz_dz;
            dBr_dz = dBz_dz;
            dBr_dr = dBz_dz;                        
            if h.I == 0    % Skip null-current loop
                return;
            end            
            % Multiplying constant
            C = h.I*h.const.mu0/(2*pi);             
            % Intermediate variables
            z = z - h.ZL; % z only appears like this
            Rp = abs(r) + h.RL;
            Rm = abs(r) - h.RL;
            R2 = r.^2;
            Z2 = z.^2;
            MP = (z.^2 + Rp.^2); % minus for z, plus for r (SQRT^2)
            MM = (z.^2 + Rm.^2); % minus for z, minus for r
            TEMP = (z.^2 + r.^2 + h.RL.^2).*(r.^2 + h.RL.^2);
            k2 = 4.*h.RL.*abs(r)./MP;
            if h.interpKE % test for faster computation of KE
                [K,E] = h.ellipke(k2);
            else
                [K,E] = ellipke(k2);
            end
            E = E./MM;             
            % Calculate derivatives. Apply rules when possible
            dBz_dz = - C .* z ./ (MM.*MP.^(3/2)) .* ...
                ( (Z2 + R2 - h.RL.^2).*K - ...
                ( (Z2 + R2 - 3.*h.RL.^2).^2 -16.*h.RL.^4 + 12.*R2.*h.RL.^2).*E);
            dBz_dr = - C  ./ (r.*MM.*MP.^(3/2)) .* ...
                ( (TEMP - 4.*R2.*h.RL.^2).*K - ...
                (TEMP.*(Z2 + R2 + h.RL.^2) - 4.*R2.*h.RL.^2.*(4.*Z2 + R2 + h.RL.^2)).*E);
            dBr_dz = dBz_dr;
            dBr_dr = - dBz_dz + (C .* z ./ (r.*sqrt(MP)) .* (K- (R2+h.RL.^2+Z2).*E))./r; 
            % Remove that contribution if (z,r) == (ZL,RL), and remove NaNs (point at the axis)
            bad = (z==h.ZL & r==h.RL);
            dBz_dz(bad) = 0;
            dBz_dr(bad) = 0;
            dBr_dz(bad) = 0;
            dBr_dr(bad) = 0;
            bad = r==0;
            dBz_dr(bad) = 0;
            dBr_dz(bad) = 0;
            dBr_dr(bad) = -dBz_dz(bad)/2; 
        end 
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_2d(h,varargin) 
            % Plot                       
            h_line = line('XData',[h.ZL,h.ZL],'YData',[h.RL,-h.RL],'DisplayName',inputname(1));
            set(h_line,'Marker','s','linestyle','none',varargin{:});
        end
    end   
%----------------------------------------------------------------------
    methods (Static = true)
        function [k,e] = ellipke(m)
        % Fast interpolation of elliptic integrals K,E with persistent
        % table
            persistent m_ k_ e_ n_library_1
            if isempty(m_)
                n_library_1 = 9999;
                m_ = linspace(0,1,n_library_1+1);
                [k_,e_] = ellipke(m_);
                k_(end) = k_(end-1); % avoid the infinite within the library
                k_ = k_.';
                e_ = e_.';
            end
            k = m; e = m; % Allocate and give shape
            % Calculate m in index space. Make it column
            m = m(:)*(n_library_1)+1;            
            % Linear interpolation method
            fm = floor(m);          % floor
            cm = ceil(m);           % ceiling (needs both for when it is on the border)
            m = m-fm;              % increment (saved on xi to save space)
            % interpolate with libraries. NO INPUT CHECKING AT ALL!
            k(:) = k_(fm).*(1-m) + k_(cm).*m;
            e(:) = e_(fm).*(1-m) + e_(cm).*m;
        end
    end
%----------------------------------------------------------------------    
end
