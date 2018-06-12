%{ 
This class interpolates the magnetic field from library arrays. NaNs are
returned if you try to extrapolate

IMPORTANT:
The grid, ZZ,RR, needs to be perfectly rectangular and as produced by
meshgrid (i.e., ZZ(i,j) grows with j). 

You must provide the coordinates of the field generator if you want to
plot it, otherwise plot_generator_2d returns empty.

MMM20121015
%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121015
%----------------------------------------------------------------------
%}
classdef library_2d < magnetic_field.element_2d
%----------------------------------------------------------------------
    properties % interpolation arrays
        ZZ % grid
        RR 
        PSI % field
        BZ
        BR
        dBz_dz % derivatives
        dBz_dr
        dBr_dr
        is_axi
        Zgenerator % cell arrays of numeric arrays for plotting generator
        Rgenerator
    end
%----------------------------------------------------------------------    
    properties (Dependent = true)
        dBr_dz % dBr_dz = dBz_dr (ie., works as an alias)
        grid_size
    end
%----------------------------------------------------------------------    
    methods % Constructor, set/get
        function h = library_2d(varargin) 
            % Generate defaults
            [ZZ,RR] = meshgrid(linspace(0,10,100),linspace(0,8,120));
            loop = magnetic_field.loop_2d;
            [PSI,BZ,BR] = loop.field_2d(ZZ,RR);
            [dBz_dz,dBz_dr,dBr_dr] = loop.derivatives_2d(ZZ,RR);
            % Validate input
            p = inputParser;         
            p.addParameter('ZZ',ZZ,@isnumeric);            
            p.addParameter('RR',RR,@isnumeric);
            p.addParameter('PSI',PSI,@isnumeric);
            p.addParameter('BZ',BZ,@isnumeric);
            p.addParameter('BR',BR,@isnumeric);
            p.addParameter('dBz_dz',dBz_dz,@isnumeric);
            p.addParameter('dBz_dr',dBz_dr,@isnumeric);
            p.addParameter('dBr_dr',dBr_dr,@isnumeric);
            p.addParameter('Zgenerator',{},@iscell);
            p.addParameter('Rgenerator',{},@iscell);
            p.addParameter('is_axi',1,@isnumeric);
            % Parse and assign input
            p.parse(varargin{:});            
            h.ZZ = p.Results.ZZ;
            h.RR = p.Results.RR;
            h.PSI = p.Results.PSI;
            h.BZ = p.Results.BZ;
            h.BR = p.Results.BR;
            h.dBz_dz= p.Results.dBz_dz;
            h.dBz_dr = p.Results.dBz_dr;
            h.dBr_dr = p.Results.dBr_dr;
            h.Zgenerator = p.Results.Zgenerator;
            h.Rgenerator= p.Results.Rgenerator;
            h.is_axi = p.Results.is_axi;
        end
        function v = get.grid_size(h)
            v = size(h.ZZ);
        end
        function v = get.dBr_dz(h)
            v = h.dBz_dr;
        end
        function set.dBr_dz(h,v)
            h.dBz_dr = v;
        end
    end
%----------------------------------------------------------------------    
    methods % Calculation methods
        function h = set_B0(h,B0)
            if ~exist('B0','var')
                B0 = 1; % if nothing provided, normalize to 1
            end  
            oldB0 = h.B_2d(0,0);
            h.PSI = h.PSI * B0/oldB0;
            h.BZ = h.BZ * B0/oldB0;
            h.BR = h.BR * B0/oldB0;
            h.dBz_dz = h.dBz_dz * B0/oldB0;
            h.dBz_dr = h.dBz_dr * B0/oldB0;
            h.dBr_dr = h.dBr_dr * B0/oldB0;
        end
        function [psi,Bz,Br] = field_2d(h,z,r)
            % Get interpolation indices
            [lz,hz,sz,iz] = utilities.qinterpi(h.ZZ(1,:),z);            
            [lr,hr,sr,ir] = utilities.qinterpi(h.RR(:,1),r);            
            % Allocate
            psi = z.*NaN; 
            Bz = psi;
            Br = psi;     
            % Compute
            for i=1:numel(z)
                if iz(i) && ir(i)
                    psi(i) = h.PSI(lr(i),lz(i)).*(1-sr(i)).*(1-sz(i)) + h.PSI(hr(i),lz(i)).*sr(i).*(1-sz(i)) + ...
                          h.PSI(lr(i),hz(i)).*(1-sr(i)).*sz(i) + h.PSI(hr(i),hz(i)).*sr(i).*sz(i);
                    Bz(i)  = h.BZ(lr(i),lz(i)).*(1-sr(i)).*(1-sz(i)) + h.BZ(hr(i),lz(i)).*sr(i).*(1-sz(i)) + ...
                          h.BZ(lr(i),hz(i)).*(1-sr(i)).*sz(i) + h.BZ(hr(i),hz(i)).*sr(i).*sz(i);
                    Br(i)  = h.BR(lr(i),lz(i)).*(1-sr(i)).*(1-sz(i)) + h.BR(hr(i),lz(i)).*sr(i).*(1-sz(i)) + ...
                          h.BR(lr(i),hz(i)).*(1-sr(i)).*sz(i) + h.BR(hr(i),hz(i)).*sr(i).*sz(i);
                end
            end
        end 
        function [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = derivatives_2d(h,z,r)
            % Get interpolation indices            
            [lz,hz,sz,iz] = utilities.qinterpi(h.ZZ(1,:),z);            
            [lr,hr,sr,ir] = utilities.qinterpi(h.RR(:,1),r);  
            % Allocate
            dBz_dz = z.*NaN; 
            dBz_dr = dBz_dz;
            dBr_dr = dBz_dz;  
            % Compute
            for i=1:numel(z)                     
                if iz(i) && ir(i)
                    dBz_dz(i) = h.dBz_dz(lr(i),lz(i)).*(1-sr(i)).*(1-sz(i)) + h.dBz_dz(hr(i),lz(i)).*sr(i).*(1-sz(i)) + ...
                          h.dBz_dz(lr(i),hz(i)).*(1-sr(i)).*sz(i) + h.dBz_dz(hr(i),hz(i)).*sr(i).*sz(i);
                    dBz_dr(i) = h.dBz_dr(lr(i),lz(i)).*(1-sr(i)).*(1-sz(i)) + h.dBz_dr(hr(i),lz(i)).*sr(i).*(1-sz(i)) + ...
                          h.dBz_dr(lr(i),hz(i)).*(1-sr(i)).*sz(i) + h.dBz_dr(hr(i),hz(i)).*sr(i).*sz(i);

                    dBr_dr(i) = h.dBr_dr(lr(i),lz(i)).*(1-sr(i)).*(1-sz(i)) + h.dBr_dr(hr(i),lz(i)).*sr(i).*(1-sz(i)) + ...
                          h.dBr_dr(lr(i),hz(i)).*(1-sr(i)).*sz(i) + h.dBr_dr(hr(i),hz(i)).*sr(i).*sz(i);
                end
            end
            dBr_dz = dBz_dr;
        end 
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_2d(h,varargin) 
            h_line = []; % Empty to start with
            for i = 1:length(h.Zgenerator)
                h_line(i) = line('XData',h.Zgenerator{i},'YData',h.Rgenerator{i},'DisplayName',inputname(1));
                set(h_line,'Marker','s','linestyle','none',varargin{:});    
            end
        end
    end    
%----------------------------------------------------------------------
end
