%{ 
This class interpolates the magnetic field from library arrays (3d
case). NaNs are returned if you try to extrapolate

IMPORTANT:
The grid, XX,YY,ZZ, needs to be perfectly rectangular and as produced by
meshgrid (i.e., XX(i,j,k) grows with j and YY(i,j,k) grows with i). 

You must provide the coordinates of the field generator if you want to
plot it, otherwise plot_generator_2d returns empty.

MMM20121015
%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121015
%----------------------------------------------------------------------
%}
classdef library_3d < magnetic_field.element_3d 
%----------------------------------------------------------------------
    properties
        XX = []; % grid
        YY = []; 
        ZZ = []; % FIELD
        BX = []; 
        BY = []; 
        BZ = []; 
        dBx_dx = []; % DERIVATIVES
        dBx_dy = []; 
        dBx_dz = []; 
        dBy_dy = [];
        dBy_dz = [];
        dBz_dz = [];
        Xgenerator % cell arrays of numeric arrays for plotting generator
        Ygenerator
        Zgenerator
    end 
%----------------------------------------------------------------------    
    properties (Dependent = true)
        dBy_dx; % dBy_dx = dBx_dy (ie., works as an alias)
        dBz_dx;
        dBz_dy;
        grid_size
    end  
%----------------------------------------------------------------------
    methods % Constructor, set/get
        function h = library_3d(varargin)
            % Generate defaults
            [XX,YY,ZZ] = meshgrid(linspace(0,10,100),linspace(0,8,120),linspace(0,5,60));
            % Validate input
            p = inputParser;         
            p.addParameter('XX',XX,@isnumeric);
            p.addParameter('YY',YY,@isnumeric);
            p.addParameter('ZZ',ZZ,@isnumeric);   
            p.addParameter('BX',XX*0,@isnumeric);
            p.addParameter('BY',XX*0,@isnumeric);
            p.addParameter('BZ',XX*0,@isnumeric);            
            p.addParameter('dBx_dx',XX*0,@isnumeric);
            p.addParameter('dBx_dy',XX*0,@isnumeric);
            p.addParameter('dBx_dz',XX*0,@isnumeric);         
            p.addParameter('dBy_dy',XX*0,@isnumeric);
            p.addParameter('dBy_dz',XX*0,@isnumeric);
            p.addParameter('dBz_dz',XX*0,@isnumeric);
            p.addParameter('Xgenerator',{},@iscell);
            p.addParameter('Ygenerator',{},@iscell);
            p.addParameter('Zgenerator',{},@iscell); 
            % Parse and assign input
            p.parse(varargin{:});            
            h.XX = p.Results.XX;
            h.YY = p.Results.YY;
            h.ZZ = p.Results.ZZ;            
            h.BX= p.Results.BX;
            h.BY= p.Results.BY;
            h.BZ = p.Results.BZ;            
            h.dBx_dx= p.Results.dBx_dx;
            h.dBx_dy= p.Results.dBx_dy;
            h.dBx_dz= p.Results.dBx_dz;
            h.dBy_dy= p.Results.dBy_dy;
            h.dBy_dz= p.Results.dBy_dz;            
            h.dBz_dz= p.Results.dBz_dz;            
            h.Xgenerator = p.Results.Xgenerator;
            h.Ygenerator = p.Results.Ygenerator;
            h.Zgenerator = p.Results.Zgenerator;
        end
        function v = get.grid_size(h)
            v = size(h.ZZ);
        end 
        function v = get.dBy_dx(h)
            v = h.dBx_dy;
        end
        function h = set.dBy_dx(h,v)
            h.dBx_dy = v;
        end
        function v = get.dBz_dx(h)
            v = h.dBx_dz;
        end
        function h = set.dBz_dx(h,v)
            h.dBx_dz = v;
        end
        function v = get.dBz_dy(h)
            v = h.dBy_dz;
        end
        function h = set.dBz_dy(h,v)
            h.dBy_dz = v;
        end
    end
%----------------------------------------------------------------------    
    methods % Calculation methods
        function h = set_B0(h,B0)
            if ~exist('B0','var')
                B0 = 1; % if nothing provided, normalize to 1
            end      
            h.BX = h.BX * B0/h.B_3d(0,0,0);
            h.BY = h.BY * B0/h.B_3d(0,0,0);
            h.BZ = h.BZ * B0/h.B_3d(0,0,0);
            h.dBx_dx = h.dBx_dx * B0/h.B_3d(0,0,0); % DERIVATIVES
            h.dBx_dy = h.dBx_dy * B0/h.B_3d(0,0,0); 
            h.dBx_dz = h.dBx_dz * B0/h.B_3d(0,0,0);
            h.dBy_dy = h.dBy_dy * B0/h.B_3d(0,0,0);
            h.dBy_dz = h.dBy_dz * B0/h.B_3d(0,0,0);
            h.dBz_dz = h.dBz_dz * B0/h.B_3d(0,0,0);
        end
        function [Bx,By,Bz] = field_3d(h,x,y,z) 
            % Get interpolation indices
            [lx,hx,sx,ix] = utilities.qinterpi(h.XX(1,:,1),x);
            [ly,hy,sy,iy] = utilities.qinterpi(h.YY(:,1,1),y);
            [lz,hz,sz,iz] = utilities.qinterpi(h.ZZ(1,1,:),z);
            % Allocate
            Bx = x.*NaN;
            By = Bx;
            Bz = Bx;
            % Compute
            for i=1:numel(x)       
                if ix(i) && iy(i) && iz(i)
                    Bx(i)  = h.BX(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                           + h.BX(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                           + h.BX(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                           + h.BX(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                           + h.BX(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                           + h.BX(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                           + h.BX(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                           + h.BX(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                    By(i)  = h.BY(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                           + h.BY(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                           + h.BY(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                           + h.BY(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                           + h.BY(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                           + h.BY(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                           + h.BY(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                           + h.BY(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                    Bz(i)  = h.BZ(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                           + h.BZ(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                           + h.BZ(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                           + h.BZ(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                           + h.BZ(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                           + h.BZ(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                           + h.BZ(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                           + h.BZ(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                end
            end
        end 
        function [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = derivatives_3d(h,x,y,z) % Calculates dervatives of Bz, Br loop by loop on z,r
            % !!! IMPLEMENTED BUT NOT CHECKED/TESTED PROPERLY 20170314.
            % Get interpolation indices    
            [lx,hx,sx,ix] = utilities.qinterpi(h.XX(1,:,1),x);
            [ly,hy,sy,iy] = utilities.qinterpi(h.YY(:,1,1),y);
            [lz,hz,sz,iz] = utilities.qinterpi(h.ZZ(1,1,:),z);
            % Allocate           
            dBx_dx = x.*NaN; 
            dBx_dy = dBx_dx;
            dBx_dz = dBx_dx;
            dBy_dy = dBx_dx;
            dBy_dz = dBx_dx;
            dBz_dz = dBx_dx;
            % Compute
            for i=1:numel(x)     
                if ix(i) && iy(i) && iz(i)
                    dBx_dx(i)  = h.dBx_dx(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBx_dx(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBx_dx(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBx_dx(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                               + h.dBx_dx(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBx_dx(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBx_dx(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBx_dx(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                    dBx_dy(i)  = h.dBx_dy(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBx_dy(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBx_dy(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBx_dy(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                               + h.dBx_dy(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBx_dy(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBx_dy(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBx_dy(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                    dBx_dz(i)  = h.dBx_dz(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBx_dz(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBx_dz(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBx_dz(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                               + h.dBx_dz(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBx_dz(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBx_dz(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBx_dz(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                    dBy_dy(i)  = h.dBy_dy(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBy_dy(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBy_dy(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBy_dy(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                               + h.dBy_dy(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBy_dy(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBy_dy(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBy_dy(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                    dBy_dz(i)  = h.dBy_dz(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBy_dz(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBy_dz(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBy_dz(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                               + h.dBy_dz(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBy_dz(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBy_dz(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBy_dz(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                    dBz_dz(i)  = h.dBz_dz(lx(i),ly(i),lz(i)).*(1-sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBz_dz(lx(i),ly(i),hz(i)).*(1-sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBz_dz(lx(i),hy(i),lz(i)).*(1-sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBz_dz(lx(i),hy(i),hz(i)).*(1-sx(i)).*(sy(i)).*(sz(i)) ...
                               + h.dBz_dz(hx(i),ly(i),lz(i)).*(sx(i)).*(1-sy(i)).*(1-sz(i)) ...
                               + h.dBz_dz(hx(i),ly(i),hz(i)).*(sx(i)).*(1-sy(i)).*(sz(i)) ...
                               + h.dBz_dz(hx(i),hy(i),lz(i)).*(sx(i)).*(sy(i)).*(1-sz(i)) ...
                               + h.dBz_dz(hx(i),hy(i),hz(i)).*(sx(i)).*(sy(i)).*(sz(i));
                end                       
            end
            dBy_dx = dBx_dy;
            dBz_dx = dBx_dz;
            dBz_dy = dBy_dz;
        end 
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_line = plot_generator_3d(h,varargin) 
            h_line = []; % Empty to start with
            for i = 1:length(h.Xgenerator)
                h_line(i) = line('XData',h.Xgenerator{i},'YData',h.Ygenerator{i},'ZData',h.Zgenerator{i},'DisplayName',inputname(1));
                set(h_line,varargin{:});    
            end
        end
    end     
%----------------------------------------------------------------------
end
