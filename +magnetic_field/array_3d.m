%{
This class serves to join several 3d generators of different types in
one array. 

Note that all generators must be referred to the same coordinate system.

Methods for calculating, plotting, etc the whole set of generators are
included. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121026
%----------------------------------------------------------------------
%} 
classdef array_3d < magnetic_field.array & magnetic_field.element_3d
%----------------------------------------------------------------------
    properties
        generators % cell array where the element generators are stored 
    end
%----------------------------------------------------------------------
     methods % Constructor, get/set
        function h = array_3d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('generators',{magnetic_field.loop_3d},@iscell);            
            % Parse and assign input
            p.parse(varargin{:});            
            h.generators = p.Results.generators;  
        end 
     end
%----------------------------------------------------------------------
    methods % Calculation methods        
        function set_B0(h,B0)
            if ~exist('B0','var')
                B0 = 1; % if nothing provided, normalize to 1
            end
            factor = B0/h.B_3d(0,0,0);
            for i=1:h.n_generators
                B0elem = h.generators{i}.B_3d(0,0,0);
                h.generators{i}.set_B0(factor*B0elem); % divide by n to make actual B0
            end
        end
        function [Bx,By,Bz] = field_3d(h,x,y,z)
        % simply calls each element and adds the results
            [Bx,By,Bz] = h.generators{1}.field_3d(x,y,z);
            for i = 2:h.n_generators
                [Bx_,By_,Bz_] = h.generators{i}.field_3d(x,y,z);
                Bx = Bx + Bx_;
                By = By + By_;
                Bz = Bz + Bz_;
            end
        end
        function [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = derivatives_3d(h,x,y,z) % Calculates dervatives of Bz, Br loop by loop on z,r
        % simply calls each element and adds the results
            [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = ...
                                            h.generators{1}.derivatives_3d(x,y,z);
            for i = 2:h.n_generators
                [dBx_dx_,dBx_dy_,dBx_dz_,dBy_dx_,dBy_dy_,dBy_dz_,dBz_dx_,dBz_dy_,dBz_dz_] = ...
                                                        h.generators{i}.derivatives_3d(x,y,z);
                dBx_dx = dBx_dx + dBx_dx_;
                dBx_dy = dBx_dy + dBx_dy_;
                dBx_dz = dBx_dz + dBx_dz_;
                
                dBy_dx = dBy_dx + dBy_dx_;
                dBy_dy = dBy_dy + dBy_dy_;
                dBy_dz = dBy_dz + dBy_dz_;
                
                dBz_dx = dBz_dx + dBz_dx_;
                dBz_dy = dBz_dy + dBz_dy_;
                dBz_dz = dBz_dz + dBz_dz_;
            end
        end
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_lines = plot_generator_3d(h,varargin)
        % simply calls each element and concatenates the results
        h_lines = []; % start with empty array
            for i = 1:h.n_generators
                temp = h.generators{i}.plot_generator_3d(varargin{:});
                h_lines = [h_lines,temp]; % concatenate
            end
        end
    end 
%----------------------------------------------------------------------    
end
