%{
This class serves to join several 2d generators of different types in
one array. 

Note that all generators must be referred to the same coordinate system.

Methods for calculating, plotting, etc the whole set of generators are
included. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121026
%----------------------------------------------------------------------
%}  
classdef array_2d < magnetic_field.array & magnetic_field.element_2d 
%----------------------------------------------------------------------
    properties
        generators % cell array where the element generators are stored
        is_axi % whether axisimetric (1) or planar (0) geometry is used
    end
%----------------------------------------------------------------------    
     methods % Constructor, get/set
        function h = array_2d(varargin)
            % Validate input
            p = inputParser;               
            p.addParameter('generators',{magnetic_field.loop_2d},@iscell);            
            % Parse and assign input
            p.parse(varargin{:});            
            h.generators = p.Results.generators;  
        end
        function set.generators(h,v)
            h.generators = v;    
            h.is_axi = []; % trigger an is_axi computation and check
        end
        function set.is_axi(h,~)
            % Check axisymmetric or planar, but not mix
            for i = 1:h.n_generators
                is_axi(i) = h.generators{i}.is_axi;
            end
            if all(is_axi == 1) 
                h.is_axi = 1; % everything is axisymmetric
            elseif all(is_axi == 0)
                h.is_axi = 0; % everything is planar
            else
                h.is_axi = NaN; % mixing is not possible
                error('You cannot mix axisymmetric and planar field generators in object %s',inputname(1));
            end          
        end
    end        
%----------------------------------------------------------------------
    methods % Calculation methods
        function set_B0(h,B0)
            if ~exist('B0','var')
                B0 = 1; % if nothing provided, normalize to 1
            end
            factor = B0/h.B_2d(0,0);
            for i=1:h.n_generators
                B0elem = h.generators{i}.B_2d(0,0);
                h.generators{i}.set_B0(factor*B0elem); % divide by n to make actual B0
            end
        end
        function [psi,Bz,Br] = field_2d(h,z,r)
        % simply calls each element and adds the results
            [psi,Bz,Br] = h.generators{1}.field_2d(z,r);
            for i = 2:h.n_generators
                [psi_,Bz_,Br_] = h.generators{i}.field_2d(z,r);
                Bz = Bz + Bz_;
                Br = Br + Br_;
                psi = psi + psi_;
            end
        end
        function [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = derivatives_2d(h,z,r) 
        % simply calls each element and adds the results
            [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = h.generators{1}.derivatives_2d(z,r);    
            for i = 2:h.n_generators
                [dBz_dz_,dBz_dr_,dBr_dz_,dBr_dr_] = h.generators{i}.derivatives_2d(z,r);  
                dBz_dz = dBz_dz + dBz_dz_;
                dBz_dr = dBz_dr + dBz_dr_;
                dBr_dz = dBr_dz + dBr_dz_;
                dBr_dr = dBr_dr + dBr_dr_;
            end
        end
    end
%----------------------------------------------------------------------
    methods % Plotting methods
        function h_lines = plot_generator_2d(h,varargin)
        % simply calls each element and concatenates the results
        h_lines = []; % start with empty array
            for i = 1:h.n_generators
                temp = h.generators{i}.plot_generator_2d(varargin{:});
                h_lines = [h_lines,temp]; % concatenate
            end
        end
    end
%---------------------------------------------------------------------- 
end
