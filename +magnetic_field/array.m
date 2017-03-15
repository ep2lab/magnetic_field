%{
This abstract class serves to join several generators of different types
in one array. Methods for calculating, plotting, etc the whole set of generators 
are established in children classes

The actual classes that the user can instantiate are array_2d and
array_3d.

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121026
%----------------------------------------------------------------------
%}
classdef array < magnetic_field.element
%----------------------------------------------------------------------
    properties (Abstract = true)
        generators % cell array where the element generators are stored
    end 
%----------------------------------------------------------------------
    properties (Dependent = true)
        n_generators
    end 
%----------------------------------------------------------------------
    methods % Constructor, set/get
        function v = get.n_generators(h)
            v = numel(h.generators);
        end 
    end 
%----------------------------------------------------------------------    
end
