%{
Parent class of all magnetic field generators

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20121011
%----------------------------------------------------------------------
%}
classdef element < handle
%----------------------------------------------------------------------
    properties
        const = constants_and_units.constants;
    end 
%----------------------------------------------------------------------    
    methods (Abstract = true)
        set_B0(h,B0,mode); % sets currents so that B = B0 at the origin
    end
%----------------------------------------------------------------------
end

