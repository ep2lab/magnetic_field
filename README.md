MAGNETIC_FIELD
==============
  
Magnetic field computations based on the analytical solution for coils and
wire segments, in 2D and 3D.

The code is structured into several Matlab classes as follows:

- element,  element_2d, element_3d (abstract classes): all magnetic field
generator elements inherit from these parent classes. Used to set some
standard interfaces; not intended for direct usage by user. They also provide
some commodity functions to calculate the field on a plane and to plot a field
component on a plane.

* loop_2d, loop_3d: analytical magnetic field of a single current loop of zero
thickness, after solving the singularity for Br at the axis. The 2d version
defines the magnetic streamfunction psi. Inherit from element_2d and
element_3d, respectively. Loop_3d, in addition, inherits from loop_2d.

* library_2d, library_3d: allows to define some interpolation libraries for
the fields and interpolate there. Inherit from element_2d and element_3d,
respectively.

* path_3d: the magnetic field of a wire defined by a set of nodes. Inherits
from element_3d.

* uniform_2d, uniform_3d: a constant, uniform magnetic field in a given
direction. Inherit from element_2d and element_3d, respectively.

* array (abstract class): sets the basic interfaces for arrays of generators 
(i.e. collections of generator elements). Inherits from element.

* array_2d, array_3d: these and 3d are classes that collect many generator
elements in an array, and redefine some methods to compute field of all
elements simultaneously. Inherit from array.

Testing
-------

Unit tests are found in the /test subdirectory. After adding the package to
your Matlab path, you can run all tests by executing 'runtests' from this 
subdirectory.

License
-------

Copyright (c) 2017 Mario Merino. All rights reserved
