MAGNETIC_FIELD
==============

[![DOI](https://zenodo.org/badge/85017158.svg)](https://zenodo.org/badge/latestdoi/85017158)

Magnetic field classes based on the analytical solution for coils and
wire segments, in 2D and 3D.

## Installation

Installation requires simply that you 
[download MAGNETIC_FIELD](https://github.com/mariomerinomartinez/magnetic_field/archive/master.zip) 
and add the base directory (the one that contains the `+magnetic_field` directory) to 
your Matlab path.

### Dependencies

A recent version of Matlab is needed to run the code. 
The code has been developed in Matlab 2016a Academic version. 

MAGNETIC_FIELD 
depends on other Matlab packages that you can download from my GitHub
account:
[utilities](https://github.com/mariomerinomartinez/utilities)
and
[constants_and_units](https://github.com/mariomerinomartinez/constants_and_units).
These packages must be installed and added to your Matlab path beforehand.

## Usage

The code is structured into several Matlab classes as follows:

* `element`,  `element_2d`, `element_3d` (abstract classes): all magnetic field
generator elements inherit from these parent classes. Used to set some
standard interfaces; not intended for direct usage by user. They also provide
some commodity functions to calculate the field on a plane and to plot a field
component on a plane.
* `loop_2d`, `loop_3d`: analytical magnetic field of a single current loop of
zero thickness, after solving the singularity for Br at the axis. The 2d
version defines the magnetic streamfunction psi. The most thoroughly tested 
classes.
* `library_2d`, `library_3d`: allow to define some interpolation libraries for
the fields and interpolate there. 
* `path_3d`: the magnetic field of a wire defined by a set of nodes, computed
in a step-wise manner with Biot-Savart's Law.  Not extensively tested yet;
could contain errors.
* `wire_2d`, `wire_3d`: the analytic magnetic field of a straight infinite 
conductor. Not extensively tested yet; could contain errors.
* `uniform_2d`, `uniform_3d`: a constant, uniform magnetic field in a given
direction.
* `array` (abstract class): sets the basic interfaces for arrays of generators 
(i.e. collections of generator elements).
* `array_2d`, `array_3d`: these are classes that allow collecting many
generator elements in an array, and redefine some methods to compute field of
all elements simultaneously. 

Usage is straight forward. Start by creating an object of the desired class and setting its parameters, then call one of its methods:

```Matlab
mf = magnetic_field.loop_2d('RL',5,'I',10); % Create loop_2d object
[psi,Bz,Br] = mf.field_2d(0,0); % query the field at z=0, r=0
```

### Testing

Unit tests are found in the `/test` subdirectory. After adding the package to
your Matlab path, you can run all tests by executing `runtests` from this 
subdirectory.

## Contributing

If you have any comments for improvement or 
are interested in contributing to the continued 
development of this or any of my other codes, you can contact us
through our [website](http://ep2.uc3m.es/). 

For updates and news, follow us on Twitter: [@ep2lab.](https://twitter.com/ep2lab).

## Acknowledging 

This program is released as open source in the hope that it will be useful to
other people. 

If you find it useful and/or use it in any of your works, we kindly ask you
to acknowledge it by citing the code directly as: 

> Mario Merino (2017). mariomerinomartinez/magnetic_field: First release DOI:10.5281/zenodo.496131
 
## License

Copyright (c) 2017 Mario Merino. The software is released as open 
source with the [MIT License](LICENSE.md).

 
