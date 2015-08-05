# Lir #

Lir is a set of tools developed in the area of computational electromagnetics, by the RF Modelling and Simulation Group in the RINCE Institute, Dublin City University and licensed under the Apache License.

The initial version of Lir has been developed by Ian Kavanagh as part of his PhD under the guidance of Dr. Conor Brennan. It consists of a full wave 2D indoor propagation modelling tool. It is planned to expand Lir to include tools for producing time domain information as well as modelling urban and rural environments.

## Table of Contents ##

- [Compiling and Executing](#compiling-and-executing)
    - [Prerequisites](#prerequisites)
    - [2D](#2d)
        - [Sample](#2d)
    - [Plotting](#plotting)
- [Contributing](#contributing)
- [Documentation](#documentation)
- [Authors](#authors)
- [Copyright and Licensing Information](#copyright-and-licensing-information)

## Compiling and Executing ##

### Prerequisites ###

Lir is built on top of the [Intel Math Kernel Library (MKL)][mkl] and thus requires an installation of it before it can be compiled or run. Free versions of [Intel MKL][mkl] are available for [students and academic researchers][intel-academic].

The parallel version of Lir also requires the use of the [Intel C/C++ Compiler][icc], `icc` as there is a bug present when compiling with `gcc`. A free version of the [Intel C/C++ Compiler][icc] is available for [students and academic researchers][intel-academic].

### 2D ###

A `makefile` is provided which will compile the 2D version of the tool. You can use

```
make 2D
```

or for the initial version

```
make
```

to compile the 2D version. This will produce an executable named `2D` in the top-level directory which can be executed with

```
./2D <frequency> <antenna x location> <antenna y location> [input building file]
```

3 variables must be passed to the program on the command line:
- `<frequency>` is the frequency of the antenna used for computing the fields within the building. It should be a floating point value.
- `<antenna x location> <antenna y location>` are the **x** and **y** locations of the antenna within the building.
- `[input building file]` is an optional argument pointing to the location of a file which can be used as the input building. If one is not present the file at `input/building.txt` is used.

#### Sample ####

Executing the 2D tool with

```
./2D 700e6 -0.8 -0.8 input/building.txt
```

is known to work and produce valid results. Sample output files for the different elements used in modelling the propagation of the building in `input/building.txt` with a Hertzian dipole antenna radiating at 700MHz positioned at (-0.8,-0.8) can be found in `samples/`.

### Plotting ###

[MATLAB](https://uk.mathworks.com/products/matlab/)/[GNU Octave](https://www.gnu.org/software/octave/) files which will plot the output of the electric field overlayed on top of the building can be found in `plot/`. Executing `plot_fields` will plot the electric field found in `output/E.txt` on a surface plot for the positions in `output/position.txt` and overlay the building defined in `output/shape.txt` on top of it.

## Contributing ##

All bug reports, feature requests and feedback should be provided on the [Issues](https://github.com/IKavanagh/Lir/issues) tab of [Github](https://github.com/IKavanagh/Lir/issues).

There are a large number of items which need to be fixed, rewritten or completed which are documented inline (and in my head) but will be moved to the [Issues](https://github.com/IKavanagh/Lir/issues) tab shortly.

Contact details for authors, collaborators and contributors can be found in [Authors](#authors)

## Documentation ##

Most documentation currently resides within the header files in `include/`. It is intended to move all this to the `doc/` folder which houses the base for a [Sphinx](http://sphinx-doc.org/) documentation set which will either be hosted on [Read The Docs](https://readthedocs.org/) or [my personal website](https://iankavanagh.me).

## Authors ##

### Owner ###
Owner: Ian Kavanagh  
Website: [IanKavanagh.me](https://iankavanagh.me)

## Copyright and Licensing Information ##

Copyright (c) 2015 Ian Kavanagh  
This file is part of Lir.

Licensed under the Apache License, Version 2.0 (the "License");  
you may not use this file except in compliance with the License.  
You may obtain a copy of the License at  

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software  
distributed under the License is distributed on an "AS IS" BASIS,  
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  
See the License for the specific language governing permissions and  
limitations under the License.  

[mkl]: https://software.intel.com/en-us/intel-mkl "Intel MKL"
[icc]: https://software.intel.com/en-us/c-compilers "Intel C/C++ Compiler"
[intel-academic]: https://software.intel.com/en-us/qualify-for-free-software
