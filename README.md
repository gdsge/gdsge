# GDSGE: A Toolbox for Solving Global DSGE Models

## The GDSGE Toolbox

GDSGE is a toolbox that solves nonlinear Dynamic Stochastic General Equilibrium (DSGE) models with a global method based on the Simultaneous Transition and Policy Function Iteration (STPFI) algorithm introduced in [Cao, Luo and Nie (2023)]( https://www.sciencedirect.com/science/article/pii/S1094202523000017). It allows users to define economic models in compact and intuitive scripts, called gmod files (gmod stands for global model). It parses the scripts into dynamic libraries which implement the actual computations (policy function iterations and Monte Carlo simulations) efficiently in C++, and provides a convenient MATLAB interface to researchers.

The toolbox can be used to solve models in macroeconomics, international finance, asset pricing, and related fields.

See the [toolbox website for examples and documentation](http://www.gdsge.com/).

## Run on MATLAB Online

First, log into your [MATLAB Online](https://matlab.mathworks.com/)

Then, download the toolbox release and unzip in the MATLAB command window:

```matlab
websave('gdsge.zip','https://github.com/gdsge/gdsge/archive/refs/tags/v0.1.5.zip')
unzip gdsge.zip
```

Finally, change directory to the unzipped folder, set up the default mex compiler, and run tests in the MATLAB command window:

```matlab
mex -setup c++
cd gdsge-0.1.5/tests
runtests
```

This produces all results in the companion paper [Cao, Luo and Nie (2023)]( https://www.sciencedirect.com/science/article/pii/S1094202523000017)

## Requirements for the local compiler

* Windows
  - MATLAB 2017b+. 
  - MATLAB toolboxes: Symbolic Math
  - Compiler: MinGW64 C++ (installed via MATLAB Add-Ons); Visual Studio C++ Compiler 2019, 2022 (community version is fine); Intel C++ Compiler 2017 or newer
* macOS (silicon processor): 
  * MATLAB silicon version (R2023b+).
  * MATLAB toolboxes: Symbolic Math
  * Compiler: Xcode

* macOS (Intel processor):
  * MATLAB 2017b+. 
  * MATLAB toolboxes: Symbolic Math
  * Compiler: g++ 8.5; see [the instruction for how to setup the g++8.5 compiler](README_compiler_macOS.md)


## Installation of the local compiler

First, Configure your mex C++ compiler by running in MATLAB

  ```matlab
  mex -setup c++
  ```

Then, acquire the source code by cloning the git repository (*the local folder name should not contain spaces*):

```git
git clone https://github.com/gdsge/gdsge
```

Next, in MATLAB, change directory to gdsge/tests, run

```matlab
runtests
```

which runs all the tests and produce all results in the companion paper [Cao, Luo, and Nie (2023)](https://www.sciencedirect.com/science/article/pii/S1094202523000017).

To compile a gmod file, add folder "source" to MATLAB search path and run *gdsge_codegen* after changing the working directory to the one that contains the gmod file. For example, suppose you have located tests/HeatonLucas1996 with HL1996.gmod in the working directory, then simply run 

```matlab
gdsge_codegen('HL1996')
```

which will generate all the source codes and call the C++ compiler to compile the mex files.

## License

GDSGE is released under the Apache License, Version 2.0,  which is available at http://www.apache.org/licenses/LICENSE-2.0. In short, this license allows you to use, compose and distribute the GDSGE compiler or generated codes freely. However, it is requested that the companion paper be cited:

**Dan Cao, Wenlan Luo, and Guangyu Nie (2023). Global DSGE models. Review of Economic Dynamics, Volume 51, December 2023. Available at: https://www.sciencedirect.com/science/article/pii/S1094202523000017**

GDSGE relies on the following external libraries, with their licenses described below and attached under folder licenses/:

* Adept: A combined automatic differentiation and array library for C++.

  Licensed under the Apache License, Version 2.0. Citation to the academic paper:

  * Hogan, R. J., 2014: Fast reverse-mode automatic differentiation using expression templates in C++. ACM Trans.
    Math. Softw., 40, 26:1-26:16.

* CoDoSol: a bound-constrained nonlinear equations solver.

  Citation to the academic paper:

  * Bellavia, S., M. Macconi, and S. Pieraccini (2012). Constrained dogleg methods for nonlinear systems with simple bounds. Computational Optimization and Applications 53(3), 771–794.

* myppual: Construct and Evaluate splines in ppform at flexible vector-valued dimensions, table look-up index, and spline dimension reduction in both vectorized pure MATLAB code and CMEX implementation.

  Copyright (c) 2014 Jinhui Bai (jinhui.bai@gmail.com) and Wenlan Luo (luowenlan@gmail.com)

* v2struct: Pack/Unpack Variables to/from a scalar structure.

  Copyright (c) 2014, Adi Navve, released under the MATLAB File Exchange License (BSD License)

* flat_hash_map: Copyright Malte Skarupke 2017.

  Distributed under the Boost Software License, Version 1.0 (http://www.boost.org/LICENSE_1_0.txt)

* Intel Math Kernel Library.

  Licensed under the Intel Simplified Software License (Version February 2020)

* Eigen: a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
  Starting from the 3.1.1 version, it is licensed under the [MPL2](https://www.mozilla.org/en-US/MPL/2.0/)
