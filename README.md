# GDSGE: A Toolbox for Solving Global DSGE Models

## The GDSGE Toolbox

GDSGE is a toolbox that solves nonlinear Dynamic Stochastic General Equilibrium (DSGE) models with a global method based on the Simultaneous Transition and Policy Function Iteration (STPFI) algorithm introduced in [Cao, Luo, and Nie (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3569013). It allows users to define economic models in compact and intuitive scripts, called gmod files (gmod stands for global model). It parses the script into dynamic libraries which implement the actual computations (policy function iterations and Monte Carlo simulations) efficiently in C++, and provides a convenient MATLAB interface to researchers.

The toolbox can be used to solve models in macroeconomics, international finance, asset pricing, and related fields.

## Requirements for running remote compiled code

* MATLAB ver>=2017b. MATLAB toolbox: curve fitting
* Upload your gmod file following the instruction here: [GDSGE: A Toolbox for Solving DSGE Models with Global Methods — GDSGE Homepage](http://www.gdsge.com/)

## Requirements for local compiler

* MATLAB ver>=2017b. MATLAB toolboxes: symbolic math, curve fitting

* Visual C++ 2019 / Intel C++ Compiler 2017 or other MATLAB-compatible compilers. 

## Installation of local compiler

First, Configure your mex C++ compiler by running in MATLAB

  ```matlab
  mex -setup c++
  ```

Then, acquire the source code by cloning the git repository:

```
git clone https://github.com/gdsge/gdsge
```

Finally, in MATLAB, change directory to gdsge/tests, run

```matlab
runtests
```

which runs all the tests and produce all results in the companion paper [Cao, Luo, and Nie (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3569013).

## License

GDSGE is released under the Apache License, Version 2.0,  which is available at http://www.apache.org/licenses/LICENSE-2.0. In short, this license allows you to use, compose and distribute the GDSGE compiler or generated codes freely. However, it is requested that the companion paper be cited:

​	Cao, Dan and Luo, Wenlan and Nie, Guangyu, Global DSGE Models (April 1, 2020). Available at SSRN: https://ssrn.com/abstract=3569013 or http://dx.doi.org/10.2139/ssrn.3569013

GDSGE relies on the following external libraries, with their licenses described below and attached under folder licenses/:

* Adept: A combined automatic differentiation and array library for C++.

  Licensed under the Apache License, Version 2.0. Citation to the academic paper:

  * Hogan, R. J., 2014: Fast reverse-mode automatic differentiation using expression templates in C++. ACM Trans.
    Math. Softw., 40, 26:1-26:16.

* CoDoSol: a bound-constrained nonlinear equations solver.

  Citation to the academic paper:

  * Bellavia, S., M. Macconi, and S. Pieraccini (2012). Constrained dogleg methods for nonlinear systems with simple bounds. Computational Optimization and Applications 53(3), 771–794.

* myppual: Construct and Evaluate splines in ppform at flexible vector-valued dimensions, table look-up index, and spline dimension reduction in both vectorized pure MATLAB code and CMEX implementation.

  Copyright Jinhui Bai (jinhui.bai@gmail.com) and Wenlan Luo (luowenlan@gmail.com)

* v2struct: Pack/Unpack Variables to/from a scalar structure.

  Copyright Adi Navve, released under the MATLAB File Exchange License (BSD License)

* flat_hash_map: Copyright Malte Skarupke 2017.

  Distributed under the Boost Software License, Version 1.0 (http://www.boost.org/LICENSE_1_0.txt)

* Intel Math Kernel Library.

  Licensed under the Intel Simplified Software License (Version February 2020)
