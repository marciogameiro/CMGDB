# CMGDB
Conley Morse Graph Database

## Overview

This project uses combinatorial and topological methods to compute dynamics of discrete dynamical systems.

## Installation

Install the latest tagged version:

	pip install CMGDB

To uninstall:

	pip uninstall CMGDB

## Documentation and examples

To get started on how to run the code see the examples in the Jupyter notebooks in the [examples](examples) folder.

In particular the notebooks [Examples.ipynb](examples/Examples.ipynb), [Gaussian\_Process\_Example.ipynb](examples/Gaussian_Process_Example.ipynb), and [Conley\_Index\_Examples.ipynb](examples/Conley_Index_Examples.ipynb) present basic examples on how to run the code and are a good starting point.

Here is an old [survey](http://chomp.rutgers.edu/Projects/survey/cmdbSurvey.pdf) and a
[talk](http://chomp.rutgers.edu/Projects/Databases_for_the_Global_Dynamics/software/LorentzCenterAugust2014.pdf) that might be useful.

## Installing from source and dependencies

To install from source you need a C++ compiler and the following dependencies installed: [Boost](https://www.boost.org/), [GMP](https://gmplib.org/), and the [Succinct Data Structure Library (SDSL)](https://github.com/simongog/sdsl-lite). Assuming you have these dependencies installed in your system, you can install from source with the command:

	pip install --force-reinstall --no-deps --no-cache-dir git+https://github.com/marciogameiro/CMGDB.git

Alternatively, you can clone the GitHub repository and install with:

	git clone https://github.com/marciogameiro/CMGDB.git
	cd CMGDB
	./install.sh
