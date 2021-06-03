# CMGDB
Conley Morse Graph Database

## Overview

This project uses combinatorial and topological methods to compute dynamics of discrete dynamical systems.

## Dependencies and Installation

The dependencies and install instructions on a Mac using [Homebrew](https://brew.sh/) are listed below.

### Boost and GMP

Install [Boost](https://www.boost.org/) and [GMP](https://gmplib.org/)

	brew install boost
	brew install gmp

### SDSL

Install the [Succinct Data Structure Library (SDSL)](https://github.com/simongog/sdsl-lite)

	git clone https://github.com/simongog/sdsl-lite.git
	cd sdsl-lite
	./install.sh /usr/local/

It is also possible to intall the `sdsl-lite` library with Homebrew using this [repository](https://repology.org/project/sdsl-lite/versions) with the commands

	brew tap Brewsci/bio
	brew install sdsl-lite

### Eigen3

Install Eigen3

	brew install eigen

### graphviz

Install graphviz

	brew install graphviz

### CMake

Install CMake

	brew install cmake

### jupyter and graphviz

Install the jupyter and graphviz Python packages

	python -m pip install jupyter graphviz

### Install CMGDB

	git clone https://github.com/marciogameiro/CMGDB.git
	cd CMGDB
	./install.sh

# Documentation and Examples

To get started on how to run the code see the examples in the Jupyter notebooks in the examples folder.

In particular the notebooks [Examples.ipynb](examples/Examples.ipynb), [Gaussian\_Process\_Example.ipynb](examples/Gaussian_Process_Example.ipynb), and [Conley\_Index\_Examples.ipynb](examples/Conley_Index_Examples.ipynb) present basic examples on how to run the code and are a good starting point.

Here is an old [survey](http://chomp.rutgers.edu/Projects/survey/cmdbSurvey.pdf) and a
[talk](http://chomp.rutgers.edu/Projects/Databases_for_the_Global_Dynamics/software/LorentzCenterAugust2014.pdf) that might be useful.
