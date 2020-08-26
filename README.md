# CMGDB
Conley Morse Graph Database

Install [Boost](https://www.boost.org/) and [GMP](https://gmplib.org/)

	brew install boost
	brew install gmp

Install the [Succinct Data Structure Library (SDSL)](https://github.com/simongog/sdsl-lite)

	git clone https://github.com/simongog/sdsl-lite.git
	cd sdsl-lite
	./install.sh /usr/local/

It is also possible to intall the `sdsl-lite` library with Homebrew using the [repository](https://repology.org/project/sdsl-lite/versions) with the commands

	brew tap Brewsci/bio
	brew install sdsl-lite

Install CMGDB

	git clone https://github.com/marciogameiro/CMGDB.git
	cd CMGDB
	./install.sh

To get started on how to run the code see the examples in the Jupyter notebooks in the examples folder.

Here is an old [survey](http://chomp.rutgers.edu/Projects/survey/cmdbSurvey.pdf) and a
[talk](http://chomp.rutgers.edu/Projects/Databases_for_the_Global_Dynamics/software/LorentzCenterAugust2014.pdf) that might be useful.

