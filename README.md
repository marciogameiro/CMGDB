# CMGDB
Conley Morse Graph Data Base

Install [Boost](https://www.boost.org/)

	brew install boost

Install the [Succinct Data Structure Library (SDSL)](https://github.com/simongog/sdsl-lite)

	git clone https://github.com/simongog/sdsl-lite.git
	cd sdsl-lite
	./install.sh /usr/local/

Install CMGDB

	git clone https://github.com/marciogameiro/CMGDB.git
	cd CMGDB
	./install.sh

This will compile the code with the map defined in the file `src/ModelMap.h` (currently the Leslie model) and will put the executable in the folder `bin`.

To run the code in the `bin` directory with the config parameters defined in `config.xml` type

	cd bin
	./CMGDB . 19.0 20.0

This will run the code in the current directory (the `.` option) with parameter values `19.0` and `20.0` and save the Morse graph to `morsegraph.gv` and the results to `data.mg`.

Here is an old [survey](http://chomp.rutgers.edu/Projects/survey/cmdbSurvey.pdf) and a
[talk](http://chomp.rutgers.edu/Projects/Databases_for_the_Global_Dynamics/software/LorentzCenterAugust2014.pdf) that might be useful.

