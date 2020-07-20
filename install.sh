
rm -rf build
rm -rf bin

mkdir build
cd build
cmake ..
make

mkdir ../bin
mv CMGDB ../bin
cp ../src/config.xml ../bin
