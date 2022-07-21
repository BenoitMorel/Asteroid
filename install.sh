rm -rf build
mkdir build
cd build
cmake ..
make -j 4


echo "Now running unit tests..."
cd ..

run_tests.sh



