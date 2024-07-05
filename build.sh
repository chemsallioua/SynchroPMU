#!/bin/bash

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="OTHER"
esac
echo Machine: ${machine}, Name: ${unameOut}

if [ "$machine" == "OTHER" ]; then
    echo "Unsupported platform"
    exit 1
fi

flags="-DNUM_CHANLS=1"
num_chanls=1
while getopts "D:N:" opt; do
  case $opt in
    D)
      flags="${flags} -DLOGGING_LEVEL=$OPTARG"
      ;;
    N)
      flags="${flags} -DNUM_CHANLS=$OPTARG"
      num_chanls=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# create the cmake_build folder
rm -rf cmake_build
mkdir cmake_build
cd cmake_build

# run cmake and build the shared library
if [ "${machine}" == "MinGw" ]; then
  cmake .. -G "MinGW Makefiles" $flags
elif [ "${machine}" == "Cygwin" ]; then
  cmake .. -G "Unix Makefiles" $flags
elif [ "${machine}" == "Linux" ]; then
  cmake .. -G "Unix Makefiles" $flags
fi

# Run the Python script to modify pmu_estimator.h
python3 ../scripts/process_header.py ../src/pmu_estimator.h $num_chanls

make PmuEstimatorStatic
make PmuEstimatorShared

# installing the library
cmake --install .

# create the build folder and copy the necessary files
rm -rf ../build
mkdir ../build
mkdir ../build/config
mkdir ../build/headers

if [ "${machine}" == "MinGw" ]; then
  cp libpmu_estimator.a ../build
  cp libpmu_estimator.dll ../build
elif [ "${machine}" == "Cygwin" ]; then
  cp libpmu_estimator.a ../build
  cp libpmu_estimator.dll ../build
elif [ "${machine}" == "Linux" ]; then
  cp libpmu_estimator.a ../build
  cp libpmu_estimator.so ../build
fi

cp ../config/config.ini ../build/config
cp ../config/m_class_config.ini ../build/config
cp ../config/p_class_config.ini ../build/config
cp ../src/pmu_estimator.h ../build/headers
cp ../src/func_stubs.h ../build/headers

# Run the Python script to modify pmu_estimator.h back to it's original value
python3 ../scripts/process_after_build.py ../src/pmu_estimator.h

# remove the cmake_build folder
cd ..
rm -rf cmake_build


