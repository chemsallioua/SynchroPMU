
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

flags="-DDEBUG=OFF -DNUM_CHANLS=1"
while getopts "DN:" opt; do
  case $opt in
    D)
      flags="-DDEBUG=ON ${flags}"
      ;;
    N)
      flags="${flags} -DNUM_CHANLS=$OPTARG"
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
if [ "${machine}" == "Linux" ]; then
  cmake .. -G "Unix Makefiles" $flags
elif [ "${machine}" == "MinGw" ]; then
  cmake .. -G "MinGW Makefiles" $flags
elif [ "${machine}" == "Cygwin" ]; then
  cmake .. -G "Unix Makefiles" $flags
fi

make PmuEstimatorStatic
make PmuEstimatorShared

# installing the library
cmake --install .

# create the build folder and copy the necessary files
rm -rf ../build
mkdir ../build
mkdir ../build/config
mkdir ../build/headers

if [ "${machine}" == "Linux" ]; then
  cp libpmu_estimator.a ../build
  cp libpmu_estimator.so ../build
elif [ "${machine}" == "MinGw" ]; then
  cp libpmu_estimator.a ../build
  cp libpmu_estimator.dll ../build
elif [ "${machine}" == "Cygwin" ]; then
  cp libpmu_estimator.a ../build
  cp libpmu_estimator.dll ../build
fi

cp config.ini ../build/config
cp m_class_config.ini ../build/config
cp p_class_config.ini ../build/config
cp pmu_estimator.h ../build/headers
cp func_stubs.h ../build/headers

# remove the cmake_build folder
cd ..
rm -rf cmake_build