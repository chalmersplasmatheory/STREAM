#!/usr/bin/bash

if [ -z "$STREAMPATH" ]; then
	export STREAMPATH=$(realpath "$(dirname $0)"/..)
fi

module purge

module load python-3.10-base-toolkit-gcc-11.2.1
module load hdf5-1.12.1-gcc-11.2.1-java-osdk-17-serial
module load openmpi-4.1.2-gcc-11.2.1-java-osdk-17-cuda-11.6-hwloc-2.7.0-pmix
module load intel-oneapi-2022.1.2.146/mkl/2022.0.2
module load cmake-3.22.2-gcc-11.2.1 gsl-2.7.1
module load petsc-3.16.5-real

