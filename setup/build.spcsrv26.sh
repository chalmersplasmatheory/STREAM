#!/usr/bin/bash
#

source environment.spcsrv26.sh

function install_stream {
	cd "$STREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_EXECUTABLE_RUNS=YES -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH &&
	make -j8
}

install_stream

