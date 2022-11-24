#!/usr/bin/bash

USERHOME=$(eval echo ~$USER)

if [ -z "$PETSC_DIR" ]; then
	export PETSC_DIR="$USERHOME/petsc"
fi
if [ -z "$PETSC_ARCH" ]; then
	export PETSC_ARCH=arch-linux-c-opt
fi
if [ -z "$STREAM_DIR" ]; then
	export STREAM_DIR="$USERHOME/STREAM"
fi

alias streamviz="python -i $STREAM_DIR/py/cli/cli.py"

