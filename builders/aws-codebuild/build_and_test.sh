#!/bin/bash

set -x
export PATH=$HOME/mc/bin:$PATH
source activate test

export SRC=$PWD
export EXPORT_ROOT=$CONDA_PREFIX
export MCVINE_PHONON_BLD_ROOT=$SRC/build
mkdir -p $MCVINE_PHONON_BLD_ROOT && cd $MCVINE_PHONON_BLD_ROOT 
cmake \
    -DCMAKE_INSTALL_PREFIX=$EXPORT_ROOT \
    -DPYTHON_LIBRARY=$CONDA_PREFIX/lib/libpython${PYTHON_VERSION}.so \
    -DPYTHON_INCLUDE_DIR=$CONDA_PREFIX/include/python${PYTHON_VERSION} \
    $SRC
make install
cd -
cd tests && py.test -s
