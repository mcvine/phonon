#!/bin/bash

set -x
export PATH=$HOME/mc/bin:$PATH
conda create -n test python=$PYTHON_VERSION
source activate test
conda install pytest awscli
conda install -c mcvine/label/unstable mcvine
conda list mcvine
conda list mcvine-core
conda install tqdm multiphonon periodictable
