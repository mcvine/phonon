#!/usr/bin/env python

import numpy as np, os, glob, shutil

# create work dir
workdir = "work.IDF"
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.makedirs(workdir)

# copy files
datasrc = "../tests/data/graphite"
files = ['POSCAR', 'SPOSCAR', 'FORCE_CONSTANTS']
for f in files:
    shutil.copyfile(os.path.join(datasrc, f), os.path.join(workdir, f))

# move to work dir
os.chdir(workdir)

#
from mcvine.phonon.from_phonopy import idf, make_all
make_all(
    qgrid_dims=[31,31,31],
    fix_pols_phase=True,
    force_constants='FORCE_CONSTANTS', poscar='POSCAR', sposcar='SPOSCAR'
)

# SQE
import mcvine.phonon.powderSQE.IDF as psidf
disp = psidf.disp_from_datadir('.')
from mccomponents.sample import phonon as mcphonon
doshist = mcphonon.read_dos.dos_fromidf('DOS').doshist
IQE = psidf.from_data_dir(
    datadir='.',
    disp=disp, 
    N = int(1e6),
    Q_bins=np.arange(0, 4, 0.04), E_bins=np.arange(0,30,.2),
    doshist=doshist,
    T=300., Ei=30., max_det_angle=60.,
    include_multiphonon=False,
)
