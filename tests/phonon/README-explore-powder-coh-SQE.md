There are several notebooks in this directory with filename `explore*.ipynb`. They were used to explore calculation of 
coherent powder SQE from phonon data.

The best one is `explore-random-sampling-big-Qspace-box-no-importance-sampleing-along-Qmag.ipynb`.

The algorithm there is
* Randomly choose Q points in a box in reciprocal space
* Calculate phonon data using phonopy (first convert Q points to reciprocal units)
* Calculate scattering powder for all phonon modes at these Q points, and create a histogram of I(Q,E)
* Apply corrections to get S(Q,E)

When comparing the S(Q,E) so calculated with the single phonon S(Q,E) from multiphonon code, we found a scaling factor near 8$\pi$.
It still is unclear to me where it comes from. I can think of 4$\pi$ from solid angle. The addional 2 is mysterious.

After adjust for that factor, the comparison to experiment looks quite good.

![graphtie-iqe-compare](https://user-images.githubusercontent.com/1796155/43750400-403f4b36-99c6-11e8-875f-4cbfb68c7aed.png)


A few things to note:
* The polarization vectors from phonopy should be normalized
* Make sure to sample the reciprocal space correctly (either evenly as in `explore-random-sampling-big-Qspace-box-no-importance-sampleing-along-Qmag.ipynb`, or with important sampling as in `explore-random-sampling-big-Qspace-box.ipynb` which needs additional weight adjustment Q^2) 
* The reciprocal base vectors from diffpy.Structure is not what I expected. It seems struct.lattice.recbase is 3-col vectors. Struct.lattice.base is 3-row vectors.
* Also there is the long-standing phase-factor problem. See https://github.com/mcvine/phonon/blob/5b4224548d72ca25afd00d516c5a13298e06a261/tests/phonon/phase-factor.ipynb 
