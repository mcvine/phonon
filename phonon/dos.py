#!/usr/bin/env python

from mccomponents.sample.idf import Omega2, DOS, units
import numpy as np

def fromOmaga2():
    """create DOS in IDF format from Omega2 in IDF format"""
    (filetype, version, comment), omega2 = Omega2.read()
    omega2[ omega2<0 ] = 0
    energies = np.sqrt(omega2) * units.hertz2mev

    max = np.max(energies)
    print max

    max = (max+10)//10*10.
    step = max/500; step = 10**np.round(np.log10(step))
    boundaries = np.arange(0, max, step)

    hist, edges = np.histogram(energies, boundaries)

    x = (edges[1:] + edges[:-1])/2
    #import pylab
    #pylab.plot(x, hist)
    #pylab.show()

    DOS.write(x, hist, E_unit="meV")
    return
