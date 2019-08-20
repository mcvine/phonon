# -*- Python -*-

import numpy as np
import periodictable as pt


class AtomicScattering:
    
    """This class gather methods related to calculation of scattering
    property of an element.
    """

    def __init__(self, element, isotope=None):
        """AtomicScattering("H", 2)
        AtomicScattering("H")
        """
        element = ''.join(c for c in element if c.isalpha()) # remove numbers
        self.element = element
        ptel = getattr(pt, element)
        if isotope is not None:
            self._e = ptel[isotope]
        else:
            self._e = ptel
        self.mass = self._e.mass
        return


    def sigma(self):
        "total cross section. barn"
        return self._e.neutron.total


    def b(self):
        "bound scattering length. fm"
        return self._e.neutron.b_c


    def sigma_inc(self):
        "incoherent scattering cross section. barn"
        return self._e.neutron.incoherent


    def sigma_abs(self):
        "absorption cross section. barn"
        return self._e.neutron.absorption


# End of file
