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
            ns = ptel[isotope].neutron
        else:
            ns = ptel.neutron
        self.ns = ns
        return


    def sigma(self):
        "total cross section. barn"
        return self.ns.total


    def b(self):
        "bound scattering length. fm"
        return self.ns.b_c


    def sigma_inc(self):
        "incoherent scattering cross section. barn"
        return self.ns.incoherent


    def sigma_abs(self):
        "absorption cross section. barn"
        return self.ns.absorption


# End of file
