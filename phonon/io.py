
def readQgridinfo(path):
    # read reciprocal basis vectors from Qgridinfo
    import os
    lines = open(path).readlines(6)
    d = {}
    for line in lines[:6]:
        exec(line, d)
        continue
    reci_basis = [d['b1'], d['b2'], d['b3']]
    shape = d['n1'], d['n2'], d['n3']
    return reci_basis, shape


