
def readQgridinfo(path):
    from mccomponents.sample.idf.Qgridinfo import read
    reci_basis, shape = read(path)
    return reci_basis, shape
