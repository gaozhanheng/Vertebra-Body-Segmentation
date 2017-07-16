import numpy as np

import meshio

filename = '../Data/registered_reg1.off'

varray,farray = meshio.readOff(filename)

varray = np.array(varray)
ra = 0.
for v in varray:
    r = np.sqrt(((v - varray[0])**2).sum())
    ra = np.maximum(r,ra)

print 'radius of',filename,'is ',ra