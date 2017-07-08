import meshio
import numpy as np



filename = './Data/reg1.off'

verts,faces = meshio.readOff(filename)
verts = np.array(verts)
faces = np.array(faces)

verts4 = np.empty(shape=(verts.shape[0],4),dtype=np.float32)
k = 0
for v in verts:
    verts4[k][0:3] = v
    verts4[k][3] = 1.
    k += 1

verts4 = np.matrix(verts4)
print 'before transform'
print verts4
verts4 = verts4.transpose()

transformMatrix = np.matrix('1 0 0 0;0 1 0 0 ;0 0 1 0; 0 0 0 1')

translateZ_m10 = np.matrix('1 0 0 0; 0 1 0 0; 0 0 1 -10; 0 0 0 1')

transformMatrix *= translateZ_m10

verts4 = transformMatrix * verts4
verts4 = verts4.transpose()
verts4 = np.array(verts4)

print 'after transform'
print verts4

k = 0
for v in verts4:
    verts[k][0:3] = verts4[k][0:3] / verts4[k][3]
    k += 1

print 'real output'
print verts




