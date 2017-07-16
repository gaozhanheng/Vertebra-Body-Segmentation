import meshio
import numpy as np

'''
We can also use functions in module vispy.util.transform to generate the transformation matrices easily
'''



def transform(verts,translate=None,scale = None ,translateFirst=False,quiet=True):
    if not translate:
        translate = np.array([0,0,0],dtype=np.float32)
    if not scale:
        scale = np.array([1.,1.,1.],dtype=np.float32)
    verts4 = np.empty(shape=(verts.shape[0],4),dtype=np.float32)
    k = 0
    for v in verts:
        verts4[k][0:3] = v
        verts4[k][3] = 1.
        k += 1

    verts4 = np.matrix(verts4)
    if not quiet:
        print 'before transform'
        print verts4
    verts4 = verts4.transpose()

    transformMatrix = np.matrix('1 0 0 0;0 1 0 0 ;0 0 1 0; 0 0 0 1')

    translateZ_m10 = np.matrix('1 0 0 %s; 0 1 0 %s; 0 0 1 %s; 0 0 0 1' % tuple(translate))

    scale_0p1 = np.matrix('%s 0 0 0;0 %s 0 0; 0 0 %s 0; 0 0 0 1' % tuple(scale))

    if not translateFirst:
        # this way,  the scale is firstly applied and the translation is applied secondly,
        # both are applied TO THE POINTS  and according to the ORIGINAL coordinate system.
        transformMatrix = transformMatrix * translateZ_m10 # note: it's not right if use *=
        transformMatrix = transformMatrix * scale_0p1
    else:
        # this way,  the scale is secondly applied and the translation is applied firstly,
        # both are applied TO THE POINTS  and according to the ORIGINAL coordinate system.
        # and thus is the inverse transform to the above one
        transformMatrix = transformMatrix * scale_0p1   # note: it's not right if use *=
        transformMatrix = transformMatrix * translateZ_m10

    verts4 = transformMatrix * verts4
    verts4 = verts4.transpose()
    verts4 = np.array(verts4)

    if not quiet:
        print 'after transform'
        print verts4

    k = 0
    for v in verts4:
        verts[k][0:3] = verts4[k][0:3] / verts4[k][3]
        k += 1

    if not quiet:
        print 'real output'
        print verts

if __name__ == '__main__':
    filename = './Data/transformed_reg1_simplifiedTo2000.off'

    verts, faces = meshio.readOff(filename)
    verts = np.array(verts)
    faces = np.array(faces)
    trans = [-5,0,4] # [50,50,50] is ok, try [44,44,44],[30,30,30],[25,25,25], [20,20,20], [15,15,15],[10,10,10],[5,5,5], quite funny,
    # and [-50,-50,-50]
    scale = [1.,1.,1.]
    transform(verts,trans,scale)
    meshio.writeOff('./Data/tmp.off',verts,faces)




