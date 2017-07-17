import meshio
import numpy as np
import math

'''
We can also use functions in module vispy.util.transform to generate the transformation matrices easily
'''



def transform(verts,translate=None,scale = None ,rx=None,ry=None,rz=None,translateFirst=False,quiet=True):
    if not translate:
        translate = np.array([0,0,0],dtype=np.float32)
    if not scale:
        scale = np.array([1.,1.,1.],dtype=np.float32)
    if not rx:
        rx = 0
    if not ry:
        ry = 0
    if not rz:
        rz = 0

    verts4 = np.empty(shape=(verts.shape[0],4),dtype=np.float32)
    k = 0
    for v in verts:
        verts4[k][0:3] = v
        verts4[k][3] = 1.
        k += 1

    barycenter = verts4.sum(axis = 0) / verts4.shape[0]
    print 'center of model'
    print barycenter

    verts4 = np.matrix(verts4)
    if not quiet:
        print 'before transform'
        print verts4
    verts4 = verts4.transpose()

    # first, translate so that the center of the mesh is the origin
    # then, do rotation. the rotation order is: z,y,x
    # then translate back to the center of the mesh
    # after that, do translate and scale or scale and translate
    rotationMatrix = np.matrix('1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1')
    translateToBarycenter = np.matrix('1 0 0 %s;0 1 0 %s; 0 0 1 %s; 0 0 0 1' % (barycenter[0],barycenter[1],barycenter[2]))
    rotationMatrix = rotationMatrix * translateToBarycenter
    c = math.cos(rx)
    s = math.sin(rx)
    rxMatrix = np.matrix('1 0 0 0;0 %s %s 0; 0 %s %s 0; 0 0 0 1' %(c,s,-s,c))
    rotationMatrix = rotationMatrix * rxMatrix
    c = math.cos(ry)
    s = math.sin(ry)
    ryMatrix = np.matrix('%s 0 %s 0;0 1 0 0;%s 0 %s 0;0 0 0 1' %(c,s,-s,c))
    rotationMatrix = rotationMatrix * ryMatrix
    c = math.cos(rz)
    s = math.sin(rz)
    rzMatrix = np.matrix('%s %s 0 0;%s %s 0 0;0 0 1 0;0 0 0 1' % (c,s,-s,c))
    rotationMatrix = rotationMatrix * rzMatrix
    translateToOrigin = np.matrix('1 0 0 %s;0 1 0 %s; 0 0 1 %s; 0 0 0 1' % (-barycenter[0],-barycenter[1],-barycenter[2]))
    rotationMatrix = rotationMatrix * translateToOrigin

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

    transformMatrix = transformMatrix * rotationMatrix

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
    filename = './Data/.off'

    verts, faces = meshio.readOff(filename)
    verts = np.array(verts)
    faces = np.array(faces)
    trans = [0,22,5] # [50,50,50] is ok, try [44,44,44],[30,30,30],[25,25,25], [20,20,20], [15,15,15],[10,10,10],[5,5,5], quite funny,
    # and [-50,-50,-50]
    scale = [1.,1.,1.]
    rx = 0.
    ry = 0.
    rz = 0.
    transform(verts,trans,scale,rx,ry,rz)
    meshio.writeOff('./Data/tmp.off',verts,faces)




