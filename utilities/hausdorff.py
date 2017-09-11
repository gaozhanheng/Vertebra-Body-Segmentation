'''
To compute the Hausdoeff distance between two surface meshes

input: two surface mesh in .off format
output: the Huasdorff distance defined by max(D(A,B),D(B,A)), in which D(X,Y) is the distance of point set X to point
set Y. D(X,Y) = max_{x\inX}min_{y\inY}d(x,y), where d(.,.) could be any metric of the space.
'''

import numpy as np
import os
import meshio
from utilities import fn_timer

def diagonalBBX(pointSet):
    '''
    compute the diagonal length of the bounding box of the given point set
    :param pointSet:
    :return:
    '''
    min = np.min(pointSet,axis=0)
    max = np.max(pointSet,axis=0)
    diag = max - min
    return np.sqrt(diag.dot(diag))


def distance(A,B):
    '''
    compute the distance from A to B
    :param A: point set
    :param B: point set
    :return: distance from A to B
    '''
    def pd(a):
        '''
        point distance from a to B
        :param a: a point from A
        :return:
        '''
        c = a - B
        c2 = c * c
        csum = np.sum(c2,axis=1)
        minV = min(csum)
        return minV

    pds = [pd(a) for a in A]
    maxV = max(pds)



    return maxV

@fn_timer
def Hausdorff(surf1,surf2):
    '''
    :param surf1: the file name of surface mesh 1, in .off format
    :param surf2: the file name of surface mesh 2, in .off format
    :return: the Huasforff distance between surface mesh 1 and 2
    '''
    [verts1,faces1] = meshio.readOff(surf1)
    [verts2,faces2] = meshio.readOff(surf2)

    verts1 = np.array(verts1)
    verts2 = np.array(verts2)

    diag1 = diagonalBBX(verts1)
    print 'diagonal length of surface mesh 1:',diag1
    diag2 = diagonalBBX(verts2)
    print 'diagonal length of surface mesh 2:',diag2
    diag = (diag1 + diag2) / 2.
    print 'diagonal length used is:',diag

    dAB = distance(verts1,verts2)
    print 'distance from A to B:',dAB
    dBA = distance(verts2,verts1)
    print 'distance from B to A:',dBA

    Hausdorff = max(dAB,dBA) / diag

    print 'Relative Hausdorff diatance:',Hausdorff * 100.,'\%, relative the average diagonal length of the two bounding' \
                                                          'boxes of the meshes.'

if __name__ == '__main__':
    surf1 = './Data/reg1.off'
    surf2 = './Data/reg1.off'

    Hausdorff(surf1,surf2)