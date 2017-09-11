#!/usr/bin/env python
# encoding: utf-8

'''
input 3D Dicom images containing spine, this code will segmente all the vertebra bodies
'''
from utilities import *
from utilities.utilities import *

from algorithms import *


def main():
    ''' Segment a 3D spine Dicom image with pre-trained mesh model.'''
    '''
    for EI, do the initial process
    for SSM, do the right icp
    for SB, put the control grid in the near neighbor of mesh
    '''

    # test transformmesh
    # verts = np.array([[2,1,0]],dtype=np.float32)
    # transformmesh.transform(verts,np.array([-1,-2,-3]),np.array([2,3,4]))
    # transformmesh.transform(verts,np.array([1,2,3]),np.array([1./2,1./3,1./4]),inverse=True)
    # return


    # First, load the mean model and the 3D Dicom data
    meshinfo = {'filename':'./Data/transformed_reg1_simplifiedTo4000.off','centerLowerEndplate':[0,0,0],'centerUpperEndplate':[0,0,1]}
    meanMesh = loadMesh('meanmesh',meshinfo)
    dicominfo = {'dirname':'/Users/mac/Downloads/Dicom Data/2585255-profYu','centerLowerEndplate':[0,0,0],'centerUpperEndplate':[0,0,1]}
    dicomImage = loadDICOM(dicominfo)

    # to convert mesh to fit the volume
    shapeForComp = dicomImage.vol_data_forComp.shape # to get the shape of the dicom data for computation
    spacing = [dicomImage.spacingX,dicomImage.spacingY,dicomImage.spacingSlice] # to get the spacing of the raw dicom images
    firstPos = dicomImage.firstLocation # the first pos is one end of the dicom images
    # print 'firstpos',firstPos
    imagePos = dicomImage.imagePosition
    # print 'impos',imagePos
    lastPos = imagePos[2]
    # print 'lastpos',lastPos
    width = (shapeForComp[0] - 1.) * spacing[0] # to compute the width, height and thickness of dicom images
    height = (shapeForComp[1] - 1.) * spacing[1]
    thickness = np.abs(firstPos - lastPos)
    # print width,height,thickness
    trans = [width - imagePos[0],height - imagePos[1],-lastPos] # to transform the data to fit integral space
    print 'trans:',trans
    scale = [shapeForComp[0]/width,shapeForComp[1]/height,shapeForComp[2]/thickness]
    print 'scale:',scale
    transformmesh.transform(meanMesh.verts, trans, scale,translateFirst=True)
    # meanMesh.updateMesh()

    # show the mesh and 3D image in the same window
    # meanMesh.show()
    # ct_aaa = getColorMaps('ct_aaa')
    # dicomImage.show(cmap = ct_aaa)

    # to test volume visualization, but no use. cmap is the only way to define the colors
    # vol_data = np.ones(100*100*100*4).reshape(10,100,1000,4) # * 0.7
    # for x in range(10):
    #     for y in range(100):
    #         for z in range(1000):
    #             vol_data[x][y][z][0] = 1.#np.abs(x - 4.5) / 4.5
    #             vol_data[x][y][z][1] = .5#np.abs(y - 49.5) / 49.5
    #             vol_data[x][y][z][2] = 1. #np.abs(z - 499.5) / 499.5
    #             vol_data[x][y][z][3] = 1. #random.random()
    # volume = visuals.Volume(vol_data, clim=np.array([0,1]), method='iso',threshold = 0.1, relative_step_size=0.8, emulate_texture=False)
    # view.add(volume)

    # v, f = openmesh.readOff('./Data/sphere.off', quiet=True)
    # sphere = visuals.Mesh(vertices=np.array(v),faces=np.array(f),color=np.array([[1.0, 1.0, 1.0,0.9]]),shading='flat')
    # view.add(sphere)
    # v, f = openmesh.readOff('./Data/x.off', quiet=True)
    # x = visuals.Mesh(vertices=np.array(v), faces=np.array(f), color=np.array([[1.0, 0.0, 0.0, 0.9]]),
    #                       shading='flat')
    # view.add(x)
    # v, f = openmesh.readOff('./Data/y.off', quiet=True)
    # y = visuals.Mesh(vertices=np.array(v), faces=np.array(f), color=np.array([[1.0, 1.0, 0.0, 0.9]]),
    #                       shading='flat')
    # view.add(y)
    # v, f = openmesh.readOff('./Data/z.off', quiet=True)
    # z = visuals.Mesh(vertices=np.array(v), faces=np.array(f), color=np.array([[0.0, 0.0, 1.0, 0.9]]),
    #                       shading='flat')
    # view.add(z)
    # view.camera.set_range()

    # view.add(visuals.Cube())
    # meanMesh.hide()

    # Second, initialize the mesh to the approximate position in the image
    # now, this function does nothing
    initMeshInDicom(meanMesh,dicomImage)

    # Third, Optimize the mesh

    model = SIM.SIM(meanMesh,dicomImage)
    # model = LevelSet.LevelSet(meanMesh,dicomImage)
    S = model.Optimize()

    # OptimizeMesh(meanMesh,Volume())

    # S.show()
    itrans = [-trans[0],-trans[1],-trans[2]]
    iscale = [1./scale[0],1./scale[1],1./scale[2]]
    if isinstance(S[0],np.ndarray) and S[0].all():
        transformmesh.transform(S[0],itrans,iscale)
        openmesh.writeOff('./Data/registered_reg1_SI.off', S[0], meanMesh.faces)
    if isinstance(S[1],np.ndarray) and S[1].all():
        transformmesh.transform(S[1], itrans, iscale)
        openmesh.writeOff('./Data/registered_reg1_SSM.off', S[1], meanMesh.faces)
    if isinstance(S[2],np.ndarray) and S[2].all():
        # pass
        transformmesh.transform(S[2], itrans, iscale)
        openmesh.writeOff('./Data/registered_reg1_SB.off', S[2], meanMesh.faces)

    # vispy.app.run()

if __name__ == "__main__":
    main()
