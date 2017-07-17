#!/usr/bin/env python
# encoding: utf-8

'''
input 3D Dicom images containing spine, this code will segmente all the vertebra bodies
'''
import sys,os,shutil,getopt,time,random
# import vtk
import copy
from PIL import Image
import math
import numpy as np
import scipy
from scipy import optimize as opt,io,linalg,ndimage,misc,sparse,interpolate
from scipy.interpolate import UnivariateSpline,BivariateSpline,BSpline
from scipy.sparse import linalg as splinalg
import skimage
from skimage import color,io,data,segmentation,img_as_ubyte,img_as_float,feature,\
    filters,measure,graph,transform,draw,external,restoration
import matplotlib.pyplot as pl

import dicom # this is the module name of pydicom

import numpy as np
import vispy.scene
from vispy.color import *
from vispy import io,gloo
from vispy.scene import visuals
from vispy.gloo import wrappers


from icp import basicICP


# from vispy import visuals

# try:
#     import openmesh
# except:
#     print 'no openmesh module. using meshio to save mesh (only support .off files)'
import meshio as openmesh
import transformmesh

from utilities import fn_timer
import re

# print wrappers.get_current_canvas()
# print wrappers.get_gl_configuration()
# print wrappers.get_state_presets()

# wrappers.global_gloo_functions.set_polygon_offset(1,1)

#
    # Make a canvas and add simple view
    #
canvas = vispy.scene.SceneCanvas(keys='interactive', show=True,config={'depth_size':24},bgcolor='white')
view = canvas.central_widget.add_view()
GLContext = canvas.context
GLContext.set_polygon_offset(0,0)


view.camera = 'arcball'  # or try 'arcball' . input an invalid and see the error for more options
view.camera.fov = 45 # important

# add a colored 3D axis for orientation
# axis = visuals.XYZAxis(parent=view.scene)
# visuals.GridLines(parent = view.scene)


def addPointCloud():
    ### this is the example code of vispy
    # first, show a point cloud given the coordinates
    pos = np.random.normal(size=(100000, 3), scale=0.2)
    # one could stop here for the data generation, the rest is just to make the
    # data look more interesting. Copied over from magnify.py
    centers = np.random.normal(size=(50, 3))
    indexes = np.random.normal(size=100000, loc=centers.shape[0] / 2.,
                               scale=centers.shape[0] / 3.)
    indexes = np.clip(indexes, 0, centers.shape[0] - 1).astype(int)
    scales = 10 ** (np.linspace(-2, 0.5, centers.shape[0]))[indexes][:, np.newaxis]
    pos *= scales
    pos += centers[indexes]
    print 'type of pos:', type(pos)
    print 'shape of pos:', pos.shape

    # create scatter object and fill in the data
    scatter = visuals.Markers()
    scatter.set_data(pos, symbol='o', edge_color=None, face_color=(1, 1, 1, .5), size=5)
    view.add(scatter)

    # vispy.app.run()

def getColorMaps(mapName,withgradientopacity=True):
    import colormaps

    colors = colormaps.cmaps[mapName]['colors']
    # controls = colormaps.cmaps[mapName]['controls']
    # mincontrol = min(controls)
    # maxcontrol = max(controls)
    #
    # controls = (controls - mincontrol) / (maxcontrol - mincontrol)
    # controls[0] = 0.
    # controls[-1] = 1.
    colors[0][3] = 0.
    colors[1][3] = .0 #4./255
    colors[2][3] = 0.
    colors[3][3] = .0 #1./255
    colors[4][3] = 1.
    colors[5][3] = .0 #1./255

    cmap = Colormap(colors=colors)

    return cmap

'''
to count and print the time of running function automatically.
:param function: any function
:return: None
'''
# print the information to console and file simultaneously
def myPrint(filename,content,mode = 'a'):
    print content
    with open(filename,mode) as f:
        f.writelines(content)
        f.writelines('\n')

class Mesh(object):
    ''' Class to represet a 3D triangular mesh. nv and nf are the number of vertices and faces of the mesh, verts and
    faces are the numpy.array objects to store the coordinates of the vertices and the indices of the faces.'''
    def __init__(self):
        self._nv = 0
        self._nf = 0
        self._verts = np.empty(shape = (self._nv,3),dtype=np.float32)
        self._faces = np.empty(shape = (self._nf,3),dtype=np.int32)

        self._mesh = None # to store the current visual mesh so that we only add it into the scene once
        self._outline = None

        # the centers of the low and upper endplates of the vb mesh
        self._centerLowerEndplate = np.zeros((1,3),dtype=np.float32)
        self._centerUpperEndplate = np.zeros((1,3),dtype=np.float32)

        self._phim = None # the \phi m in SSM model. a list, each element is an nv * 3 ndarray

    def __copy__(self):
        mesh = Mesh()
        mesh.verts = self._verts.copy()
        mesh.faces = self._faces.copy()
        mesh.phim = copy.deepcopy(self._phim)
        return mesh

    @property
    def nv(self):
        return self._nv

    @property
    def nf(self):
        return self._nf

    @property
    def verts(self):
        return self._verts
    @verts.setter
    def verts(self,verts):
        if not verts.shape[1] == 3:
            return ValueError('Invalid input verts: it must be of shape [x][3]')
        self._verts = verts
        self._nv = verts.shape[0]
        self._updateMesh()

    @property
    def faces(self):
        return self._faces
    @faces.setter
    def faces(self,faces):
        if not faces.shape[1] == 3:
            return ValueError('Invalid input faces: it must be of shape [x][3]')
        self._faces = faces
        self._nf = faces.shape[0]
        self._updateMesh()

    @property
    def phim(self):
        return self._phim
    @phim.setter
    def phim(self,phim):
        self._phim = phim

    @property
    def centerLowerEndplate(self):
        return self._centerLowerEndplate
    @centerLowerEndplate.setter
    def centerLowerEndplate(self,c):
        self._centerLowerEndplate = c

    @property
    def centerUpperEndplate(self):
        return self._centerUpperEndplate
    @centerUpperEndplate.setter
    def centerUpperEndplate(self,c):
        self._centerUpperEndplate = c

    def _updateMesh(self):
        # add mesh data into the scene
        # Vertices positions
        # p = np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1],
        #               [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, -1]])
        # # p = p * 1000.0
        # faces = np.array([[0, 2, 1],
        #                   [0, 3, 2],
        #                   [0, 4, 3],
        #                   [0, 5, 4],
        #                   [0, 6, 5],
        #                   [0, 1, 6],
        #                   [1, 7, 6],
        #                   [1, 2, 7],
        #                   [7, 3, 4],
        #                   [7, 2, 3],
        #                   [4, 6, 7],
        #                   [4, 5, 6]
        #                   ])
        p = self._verts.copy()
        faces = []
        for face in self._faces:
            faces.append([face[2],face[1],face[0]])
        faces = np.array(faces) # otherwise, it is dark.
        # faces = self._faces


        # for more deatils on color, please refer to vispy/colors/colorarray.py
        # colors_p = np.array(
        #     [[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1],
        #      [1, 1, 1, 1]])
        # colors_f = np.array([(1.0, 0.1, 0.0, 0.3) for k in range(faces.shape[0])])
        color_mesh = np.array([[1.0, 1.0, 0.3,0.9]])
        # for more details see vispy/visuals/mesh.py
        mesh = visuals.Mesh(vertices=p, faces=faces,
                            # vertex_colors = colors_p,
                            # face_colors = colors_f,
                            color=color_mesh,
                            shading='flat',
                            mode='triangles')

        self._mesh = mesh

        # outline_p = []
        # for face in faces:
        #     outline_p.append(p[face[0]])
        #     outline_p.append(p[face[1]])
        #     outline_p.append(p[face[1]])
        #     outline_p.append(p[face[2]])
        #     outline_p.append(p[face[2]])
        #     outline_p.append(p[face[0]])
        # outline_p = np.array(outline_p)
        #
        # outline = visuals.Line(pos=outline_p, color=(0.1, 0.1, 0.1, 1.), width=1,
        #                        connect='segments', method='gl', antialias=False)
        # self._outline = outline

    def updateMesh(self):
        self._updateMesh()

    def show(self,withOutline=False):
        view.add(self._mesh)
        if withOutline:
            GLContext.set_polygon_offset(1, 1)
            self._mesh.update_gl_state(**{'polygon_offset_fill': True})

            view.add(self._outline)

        # cube = visuals.Cube(edge_color = [1.0,0,0,1])
        # view.add(cube)

        view.camera.set_range()  # similar to viewAll in Coin3D
        # vispy.app.run

    def hide(self):
        self._mesh.remove_parent(view.scene)
        if self._outline in view.scene.children:
            self._outline.remove_parent(view.scene)

class Volume(object):
    ''' Class to represent a 3D Dicom object. numX, numY and numZ are the number of pixels in each directions.
    pixel_array stores the pixel value of the Dicom object slice by slice. Note that the Z direction is Axis for
    CT data and spinal for MR image.'''
    def __init__(self):

        self._spacingX = 0.
        self._spacingY = 0.
        self._spacingSlice = 0. # the spacing along each direction

        self._pixel_array = None # the copy of pixel_array of dicom files read by pydicom
        self._gradient_array = None # the gradient of the pixel array for computing E_I
        self._hessian_array = None # this is used for the same as gradient array

    @property
    def spacingX(self):
        return self._spacingX
    @spacingX.setter
    def spacingX(self,s):
        self._spacingX = s
        # self._updateVolume()
    @property
    def spacingY(self):
        return self._spacingY
    @spacingY.setter
    def spacingY(self,s):
        self._spacingY = s
        # self._updateVolume()
    @property
    def spacingSlice(self):
        return self._spacingSlice
    @spacingSlice.setter
    def spacingSlice(self,s):
        self._spacingSlice = s
        # self._updateVolume()

    @property
    def pixelArray(self):
        return self._pixel_array
    @pixelArray.setter
    def pixelArray(self,voldata):
        self._pixel_array = voldata
        self._updateVolume()

    @property
    def gradientArray(self):
        return self._gradient_array
    @property
    def hessianArray(self):
        return self._hessian_array

    # the first dim must be the slice, while the second and third dim are the x and y pixels
    @property
    def dataShape(self):
        return self._pixel_array.shape
    @dataShape.setter
    def dataShape(self,shape):
        self._pixel_array = np.zeros(shape=shape)

    def _updateVolume(self):
        '''
        First, compute the gradient of the pixel volume
        :return:
        '''
        voldata = self.vol_data_forComp
        # print 'data',voldata
        gradientArray = np.gradient(voldata)
        self._gradient_array = gradientArray
        # print 'gradient'
        # print gradientArray[0],gradientArray[1],gradientArray[2]
        # print '100,0,0',gradientArray[0][100][0][0],gradientArray[1][100][0][0],gradientArray[2][100][0][0]
        # print '0,100,0',gradientArray[0][0][100][0],gradientArray[1][0][100][0],gradientArray[2][0][100][0]
        # print '100,100,0',gradientArray[0][100][100][0],gradientArray[1][100][100][0],gradientArray[2][100][100][0]

        hessianX = np.gradient(gradientArray[0]) # \frac{\partial^2 f,\partial x^2},\frac{\partial f,
                                                #  \partial x \partial y},\frac{\partial f, \partial x \partial z}
        # print 'length',len(hessianX)
        # print 'shape',hessianX[0].shape
        # print 'hx',hessianX[0]
        hessianY = np.gradient(gradientArray[1])
        hessianZ = np.gradient(gradientArray[2])
        # print 'hx',hessianX[0]

        # a = hessianX[1] - hessianY[0]
        # print 'xy:',a
        # a = hessianX[2] - hessianZ[0]
        # print 'xz:', a
        # a = hessianY[2] - hessianZ[1]
        # print 'yz:', a
        self._hessian_array = [hessianX,hessianY,hessianZ]

    def show(self,clim=None,cmap = None):
        # add a volume into the scene
        # vol_data = np.load(io.load_data_file('brain/mri.npz'))['data']
        # print type(vol_data)
        # print vol_data.shape
        # vol_data = np.flipud(np.rollaxis(vol_data, 1))
        # print type(vol_data)
        # print vol_data.shape

        vol_data = self._pixel_array.copy()
        # vol_data = np.gradient(self._pixel_array)
        # vol_data = np.sqrt(vol_data[0] ** 2 + vol_data[1] ** 2 + vol_data[2] ** 2)
        # print vol_data.shape
        clim = clim or [-3024, 3071]
        cmap = cmap or 'grays'
        volume = visuals.Volume(vol_data, clim=clim,cmap = cmap,method = 'mip',relative_step_size=0.8,emulate_texture=False)

        view.add(volume)

        view.camera.set_range()  # similar to viewAll in Coin3D

@fn_timer
def loadMesh(meshName,meshInfo):
    ''' Load a mesh from meshName and return a Mesh object'''

    v,f = openmesh.readOff(meshInfo['filename'],quiet=True)

    mesh = Mesh()
    mesh.verts = np.array(v,dtype=np.float32)
    mesh.faces = np.array(f,dtype=np.int32)

    # for now, the centers of lower and upper end plate are manually input
    mesh.centerUpperEndplate = np.array(meshInfo['centerUpperEndplate'])
    mesh.centerLowerEndplate = np.array(meshInfo['centerLowerEndplate'])

    phim = []
    for k in range(1):
        phim.append(np.zeros(shape=(mesh.nv * 3,)).reshape(mesh.nv,3))
    mesh.phim = phim

    return mesh

@fn_timer
def loadDICOM(dicominfo):
    ''' Load all dicom files from dicomDir. Integrate and return the data with a whole Dicom3D object.'''
    dicomDir = dicominfo['dirname']
    if not os.path.isdir(dicomDir):
        return ValueError('Input %s must be a dir.' % dicomDir)

    allfiles = os.listdir(dicomDir)
    allDicoms = [int(os.path.splitext(file)[0]) for file in allfiles if os.path.splitext(file)[1].lower() == '.dcm']
    allDicoms = np.array(allDicoms)
    allDicoms.sort()
    if len(allDicoms) == 0:
        return ValueError('No dcm files detected in the dir: %s' % dicomDir)

    image = Volume()
    vol_data = []
    k = 0
    lastSliceLocation = 0
    for dcm in allDicoms:
        filename = dicomDir + '/%s' % dcm + '.dcm'
        if not os.path.isfile(filename):
            print filename,'is not a file'
        ds = dicom.read_file(filename)

        # print '='*20
        # print dcm
        if not ds.SeriesNumber == 5:
            # print 'series number',ds.SeriesNumber
            # print 'skip dicom file: %s' % dcm
            continue

        # print 'series number:', ds.SeriesNumber
        # print 'rows and cols:',ds.Rows,ds.Columns
        # print 'slice location:',ds.SliceLocation
        # print 'single collimation width:',ds.SingleCollimationWidth
        # print 'data collectin diameter:',ds.DataCollectionDiameter
        # print 'reconstruction diameter:',ds.ReconstructionDiameter
        # print 'Window center:',ds.WindowCenter
        # print 'window width:',ds.WindowWidth
        # print 'window center width explanation:',ds.WindowCenterWidthExplanation
        # print 'image type:',ds.ImageType
        # print 'image position:',ds.ImagePositionPatient

        lastSliceLocation = ds.SliceLocation

        if k == 0:
            image.spacingX = ds.PixelSpacing[0]
            image.spacingY = ds.PixelSpacing[1]  # unit: mm

            image.imagePosition = ds.ImagePositionPatient # dynamic attribute
            image.firstLocation = ds.SliceLocation # dynamic attribute
            image.rescaleIntercept = ds.RescaleIntercept # this and the following one are
            #  used for computing the real value of ct
            image.rescaleSlope = ds.RescaleSlope
            print 'ds0'
            print ds
        if k == 1:
            print 'ds1'
            print ds
            z1 = ds.SliceLocation

            image.spacingSlice = np.abs(image.firstLocation - z1)

        # if k == 14:
        #     print 'ds14',ds.pixel_array
        k = k + 1
        vol_data.append(np.flip(np.flip(ds.pixel_array,axis=0),axis=1).astype(np.float32))

    ds.save_as()

    vol_data = np.flip(vol_data, axis=0)
    # print 'vol data',vol_data
    # vol_data = ndimage.median_filter(vol_data,size=3)
    vol_data = ndimage.gaussian_filter(np.array(vol_data),sigma=2,order=0)
    vol_data_forComp = np.rollaxis(np.rollaxis(vol_data,axis=1),axis=2)# we don't change vol_data_forComp so no need to copy vol_data

    print 'shape of vol_data:',vol_data.shape
    print 'shape of vol_data for computation:',vol_data_forComp.shape


    if True:
        vol = vol_data.copy()
        vol = np.flip(np.flip(np.flip(vol,axis=0),axis=1),axis=2)
        print 'shape of vol for observing',vol.shape
        outputDicomDir = '/Users/mac/Downloads/Dicom Data/filteredVol'
        if os.path.exists(outputDicomDir):
            shutil.rmtree(outputDicomDir)
        os.mkdir(outputDicomDir)

        k = 0
        for dcm in allDicoms:
            filenameIn = dicomDir + '/%s' % dcm + '.dcm'
            filenameOut = outputDicomDir + '/%s' % dcm + '.dcm'
            if not os.path.isfile(filenameIn):
                print filenameIn,'is not a file'
            ds = dicom.read_file(filenameIn)

            if not ds.SeriesNumber == 5:
                continue

            ds.StudyDate = time.strftime('%Y%m%d',time.localtime())
            ds.SeriesDate = time.strftime('%Y%m%d',time.localtime())
            pixelarray = vol[k].astype(np.uint16)
            # for n, val in enumerate(pixelarray.flat):
            #     if val < 300:
            #         pixelarray.flat[n] = 0
            ds.PixelData = pixelarray.tostring()
            ds.save_as(filenameOut)
            k += 1
        print 'finished writing dicom data for obversing'

    # for x in range(vol_data.shape[0]):
    #     for y in range(vol_data.shape[1]):
    #         for z in range(vol_data.shape[2]):
    #             if not vol_data[x][y][z] == vol_data_forComp[z][y][x]:
    #                 print 'roll axis error:',x,y,z

    # # to test if the optimization code works well on simple volumes
    # vol_data = np.empty(shape=(101,101,101),dtype=np.float32)
    # c = np.array([50,50,50],dtype=np.float32) # center of the sphere
    # for x in range(101):
    #     for y in range(101):
    #         for z in range(101):
    #             p = np.array([x,y,z],dtype=np.float32)
    #             ra = np.array([2,2,1],dtype=np.float32) # when this is [1,1,1], the objective function is a sphere,
    #             # eclipsoid otherwise
    #             dist = np.sqrt(np.sum((p - c) * (p - c) / ra))
    #             vol_data[x][y][z] = -np.square(dist - 25.) # the data is also critical for the cnvergence. try np.abs
    # vol_data_forComp = vol_data
    # print 'vol data for test',vol_data

    image.vol_data_forComp = vol_data_forComp  # this is a dynamic attribute

    image.pixelArray = vol_data

    image.imagePosition = np.array([0,0,0]) - image.imagePosition
    image.imagePosition[2] = lastSliceLocation
    print '='*20
    print '%s slice loaded' % k
    print 'min and max value in the volume:',np.min(vol_data),np.max(vol_data)
    print 'shape of image:',image.dataShape
    print 'image position:',image.imagePosition
    print 'image spacing:',image.spacingX,image.spacingY,image.spacingSlice

    return image

@fn_timer
def initMeshInDicom(mesh,volume):
    ''' Put the mesh to a proper position in the 3D Dicom image. Coordinates of vertices are modified.'''
    pass

@fn_timer
def OptimizeMesh(mesh,volume):
    ''' Optimize the mesh.'''
    S_I = copy.copy(mesh)
    S_B = copy.copy(mesh)
    S_SSM = copy.copy(mesh)

    global km
    km = 1
    dirname = './Data/tmp'
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.mkdir(dirname)

    def OptimizeSI(): # more details about scipy.optimization please refer to the optimize __init__.py
        nv = S_I.nv # the number of vertices used to test the optimization efficiency
        # print 'number of vertices used:',nv

        def showMeshAndNormal():
            '''
            To test if the normals is pointing from inside to outlide and if the ordering of the normals is right.
            :return:
            '''
            m = visuals.Mesh(S_B.verts,mesh.faces)
            l = visuals.Line()
            view.add(m)
            outline_p = []
            k = 0
            for v in S_B.verts:
                outline_p.append(v)
                outline_p.append(v + normals[k])
                k += 1

            outline_p = np.array(outline_p)

            outline = visuals.Line(pos=outline_p, color=(0.1, 0.1, 0.1, 1.), width=1,
                                   connect='segments', method='gl', antialias=False)
            view.add(outline)
            vispy.app.run()
        # showMeshAndNormal()

        def normalizeVec(vec):
            '''

            :param vec:
            :return:
            '''
            mag = vec.dot(vec)
            vec /= mag

        # these two parameters really affect the convergence and final result of the optimization process. try (1.,0.)
        # (.8,.2), and (.2,.8) for the eclipsoid case. Some results are funny.
        omega1 = .5
        omega2 = .3
        signGradient = 1.
        signPixel = 1. # put these two here to convinently modify. (-1.,-1.) for eclipsoid case. signPixel is critical
        # for success
        # if the function is defined well, it can also converge to the right result without help of gradient
        # with gradient help, the convergence can be speeded up a lot
        # the initial guess is QUITE important for the final result. See transform.py
        # right fun + wrong gradient can not guarantee the right result
        # fun and gradient must be concide, wired otherwise
        # of course you can only give fun without gradient, but, please make sure the given gradient is as precise as possible,
        # if you prefer to give.

        def E_I(SIVerts):
            SIVerts = SIVerts.reshape((nv,3))

            meshData = vispy.geometry.MeshData(vertices=SIVerts, faces=mesh.faces)  # the faces are all the same
            normals = meshData.get_vertex_normals()

            ei = 0.
            pixelarray = volume.vol_data_forComp
            gradientArray = volume.gradientArray
            k = 0
            # sum all the values
            gradientX = gradientArray[0]
            gradientY = gradientArray[1]
            gradientZ = gradientArray[2]
            dataShape = pixelarray.shape
            for v in SIVerts:
                def trilinearInterp():
                    x = int(v[0])
                    x2 = int(v[0] + 1.)
                    y = int(v[1])
                    y2 = int(v[1] + 1.)
                    z = int(v[2])
                    z2 = int(v[2] + 1.)

                    if x < 0 or y < 0 or z < 0 or x >= dataShape[0] or y >= dataShape[1] or z >= dataShape[2]:
                        print 'vert %s is outside of volume data in E_I:' % k, x, y, z, 'vert coordinate:', v
                        return 0,0 # pixel and then gradient

                    # use trilinear interpolation, refer to wiki
                    xd = v[0] - x
                    yd = v[1] - y
                    zd = v[2] - z

                    X = np.array([x,x2])
                    Y = np.array([y,y2])
                    Z = np.array([z,z2])
                    vPixel = np.empty((2,2,2),dtype=np.float32)
                    vGradient = np.empty((2,2,2),dtype=np.float32)
                    for i in range(2):
                        for j in range(2):
                            for kk in range(2):
                                gradient = np.array([gradientX[X[i]][Y[j]][Z[kk]],
                                                     gradientY[X[i]][Y[j]][Z[kk]],
                                                     gradientZ[X[i]][Y[j]][Z[kk]]])
                                vGradient[i][j][kk] = omega2 * (normals[k].dot(gradient))
                                vPixel[i][j][kk] = omega1 * pixelarray[X[i]][Y[j]][Z[kk]]
                    def triInt(values):
                        c00 = values[0][0][0] * (1 - xd) + values[1][0][0] * xd
                        c01 = values[0][0][1] * (1 - xd) + values[1][0][1] * xd
                        c10 = values[0][1][0] * (1 - xd) + values[1][1][0] * xd
                        c11 = values[0][1][1] * (1 - xd) + values[1][1][1] * xd
                        c0 = c00 * (1 - yd) + c10 * yd
                        c1 = c01 * (1 - yd) + c11 * yd
                        c = c0 * (1 - zd) + c1 * zd
                        return c
                    vpixel = triInt(vPixel)
                    vgradient= triInt(vGradient)

                    return vpixel,vgradient
                def noInterp():
                    x = int(v[0] + .5)
                    y = int(v[1] + .5)
                    z = int(v[2] + .5)

                    if x < 0 or y < 0 or z < 0 or x >= dataShape[0] or y >= dataShape[1] or z >= dataShape[2]:
                        print 'vert %s is outside of volume data in E_I:' % k, x, y, z, 'vert coordinate:', v
                        return 0,0 # pixel and then gradient

                    gradient = np.array([gradientX[x][y][z],
                                         gradientY[x][y][z],
                                         gradientZ[x][y][z]])
                    vgradient = omega2 * (normals[k].dot(gradient))
                    vpixel = omega1 * pixelarray[x][y][z]
                    return vpixel,vgradient

                vpixel,vgradient = trilinearInterp()
                # vpixel,vgradient = noInterp()

                # val = np.square(np.sqrt(((v - np.array([50.,50.,50.]))**2).sum()) - 25.)
                val = signGradient * vgradient - signPixel * vpixel
                ei += val
                k += 1

            return ei
        def gradient_EI(SIVerts):
            SIVerts = SIVerts.reshape((nv,3))

            meshData = vispy.geometry.MeshData(vertices=SIVerts, faces=mesh.faces)  # the faces are all the same
            normals = meshData.get_vertex_normals()

            pixelarray = volume.vol_data_forComp
            gradientArray = volume.gradientArray
            hessianArray = volume.hessianArray

            gradientX = gradientArray[0]
            gradientY = gradientArray[1]
            gradientZ = gradientArray[2]

            hessianX = hessianArray[0]
            hessianY = hessianArray[1]
            hessianZ = hessianArray[2]

            hessianXx = hessianX[0]
            hessianXy = hessianX[1]
            hessianXz = hessianX[2]

            hessianYx = hessianY[0]
            hessianYy = hessianY[1]
            hessianYz = hessianY[2]

            hessianZx = hessianZ[0]
            hessianZy = hessianZ[1]
            hessianZz = hessianZ[2]

            dataShape = pixelarray.shape
            k = 0
            g = []
            for v in SIVerts:
                def trilinearInterp():
                    x = int(v[0])
                    x2 = int(v[0] + 1.)
                    y = int(v[1])
                    y2 = int(v[1] + 1.)
                    z = int(v[2])
                    z2 = int(v[2] + 1.)

                    if x < 0 or y < 0 or z < 0 or x >= dataShape[0] or y >= dataShape[1] or z >= dataShape[2]:
                        print 'vert %s is outside of volume data in gradient E_I:' % k, x, y, z, 'vert coordinate:', v
                        return 0,0,0,0,0,0 # pixel and then gradient

                    # use trilinear interpolation, refer to wiki
                    xd = v[0] - x
                    yd = v[1] - y
                    zd = v[2] - z

                    X = np.array([x,x2])
                    Y = np.array([y,y2])
                    Z = np.array([z,z2])
                    gXgradient = np.empty((2,2,2),dtype=np.float32)
                    gXpixel = np.empty((2,2,2),dtype=np.float32)
                    gYgradient = np.empty((2, 2, 2), dtype=np.float32)
                    gYpixel = np.empty((2, 2, 2), dtype=np.float32)
                    gZgradient = np.empty((2, 2, 2), dtype=np.float32)
                    gZpixel = np.empty((2, 2, 2), dtype=np.float32)
                    for i in range(2):
                        for j in range(2):
                            for kk in range(2):
                                hessian = np.array([hessianXx[X[i]][Y[j]][Z[kk]],
                                                    hessianXy[X[i]][Y[j]][Z[kk]],
                                                    hessianXz[X[i]][Y[j]][Z[kk]]])
                                gXgradient[i][j][kk] = omega2 * (normals[k].dot(hessian))
                                gXpixel[i][j][kk] = omega1 * gradientX[X[i]][Y[j]][Z[kk]]

                                hessian = np.array([hessianYx[X[i]][Y[j]][Z[kk]],
                                                    hessianYy[X[i]][Y[j]][Z[kk]],
                                                    hessianYz[X[i]][Y[j]][Z[kk]]])
                                gYgradient[i][j][kk] = omega2 * (normals[k].dot(hessian))
                                gYpixel[i][j][kk] = omega1 * gradientY[X[i]][Y[j]][Z[kk]]

                                hessian = np.array([hessianZx[X[i]][Y[j]][Z[kk]],
                                                    hessianZy[X[i]][Y[j]][Z[kk]],
                                                    hessianZz[X[i]][Y[j]][Z[kk]]])
                                gZgradient[i][j][kk] = omega2 * (normals[k].dot(hessian))
                                gZpixel[i][j][kk] = omega1 * gradientZ[X[i]][Y[j]][Z[kk]]
                    def triInt(values):
                        c00 = values[0][0][0] * (1 - xd) + values[1][0][0] * xd
                        c01 = values[0][0][1] * (1 - xd) + values[1][0][1] * xd
                        c10 = values[0][1][0] * (1 - xd) + values[1][1][0] * xd
                        c11 = values[0][1][1] * (1 - xd) + values[1][1][1] * xd
                        c0 = c00 * (1 - yd) + c10 * yd
                        c1 = c01 * (1 - yd) + c11 * yd
                        c = c0 * (1 - zd) + c1 * zd
                        return c
                    gxgradient = triInt(gXgradient)
                    gygradient = triInt(gYgradient)
                    gzgradient = triInt(gZgradient)
                    gxpixel = triInt(gXpixel)
                    gypixel = triInt(gYpixel)
                    gzpixel = triInt(gZpixel)

                    return gxpixel,gypixel,gzpixel,gxgradient,gygradient,gzgradient
                def noInterp():
                    x = int(v[0] + .5)
                    y = int(v[1] + .5)
                    z = int(v[2] + .5)

                    if x < 0 or y < 0 or z < 0 or x >= dataShape[0] or y >= dataShape[1] or z >= dataShape[2]:
                        gx = 9999  # the same as fun
                        gy = 9999
                        gz = 9999

                        print 'vert %s is outside of volume data in gradient E_I:' % k, x, y, z, 'vert coordinate:', v
                        return 0,0,0,0,0,0 # pixel and then gradient

                    hessian = np.array([hessianXx[x][y][z],
                                        hessianXy[x][y][z],
                                        hessianXz[x][y][z]])
                    gxgradient = omega2 * (normals[k].dot(hessian))
                    gxpixel = omega1 * gradientX[x][y][z]

                    hessian = np.array([hessianYx[x][y][z],
                                        hessianYy[x][y][z],
                                        hessianYz[x][y][z]])
                    gygradient = omega2 * (normals[k].dot(hessian))
                    gypixel = omega1 * gradientY[x][y][z]

                    hessian = np.array([hessianZx[x][y][z],
                                        hessianZy[x][y][z],
                                        hessianZz[x][y][z]])
                    gzgradient = omega2 * (normals[k].dot(hessian))
                    gzpixel = omega1 * gradientZ[x][y][z]
                    return gxpixel,gypixel,gzpixel,gxgradient,gygradient,gzgradient

                gxpixel,gypixel,gzpixel,gxgradient,gygradient,gzgradient = trilinearInterp()
                # gxpixel,gypixel,gzpixel,gxgradient,gygradient,gzgradient = noInterp()

                gx = signGradient * gxgradient - signPixel * gxpixel
                gy = signGradient * gygradient - signPixel * gypixel
                gz = signGradient * gzgradient - signPixel * gzpixel
                gx2 = gx
                gy2 = gy
                gz2 = gz

                # ff = np.sqrt(((v - np.array([50., 50., 50.])) ** 2).sum())
                # f = ff - 25.
                # gx = 2.0 * (v[0] - 50.) * f / ff
                # gy = 2.0 * (v[1] - 50.) * f / ff
                # gz = 2.0 * (v[2] - 50.) * f / ff
                # print 'gx,gy,gz',gx,gy,gz
                # print 'diff of gradient',gx - gx2,gy - gy2,gz - gz2

                g.append(gx)
                g.append(gy)
                g.append(gz)
                k += 1
            g = np.array(g).reshape((nv * 3,))
            return g

        checkGrad = False
        # check the gradient computation using optimize.check_grad
        if checkGrad:
            print 'check gradient for EI in SI optimization....'
            graderr = opt.check_grad(E_I,gradient_EI,S_I.verts.copy().reshape((nv * 3,)))
            print 'graderr for EI when optimizing SI is:',graderr

        refMesh = S_SSM  # original is S_B
        def E_I_B(SIVerts):
            SIVerts = SIVerts.reshape((nv,3))

            verts_SB = refMesh.verts[0:nv]
            verts_SI = SIVerts

            verts_BI = verts_SB - verts_SI
            # return np.sqrt(np.sum(verts_BI * verts_BI, axis=1))
            return np.sum(verts_BI * verts_BI) / 2.
        def gradient_EIB(SIVerts):
            SIVerts = SIVerts.reshape((nv,3))

            verts_SB = refMesh.verts[0:nv]
            verts_SI = SIVerts
            verts_BI = verts_SI - verts_SB
            return verts_BI.reshape((nv * 3,))

        if checkGrad:
            print 'check gradient for EIB in SI optimization....'
            graderr = opt.check_grad(E_I_B, gradient_EIB, S_I.verts.copy().reshape((nv * 3,)))
            print 'graderr for EIB when optimizing SI is:', graderr

        # weights really affect critically
        omega_EI = .5
        omega_BI = .5
        def fun(SIVerts):
            verts = SIVerts.reshape((nv,3))

            e_BI = omega_BI * E_I_B(verts)
            e_EI = omega_EI * E_I(verts)
            e = e_EI + e_BI
            return e
        def gradient(SIVerts):
            verts = SIVerts.reshape((nv, 3))

            g_BI = omega_BI * gradient_EIB(verts)
            g_EI = omega_EI * gradient_EI(verts)
            g = g_EI + g_BI
            g = g.reshape((nv * 3,))
            return g

        if checkGrad:
            print 'check gradient in SI optimization....'
            graderr = opt.check_grad(fun, gradient, S_I.verts.copy().reshape((nv * 3,)))
            print 'graderr when optimizing SI is:', graderr

        def callback(verts):
            verts = verts.reshape((nv,3))
            global km
            openmesh.writeOff(dirname + '/%s_registered_SI.off' % km, verts, mesh.faces)
            km += 1

        # l-bfgs-b,SLSQP, are OK. fastest is l-bfgs-b,
        #  with gradient (or even hessian) given, the computation efficiency
        # and accuracy are highly improved.
        # without any derivative information, we can choose: l-bfgs-b and SLSQP
        # with gradient, we can choose:CG,Newton-CG,l-bfgs-b,TNC,
        # with gradient and hessian, we can choose: dogleg and trust-cg (both not tested)

        res = opt.minimize(fun, refMesh.verts.reshape((nv * 3,)), args=(),method='l-bfgs-b',jac=gradient,
                           callback=callback,tol=1e-6,options={'maxiter': 100, 'disp': True,'gtol':1e-6})
        x = res['x'].reshape((nv,3))
        print 'x:',x
        # print 'mesh',mesh.verts
        print 'success?',res['success']
        print 'message:',res['message']

        S_I.verts = x.copy()

        return x

    def OptimizeSSM():
        '''
        repeat optimize transformation and SSM parameters until converge on S_SSM
        :return:
        '''
        phim = mesh.phim
        refMesh = S_I # original is S_B
        def transformSSM(a,R,c,b):
            verts_SSM = mesh.verts.copy() # mesh is the mean mesh.

            if not len(phim) == len(b):
                raise ValueError('Length of b(%s) must equal lenght of $\phi$(%s):', len(b), len(phim))

            k = 0
            for bm in b:
                verts_SSM += phim[k] * bm
                k += 1

            verts_SSM = verts_SSM.transpose()
            verts_SSM = R.dot(verts_SSM) * a + c.reshape(3,1)  # mesh is the mean mesh
            return np.array(verts_SSM.transpose())
        def E_SSM(a, R, c,b):  # translation, rotation, scale, and SSM parameters
            verts_ref = refMesh.verts
            verts_SSM = np.array(transformSSM(a,R,c,b))

            verts_BSSM = verts_SSM - verts_ref
            return np.sum(verts_BSSM * verts_ref) / 2.

        def fun_trans(vs,b=np.zeros(len(phim))):
            a = vs[0]
            R = np.matrix(vs[1:10].reshape((3,3)))
            c = np.array(vs[10:13])
            return E_SSM(a,R,c,b)
        def rotationConstraint(vs):
            R = np.matrix(vs[1:10].reshape((3, 3)))
            return R.transpose() * R - np.identity(3,dtype=np.float32)

        def fun_SSMP(b,a=0,R=np.matrix('1 0 0;0 1 0; 0 0 1'),c=np.array([0,0,0])):
            return E_SSM(a,R,c,b)

        nv = mesh.verts.shape[0]
        def callbackTrans(x):
            a = x[0]
            R = np.matrix(x[1:10].reshape((3, 3)))
            c = np.array(x[10:13])
            global km
            verts = transformSSM(a,R,c,np.zeros(len(phim)))
            openmesh.writeOff(dirname + '/%s_registered_SSM.off' % km, verts, mesh.faces)
            km += 1
        def callbackSSMP(b):
            pass

        b = np.ones(len(phim)) # initial guess
        a = 1.
        R = np.matrix('1. 0. 0.;0. 1. 0.;0. 0. 1.')
        c = np.array([0.,0.,0.],dtype=np.float32)
        verts_SSM = None
        for k in range(1): # only iterate 5 times since not confident about how to evaluate the convergence of a,R,c and b
            initvals = np.zeros((13,))
            initvals[0] = a
            initvals[1:10] = R.reshape((9,))
            initvals[10:13] = c
            res = opt.minimize(fun_trans,initvals, args=(b), method='slsqp', jac=False,
                               callback=callbackTrans, options={'maxiter': 100, 'disp': True})
            x = res['x']
            a = x[0]
            R = np.matrix(x[1:10].reshape((3,3)))
            c = np.array(x[10:13])
            print 'x:', x
            print 'success?', res['success']
            print 'message:', res['message']

            # res = opt.minimize(fun_SSMP, b, args=(a,R,c),method='l-bfgs-b', jac=False,
            #                    callback=callbackSSMP, options={'maxiter': 100, 'disp': True})
            # x = res['x']
            # b = x
            # print 'x:', x
            # print 'success?', res['success']
            # print 'message:', res['message']

            verts_SSM = transformSSM(a,R,c,b)
            # here, we can terminate the loop by checking whether the change of verts_SSM is small enough
        S_SSM.verts = verts_SSM.copy()

        return verts_SSM.copy()

    def OptimizeSSM_ICP():
        verts_SSM = mesh.verts.copy()  # mesh is the mean mesh.
        verts_ref = S_I.verts.copy()

        initial = np.array([[0.01], [0.05], [0.01], [0.001], [0.001], [0.001]])
        params = np.empty((6,))
        basicICP.icp_point_to_point_lm(verts_SSM,verts_ref,initial=initial,loop = 0,params = params)
        print 'params',params

        S_SSM.verts = mesh.verts.copy()
        # print 'verts before',S_SSM.verts
        S_SSM.verts = basicICP.transformpoints(S_SSM._verts,params=params)
        # print 'verts after', S_SSM.verts

        global km
        openmesh.writeOff(dirname + '/%s_registered_SSM.off' % km, S_SSM._verts, mesh.faces)
        km += 1

        return S_SSM.verts.copy()

    def OptimizeSB():
        gridDim = (5, 5, 10, 3)
        numcoords = gridDim[0] * gridDim[1] * gridDim[2] * gridDim[3]
        tbounds = np.array(volume.pixelArray.shape)
        nv = mesh.verts.shape[0]

        def E_I_B(SBVerts):
            verts_SB = SBVerts
            verts_SI = S_I.verts #

            verts_BI = verts_SB - verts_SI
            e = np.sum(verts_BI * verts_BI) / 2.

            return e

        def E_SSM(SBVerts):
            verts_SB = SBVerts
            verts_SSM = S_SSM.verts #

            verts_BSSM = verts_SB - verts_SSM
            e = np.sum(verts_BSSM * verts_BSSM) / 2.
            return e

        def E_SIM(SBVerts):
            verts_SB = SBVerts

            e = 0
            return e
        def TrivariateBSpline(tbounds,z,s,k=3):# z.ndim should be 4. the first 3 are the shape of the control vertices,
            # ,and shape[3]  should be 3 to contain the coordinate of the vertices. Note, if the number of voxel is x,
            # then the number of control points should be x + k - 1
            shape = z.shape[0:3] # the number of control points
            shape = np.array(shape)
            shape -= [k - 1,k - 1,k - 1] # the number of voxels in each dimension
            tx = np.zeros(shape=(shape[0] + 2 * k,))
            tx[k:-k] = np.linspace(0.,tbounds[0],shape[0])
            tx[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[0])
            # print 'tx:',tx
            ty = np.zeros(shape=(shape[1] + 2 * k,))
            ty[k:-k] = np.linspace(0.,tbounds[1],shape[1])
            ty[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[1])
            # print 'ty:',ty
            tz = np.zeros(shape=(shape[2] + 2 * k,))
            tz[k:-k] = np.linspace(0.,tbounds[2],shape[2])
            tz[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[2])
            # print 'tz:',tz
            spl = None
            shape = z.shape[0:3]
            vv = []
            s = s.transpose()
            for ii in range(shape[0]):
                v = []
                for jj in range(shape[1]):
                    spl = BSpline(tz,z[ii][jj],k)
                    val = spl(s[2])
                    v.append(val)
                v = np.array(v)
                spl = BSpline(ty,v,k)
                val = spl(s[1])
                vv.append(spl(s[1]))
            vv = np.array(vv)
            spl = BSpline(tx,vv,k)
            spl = spl(s[0])
            return spl
        # z = np.arange(0,5 * 100 * 3)
        # z = z.reshape((5,10,10,3))
        # tbounds = [1024,512,512]
        # spl = TrivariageBSpline(tbounds,z,[2,7,0])
        # print 'spline %s' % (spl)
        @fn_timer
        def genetateBasicFunc(tbounds,z,verts):
            z = z.reshape(gridDim)

            SBverts = verts.copy()
            k = 3
            shape = z.shape[0:3]  # the number of control points
            shape = np.array(shape)
            shape -= [k - 1, k - 1, k - 1]  # the number of voxels in each dimension
            tx = np.zeros(shape=(shape[0] + 2 * k,))
            tx[k:-k] = np.linspace(0., tbounds[0], shape[0])
            tx[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[0])
            # print 'tx:',tx
            ty = np.zeros(shape=(shape[1] + 2 * k,))
            ty[k:-k] = np.linspace(0., tbounds[1], shape[1])
            ty[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[1])
            # print 'ty:',ty
            tz = np.zeros(shape=(shape[2] + 2 * k,))
            tz[k:-k] = np.linspace(0., tbounds[2], shape[2])
            tz[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[2])

            shape = z.shape[0:3]  # the number of control points
            Bx = []
            for ii in range(shape[0]):
                c = np.zeros((shape[0],), dtype=np.float32)
                c[ii] = 1.
                spl = BSpline(tx, c, k)
                kp = 0
                bvals = np.empty((SBverts.shape[0],), dtype=np.float32)
                for v in SBverts:
                    bvals[kp] = spl(v[0])
                    kp += 1
                Bx.append(bvals)
            Bx = np.array(Bx)
            By = []
            for ii in range(shape[1]):
                c = np.zeros((shape[1],), dtype=np.float32)
                c[ii] = 1.
                spl = BSpline(ty, c, k)
                kp = 0
                bvals = np.empty((SBverts.shape[0],), dtype=np.float32)
                for v in SBverts:
                    bvals[kp] = spl(v[1])
                    kp += 1
                By.append(bvals)
            By = np.array(By)
            Bz = []
            for ii in range(shape[2]):
                c = np.zeros((shape[2],), dtype=np.float32)
                c[ii] = 1.
                spl = BSpline(tz, c, k)
                kp = 0
                bvals = np.empty((SBverts.shape[0],), dtype=np.float32)
                for v in SBverts:
                    bvals[kp] = spl(v[2])
                    kp += 1
                Bz.append(bvals)
            Bz = np.array(Bz)
            # print Bx.shape,By.shape,Bz.shape

            return Bx,By,Bz

        def updateSB(tbounds,z,verts):
            '''
            to compute s_b according to the equation (4). verts is s_0
            this is a time consuming operation if len(verts) is too large
            Note: this will change verts
            :param tbounds:
            :param z:
            :param verts:
            :return:
            '''
            Bx,By,Bz = genetateBasicFunc(tbounds,z,verts)

            vs = []
            kp = 0
            for v in verts:
                # v += TrivariateBSpline(tbounds,z,v)
                for ii in range(gridDim[0]):
                    for jj in range(gridDim[1]):
                        for kk in range(gridDim[2]):
                            v += z[ii][jj][kk] * Bx[ii][kp] * By[jj][kp] * Bz[kk][kp]
                kp += 1
                vs.append(v)
            return np.array(vs)

        omega_IB = .5
        omega_SSM = .5
        omega_SIM = .0

        def fun(z,tbounds,SBVerts):
            verts = SBVerts.reshape((nv, 3)).copy()
            z = z.reshape(gridDim)
            verts = updateSB(tbounds,z,verts)

            eIB = omega_IB * E_I_B(verts)
            eSSM = omega_SSM * E_SSM(verts)
            eSIM = omega_SIM * E_SIM(verts)
            e = eIB + eSSM + eSIM

            return e

        def G_I_B_n_SSM(z,tbounds, verts,refverts):
            SIverts = refverts.copy()

            Bx,By,Bz = genetateBasicFunc(tbounds,z,verts)

            verts = updateSB(tbounds, z, verts)  # S_B
            verts = verts - SIverts
            g = []
            for ii in range(gridDim[0]):
                for jj in range(gridDim[1]):
                    for kk in range(gridDim[2]):
                        gx = 0.
                        gy = 0.
                        gz = 0.
                        kp = 0
                        for v in verts:
                            kg = v * Bx[ii][kp] * By[jj][kp] * Bz[kk][kp]
                            gx += kg[0]
                            gy += kg[1]
                            gz += kg[2]
                            kp += 1

                        g.append(gx)
                        g.append(gy)
                        g.append(gz)
            return np.array(g).reshape((numcoords,))

        def G_SIM(z,tbounds,verts):
            return np.zeros(shape=(z.shape[0] * z.shape[1] * z.shape[2] * 3,))

        def gradient(z,tbounds,SBVerts):
            verts = SBVerts.reshape((nv, 3)).copy()
            z = z.reshape(gridDim)

            gIB = omega_IB * G_I_B_n_SSM(z,tbounds,verts.copy(),S_I.verts)
            gSSM = omega_SSM * G_I_B_n_SSM(z,tbounds,verts.copy(),S_SSM.verts)
            gSIM = omega_SIM * G_SIM(z,tbounds,verts.copy())
            g = gIB + gSSM + gSIM

            print g
            return g

        if True:
            print 'check gradient in SB optimization....'
            graderr = opt.check_grad(fun, gradient,np.ones(shape=(numcoords,)),*[tbounds,S_B.verts.copy()])
            print 'graderr when optimizing SB is:', graderr

        def callback(z,tbounds,verts):
            verts = verts.reshape((nv, 3)).copy()
            z = z.reshape(gridDim)
            verts = updateSB(tbounds, z, verts)
            global km
            openmesh.writeOff(dirname + '/%s_registered_SB.off' % km, verts, mesh.faces)
            km += 1

        res = opt.minimize(fun, np.zeros(shape = (numcoords,)),\
                           args=(tbounds,S_B.verts.copy()),method='l-bfgs-b',jac=gradient,
                           callback=callback,options={'maxiter': 100, 'disp': True})

        S_B.verts = updateSB(tbounds,res['x'].reshape(gridDim),S_B.verts)
        return S_B.verts.copy()

    #repeat optimization until converge
    for k in range(10):
        SIverts = None
        SSMverts = None
        SBverts = None

        # print '------------> Optimizing SI <-------------'
        # SIverts = OptimizeSI() # will update S_I.verts
        # print '------------> Optimizing SSM <-------------'
        # SSMverts = OptimizeSSM_ICP() # will update S_SSM.verts
        print '------------> Optimizing SB <-------------'
        SBverts = OptimizeSB() # will update S_B.verts

    return SIverts, SSMverts,SBverts

def main():
    ''' Segment a 3D spine Dicom image with pre-trained mesh model.'''

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
    shapeForComp = dicomImage.vol_data_forComp.shape
    spacing = [dicomImage.spacingX,dicomImage.spacingY,dicomImage.spacingSlice]
    firstPos = dicomImage.firstLocation
    # print 'firstpos',firstPos
    imagePos = dicomImage.imagePosition
    # print 'impos',imagePos
    lastPos = imagePos[2]
    # print 'lastpos',lastPos
    width = (shapeForComp[0] - 1.) * spacing[0]
    height = (shapeForComp[1] - 1.) * spacing[1]
    thickness = np.abs(firstPos - lastPos)
    # print width,height,thickness
    trans = [width - imagePos[0],height - imagePos[1],-lastPos]
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
    initMeshInDicom(meanMesh,dicomImage)

    # Third, Optimize the mesh

    S = OptimizeMesh(meanMesh,dicomImage)
    # OptimizeMesh(meanMesh,Volume())

    # S.show()
    itrans = [-trans[0],-trans[1],-trans[2]]
    iscale = [1./scale[0],1./scale[1],1./scale[2]]
    if S[0].all():
        transformmesh.transform(S[0],itrans,iscale)
        openmesh.writeOff('./Data/registered_reg1_SI.off', S[0], meanMesh.faces)
    if S[1].all():
        transformmesh.transform(S[1], itrans, iscale)
        openmesh.writeOff('./Data/registered_reg1_SSM.off', S[1], meanMesh.faces)
    if S[2].all():
        # pass
        transformmesh.transform(S[2], itrans, iscale)
        openmesh.writeOff('./Data/registered_reg1_SB.off', S[2], meanMesh.faces)

    # vispy.app.run()

if __name__ == "__main__":
    main()
