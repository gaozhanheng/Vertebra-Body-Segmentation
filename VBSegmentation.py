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
        p = self._verts
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
        # print 'gradient',gradientArray[0],gradientArray[1],gradientArray[2]
        print 'g',gradientArray

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

        vol_data = self._pixel_array
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
        phim.append(np.ones(shape=(mesh.nv * 3)).reshape(mesh.nv,3))
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
        k = k + 1
        vol_data.append(np.flip(np.flip(ds.pixel_array,axis=0),axis=1).astype(np.float32))

    vol_data = np.flip(vol_data, axis=0)
    vol_data = ndimage.gaussian_filter(np.array(vol_data),sigma=3,order=1)
    vol_data_forComp = np.rollaxis(np.rollaxis(vol_data.copy(),1),2)
    # vol_data = vol_data_forComp

    print 'shape of vol_data:',vol_data.shape
    print 'shape of vol_data for computation:',vol_data_forComp.shape

    # for x in range(vol_data.shape[0]):
    #     for y in range(vol_data.shape[1]):
    #         for z in range(vol_data.shape[2]):
    #             if not vol_data[x][y][z] == vol_data_forComp[z][y][x]:
    #                 print 'roll axis error:',x,y,z


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

    def OptimizeSI(): # more details about scipy.optimization please refer to the optimize __init__.py
        nv = S_I.nv # the number of vertices used to test the optimization efficiency
        # print 'number of vertices used:',nv

        meshData = vispy.geometry.MeshData(vertices=S_B.verts, faces=mesh.faces)  # the faces are all the same
        normals = meshData.get_vertex_normals()


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

        omega1 = .1
        omega2 = .9
        def E_I(SIVerts):
            SIVerts = SIVerts.reshape((nv,3))

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
                x = int(v[0])
                y = int(v[1])
                z = int(v[2])
                if x < 0 or y < 0 or z < 0 or x >= dataShape[0] or y >= dataShape[1] or z >= dataShape[2]:
                    ei += 9999
                    print 'outside of volume data in E_I:',x,y,z
                else:
                    vpixel = omega1 * pixelarray[x][y][z]
                    gradient = np.array([gradientX[x][y][z],gradientY[x][y][z],gradientZ[x][y][z]])
                    # normalizeVec(gradient)
                    vgradient = omega2 * (normals[k].dot(gradient))

                    v = vgradient  + vpixel
                    ei += v
                k += 1

            # print ei
            return ei
        def gradient_EI(SIVerts):
            SIVerts = SIVerts.reshape((nv,3))

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
                x = int(v[0])
                y = int(v[1])
                z = int(v[2])
                if x < 0 or y < 0 or z < 0 or x >= dataShape[0] or y >= dataShape[1] or z >= dataShape[2]:
                    gx = 9999
                    gy = 9999
                    gz = 9999

                    print 'outside of volume data in gradient_EI:', x, y, z
                else:
                    gxpixel = omega1 * gradientX[x][y][z]
                    hessian = np.array([hessianXx[x][y][z], hessianXy[x][y][z],hessianXz[x][y][z]])
                    # normalizeVec(hessian)
                    gxgradient = omega2 * (normals[k].dot(hessian))
                    gx = gxgradient + gxpixel

                    gypixel = omega1 * gradientY[x][y][z]
                    hessian = np.array([hessianYx[x][y][z], hessianYy[x][y][z],hessianYz[x][y][z]])
                    # normalizeVec(hessian)
                    gygradient = omega2 * (normals[k].dot(hessian))
                    gy = gygradient  + gypixel

                    gzpixel = omega1 * gradientZ[x][y][z]
                    hessian = np.array([hessianZx[x][y][z], hessianZy[x][y][z],hessianZz[x][y][z]])
                    # normalizeVec(hessian)
                    gzgradient = omega2 * (normals[k].dot(hessian))
                    gz = gzgradient  + gzpixel

                g.append(gx)
                g.append(gy)
                g.append(gz)
                k += 1

            return np.array(g).reshape((nv * 3,))

        checkGrad = False
        # check the gradient computation using optimize.check_grad
        if checkGrad:
            print 'check gradient for EI in SI optimization....'
            graderr = opt.check_grad(E_I,gradient_EI,S_I.verts.copy().reshape((nv * 3,)))
            print 'graderr for EI when optimizing SI is:',graderr

        def E_I_B(SIVerts):
            SIVerts = SIVerts.reshape((nv,3))

            verts_SB = S_B.verts[0:nv]
            verts_SI = SIVerts

            verts_BI = verts_SB - verts_SI
            # return np.sqrt(np.sum(verts_BI * verts_BI, axis=1))
            return np.sum(verts_BI * verts_BI) / 2.
        def gradient_EIB(SIVerts):
            SIVerts = SIVerts.reshape((nv,3))

            verts_SB = S_B.verts[0:nv]
            verts_SI = SIVerts
            verts_BI = verts_SI - verts_SB
            return verts_BI.reshape((nv * 3,))

        if checkGrad:
            print 'check gradient for EIB in SI optimization....'
            graderr = opt.check_grad(E_I_B, gradient_EIB, S_I.verts.copy().reshape((nv * 3,)))
            print 'graderr for EIB when optimizing SI is:', graderr

        omega_EI = .8
        omega_BI = .2
        def fun(SIVerts):
            verts = SIVerts.reshape((nv,3))

            e_BI = omega_BI * E_I_B(verts)
            e_EI = omega_EI * E_I(verts)
            # print 'BI:',e_BI
            # print 'EI:',e_EI
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

        def callback(xk):
            verts = xk - S_B.verts[0:nv].reshape((nv * 3,))
            verts = verts * verts

        # res = opt.minimize(fun, S_B.verts.reshape((nv * 3,)), method='l-bfgs-b',jac=gradient,
        #                    callback=callback,options={'maxiter': 100, 'disp': True})
        # l-bfgs-b,SLSQP, are OK. fastest is l-bfgs-b,
        #  with gradient (or even hessian) given, the computation efficiency
        # and accuracy are highly improved.
        # without any derivative information, we can choose: l-bfgs-b and SLSQP
        # with gradient, we can choose:CG,Newton-CG,l-bfgs-b,TNC,
        # with gradient and hessian, we can choose: dogleg and trust-cg (both not tested)
        res = opt.minimize(fun, S_B.verts.reshape((nv * 3,)), args=(),method='l-bfgs-b',jac=gradient,
                           callback=callback,options={'maxiter': 100, 'disp': True})
        x = res['x'].reshape((nv,3))
        print 'x:',x
        print mesh.verts
        print 'success?',res['success']
        print 'message:',res['message']

        S_I.verts = x

        return x

    def OptimizeSSM():
        '''
        repeat optimize transformation and SSM parameters until converge on S_SSM
        :return:
        '''
        phim = mesh.phim
        def updateSSM(a,R,c,b):
            verts_SSM = mesh.verts.copy() # mesh is the mean mesh

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
            verts_SB = S_B.verts
            verts_SSM = np.array(updateSSM(a,R,c,b))

            verts_BSSM = verts_SSM - verts_SB
            return np.sum(verts_BSSM * verts_BSSM) / 2.

        def fun_trans(vs,b=np.zeros(len(phim))):
            a = vs[0]
            R = np.matrix(vs[1:10].reshape((3,3)))
            c = np.array(vs[10:13])
            return E_SSM(a,R,c,b)
        def fun_SSMP(b,a=0,R=np.matrix('1 0 0;0 1 0; 0 0 1'),c=np.array([0,0,0])):
            return E_SSM(a,R,c,b)

        b = np.zeros(len(phim))
        verts_SSM = None
        numiter = 1
        for k in range(1): # only iterate 5 times since not confident about how to evaluate the convergence of a,R,c and b
            res = opt.minimize(fun_trans, np.ones(shape=(13,)), args=(b), method='l-bfgs-b', jac=False,
                               callback=None, options={'maxiter': 100, 'disp': True})
            x = res['x']
            a = x[0]
            R = np.matrix(x[1:10].reshape((3,3)))
            c = np.array(x[10:13])
            print 'x:', x
            print 'success?', res['success']
            print 'message:', res['message']

            res = opt.minimize(fun_SSMP, np.ones(shape=(len(phim),)), args=(a,R,c),method='l-bfgs-b', jac=False,
                               callback=None, options={'maxiter': 100, 'disp': True})
            x = res['x']
            b = x
            print 'x:', x
            print 'success?', res['success']
            print 'message:', res['message']

            verts_SSM = updateSSM(a,R,c,b)
            # here, we can terminate the loop by checking whether the change of verts_SSM is small enough
        S_SSM.verts = verts_SSM

        return verts_SSM

    def OptimizeSB():
        nv = S_B.verts.shape[0]
        def E_I_B(SBVerts):
            verts_SB = SBVerts
            verts_SI = S_I.verts

            verts_BI = verts_SB - verts_SI
            e = np.sum(verts_BI * verts_BI) / 2.

            return e

        def E_SSM(SBVerts):
            verts_SB = SBVerts
            verts_SSM = S_SSM.verts

            verts_BSSM = verts_SB - verts_SSM
            e = np.sum(verts_BSSM * verts_BSSM) / 2.
            return e

        def E_SIM(SBVerts):
            verts_SB = SBVerts

            e = 0
            return e
        def TrivariageBSpline(tbounds,z,s,k=3):# z.ndim should be 4. the first 3 are the shape of the control vertices,
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
        def updateSB(tbounds,z,verts):
            '''
            this is a time consuming operation if len(verts) is too large
            Note: this will change verts
            :param tbounds:
            :param z:
            :param verts:
            :return:
            '''
            for v in verts:
                v += TrivariageBSpline(tbounds,z,v)

        gridDim = (5,5,10,3)
        def fun(z,tbounds,SBVerts):
            verts = SBVerts.reshape((nv, 3)).copy()
            z = z.reshape(gridDim)
            updateSB(tbounds,z,verts)
            omega_IB = 1.
            omega_SSM = 1.
            omega_SIM = 1.
            eIB = omega_IB * E_I_B(verts)
            eSSM = omega_SSM * E_SSM(verts)
            eSIM = omega_SIM * E_SIM(verts)
            e = eIB + eSSM + eSIM

            return e

        def gradient(z,tbounds,SBVerts):
            pass
        tbounds = np.array(volume.pixelArray.shape)
        res = opt.minimize(fun, np.ones(shape = gridDim), args=(tbounds,S_B.verts),method='l-bfgs-b',jac=False,
                           callback=None,options={'maxiter': 5, 'disp': True})

        updateSB(tbounds,res['x'],S_B.verts)
        return S_B.verts

    #repeat optimization until converge
    for k in range(1):
        verts = OptimizeSI()
        # verts = OptimizeSSM()
        # verts = OptimizeSB()

        S_B.verts = verts

    return S_B

def main():
    ''' Segment a 3D spine Dicom image with pre-trained mesh model.'''

    # test transformmesh
    # verts = np.array([[2,1,0]],dtype=np.float32)
    # transformmesh.transform(verts,np.array([-1,-2,-3]),np.array([2,3,4]))
    # transformmesh.transform(verts,np.array([1,2,3]),np.array([1./2,1./3,1./4]),inverse=True)
    # return


    # First, load the mean model and the 3D Dicom data
    meshinfo = {'filename':'./Data/transformed_reg1_simplifiedTo2000.off','centerLowerEndplate':[0,0,0],'centerUpperEndplate':[0,0,1]}
    meanMesh = loadMesh('meanmesh',meshinfo)
    dicominfo = {'dirname':'/Users/mac/Downloads/Dicom Data/2585255-profYu','centerLowerEndplate':[0,0,0],'centerUpperEndplate':[0,0,1]}
    dicomImage = loadDICOM(dicominfo)

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
    meanMesh.show()
    ct_aaa = getColorMaps('ct_aaa')
    dicomImage.show(cmap = ct_aaa)

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

    v, f = openmesh.readOff('./Data/sphere.off', quiet=True)
    sphere = visuals.Mesh(vertices=np.array(v),faces=np.array(f),color=np.array([[1.0, 1.0, 1.0,0.9]]),shading='flat')
    view.add(sphere)
    v, f = openmesh.readOff('./Data/x.off', quiet=True)
    x = visuals.Mesh(vertices=np.array(v), faces=np.array(f), color=np.array([[1.0, 0.0, 0.0, 0.9]]),
                          shading='flat')
    view.add(x)
    v, f = openmesh.readOff('./Data/y.off', quiet=True)
    y = visuals.Mesh(vertices=np.array(v), faces=np.array(f), color=np.array([[1.0, 1.0, 0.0, 0.9]]),
                          shading='flat')
    view.add(y)
    v, f = openmesh.readOff('./Data/z.off', quiet=True)
    z = visuals.Mesh(vertices=np.array(v), faces=np.array(f), color=np.array([[0.0, 0.0, 1.0, 0.9]]),
                          shading='flat')
    view.add(z)
    view.camera.set_range()


    # view.add(visuals.Cube())
    meanMesh.hide()

    # Second, initialize the mesh to the approximate position in the image
    initMeshInDicom(meanMesh,dicomImage)

    # Third, Optimize the mesh
    S = OptimizeMesh(meanMesh,dicomImage)
    # OptimizeMesh(meanMesh,Volume())

    S.show()
    itrans = [-trans[0],-trans[1],-trans[2]]
    iscale = [1./scale[0],1./scale[1],1./scale[2]]
    transformmesh.transform(S.verts,itrans,iscale)
    openmesh.writeOff('./Data/registered_reg1.off', S.verts, S.faces)

    # vispy.app.run()

if __name__ == "__main__":
    main()
