import time
from functools import wraps

import numpy as np

from scipy import optimize as opt,io,linalg,ndimage,misc,sparse,interpolate

import copy

from vispy.scene import visuals

import sys,os,shutil,getopt,time,random

import meshio as openmesh
import dicom # this is the module name of pydicom

import sys,os,shutil,getopt,time,random
# import vtk
from PIL import Image
import math

import scipy
from scipy.interpolate import UnivariateSpline,BivariateSpline,BSpline
from scipy.sparse import linalg as splinalg
import skimage
from skimage import color,io,data,segmentation,img_as_ubyte,img_as_float,feature,\
    filters,measure,graph,transform,draw,external,restoration
import matplotlib.pyplot as plt

import vispy.scene
from vispy.color import *
from vispy import io,gloo
from vispy.gloo import wrappers


from icp import basicICP

import re

from utilities import *



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

def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        totaltime = t1 - t0
        if totaltime < 60:
            print ("Total time running \'%s\': %s seconds" %
                (function.func_name, str(totaltime))
                )
        elif totaltime < 3600:
            print ("Total time running \'%s\': %s minutes and %s seconds" %
               (function.func_name, str(int(totaltime) / 60), str((totaltime % 60)))
               )
        elif totaltime < 3600 * 24:
            print ("Total time running \'%': %s hours and %s minutes" %
               (function.func_name, str(int(totaltime) / 3600),str((totaltime % 3600) / 60.))
               )
        else:
            print ("Total time running \'%s\': %s days,%s hours, and %s minutes" %
               (function.func_name, str(int(totaltime) / (3600 * 24)),str(int(totaltime % (3600 * 24)) / 3600)
                ,str((totaltime % 3600) / 60.))
               )
        return result
    return function_timer

def initMeshInDicom(mesh,volume):
    ''' Put the mesh to a proper position in the 3D Dicom image. Coordinates of vertices are modified.'''
    pass



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
        # self._computeDerivatives()
    @property
    def spacingY(self):
        return self._spacingY
    @spacingY.setter
    def spacingY(self,s):
        self._spacingY = s
        # self._computeDerivatives()
    @property
    def spacingSlice(self):
        return self._spacingSlice
    @spacingSlice.setter
    def spacingSlice(self,s):
        self._spacingSlice = s
        # self._computeDerivatives()

    @property
    def pixelArray(self):
        return self._pixel_array
    @pixelArray.setter
    def pixelArray(self,voldata):
        self._pixel_array = voldata
        self._computeDerivatives()

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

    def _computeDerivatives(self):
        '''
        First, compute the gradient of the pixel volume
        :return:
        '''
        return

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
            # print 'ds0'
            # print ds
        if k == 1:
            # print 'ds1'
            # print ds
            z1 = ds.SliceLocation

            image.spacingSlice = np.abs(image.firstLocation - z1)

        # if k == 14:
        #     print 'ds14',ds.pixel_array
        k = k + 1
        vol_data.append(np.flip(np.flip(ds.pixel_array,axis=0),axis=1).astype(np.float32))

    vol_data = np.flip(vol_data, axis=0)
    # print 'vol data',vol_data
    # vol_data = ndimage.median_filter(vol_data,size=3)
    vol_data = ndimage.gaussian_filter(np.array(vol_data),sigma=2,order=0)
    vol_data_forComp = np.rollaxis(np.rollaxis(vol_data,axis=1),axis=2)# we don't change vol_data_forComp so no need to copy vol_data

    print 'shape of vol_data:',vol_data.shape
    print 'shape of vol_data for computation:',vol_data_forComp.shape


    if False:
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