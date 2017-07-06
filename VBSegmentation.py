#!/usr/bin/env python
# encoding: utf-8

'''
input 3D Dicom images containing spine, this code will segmente all the vertebra bodies
'''
import sys,os,shutil,getopt
from PIL import Image
import math
import numpy as np
import scipy
from scipy import optimize as opt,io,linalg,ndimage,misc,sparse
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

import time
from functools import wraps
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
            print ("Total time running %s: %s seconds" %
                (function.func_name, str(totaltime))
                )
        elif totaltime < 3600:
            print ("Total time running %s: %s minutes" %
               (function.func_name, str(totaltime / 60.))
               )
        elif totaltime < 3600 * 24:
            print ("Total time running %s: %s hours and %s minutes" %
               (function.func_name, str(totaltime / 3600),str((totaltime % 3600) / 60.))
               )
        else:
            print ("Total time running %s: %s days,%s hours, and %s minutes" %
               (function.func_name, str(totaltime / (3600 * 24)),str((totaltime % (3600 * 24)) / 3600),str((totaltime % 3600) / 60.))
               )
        return result
    return function_timer

class Mesh(object):
    ''' Class to represet a 3D triangular mesh. nv and nf are the number of vertices and faces of the mesh, verts and
    faces are the numpy.array objects to store the coordinates of the vertices and the indices of the faces.'''
    def __init__(self):
        self._nv = 0
        self._nf = 0
        self._verts = np.empty(shape = (self._nv,3),dtype=np.float32)
        self._faces = np.empty(shape = (self._nf,3),dtype=np.int32)

        # the centers of the low and upper endplates of the vb mesh
        self._centerLowerEndplate = np.zeros((1,3),dtype=np.float32)
        self._centerUpperEndplate = np.zeros((1,3),dtype=np.float32)

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

    @property
    def faces(self):
        return self._faces
    @faces.setter
    def faces(self,faces):
        if not faces.shape[1] == 3:
            return ValueError('Invalid input faces: it must be of shape [x][3]')
        self._faces = faces
        self._nf = faces.shape[0]

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

    def show(self,withoutline=True):
        '''Show current mesh and the centers of lower and upper endplate in current view'''
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
        view.add(mesh)

        if withoutline:
            GLContext.set_polygon_offset(1, 1)
            mesh.update_gl_state(**{'polygon_offset_fill': True})

            outline_p = []
            for face in faces:
                outline_p.append(p[face[0]])
                outline_p.append(p[face[1]])
                outline_p.append(p[face[1]])
                outline_p.append(p[face[2]])
                outline_p.append(p[face[2]])
                outline_p.append(p[face[0]])
            outline_p = np.array(outline_p)

            outline = visuals.Line(pos=outline_p, color=(0.1, 0.1, 0.1, 1.), width=1,
                 connect='segments', method='gl', antialias=False)
            view.add(outline)

        # cube = visuals.Cube(edge_color = [1.0,0,0,1])
        # view.add(cube)

        view.camera.set_range() # similar to viewAll in Coin3D
        # vispy.app.run()

class Volume(object):
    ''' Class to represent a 3D Dicom object. numX, numY and numZ are the number of pixels in each directions.
    pixel_array stores the pixel value of the Dicom object slice by slice. Note that the Z direction is Axis for
    CT data and spinal for MR image.'''
    def __init__(self):

        self._spacingX = 0.
        self._spacingY = 0.
        self._spacingSlice = 0. # the spacing along each direction

        self._pixel_array = None # the copy of pixel_array of dicom files read by pydicom

    @property
    def spacingX(self):
        return self._spacingX
    @spacingX.setter
    def spacingX(self,s):
        self._spacingX = s
    @property
    def spacingY(self):
        return self._spacingY
    @spacingY.setter
    def spacingY(self,s):
        self._spacingY = s
    @property
    def spacingSlice(self):
        return self._spacingSlice
    @spacingSlice.setter
    def spacingSlice(self,s):
        self._spacingSlice = s

    @property
    def pixelArray(self):
        return self._pixel_array

    # the first dim must be the slice, while the second and third dim are the x and y pixels
    @property
    def dataShape(self):
        return self._pixel_array.shape
    @dataShape.setter
    def dataShape(self,shape):
        self._pixel_array = np.zeros(shape=shape)

    def setPixelArrayData(self,iSlice,data):
        self._pixel_array[iSlice] = data
    def setPixelArrayData(self,voldata):
        self._pixel_array = voldata

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
        volume = visuals.Volume(vol_data, clim=clim,cmap = cmap,method = 'mip',relative_step_size=0.8)


        view.add(volume)

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
    for dcm in allDicoms:
        filename = dicomDir + '/%s' % dcm + '.dcm'
        if not os.path.isfile(filename):
            print filename,'is not a file'
        ds = dicom.read_file(filename)

        print '='*20
        print dcm
        if not ds.SeriesNumber == 5:
            # print 'series number',ds.SeriesNumber
            print 'skip dicom file: %s' % dcm
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
        print 'image position:',ds.ImagePositionPatient

        k = k + 1
        if k == 0:
            image.spacingX = ds.PixelSpacing[0]
            image.spacingY = ds.PixelSpacing[1]  # unit: mm
        vol_data.append(ds.pixel_array)

    vol_data = np.array(vol_data)
    # vol_data = np.flipud(np.rollaxis(vol_data, 1))
    image.setPixelArrayData(vol_data)
    print '='*20
    print '%s slice loaded' % k
    print 'shape of image:',image.dataShape



    return image

@fn_timer
def initMeshInDicom(mesh,volume):
    ''' Put the mesh to a proper position in the 3D Dicom image. Coordinates of vertices are modified.'''
    pass

@fn_timer
def OptimizeMesh(mesh,volume):
    ''' Optimize the mesh.'''
    pass

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

    vispy.app.run()

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
    colors[1][3] = 4./255
    colors[2][3] = 0.
    colors[3][3] = 1./255
    colors[4][3] = 1.
    colors[5][3] = 1./255

    cmap = Colormap(colors=colors)



    return cmap

def main():
    ''' Segment a 3D spine Dicom image with pre-trained mesh model.'''



    # First, load the mean model and the 3D Dicom data
    meshinfo = {'filename':'./Data/reg1.off','centerLowerEndplate':[0,0,0],'centerUpperEndplate':[0,0,1]}
    meanMesh = loadMesh('meanmesh',meshinfo)
    dicominfo = {'dirname':'/Users/mac/Downloads/Dicom Data/2585255-profYu','centerLowerEndplate':[0,0,0],'centerUpperEndplate':[0,0,1]}
    dicomImage = loadDICOM(dicominfo)

    # show the mesh and 3D image in the same window
    meanMesh.show()
    ct_aaa = getColorMaps('ct_aaa')
    dicomImage.show(cmap = ct_aaa)

    # Second, initialize the mesh to the approximate position in the image
    initMeshInDicom(meanMesh,dicomImage)

    # Third, Optimize the mesh
    OptimizeMesh(meanMesh,dicomImage)

    vispy.app.run()

if __name__ == "__main__":
    main()
