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
from vispy import io
from vispy.scene import visuals
# from vispy import visuals

try:
    import openmesh
except:
    print 'no openmesh module. using meshio to save mesh (only support .off files)'
    import meshio as openmesh

import time
from functools import wraps
import re


#
    # Make a canvas and add simple view
    #
canvas = vispy.scene.SceneCanvas(keys='interactive', show=True)
view = canvas.central_widget.add_view()

view.camera = 'turntable'  # or try 'arcball' . input an invalid and see the error for more options

# add a colored 3D axis for orientation
axis = visuals.XYZAxis(parent=view.scene)
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
        self._verts = np.zeros(shape = (self._nv,3),dtype=np.float32)
        self._faces = np.zeros(shape = (self._nf,3),dtype=np.int32)

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

    def show(self):
        '''Show current mesh and the centers of lower and upper endplate in current view'''
        # add mesh data into the scene
        # Vertices positions
        p = self._verts
        faces = self._faces

        # for more deatils on color, please refer to vispy/colors/colorarray.py
        # colors_p = np.array(
        #     [[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1],
        #      [1, 1, 1, 1]])
        # colors_f = np.array([(1.0, 0.1, 0.0, 0.3) for k in range(faces.shape[0])])
        color_mesh = np.array([[1.0, 0.1, 0.9, 0.5]])
        # for more details see vispy/visuals/mesh.py
        mesh = visuals.Mesh(vertices=p, faces=faces,
                            # vertex_colors = colors_p,
                            # face_colors = colors_f,
                            color=color_mesh,
                            shading='flat',
                            mode='triangles')
        view.add(mesh)

        vispy.app.run()

class Volume(object):
    ''' Class to represent a 3D Dicom object. numX, numY and numZ are the number of pixels in each directions.
    pixel_array stores the pixel value of the Dicom object slice by slice. Note that the Z direction is Axis for
    CT data and spinal for MR image.'''
    def __init__(self):
        self._numX = 0
        self._numY = 0 # numX and numY are the resolution in each slice
        self._numSlice = 0 # numZ is the number of slices
        self._spacingX = 0.
        self._spacingY = 0.
        self._spacingSlice = 0. # the spacing along each direction

        self._pixel_array = None # the copy of pixel_array of dicom files read by pydicom

    @property
    def numX(self):
        return self._numX
    @property
    def numY(self):
        return self._numY
    @property
    def numSlice(self):
        return self._numSlice
    @property
    def spacingX(self):
        return self._spacingX
    @property
    def spacingY(self):
        return self._spacingY
    @property
    def spacingSlice(self):
        return self._spacingSlice

    @property
    def pixelArray(self):
        return self._pixel_array

    def show(self):
        # add a volume into the scene
        vol_data = np.load(io.load_data_file('brain/mri.npz'))['data']
        print type(vol_data)
        print vol_data.shape
        vol_data = np.flipud(np.rollaxis(vol_data, 1))
        clim = [32, 192]
        volume = visuals.Volume(vol_data, clim=clim)

        view.add(volume)

        vispy.app.run()

@fn_timer
def loadMesh(meshName,meshInfo):
    ''' Load a mesh from meshName and return a Mesh object'''
    mesh = Mesh()

    return mesh

@fn_timer
def loadDICOM(dicomDir):
    ''' Load all dicom files from dicomDir. Integrate and return the data with a whole Dicom3D object.'''
    image = Volume()
    return image

def showVolume(volume):
    pass

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

@fn_timer
def main():
    ''' Segment a 3D spine Dicom image with pre-trained mesh model.'''

    # First, load the mean model and the 3D Dicom data
    meanMesh = loadMesh('meanmesh','meshinfo')
    dicomImage = loadDICOM('dicomdir')

    print meanMesh
    print dicomImage

    # show the mesh and 3D image in the same window
    meanMesh.show()

    # Second, initialize the mesh to the approximate position in the image
    initMeshInDicom(meanMesh,dicomImage)

    # Third, Optimize the mesh
    OptimizeMesh(meanMesh,dicomImage)


if __name__ == "__main__":
    main()
