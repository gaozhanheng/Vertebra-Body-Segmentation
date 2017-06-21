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
try:
    import openmesh
except:
    print 'no openmesh module. using meshio to save mesh (only support .off files)'
    import meshio as openmesh

import time
from functools import wraps
import re




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

def loadImg(filename):
    '''
    load the image using scipy
    :param filename: image name
    :return: numpy.ndarray image
    '''
    image = misc.imread(filename)
    return image

def generateLN(image,image_alpha,show=False,asFastAsPossible = True):
    '''
    process the image: generate the Laplacian of the input image.
    :param image: the original image.must be rgba
    :return: processed image: one is generated using method metioned in 'normal map bas relief' paper, and the other
    using scipy.ndimage.laplace function
    '''
    if not asFastAsPossible:
        print image.dtype,image.shape,image[:,:,0:3].min(),image[:,:,0:3].max()

    image_x = image[:,:,0]
    image_y = image[:,:,1]
    image_z = image[:,:,2] # !! this is true
    # image_alpha = image[:,:,3]
    if not asFastAsPossible:
        print 'image_x',image_x.dtype, image_x.shape, image_x.min(), image_x.max()
        print 'image_y',image_y.dtype, image_y.shape, image_y.min(), image_y.max()
        print 'image_z',image_z.dtype, image_z.shape, image_z.min(), image_z.max()
        print 'image_alpha',image_alpha.dtype, image_alpha.shape, image_alpha.min(), image_alpha.max()

    if not asFastAsPossible and show:
        pl.figure()
        pl.subplot(221)
        pl.imshow(image_x, pl.cm.gray)
        pl.title('normal x')

        pl.subplot(222)
        pl.imshow(image_y, pl.cm.gray)
        pl.title('normal y')

        pl.subplot(223)
        pl.imshow(image_z, pl.cm.gray)
        pl.title('normal z')

        pl.subplot(224)
        pl.imshow(image_alpha, pl.cm.gray)
        pl.title('alpha channel')
        pl.show()

    image_Nx_Nz = image_x / image_z * image_alpha
    image_Ny_Nz = image_y / image_z * image_alpha
    if not asFastAsPossible:
        print 'image_Nx_Nz:',image_Nx_Nz.dtype,image_Nx_Nz.shape,image_Nx_Nz.min(),image_Nx_Nz.max()
        print 'image_Ny_Nz:',image_Ny_Nz.dtype, image_Ny_Nz.shape, image_Ny_Nz.min(), image_Ny_Nz.max()

    # pl.figure()
    # pl.subplot(121)
    # pl.plot(image_Nx_Nz[250,:])
    # pl.title('one horizon line of Nx_Nz')
    # pl.subplot(122)
    # pl.plot(image_Nx_Nz[:,300])
    # pl.title('one vertical line of Nx_Mz')
    # pl.show(**{'block':False})

    if not asFastAsPossible and show:
        pl.figure()
        pl.subplot(121)
        pl.imshow(image_Nx_Nz, pl.cm.gray)
        pl.title('normal Nx_Nz')

        pl.subplot(122)
        pl.imshow(image_Ny_Nz, pl.cm.gray)
        pl.title('normal Ny_Nz')
        pl.show()

    # print 'number of valid pixels before erosion:',scipy.count_nonzero(image_alpha)
    # image_alpha = scipy.ndimage.binary_erosion(image_alpha,iterations=1)
    # print 'number of valid pixels after erosion:',scipy.count_nonzero(image_alpha)
    # pl.figure()
    # pl.imshow(image_alpha,pl.cm.gray)
    # pl.show()

    if True: # method 1 for computing the gradient of an image
        img_derivative_x = scipy.gradient(image_Nx_Nz,axis=1)
        img_derivative_y = -scipy.gradient(image_Ny_Nz,axis=0) # this is right, but why??
    else: #method 2 for computing the gradient of an image
        sigma = 1  # 标准差, different visual result with different value. try use 5 and enjoy
        img_derivative_x = np.zeros(image_Nx_Nz.shape)
        ndimage.gaussian_filter(image_Nx_Nz, (sigma, sigma), (0, 1), img_derivative_x)
        img_derivative_y = np.zeros(image_Ny_Nz.shape)
        ndimage.gaussian_filter(image_Ny_Nz, (sigma, sigma), (1, 0), img_derivative_y)
    LN = (img_derivative_x + img_derivative_y) * image_alpha
    if not asFastAsPossible:
        print 'LN:',LN.dtype,LN.min(),LN.max()

    # pl.figure()
    # pl.subplot(121)
    # pl.plot(img_derivative_x[250, :])
    # pl.title('one horizon line of img_derivative_x')
    # pl.subplot(122)
    # pl.plot(img_derivative_x[:, 300])
    # pl.title('one vertical line of img_derivative_x')
    # pl.show(**{'block': False})

    if not asFastAsPossible and show:
        pl.figure()
        pl.subplot(121)
        pl.imshow(img_derivative_x, pl.cm.gray)
        pl.title('normal Nx_Nz_x')

        pl.subplot(122)
        pl.imshow(img_derivative_y, pl.cm.gray)
        pl.title('normal Ny_Nz_y')
        pl.show()

    img_laplace2 = ndimage.laplace(color.rgb2gray(image))

    return LN,img_laplace2,image_alpha

def genMeshFromImg(h0,scaler = 1.):
    '''
    generate mesh form input image. take the pixel value in h0 as height field
    :param h0: a .png file
    :return:
    '''
    mesh = openmesh.TriMesh()
    # add a a couple of vertices to the mesh
    vhs = np.ndarray(h0.shape, dtype=openmesh.VertexHandle)
    for r in range(0, h0.shape[0]):
        for c in range(0, h0.shape[1]):
            cvh = mesh.add_vertex(openmesh.TriMesh.Point(c, r, h0[r][c] * scaler))
            vhs[r][c] = cvh
    # add a couple of faces to the mesh
    fhs = np.ndarray(((h0.shape[0] - 1) * 2, h0.shape[1] - 1), dtype=openmesh.FaceHandle)
    for r in range(0, h0.shape[0] - 1):
        for c in range(0, h0.shape[1] - 1):
            cfh1 = mesh.add_face(vhs[r][c], vhs[r + 1][c + 1], vhs[r][c + 1])
            cfh2 = mesh.add_face(vhs[r][c], vhs[r + 1][c], vhs[r + 1][c + 1])
            fhs[2 * r][c] = cfh1
            fhs[2 * r + 1][c] = cfh2
    return mesh

# # the reshape operation is row by row in numpy and scipy
# @fn_timer
# def generateRelief(img_laplace,func_phi,C,miu=0,theta=0.001,*args,**kwargs):
#     '''
#     the difference of this and the above function (3) are: try to speed up the processes of generating matrixL, matrixC and matrixTheta
#     :param img_laplace:
#     :param func_phi:
#     :param miu:
#     :param theta:
#     :param C:
#     :param args:
#     :param kwargs:
#     :return:
#     '''
#
#     shape = img_laplace.shape
#
#     # get the indices of nonzeros in imgAlpha
#     # nonzeropixels = np.nonzero(imgAlpha)
#     # numHeightFieldPixels = len(nonzeropixels[0])
#     # maxtrixL = np.zeros((numHeightFieldPixels,numHeightFieldPixels))
#     # phiLN = func_phi(img_laplace)[nonzeropixels]
#     # for k in range(numHeightFieldPixels):
#
#     matrixL = sparse.lil_matrix((shape[0] * shape[1],shape[0] * shape[1]),dtype=np.float64) # lil_matrix is more efficient than csc_matrix, in both construction and compution phrases
#     t0 = time.time()
#     if 'logfile' in kwargs.keys():
#         logfile = kwargs['logfile']
#     else:
#         logfile = ''
#     myPrint(logfile,'begin generating Laplacian matrix')
#     height = shape[0]
#     width = shape[1]
#     # use the weight2 as the convolve mask
#     # set the main diagonal
#     maindiag = np.ones(width * height) * (-8)
#     matrixL.setdiag(maindiag,k=0)
#
#     updiag = np.ones(width * height - 1)
#     for k in range(1,height + 1):
#         if k * width - 1 < width * height - 1:
#             updiag[k * width - 1] = 0
#     matrixL.setdiag(updiag,k = 1)
#     matrixL.setdiag(updiag,k = -1)# the lower diagonal is the same as the upper one
#
#     diagW = np.ones(width * height - width) # diag W
#     matrixL.setdiag(diagW,k = width)
#     diagWs1 = np.ones(width * height - (width - 1))# diag (w - 1)
#     for k in range(height):
#         if k * width < width * height - (width - 1):
#             diagWs1[k * width] = 0
#     matrixL.setdiag(diagWs1,k = width - 1)
#
#     diagWp1 = np.ones(width * height - (width + 1))# diag (w + 1)
#     for k in range(1,height):
#         if k * width - 1 < width * height - (width + 1):
#             diagWp1[k * width - 1] = 0
#     matrixL.setdiag(diagWp1,k = width + 1)
#
#     # L is symtric, so...
#     matrixL.setdiag(diagW,k = -width)
#     matrixL.setdiag(diagWs1,k = - (width - 1))
#     matrixL.setdiag(diagWp1,k = - (width + 1))
#     t1 = time.time()
#     totaltime = t1 - t0
#     myPrint(logfile,'time of generating Laplace matrix: {totaltime} seconds'.format(totaltime = totaltime))
#     t0 = time.time()
#     matrixLT = matrixL.transpose() # this is a copy
#     rawMatrixL = matrixLT.transpose()
#     matrixL = matrixLT.dot(matrixL) # multiply is element by element, dot is the real multiplication of matrices.
#     t1 = time.time()
#     print 'time of computing L^T * L is:', t1 - t0
#     myPrint(logfile,'number of non-zeros in matrixL: {n}'.format(n = matrixL.count_nonzero()))
#     # myPrint(logfile,'matrixL looks like: {m}'.format(m = matrixL))
#
#     t0 = time.time()
#     myPrint(logfile,'begin generating orhter matrices')
#     print 'number of pixels in C:',np.count_nonzero(C)
#     CIndices = np.nonzero(C)
#     print 'type of CIndices:',type(CIndices)
#     CFlags = sparse.lil_matrix(C.reshape(C.shape[0] * C.shape[1],1) * 1.)
#     print 'type of CFlags:',type(CFlags)
#     print 'shape of CFlags:',CFlags.shape
#     print 'dtype of CFlags:',CFlags.dtype
#     print 'nonzeros in CFlags:',CFlags.count_nonzero()
#     matrixTheta = CFlags.dot(theta) # theta is a scalar, so multiply is ok. otherwise, use dot
#     print 'type of matrixTheta:',type(matrixTheta)
#     print 'shape of matrixTheta:',matrixTheta.shape
#     print 'dtype of matrixTheta:',matrixTheta.dtype
#     print 'nonzeros in matrixTheta:',matrixTheta.count_nonzero()
#     # print 'matrixTheta looks like:',matrixTheta
#     matrixC = sparse.lil_matrix((shape[0] * shape[1],shape[0] * shape[1]))
#     matrixC.setdiag(C.reshape(C.shape[0] * C.shape[1],) * 1.)
#     print 'type of matrixC:',type(matrixC)
#     print 'shape of marixC:',matrixC.shape
#     print 'dtype of matrixC:',matrixC.dtype
#     print 'nonzeros in matrixC:',matrixC.count_nonzero()
#     t1 = time.time()
#     totaltime = t1 - t0
#     myPrint(logfile,'time of generating C and Theta matrix: {time} seconds'.format(time = totaltime))
#     matrixCT = matrixC.transpose()
#     matrixC = matrixCT.dot(matrixC)
#     # print 'matrices looks like:', matrixC,matrixTheta
#     myPrint(logfile,'nonzeros in matrix C and Theta are: {nnzC} and {nnzT}'.format(nnzC = matrixC.count_nonzero(),nnzT = matrixTheta.count_nonzero()))
#
#     print 'begin computing matrix A and B'
#     t0 = time.time()
#     matrixPhi = func_phi(img_laplace).reshape((shape[0] * shape[1],1))
#
#     print 'types of matrixLT and matrixCT:',type(matrixLT),type(matrixCT)
#
#     matrixA = matrixL + miu * miu * matrixC
#     matrixB = np.array(matrixLT.dot(matrixPhi) + miu * miu * matrixCT.dot(matrixTheta))
#     t1 = time.time()
#     myPrint(logfile,'time of computing matrix A and B are: {t} seconds'.format(t = t1 - t0))
#     # matrixB = matrixB.reshape((shape[0] * shape[1],1))
#     # matrixB = sparse.lil_matrix(matrixB,shape=(shape[0] * shape[1],1))
#     print 'types of matrix A and B:',type(matrixA),type(matrixB)
#     myPrint(logfile,'shapes of A and B: {shapeA} and {shapeB}'.format(shapeA = matrixA.shape,shapeB = matrixB.shape))
#     myPrint(logfile,'begin sovlving')
#     t0 = time.time()
#
#     rst = splinalg.lsqr(matrixA,matrixB,show=True,iter_lim=100000)
#     t1 = time.time()
#     totaltime = t1 - t0
#     myPrint(logfile,'time of solving: {time} seconds'.format(time = totaltime))
#
#     myPrint(logfile,'H:{h}'.format(h = rst[0]))
#     myPrint(logfile,'istop: {i}'.format(i = rst[1]))
#     myPrint(logfile,'itn: {i}'.format(i = rst[1]))
#     myPrint(logfile,'r1norm: {n}'.format(n = rst[3]))
#     myPrint(logfile,'r2norm: {n}'.format(n = rst[4]))
#
#     H = rst[0].reshape(shape)
#
#     # check
#     weight = np.array([[1,1,1],[1,-8,1],[1,1,1]])
#     print 'begin checking'
#     convH = ndimage.convolve(H,weight)
#     residual = convH - matrixPhi.reshape(shape[0], shape[1])
#
#     pl.subplot(221)
#     pl.imshow(H,pl.cm.gray)
#     pl.title('$H$')
#     pl.subplot(222)
#     pl.imshow(convH,pl.cm.gray)
#     pl.title('$L(H)$\nmin={min}\nmax={max}'.format(min=convH.min(),max=convH.max()))
#     pl.subplot(223)
#     pl.imshow(matrixPhi.reshape(shape[0],shape[1]),pl.cm.gray)
#     pl.title('$\phi(LN)$\nmin={min}\nmax={max}'.format(min=matrixPhi.min(),max=matrixPhi.max()))
#     pl.subplot(224)
#     pl.imshow(residual,pl.cm.gray)
#     pl.title('$L(H)-\phi(LN)$\nmin={min}\nmax={max}'.format(min=residual.min(),max=residual.max()))
#     pl.show()
#     print 'residual:',type(residual),residual.dtype,residual.shape,residual.min(),residual.max(),np.linalg.norm(residual)
#     print 'end checking'
#
#     scalar = 1.
#     if 'scalar' in kwargs.keys():
#         scalar = kwargs['scalar']
#         print 'scaling mesh with scalar',scalar
#
#
#     H = scalar * H
#     print 'min and max values of recovered height field:',H.min(),H.max()
#     return genMeshFromImg(H),H

# the reshape operation is row by row in numpy and scipy
# there is must be enough space between boundary of image and the bas-relief area

def checkRst(shape,rst,validpixels,func_phi,LN,weight,show=False,asFastAsPossible=True):
    H = np.zeros(shape)
    H[validpixels] = rst.reshape((len(validpixels[0]),))
    H = H.reshape(shape)
    if asFastAsPossible:
        return H

    # check
    # weight = np.array([[1, 1, 1], [1, -8, 1], [1, 1, 1]])
    print 'begin checking'
    convH = ndimage.convolve(H, weight)
    # convH[validpixels] = 0. # important
    matrixPhi = func_phi(LN)
    residual = convH[validpixels] - matrixPhi[validpixels] # casue some invalid pixels in convH would result in false error

    print 'residual:', type(residual), residual.dtype, residual.shape, residual.min(), residual.max(), np.linalg.norm(
        residual)

    if show:
        pl.figure()
        pl.subplot(221)
        pl.imshow(H, pl.cm.gray)
        pl.title('$H$')
        pl.subplot(222)
        pl.imshow(convH, pl.cm.gray)
        pl.title('$L(H)$\nmin={min}\nmax={max}'.format(min=convH.min(), max=convH.max()), {'fontsize': 10})
        pl.subplot(223)
        pl.imshow(matrixPhi.reshape(shape[0], shape[1]), pl.cm.gray)
        pl.title('$\phi(LN)$\nmin={min}\nmax={max}'.format(min=matrixPhi.min(), max=matrixPhi.max()), {'fontsize': 10})
        pl.subplot(224)
        pl.hist(residual, 20)
        pl.title('$L(H)-\phi(LN)$\nmin={min}\nmax={max}'.format(min=residual.min(), max=residual.max()), {'fontsize': 10})
        pl.show()

        mesh = genMeshFromImg(H)
        openmesh.write_mesh(mesh,'meshForChecking.off')
        normalmapexec = '../normalmapCoin3D/DerivedData/normalmapCoin3D/Build/Products/Release/normalmapCoin3D'
        datafolder = './'
        imageheights = ['0']
        showNormalMap.processFiles(normalmapexec,datafolder,imageheights,['meshForChecking.off'],'./',askkforstop=False)

    print 'end checking'

    return H

def saveH(H,**kwargs):
    scalar = 1.
    if 'scalar' in kwargs.keys():
        scalar = kwargs['scalar']
        # print 'scaling mesh with scalar', scalar

    H = scalar * H
    # print 'min and max values of recovered height field:', H.min(), H.max()
    return genMeshFromImg(H), H

@fn_timer
def generateRelief(LN,func_phi,image_alpha,weight,oriimg,C=None,miu=0,theta=0.001,show=False,asFastAsPossible=True,*args,**kwargs):
    '''
    the difference of this and the above function (3) are: try to speed up the processes of generating matrixL, matrixC and matrixTheta
    :param LN:
    :param func_phi:
    :param image_alpha: alpha values of image
    :param miu:
    :param theta:
    :param C: is the true-false matrix with the same dimensions as img_alpha, and the true element set must be a subset of image_alpha > 0
    :param args:
    :param kwargs:
    :return:
    '''

    if (not C.shape == LN.shape) or (not C.shape == image_alpha.shape):
        raise ValueError('Shapes of img_laplace, C and img_alpha must be equal')
    if not weight.shape ==(3,3):
        raise ValueError('Shape of weight must be (3,3)')

    shape = LN.shape
    imheight = shape[0]
    imwidth = shape[1]

    validFlags = image_alpha > 0 # generate a true-false matrix with the same shape of image_alpha, to denote which pixels not transparent
    validpixels = np.nonzero(image_alpha)  # the index of valid pixels in image_alpha.
    # by testing, we found that nonzero is row by row as well, just like reshape
    numvalidpixels = len(validpixels[0]) # number of valid pixels
    validPixelIndex = np.ones(shape,dtype=np.int32) * -1 # we arrange the index of valid pixels into 1-D array 'validpixels',
    # and this is the index of valid pixels in validpixels
    validPixelIndex[validFlags] = range(numvalidpixels)
    for k in range(numvalidpixels):
        r = validpixels[0][k]
        c = validpixels[1][k]
        if not validPixelIndex[r][c] == k:
            print k,validPixelIndex[r][c]
            raise ValueError('index error: type 1')
    print 'pass index checking: type 1'
    for r in range(imheight):
        for c in range(imwidth):
            k = validPixelIndex[r][c]
            if k >= 0 and (not r == validpixels[0][k] or not c == validpixels[1][k]):
                print k
                print r,c
                print validpixels[0][k],validpixels[1][k]
                raise ValueError('index error: type 2')
    print 'pass index checking: type 2'
    print 'number of valid pixels:',numvalidpixels
    # pl.imshow(validPixelIndex,pl.cm.jet)
    # pl.show()
    matrixL = sparse.lil_matrix((numvalidpixels, numvalidpixels),
                                dtype=np.float64)  # lil_matrix is more efficient than csc_matrix, in both construction and compution phrases
    if 'logfile' in kwargs.keys():
        logfile = kwargs['logfile']
    else:
        logfile = ''

    myPrint(logfile,'begin generating Laplacian matrix')
    t0 = time.time()
    matrixL.setdiag(np.ones((numvalidpixels,)) * weight[1][1]) # the main diagonal is always -8
    # startIDs = np.ones((imheight,),dtype=np.int32) * -1 # get the index of first valid pixel for every row
    # startIDs[validpixels[0][0]] = 0
    # currow = validpixels[0][0]
    # compute startIDs
    # for k in range(1,numvalidpixels):
    #     r = validpixels[0][k]
    #     if r > currow:
    #         currow = r
            # startIDs[r] = k

    # compute elements for every row in matrixL
    # convolve mask must be 0-2 switching
    for k in range(numvalidpixels):
        r = validpixels[0][k]
        c = validpixels[1][k]

        if k > 0 and validpixels[0][k - 1] == r and validpixels[1][k - 1] == c - 1: # if last valid pixel is [r, c -1], set it to be 1.
            matrixL[k, k - 1] = weight[1][2]
        if k < numvalidpixels - 1 and validpixels[0][k + 1] == r and validpixels[1][k + 1] == c + 1: # if next valid pixel is [r, c + 1], set it to be 1.
            matrixL[k, k + 1] = weight[1][0]
        if r - 1 >= 0:
            if c - 1 >= 0 and validFlags[r - 1][c - 1]:
                matrixL[k,validPixelIndex[r - 1][c - 1]] = weight[2][2]
            if validFlags[r - 1][c]:
                matrixL[k,validPixelIndex[r - 1][c]] = weight[2][1]
            if c + 1 < imwidth and validFlags[r - 1][c + 1]:
                matrixL[k,validPixelIndex[r - 1][c + 1]] = weight[2][0]
        if r + 1 < imheight:
            if c - 1 >= 0 and validFlags[r + 1][c - 1]:
                matrixL[k, validPixelIndex[r + 1][c - 1]] = weight[0][2]
            if validFlags[r + 1][c]:
                matrixL[k, validPixelIndex[r + 1][c]] = weight[0][1]
            if c + 1 < imwidth and validFlags[r + 1][c + 1]:
                matrixL[k, validPixelIndex[r + 1][c + 1]] = weight[0][0]
    t1 = time.time()
    totaltime = t1 - t0
    myPrint(logfile,'time of generating Laplace matrix: {totaltime} seconds'.format(totaltime = totaltime))

    # test matrix
    if not asFastAsPossible:
        print 'begin testing matrixL'
        image = oriimg * image_alpha
        randompixels = image[validpixels]

        if show:
            pl.figure()
            pl.subplot(221)
            pl.imshow(image,pl.cm.gray)
            pl.title('original image')

        conv = ndimage.convolve(image,weight)
        conv[image_alpha < 1] = 0. # cause invalid pixels around boudnary would affect the result
        if show:
            pl.subplot(222)
            pl.imshow(conv,pl.cm.jet)
            pl.title('conv with weight')

        conv2 = matrixL.dot(randompixels)
        image[validpixels] = conv2
        if show:
            pl.subplot(223)
            pl.imshow(image,pl.cm.jet)
            pl.title('Laplacian by matrixL')

        residual = conv2.reshape((numvalidpixels,)) - conv[validpixels].reshape((numvalidpixels,))
        print 'residual:',np.linalg.norm(residual),residual.min(),residual.max()
    # print np.nonzero(residual),residual[0:10]
    # print validpixels[0][0:35],validpixels[1][0:35]
    # print conv2.shape,conv.shape
    # print conv2[0],conv[52][295]
    # print matrixL[0].nonzero()
    # print matrixL[0][matrixL[0].nonzero()]
    # print image[52][295:297],image[53][294:297]
    # print -8 * image[52][295] + image[52][296] + image[53][294] + image[53][295] + image[53][296]

        image[validpixels] = residual
        if show:
            pl.subplot(224)
            pl.hist(residual[residual != 0.],20)
            pl.title('residual hist')
            pl.show()

        print 'end of testing matrixL'
        matrixLCopy = matrixL.copy()
    # end of testing matrixL

    # generate L^T * L
    t0 = time.time()
    matrixLT = matrixL.transpose() # this is a copy
    matrixL = matrixLT.dot(matrixL) # multiply is element by element, dot is the real multiplication of matrices.
    t1 = time.time()
    if not asFastAsPossible:
        print 'time of computing L^T * L is:', t1 - t0
        myPrint(logfile,'number of non-zeros in L^T * L: {n}'.format(n = matrixL.count_nonzero()))

    # generate C and theta
    t0 = time.time()
    myPrint(logfile,'begin generating orhter matrices')
    if not asFastAsPossible:
        print 'number of pixels in C:',np.count_nonzero(C)
    matrixTheta = sparse.lil_matrix((numvalidpixels,1))  # theta is a scalar, so multiply is ok. otherwise, use dot
    matrixC = sparse.eye(numvalidpixels)
    for k in range(numvalidpixels):
        r = validpixels[0][k]
        c = validpixels[1][k]
        if not C[r][c]:
            matrixC[k,k] = 0.
        else:
            matrixTheta[k] = theta
    t1 = time.time()
    totaltime = t1 - t0
    myPrint(logfile,'time of generating C and Theta matrix: {time} seconds'.format(time = totaltime))

    # generate C^T * C
    # matrixCT = matrixC.transpose()
    # matrixC = matrixCT.dot(matrixC) # C is sysmetric and eye-like, C=CT=CT*C
    if not asFastAsPossible:
        myPrint(logfile,'nonzeros in C^T * C and Theta are: {nnzC} and {nnzT}'.format(nnzC = matrixC.count_nonzero(),nnzT = matrixTheta.count_nonzero()))

    # generate matrix A and B
    print 'begin computing matrix A and B'
    t0 = time.time()
    matrixPhi = func_phi(LN)
    matrixPhi = matrixPhi[validpixels].reshape((numvalidpixels,1))

    matrixA = matrixL + miu * miu * matrixC
    matrixB = np.array(matrixLT.dot(matrixPhi) + miu * miu * matrixTheta) # since matrixC is diagonal and similar to eye, it no use to multiply C^T with \theta
    t1 = time.time()
    myPrint(logfile,'time of computing matrix A and B are: {t} seconds'.format(t = t1 - t0))
    if not asFastAsPossible:
        myPrint(logfile,'shapes of A and B: {shapeA} and {shapeB}'.format(shapeA = matrixA.shape,shapeB = matrixB.shape))

    # solving the system
    myPrint(logfile, 'begin system solving')
    t0 = time.time()

    # use different methods to solve Ax=b, see document in sparse.linalg

    # directive methods
    validMethods = ['spsolve']# ''spsolve' ,'factorized','bicg','bicgstab','cg','cgs','no-gmres','lgmres','minres','qmr','no-lsqr','lsmr']
    def callSolver(solver,*args,**kwargs):
        try:
            print '#' * 30
            print 'using %s' % solver.__name__
            tt0 = time.time()
            rst = solver(*args,**kwargs)
            rawrst = rst.copy()
            print 'time:', time.time() - tt0
            if type(rst) == tuple:
                print rst[1]
                rst = rst[0]
            residual = matrixA.dot(rst.reshape((numvalidpixels, 1))) - matrixB

            print 'equation system residual:', np.linalg.norm(residual, ord=2), residual.min(), residual.max()

            H = checkRst(shape, rst, validpixels, func_phi, LN, weight, show=show,asFastAsPossible=asFastAsPossible)
            return H,rawrst
        except:
            return None
    if 'spsolve' in validMethods:
        H,rst = callSolver(splinalg.spsolve,matrixA,matrixB)
    if 'factorized' in validMethods:
        solver = splinalg.factorized(matrixA)
        H,rst = callSolver(solver,matrixB)
    if 'bicg' in validMethods: # a little inaccurate
        H,rst = callSolver(splinalg.bicg,matrixA,matrixB,maxiter = 100000)
    if 'bicgstab' in validMethods:# a little inaccurate
        H,rst = callSolver(splinalg.bicgstab,matrixA,matrixB,maxiter=100000)
    if 'cg' in validMethods: # a little inaccurate
        H,rst = callSolver(splinalg.cg, matrixA, matrixB, maxiter=100000)
    if 'cgs' in validMethods: # a little inaccurate
        H,rst = callSolver(splinalg.cgs, matrixA, matrixB, maxiter=100000)
    if 'gmres' in validMethods: # too much time
        H,rst = callSolver(splinalg.gmres, matrixA, matrixB, restart=20,maxiter=10000)
    if 'lgmres' in validMethods:# a little inaccurate
        H,rst = callSolver(splinalg.lgmres, matrixA, matrixB, maxiter=100000)
    if 'minres' in validMethods: # a little inaccurate but pretty fast
        H,rst = callSolver(splinalg.minres, matrixA, matrixB, shift=0.0,maxiter=100000,show=True)
    if 'qmr' in validMethods:# a little inaccurate
        H,rst = callSolver(splinalg.qmr, matrixA, matrixB, maxiter=100000)
    if 'lsqr' in validMethods:# a little inaccurate and too much time
        H,rst= callSolver(splinalg.lsqr, matrixA,matrixB,damp=0.0,show=True,iter_lim=100000)
        print rst[1:5]
    if 'lsmr' in validMethods: # a little inaccurate and more time
        H, rst = callSolver(splinalg.lsmr, matrixA, matrixB, damp=0.0, show=True, iter_lim=100000)
        print rst[1:5]

    t1 = time.time()
    totaltime = t1 - t0
    myPrint(logfile, 'time of solving: {time} seconds'.format(time=totaltime))

    # test matrix
    if not asFastAsPossible:
        print 'begin testing matrixL'
        matrixL = matrixLCopy
        image = H.copy()
        randompixels = H[validpixels]

        if show:
            pl.figure()
            pl.subplot(331)
            pl.imshow(image, pl.cm.gray)
            pl.title('height field')

        conv = ndimage.convolve(image, weight)
        if show:
            pl.subplot(332)
            pl.imshow(conv, pl.cm.gray)
            pl.title('conv with weight')

        conv2 = matrixL.dot(randompixels)
        image[validpixels] = conv2
        if show:
            pl.subplot(333)
            pl.imshow(image, pl.cm.jet)
            pl.title('Laplacian by matrixL')

        residual = conv2.reshape((numvalidpixels,)) - conv[validpixels].reshape((numvalidpixels,))
        print 'residual: matrixL - conv', np.linalg.norm(residual), residual.min(), residual.max()
        image[validpixels] = residual
        if show:
            pl.subplot(334)
            pl.imshow(image, pl.cm.gray)
            pl.title('residual: matrixL - conv')

        residual = conv2.reshape((numvalidpixels,)) - matrixPhi.reshape((numvalidpixels,))
        print 'residual: matrixL - \phi LN', np.linalg.norm(residual), residual.min(), residual.max()
        image[validpixels] = residual
        if show:
            pl.subplot(335)
            pl.imshow(image, pl.cm.gray)
            pl.title('residual: matrixL - $\phi LN$')

        residual = conv[validpixels].reshape((numvalidpixels, )) - matrixPhi.reshape((numvalidpixels,))
        print 'residual: conv - \phi LN', np.linalg.norm(residual), residual.min(), residual.max()
        image[validpixels] = residual
        if show:
            pl.subplot(336)
            pl.imshow(image, pl.cm.gray)
            pl.title('residual: conv - $\phi LN$')

        validLN = func_phi(LN[validpixels]).reshape((numvalidpixels,))
        residual = conv2.reshape((numvalidpixels,)) - validLN
        print 'residual: matrixL - LN', np.linalg.norm(residual), residual.min(), residual.max()
        image[validpixels] = residual
        if show:
            pl.subplot(337)
            pl.imshow(image, pl.cm.gray)
            pl.title('residual: matrixL - LN')

        residual = conv[validpixels].reshape((numvalidpixels,)) -validLN
        print 'residual: conv - LN', np.linalg.norm(residual), residual.min(), residual.max()
        image[validpixels] = residual
        if show:
            pl.subplot(338)
            pl.imshow(image, pl.cm.gray)
            pl.title('residual: conv - $LN$')

        residual = validLN - matrixPhi.reshape((numvalidpixels,))
        print 'residual: LN - \phi LN', np.linalg.norm(residual), residual.min(), residual.max()
        image[validpixels] = residual
        if show:
            pl.subplot(339)
            pl.imshow(image, pl.cm.gray)
            pl.title('residual: LN - $LN$')
            pl.show()

        print scipy.count_nonzero(conv)
        conv[image_alpha < 1] = 0.
        print scipy.count_nonzero(conv)
        # LN[image_alpha < 1] = 0. # this should be true, without explicitly setting
        residual = conv - func_phi(LN)
        print 'conv - LN',np.linalg.norm(residual),residual.min(),residual.max()

        if show:
            pl.figure()
            pl.subplot(131)
            pl.imshow(conv,pl.cm.gray)
            pl.title('conv')
            pl.subplot(132)
            pl.imshow(LN,pl.cm.gray)
            pl.title('LN')
            pl.subplot(133)
            pl.imshow(residual,pl.cm.gray)
            pl.title('conv - LN')
            pl.show()
        print 'end of testing matrixL'
        # end of testing matrixL

    # myPrint(logfile,'H:{h}'.format(h = rst[0]))
    # myPrint(logfile,'istop: {i}'.format(i = rst[1]))
    # myPrint(logfile,'itn: {i}'.format(i = rst[1]))
    # myPrint(logfile,'r1norm: {n}'.format(n = rst[3]))
    # myPrint(logfile,'r2norm: {n}'.format(n = rst[4]))

    mesh,H = saveH(H,**kwargs)
    return mesh,H,LN

def main():
    allt0 = time.time()
    timeSavingFiles = 0
    reliefImFolder = './ReliefImages'
    reliefMeshFolder = './ReliefMeshes'
    show = False
    asFastAsPossible = False

    options,args = getopt.getopt(sys.argv[1:],'',['input=','reliefImDir=','reliefMeshDir='])

    filename = ''
    for name,value in options:
        if name in ('--input',):
            filename = value
        elif name in ('--reliefImDir',):
            reliefImFolder = value
        elif name in ('--reliefMeshDir',):
            reliefMeshFolder = value
        else:
            print 'unprocessed argument:',name, 'with value',value

    image = loadImg(filename) # 0 - 255
    # image = scipy.ndimage.gaussian_filter(image,sigma=[1,1,0])
    image[:,:,0:3] = scipy.ndimage.gaussian_filter(image[:,:,0:3], sigma=[1.5,1.5,0.]) # affect the details and singular values in LN

    # pl.figure()
    # pl.imshow(image[:,:,3], pl.cm.gray)
    # pl.show()

    image_alpha = img_as_float(image[:, :, 3])
    # convert the range of pixel values from 0 - 255 to -1 - 1, which is the correct range for normal values
    image = img_as_float(image) # 0 - 1
    oriimg = image = (image - 0.5) * 2 # this is because that the conversion between normal and color in fragment shader is: color * 0.5 + 0.5
    # image[image_alpha < 1,:] = 0. # very very important
    # image = image[image_alpha > 0] # the following two lines to test whether the sum of square of rgb is 1.0 for every pixel.
    # image = image[:,  0] ** 2 + image[:,  1] ** 2 + image[:,  2] ** 2
    if not asFastAsPossible:
        print 'original image -- type:{type}, dtype:{dtype} and shape:{shape}'.format(type=type(image), dtype=image.dtype,
                                                                                      shape=image.shape)
        print 'original image -- min and max values:', image[:,:,0:3].min(), image[:,:,0:3].max()

    filename = os.path.basename(filename)
    filebasename = filename.split('.')[0]
    # print filebasename
    m = re.search(r'^(normalmapOf)(.+)', filebasename)  # only support png image
    rst = m.groups()
    if len(rst) > 1:
        filebasename = rst[1]
        # print filebasename
        # print rst

    # print 'sample of original image:',image[0:2,0:10,0:4]

    LN,img_laplace_common,image_alpha = generateLN(image,image_alpha,show=show,asFastAsPossible=asFastAsPossible)
    # LN = scipy.ndimage.gaussian_filter(LN,sigma=1)
    if not asFastAsPossible:
        print LN.shape,scipy.count_nonzero(image_alpha)
        # LN[image_alpha < 1] = 0. # important
    basImageName = reliefImFolder + '/' + filebasename + '_LN.png'
    if not asFastAsPossible:
        pl.imsave(basImageName, LN, cmap=pl.cm.gray)

    if not asFastAsPossible and show:
        pl.figure()
        pl.subplot(121)
        pl.imshow(LN,pl.cm.gray)
        pl.title('LN')
        pl.subplot(122)
        minln = LN.min()
        maxln = LN.max()
        pl.imshow(img_as_ubyte((LN - minln) / (maxln - minln)),pl.cm.gray)
        pl.title('converted LN to ubyte')
        pl.show(**{'block':False})

    weight1 = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]])
    weight2 = np.array([[1, 1, 1], [1, -8, 1], [1, 1, 1]])
    # weight3 = np.array([[0, -1, 0], [-1, 4, -1], [0, -1, -1]]) # wrong weight
    weight4 = np.array([[-1, 1, -1], [1, 8, -1], [-1, 1, -1]])
    weight5 = np.array([[19, 27, 19], [27, -184, 27], [19, 27, 19]])
    weight = weight1 # try other weights !!

    if not asFastAsPossible:
        print 'Laplacian image from normal map -- type:{type}, dtype:{dtype} and shape:{shape}'.format(type = type(LN),dtype=LN.dtype,shape=LN.shape)
        print 'Laplacian image from normal map -- min and max values:',LN.min(),LN.max()

    if not asFastAsPossible and show:
        pl.figure()
        pl.subplot(221)
        pl.imshow(image)
        pl.title('normal values')

        pl.subplot(222)
        pl.imshow(LN, pl.cm.gray)
        pl.title('Laplacian image form normal map: $L(\mathbf{N})$')

        pl.subplot(223)
        pl.title('laplace image using ndimage.laplace')
        pl.imshow(img_laplace_common, pl.cm.gray)

        pl.show(**{'block': True})

    if not asFastAsPossible:
        print 'Laplacian of original image -- type:{type}, dtype:{dtype} and shape:{shape}'.format(
            type=type(img_laplace_common), dtype=img_laplace_common.dtype, shape=img_laplace_common.shape)
        print 'Laplacian of original image -- min and max values:', img_laplace_common.min(), img_laplace_common.max()

    imLaplaceName = reliefImFolder + '/' + filebasename + '_Laplacian.png'
    if not asFastAsPossible:
        pl.imsave(imLaplaceName,img_laplace_common,cmap=pl.cm.gray)

    # some other kind of bas relief generated from images
    image = color.rgb2gray(oriimg)
    # image = scipy.ndimage.gaussian_filter(image,1)
    if not asFastAsPossible:
        print 'generating mesh from original image. min and max values:',image.min(),image.max()
        mesh = genMeshFromImg(image,image.shape[0] / 30)
        imName = reliefImFolder + '/' + filebasename + '_OriGray.png'
        pl.imsave(imName, image, cmap=pl.cm.gray)
        filename = reliefMeshFolder + '/' + filebasename + '_BasRelief_fromOriImg.off'
        t0 = time.time()
        openmesh.write_mesh(mesh, filename)
        timeSavingFiles += time.time() - t0

    # print 'generating mesh from processed image. min and max values:', img_processed.min(), img_processed.max()
    # mesh = genMeshFromImg(img_processed)
    # filename = reliefMeshFolder + '/' + filebasename + '_BasRelief_fromImgProcessed.off'
    # openmesh.write_mesh(mesh, filename)
    #
    # print 'generating mesh from image laplacian. min and max values:', img_laplace2.min(), img_laplace2.max()
    # mesh = genMeshFromImg(img_laplace2)
    # filename = reliefMeshFolder + '/' + filebasename + '_BasRelief_fromImgLaplace.off'
    # openmesh.write_mesh(mesh, filename)
    #
    # image = pl.imread(basImageName)
    # image = color.rgb2gray(image)
    # print 'generating mesh from IMAGE of processed img. min and max values:', image.min(), image.max()
    # mesh = genMeshFromImg(image)
    # filename = reliefMeshFolder + '/' + filebasename + '_BasRelief_fromImgProcessedImg.off'
    # openmesh.write_mesh(mesh, filename)
    #
    # image = pl.imread(imLaplaceName)
    # image = color.rgb2gray(image)
    # print 'generating mesh from IMAGE of image Laplacian. min and max values:', image.min(), image.max()
    # mesh = genMeshFromImg(image)
    # filename = reliefMeshFolder + '/' + filebasename + '_BasRelief_fromImgLaplaceImg.off'
    # openmesh.write_mesh(mesh, filename)

    # mesh generated using optimization
    # img_alpha = oriimg[:,:,3]
    # mesh = generateRelief(img_processed,lambda t: 1. * t,miu = 0.001,theta=0.05,C=img_alpha > 0,imgAlpha=img_alpha)
    logfile = reliefImFolder + '/' + filebasename + '_log.txt'
    with open(logfile,'w') as f:
        pass
    # print 'samples of $LN$:',LN[1:100,0]
    phivalue = 1.0
    miuvalue = 0.001
    thetavalue = 0.05
    mesh,heightimg,LN2 = generateRelief(LN,
                                        lambda t: phivalue * t, image_alpha=image_alpha,
                                        C = image_alpha > 0,
                                        oriimg=image,weight=weight,
                                        miu=miuvalue, theta= - thetavalue,# !! why it should be -3 instead of 3??
                                        logfile = logfile,
                                        scalar = 1.,
                                        show=show,asFastAsPossible=asFastAsPossible)
    filename = reliefMeshFolder + '/' + filebasename + '_BasRelief-phi({phi})-miu({miu})-theta({theta}).off'.format(
        phi=phivalue,miu=miuvalue,theta=thetavalue)
    print 'saved mesh file name:',filename
    t0 = time.time()
    openmesh.write_mesh(mesh,filename)
    timeSavingFiles += time.time() - t0

    if not asFastAsPossible:
        imLaplaceName = reliefImFolder + '/' + filebasename + '_HeightField.png'
        print 'Recoveried height field image -- type:{type}, dtype:{dtype} and shape:{shape}'.format(
            type=type(heightimg), dtype=heightimg.dtype, shape=heightimg.shape)
        print 'Recoveried height field image -- min and max values:', heightimg.min(), heightimg.max()
        pl.imsave(imLaplaceName, heightimg, cmap=pl.cm.gray)

    allt1 = time.time()
    myPrint(logfile,'all time: {t} seconds, in which time of saving meshes is: {t2} seconds'.format(t = allt1 - allt0,t2 = timeSavingFiles))

if __name__ == "__main__":
    main()
