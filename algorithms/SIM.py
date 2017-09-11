from SegAlgo import *

class SIM(SegAlgo):
    def __init__(self,mesh, volume):
        super(SIM,self).__init__(mesh,volume)

    @fn_timer
    def Optimize(self,*args,**kwargs):
        ''' Optimize the mesh.'''
        mesh = self._mesh
        volume = self._volume

        S_I = copy.copy(mesh)
        S_B = copy.copy(mesh)
        S_SSM = copy.copy(mesh)

        global km
        km = 1
        dirname = './Data/tmp'
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.mkdir(dirname)

        def OptimizeSI():  # more details about scipy.optimization please refer to the optimize __init__.py
            nv = S_I.nv  # the number of vertices used to test the optimization efficiency

            # print 'number of vertices used:',nv

            def showMeshAndNormal():
                '''
                To test if the normals is pointing from inside to outlide and if the ordering of the normals is right.
                :return:
                '''
                m = visuals.Mesh(S_B.verts, mesh.faces)
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
            signPixel = 1.  # put these two here to convinently modify. (-1.,-1.) for eclipsoid case. signPixel is critical

            # for success
            # if the function is defined well, it can also converge to the right result without help of gradient
            # with gradient help, the convergence can be speeded up a lot
            # the initial guess is QUITE important for the final result. See transform.py
            # right fun + wrong gradient can not guarantee the right result
            # fun and gradient must be concide, wired otherwise
            # of course you can only give fun without gradient, but, please make sure the given gradient is as precise as possible,
            # if you prefer to give.

            def E_I(SIVerts):
                SIVerts = SIVerts.reshape((nv, 3))

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
                            return 0, 0  # pixel and then gradient

                        # use trilinear interpolation, refer to wiki
                        xd = v[0] - x
                        yd = v[1] - y
                        zd = v[2] - z

                        X = np.array([x, x2])
                        Y = np.array([y, y2])
                        Z = np.array([z, z2])
                        vPixel = np.empty((2, 2, 2), dtype=np.float32)
                        vGradient = np.empty((2, 2, 2), dtype=np.float32)
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
                        vgradient = triInt(vGradient)

                        return vpixel, vgradient

                    def noInterp():
                        x = int(v[0] + .5)
                        y = int(v[1] + .5)
                        z = int(v[2] + .5)

                        if x < 0 or y < 0 or z < 0 or x >= dataShape[0] or y >= dataShape[1] or z >= dataShape[2]:
                            print 'vert %s is outside of volume data in E_I:' % k, x, y, z, 'vert coordinate:', v
                            return 0, 0  # pixel and then gradient

                        gradient = np.array([gradientX[x][y][z],
                                             gradientY[x][y][z],
                                             gradientZ[x][y][z]])
                        vgradient = omega2 * (normals[k].dot(gradient))
                        vpixel = omega1 * pixelarray[x][y][z]
                        return vpixel, vgradient

                    vpixel, vgradient = trilinearInterp()
                    # vpixel,vgradient = noInterp()

                    # val = np.square(np.sqrt(((v - np.array([50.,50.,50.]))**2).sum()) - 25.)
                    val = signGradient * vgradient - signPixel * vpixel
                    ei += val
                    k += 1

                return ei

            def gradient_EI(SIVerts):
                SIVerts = SIVerts.reshape((nv, 3))

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
                            return 0, 0, 0, 0, 0, 0  # pixel and then gradient

                        # use trilinear interpolation, refer to wiki
                        xd = v[0] - x
                        yd = v[1] - y
                        zd = v[2] - z

                        X = np.array([x, x2])
                        Y = np.array([y, y2])
                        Z = np.array([z, z2])
                        gXgradient = np.empty((2, 2, 2), dtype=np.float32)
                        gXpixel = np.empty((2, 2, 2), dtype=np.float32)
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

                        return gxpixel, gypixel, gzpixel, gxgradient, gygradient, gzgradient

                    def noInterp():
                        x = int(v[0] + .5)
                        y = int(v[1] + .5)
                        z = int(v[2] + .5)

                        if x < 0 or y < 0 or z < 0 or x >= dataShape[0] or y >= dataShape[1] or z >= dataShape[2]:
                            gx = 9999  # the same as fun
                            gy = 9999
                            gz = 9999

                            print 'vert %s is outside of volume data in gradient E_I:' % k, x, y, z, 'vert coordinate:', v
                            return 0, 0, 0, 0, 0, 0  # pixel and then gradient

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
                        return gxpixel, gypixel, gzpixel, gxgradient, gygradient, gzgradient

                    gxpixel, gypixel, gzpixel, gxgradient, gygradient, gzgradient = trilinearInterp()
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
                graderr = opt.check_grad(E_I, gradient_EI, S_I.verts.copy().reshape((nv * 3,)))
                print 'graderr for EI when optimizing SI is:', graderr

            refMesh = S_B  # original is S_B

            def E_I_B(SIVerts):
                SIVerts = SIVerts.reshape((nv, 3))

                verts_SB = refMesh.verts[0:nv]
                verts_SI = SIVerts

                verts_BI = verts_SB - verts_SI
                # return np.sqrt(np.sum(verts_BI * verts_BI, axis=1))
                return np.sum(verts_BI * verts_BI) / 2.

            def gradient_EIB(SIVerts):
                SIVerts = SIVerts.reshape((nv, 3))

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
                verts = SIVerts.reshape((nv, 3))

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
                verts = verts.reshape((nv, 3))
                global km
                openmesh.writeOff(dirname + '/%s_registered_SI.off' % km, verts, mesh.faces)
                km += 1

            # l-bfgs-b,SLSQP, are OK. fastest is l-bfgs-b,
            #  with gradient (or even hessian) given, the computation efficiency
            # and accuracy are highly improved.
            # without any derivative information, we can choose: l-bfgs-b and SLSQP
            # with gradient, we can choose:CG,Newton-CG,l-bfgs-b,TNC,
            # with gradient and hessian, we can choose: dogleg and trust-cg (both not tested)

            res = opt.minimize(fun, refMesh.verts.reshape((nv * 3,)), args=(), method='l-bfgs-b', jac=gradient,
                               callback=callback, tol=1e-6, options={'maxiter': 100, 'disp': True, 'gtol': 1e-6})
            x = res['x'].reshape((nv, 3))
            print 'x:', x
            # print 'mesh',mesh.verts
            print 'success?', res['success']
            print 'message:', res['message']

            S_I.verts = x.copy()

            return x

        def OptimizeSSM():
            '''
            repeat optimize transformation and SSM parameters until converge on S_SSM
            :return:
            '''
            phim = mesh.phim
            refMesh = S_I  # original is S_B

            def transformSSM(a, R, c, b):
                verts_SSM = mesh.verts.copy()  # mesh is the mean mesh.

                if not len(phim) == len(b):
                    raise ValueError('Length of b(%s) must equal lenght of $\phi$(%s):', len(b), len(phim))

                k = 0
                for bm in b:
                    verts_SSM += phim[k] * bm
                    k += 1

                verts_SSM = verts_SSM.transpose()
                verts_SSM = R.dot(verts_SSM) * a + c.reshape(3, 1)  # mesh is the mean mesh
                return np.array(verts_SSM.transpose())

            def E_SSM(a, R, c, b):  # translation, rotation, scale, and SSM parameters
                verts_ref = refMesh.verts
                verts_SSM = np.array(transformSSM(a, R, c, b))

                verts_BSSM = verts_SSM - verts_ref
                return np.sum(verts_BSSM * verts_ref) / 2.

            def fun_trans(vs, b=np.zeros(len(phim))):
                a = vs[0]
                R = np.matrix(vs[1:10].reshape((3, 3)))
                c = np.array(vs[10:13])
                return E_SSM(a, R, c, b)

            def rotationConstraint(vs):
                R = np.matrix(vs[1:10].reshape((3, 3)))
                return R.transpose() * R - np.identity(3, dtype=np.float32)

            def fun_SSMP(b, a=0, R=np.matrix('1 0 0;0 1 0; 0 0 1'), c=np.array([0, 0, 0])):
                return E_SSM(a, R, c, b)

            nv = mesh.verts.shape[0]

            def callbackTrans(x):
                a = x[0]
                R = np.matrix(x[1:10].reshape((3, 3)))
                c = np.array(x[10:13])
                global km
                verts = transformSSM(a, R, c, np.zeros(len(phim)))
                openmesh.writeOff(dirname + '/%s_registered_SSM.off' % km, verts, mesh.faces)
                km += 1

            def callbackSSMP(b):
                pass

            b = np.ones(len(phim))  # initial guess
            a = 1.
            R = np.matrix('1. 0. 0.;0. 1. 0.;0. 0. 1.')
            c = np.array([0., 0., 0.], dtype=np.float32)
            verts_SSM = None
            for k in range(
                    1):  # only iterate 5 times since not confident about how to evaluate the convergence of a,R,c and b
                initvals = np.zeros((13,))
                initvals[0] = a
                initvals[1:10] = R.reshape((9,))
                initvals[10:13] = c
                res = opt.minimize(fun_trans, initvals, args=(b), method='slsqp', jac=False,
                                   callback=callbackTrans, options={'maxiter': 100, 'disp': True})
                x = res['x']
                a = x[0]
                R = np.matrix(x[1:10].reshape((3, 3)))
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

                verts_SSM = transformSSM(a, R, c, b)
                # here, we can terminate the loop by checking whether the change of verts_SSM is small enough
            S_SSM.verts = verts_SSM.copy()

            return verts_SSM.copy()

        def OptimizeSSM_ICP():
            verts_SSM = mesh.verts.copy()  # mesh is the mean mesh.
            verts_ref = S_I.verts.copy()

            initial = np.array([[0.01], [0.05], [0.01], [0.001], [0.001], [0.001]])
            params = np.empty((6,))
            basicICP.icp_point_to_point_lm(verts_SSM, verts_ref, initial=initial, loop=0, params=params)
            print 'params', params

            S_SSM.verts = mesh.verts.copy()
            # print 'verts before',S_SSM.verts
            S_SSM.verts = basicICP.transformpoints(S_SSM._verts, params=params)
            # print 'verts after', S_SSM.verts

            global km
            openmesh.writeOff(dirname + '/%s_registered_SSM.off' % km, S_SSM._verts, mesh.faces)
            km += 1

            return S_SSM.verts.copy()

        def OptimizeSB():
            gridDim = (5, 5, 5, 3)  # each dim must be > 5
            numcoords = gridDim[0] * gridDim[1] * gridDim[2] * gridDim[
                3]  # number of coordinates of control points when flattened
            tbounds = np.array(volume.vol_data_forComp.shape)
            nv = mesh.verts.shape[0]

            @fn_timer
            def genetateBSplineBaseVals(tbounds, gridDim, verts):
                '''
                compute the BSpline values of each BSpline basic function on all vertex in verts
                :param tbounds:
                :param z:
                :param verts:
                :return:
                '''
                SBverts = verts.copy()
                k = 3
                shape = gridDim[0:3]  # the number of control points
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

                print 'tx', tx
                print 'ty', ty
                print 'tz', tz

                shape = gridDim[0:3]  # the number of control points
                Bx = []
                BasicX = []
                for ii in range(shape[0]):
                    c = np.zeros((shape[0],), dtype=np.float32)
                    c[ii] = 1.
                    spl = BSpline(tx, c, k)
                    bvals = map(spl, SBverts.transpose()[0])
                    Bx.append(bvals)
                    BasicX.append(spl)
                Bx = np.array(Bx, dtype=np.float32)

                By = []
                BasicY = []
                for ii in range(shape[1]):
                    c = np.zeros((shape[1],), dtype=np.float32)
                    c[ii] = 1.
                    spl = BSpline(ty, c, k)
                    bvals = map(spl, SBverts.transpose()[1])
                    By.append(bvals)
                    BasicY.append(spl)
                By = np.array(By, dtype=np.float32)

                Bz = []
                BasicZ = []
                for ii in range(shape[2]):
                    c = np.zeros((shape[2],), dtype=np.float32)
                    c[ii] = 1.
                    spl = BSpline(tz, c, k)
                    bvals = map(spl, SBverts.transpose()[2])
                    Bz.append(bvals)
                    BasicZ.append(spl)
                Bz = np.array(Bz, dtype=np.float32)
                # print Bx.shape,By.shape,Bz.shape

                return tx, ty, tz, Bx, By, Bz, BasicX, BasicY, BasicZ

            # Note: the s_{0}^p is unchanged during the optimization, we can compute the Bx,By,Bz once and use them many times
            tx, ty, tz, Bx, By, Bz, BasicX, BasicY, BasicZ = genetateBSplineBaseVals(tbounds, gridDim,
                                                                                     S_B._verts.copy())

            # for ii in range(gridDim[0]):
            #     for jj in range(gridDim[1]):
            #         for kk in range(gridDim[2]):
            #             z = Bx[ii] * By[jj] * Bz[kk]
            # print 'ii,jj,kk:',ii,jj,kk
            # print 'z[0]:',z[0]
            # @fn_timer
            def updateSB(z, verts):
                '''
                to compute s_b according to the equation (4). verts is s_0
                this is a time consuming operation if len(verts) is too large
                Note: this will change verts
                :param tbounds:
                :param z:
                :param verts:
                :return:
                '''
                vv = np.zeros(shape=(3, verts.shape[0]), dtype=np.float32)
                z = z.reshape(gridDim)
                for ii in range(gridDim[0]):
                    for jj in range(gridDim[1]):
                        for kk in range(gridDim[2]):
                            vv += (z[ii][jj][kk]).reshape(3, 1) * (Bx[ii] * By[jj] * Bz[kk])
                vs = verts + vv.transpose()
                return np.array(vs)

            def E_I_B(SBVerts):
                verts_SB = SBVerts
                verts_SI = S_I.verts  #

                verts_BI = verts_SB - verts_SI
                e = np.sum(verts_BI * verts_BI) / 2.

                return e

            def E_SSM(SBVerts):
                verts_SB = SBVerts
                verts_SSM = S_SSM.verts  #

                verts_BSSM = verts_SB - verts_SSM
                e = np.sum(verts_BSSM * verts_BSSM) / 2.
                return e

            # @fn_timer
            def G_I_B_n_SSM(z, verts, refverts):
                SIverts = refverts.copy()

                verts = updateSB(z, verts)  # S_B
                verts = verts - SIverts
                g = []
                for ii in range(gridDim[0]):
                    for jj in range(gridDim[1]):
                        for kk in range(gridDim[2]):
                            vv = verts.transpose() * (Bx[ii] * By[jj] * Bz[kk])
                            vv = vv.sum(axis=1)
                            g.append(vv[0])
                            g.append(vv[1])
                            g.append(vv[2])
                            # print 'when computing gradient'
                            # print 'ii,jj,kk:', ii, jj, kk
                            # print 'vv[0]:', vv[0]
                return np.array(g).reshape((numcoords,))

            def E_SIM(SBVerts):
                verts_SB = SBVerts

                e = 0
                return e

            def G_SIM(z, verts):
                return np.zeros(shape=(z.shape[0] * z.shape[1] * z.shape[2] * 3,))

            # def TrivariateBSpline(tbounds,z,s,k=3):# z.ndim should be 4. the first 3 are the shape of the control vertices,
            #     # ,and shape[3]  should be 3 to contain the coordinate of the vertices. Note, if the number of voxel is x,
            #     # then the number of control points should be x + k - 1
            #     shape = z.shape[0:3] # the number of control points
            #     shape = np.array(shape)
            #     shape -= [k - 1,k - 1,k - 1] # the number of voxels in each dimension
            #     tx = np.zeros(shape=(shape[0] + 2 * k,))
            #     tx[k:-k] = np.linspace(0.,tbounds[0],shape[0])
            #     tx[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[0])
            #     # print 'tx:',tx
            #     ty = np.zeros(shape=(shape[1] + 2 * k,))
            #     ty[k:-k] = np.linspace(0.,tbounds[1],shape[1])
            #     ty[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[1])
            #     # print 'ty:',ty
            #     tz = np.zeros(shape=(shape[2] + 2 * k,))
            #     tz[k:-k] = np.linspace(0.,tbounds[2],shape[2])
            #     tz[-1:-(k + 1):-1] = np.ones(shape=(k,)) * (tbounds[2])
            #     # print 'tz:',tz
            #     spl = None
            #     shape = z.shape[0:3]
            #     vv = []
            #     s = s.transpose()
            #     for ii in range(shape[0]):
            #         v = []
            #         for jj in range(shape[1]):
            #             spl = BSpline(tz,z[ii][jj],k)
            #             val = spl(s[2])
            #             v.append(val)
            #         v = np.array(v)
            #         spl = BSpline(ty,v,k)
            #         val = spl(s[1])
            #         vv.append(spl(s[1]))
            #     vv = np.array(vv)
            #     spl = BSpline(tx,vv,k)
            #     spl = spl(s[0])
            #     return spl
            # z = np.arange(0,5 * 100 * 3)
            # z = z.reshape((5,10,10,3))
            # tbounds = [1024,512,512]
            # spl = TrivariageBSpline(tbounds,z,[2,7,0])
            # print 'spline %s' % (spl)

            omega_IB = .5
            omega_SSM = .5
            omega_SIM = .0

            def fun(z, SBVerts):
                verts = SBVerts.reshape((nv, 3)).copy()
                z = z.reshape(gridDim)
                verts = updateSB(z, verts)

                eIB = omega_IB * E_I_B(verts)
                # eSSM = omega_SSM * E_SSM(verts)
                # eSIM = omega_SIM * E_SIM(verts)
                # e = eIB + eSSM + eSIM
                return eIB

            def gradient(z, SBVerts):
                verts = SBVerts.reshape((nv, 3)).copy()
                z = z.reshape(gridDim)

                gIB = omega_IB * G_I_B_n_SSM(z, verts.copy(), S_I.verts)
                # gSSM = omega_SSM * G_I_B_n_SSM(z,verts.copy(),S_SSM.verts)
                # gSIM = omega_SIM * G_SIM(z,verts.copy())
                # g = gIB + gSSM + gSIM

                return gIB

            if True:
                print 'check gradient in SB optimization....'

                def compGrad(xk, f, epsilon):
                    f0 = f(xk, S_B.verts)
                    grad = np.zeros((len(xk),), float)
                    ei = np.zeros((len(xk),), float)
                    ii = 0
                    jj = 0
                    kk = 0

                    def compUpdate(z):
                        verts = S_B.verts.copy()
                        vv = np.zeros(shape=(3, verts.shape[0]), dtype=np.float32)
                        z = z.reshape(gridDim)
                        for ii in range(gridDim[0]):
                            for jj in range(gridDim[1]):
                                for kk in range(gridDim[2]):
                                    if ii == 1 and jj == 0 and kk == 0:
                                        print 'Bx*By*Bz', (Bx[ii] * By[jj] * Bz[kk])[0]
                                        print 'z[ii][jj][kk]', z[ii][jj][kk][0]
                                    vv += (z[ii][jj][kk]).reshape(3, 1) * (Bx[ii] * By[jj] * Bz[kk])
                        vs = verts + vv.transpose()
                        print 'vs', vs
                        return vs

                    for k in range(len(xk)):
                        ei[k] = 1.0
                        d = epsilon * ei
                        f1 = f(xk + d, S_B.verts)
                        grad[k] = (f1 - f0) / d[k]
                        ei[k] = 0.0

                        if k % 3 == 0:
                            print 'compGrad'
                            print 'ii,jj,kk', ii, jj, kk
                            print 'grad[k]', grad[k]
                            if ii == 1 and jj == 0 and kk == 0:
                                print 'd', d
                                print 'f0', f0
                                print 'f1', f1
                                print 'update for xk'
                                vsxk = compUpdate(xk)
                                v1 = ((vsxk - S_I.verts) ** 2).sum()
                                print 'update for xk + d'
                                vsxkd = compUpdate(xk + d)
                                v2 = ((vsxkd - S_I.verts) ** 2).sum()
                                print v1, v2, v1 - v2
                                print ((vsxkd - vsxk) ** 2).sum()
                            kk += 1
                            if kk == gridDim[2]:
                                kk = 0
                                jj += 1
                                if jj == gridDim[1]:
                                    jj = 0
                                    ii += 1

                    return grad

                epsilon = 1.  # 0.001 #np.sqrt(np.finfo(float).eps)
                zf = np.random.random(numcoords) * 1.
                # grad = compGrad(zf, fun,epsilon)
                # graderr = np.sqrt(np.sum((gradient(zf, S_B.verts) - grad) ** 2))
                graderr = opt.check_grad(fun, gradient, zf, S_B.verts, **{'epsilon': epsilon})
                print 'graderr when optimizing SB is:', graderr

                # verts = S_B.verts.transpose()
                # plt.figure(1)
                # plt.plot(verts[0],'r')
                # plt.plot(verts[1],'g')
                # plt.plot(verts[2],'b')
                # plt.legend(('$v_{kx}$','$v_{ky}$','$v_{kz}$'))
                # plt.title('mesh vertex coordinates')

                # plt.figure(2)
                # plt.title('BSpline Base Funcs along X')
                # style = ['r.','r:','r-','r-.','r--']
                # for i in range(gridDim[0]):
                #     x = map(BasicX[i],range(tbounds[0]))
                #     plt.plot(x,style[i])
                # plt.legend(('$B_0$','$B_1$','$B_2$','$B_3$','$B_4$'))

                # plt.figure(3)
                # plt.title('BSpline Base Funcs along y')
                # style = ['g.', 'g:', 'g-', 'g-.', 'g--']
                # for i in range(gridDim[1]):
                #     y = map(BasicY[i], range(tbounds[1]))
                #     plt.plot(y, style[i])
                # plt.legend(('$B_0$', '$B_1$', '$B_2$', '$B_3$', '$B_4$'))
                #
                # plt.figure(4)
                # plt.title('BSpline Base Funcs along Z')
                # style = ['b.', 'b:', 'b-', 'b-.', 'b--']
                # for i in range(gridDim[2]):
                #     x = map(BasicZ[i], range(tbounds[2]))
                #     plt.plot(x, style[i])
                # plt.legend(('$B_0$', '$B_1$', '$B_2$', '$B_3$', '$B_4$'))

                # plt.figure(5)
                # plt.subplot(321)
                # plt.plot(Bx[0],'r.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bx_{0}$')
                #
                # plt.subplot(322)
                # plt.plot(Bx[1], 'r:')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bx_{1}$')
                #
                # plt.subplot(323)
                # plt.plot(Bx[2], 'r-')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bx_{2}$')
                #
                # plt.subplot(324)
                # plt.plot(Bx[3], 'r-.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bx_{3}$')
                #
                # plt.subplot(325)
                # plt.plot(Bx[4], 'r--')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bx_{4}$')

                # plt.figure(6)
                # plt.plot(Bx[0], 'r.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bx_{0}$')

                # plt.figure(7)
                # plt.subplot(321)
                # plt.plot(By[0], 'g.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$By_{0}$')
                #
                # plt.subplot(322)
                # plt.plot(By[1], 'g:')
                # plt.title('Vals of BSpline Base func at mesh vertices:$By_{1}$')
                #
                # plt.subplot(323)
                # plt.plot(By[2], 'g-')
                # plt.title('Vals of BSpline Base func at mesh vertices:$By_{2}$')
                #
                # plt.subplot(324)
                # plt.plot(By[3], 'g-.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$By_{3}$')
                #
                # plt.subplot(325)
                # plt.plot(By[4], 'g--')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bz_{4}$')

                # plt.figure(8)
                # plt.plot(By[0], 'g.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$By_{0}$')

                # plt.figure(9)
                # plt.subplot(321)
                # plt.plot(Bz[0], 'b.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bz_{0}$')
                #
                # plt.subplot(322)
                # plt.plot(Bz[1], 'b:')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bz_{1}$')
                #
                # plt.subplot(323)
                # plt.plot(Bz[2], 'b-')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bz_{2}$')
                #
                # plt.subplot(324)
                # plt.plot(Bz[3], 'b-.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bz_{3}$')
                #
                # plt.subplot(325)
                # plt.plot(Bz[4], 'b--')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bz_{4}$')

                # plt.figure(10)
                # plt.plot(Bz[0], 'b.')
                # plt.title('Vals of BSpline Base func at mesh vertices:$Bz_{0}$')

                plt.show()
                return

            def callback(z):
                '''
                callback only accept one arguments: the variables
                :param z:
                :return:
                '''
                verts = S_B.verts.copy()
                z = z.reshape(gridDim)
                verts = updateSB(z, verts)
                global km
                openmesh.writeOff(dirname + '/%s_registered_SB.off' % km, verts, mesh.faces)
                km += 1

            res = opt.minimize(fun, np.random.random(numcoords), \
                               args=(S_B.verts.copy()), method='l-bfgs-b', jac=gradient,
                               callback=callback, options={'maxiter': 100, 'disp': True})
            print 'x', res['x'].reshape(gridDim)
            print 'success?', res['success']
            print 'message:', res['message']
            S_B.verts = updateSB(res['x'].reshape(gridDim), S_B.verts)
            return S_B.verts.copy()

        # repeat optimization until converge
        for k in range(1):
            SIverts = None
            SSMverts = None
            SBverts = None

            # print '------------> Optimizing SI <-------------'
            # SIverts = OptimizeSI() # will update S_I.verts
            # print '------------> Optimizing SSM <-------------'
            # SSMverts = OptimizeSSM_ICP() # will update S_SSM.verts
            print '------------> Optimizing SB <-------------'
            SBverts = OptimizeSB()  # will update S_B.verts

        return SIverts, SSMverts, SBverts
