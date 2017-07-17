

from . import SegAlgo

class SIM(SegAlgo):
    def __init__(self):
        pass

    @fn_timer
    def Optimize(self,mesh,volume):
        ''' Optimize the mesh.'''
        S_I = copy.copy(mesh)
        S_B = copy.copy(mesh)
        S_SSM = copy.copy(mesh)

        global km
        km = 1
        dirname = '../Data/tmp'
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.mkdir(dirname)

        def OptimizeSI():  # more details about scipy.optimization please refer to the optimize __init__.py
            nv = S_I.nv  # the number of vertices used to test the optimization efficiency

            # print 'number of vertices used:',nv

            # def showMeshAndNormal():
            #     '''
            #     To test if the normals is pointing from inside to outlide and if the ordering of the normals is right.
            #     :return:
            #     '''
            #     m = visuals.Mesh(S_B.verts, mesh.faces)
            #     l = visuals.Line()
            #     view.add(m)
            #     outline_p = []
            #     k = 0
            #     for v in S_B.verts:
            #         outline_p.append(v)
            #         outline_p.append(v + normals[k])
            #         k += 1
            #
            #     outline_p = np.array(outline_p)
            #
            #     outline = visuals.Line(pos=outline_p, color=(0.1, 0.1, 0.1, 1.), width=1,
            #                            connect='segments', method='gl', antialias=False)
            #     view.add(outline)
            #     vispy.app.run()
            #
            # # showMeshAndNormal()

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

            refMesh = S_SSM  # original is S_B

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
            nv = S_SSM.verts.shape[0]

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

            def E_SIM(SBVerts):
                verts_SB = SBVerts

                e = 0
                return e

            def TrivariateBSpline(tbounds, z, s,
                                  k=3):  # z.ndim should be 4. the first 3 are the shape of the control vertices,
                # ,and shape[3]  should be 3 to contain the coordinate of the vertices. Note, if the number of voxel is x,
                # then the number of control points should be x + k - 1
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
                # print 'tz:',tz
                spl = None
                shape = z.shape[0:3]
                vv = []
                s = s.transpose()
                for ii in range(shape[0]):
                    v = []
                    for jj in range(shape[1]):
                        spl = BSpline(tz, z[ii][jj], k)
                        val = spl(s[2])
                        v.append(val)
                    v = np.array(v)
                    spl = BSpline(ty, v, k)
                    val = spl(s[1])
                    vv.append(spl(s[1]))
                vv = np.array(vv)
                spl = BSpline(tx, vv, k)
                spl = spl(s[0])
                return spl

            # z = np.arange(0,5 * 100 * 3)
            # z = z.reshape((5,10,10,3))
            # tbounds = [1024,512,512]
            # spl = TrivariageBSpline(tbounds,z,[2,7,0])
            # print 'spline %s' % (spl)
            def updateSB(tbounds, z, verts):
                '''
                this is a time consuming operation if len(verts) is too large
                Note: this will change verts
                :param tbounds:
                :param z:
                :param verts:
                :return:
                '''
                vs = []
                for v in verts:
                    v += TrivariateBSpline(tbounds, z, v)
                    vs.append(v)
                return np.array(vs)

            gridDim = (5, 5, 10, 3)
            omega_IB = .5
            omega_SSM = .5
            omega_SIM = .0

            def fun(z, tbounds, SBVerts):
                verts = SBVerts.reshape((nv, 3)).copy()
                z = z.reshape(gridDim)
                verts = updateSB(tbounds, z, verts)

                eIB = omega_IB * E_I_B(verts)
                eSSM = omega_SSM * E_SSM(verts)
                eSIM = omega_SIM * E_SIM(verts)
                e = eIB + eSSM + eSIM

                return e

            def G_I_B_n_SSM(z, tbounds, verts, refverts):
                SBverts = verts.copy()
                SIverts = refverts.copy()

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

                verts = updateSB(tbounds, z, verts)  # S_B
                g = []
                for ii in range(shape[0]):
                    for jj in range(shape[1]):
                        for kk in range(shape[2]):
                            gx = 0.
                            gy = 0.
                            gz = 0.
                            kp = 0
                            for v in verts:
                                vb = SBverts[kp]
                                kg = (v - SIverts) * Bx[ii][kp] * By[jj][kp] * Bz[kk][kp]
                                gx += kg[0]
                                gy += kg[1]
                                gz += kg[2]
                                kp += 1

                            g.append(gx)
                            g.append(gy)
                            g.append(gz)
                return np.array(g)

            def G_SIM(z, tbounds, verts):
                return np.zeros(shape=(z.shape[0] * z.shape[1] * z.shape[2],))

            def gradient(z, tbounds, SBVerts):
                verts = SBVerts.reshape((nv, 3)).copy()
                z = z.reshape(gridDim)

                gIB = omega_IB * G_I_B_n_SSM(z, tbounds, verts.copy(), S_I.verts)
                gSSM = omega_SSM * G_I_B_n_SSM(z, tbounds, verts.copy(), S_SSM.verts)
                gSIM = omega_SIM * G_SIM(z, tbounds, verts.copy())
                g = gIB + gSSM + gSIM

                return g

            tbounds = np.array(volume.pixelArray.shape)

            if True:
                print 'check gradient in SB optimization....'
                graderr = opt.check_grad(fun, gradient, np.ones(gridDim), tbounds, S_B.verts.copy())
                print 'graderr when optimizing SB is:', graderr

            nv = mesh.verts.shape[0]

            def callback(z, tbounds, verts):
                verts = verts.reshape((nv, 3)).copy()
                z = z.reshape(gridDim)
                verts = updateSB(tbounds, z, verts)
                global km
                openmesh.writeOff(dirname + '/%s_registered_SB.off' % km, verts, mesh.faces)
                km += 1

            res = opt.minimize(fun, np.ones(shape=gridDim), args=(tbounds, S_B.verts.copy()), method='l-bfgs-b',
                               jac=gradient,
                               callback=callback, options={'maxiter': 5, 'disp': True})

            S_B.verts = updateSB(tbounds, res['x'], S_B.verts)
            return S_B.verts.copy()

        # repeat optimization until converge
        for k in range(10):
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