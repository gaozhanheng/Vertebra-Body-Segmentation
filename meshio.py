#!/usr/bin/env python
#encoding: utf-8
'''
to read and write 3d mesh file with format: off

'''

import os
import re # regular expression for processing string

retrieveNumberRule = r'-*\d+\.?\d*'
def readPly(filename,quiet=True):
    pass
def writePly(filename,quiet=True):
    pass
def readOff(filename,quiet=True):
    with open(filename,'r') as file:
        if not quiet:
            print 'reading file: %s' % filename
        allLines = file.readlines()
        if not allLines[0].strip().lower() == 'off':
            print 'current file is not an off file: %s' % filename
            file.close()
            return None
        # numberStr = re.findall(retrieveNumberRule,allLines[1].strip())
        numberStr = allLines[1].split()
        numberEles = []
        for s in numberStr:
            if len(s) == 0:
                continue
            numberEles.append(int(s))
        nv = numberEles[0]
        nf = numberEles[1]
        if not quiet:
            print 'number of vertices: %d\nnumber of faces: %d\nnumber of edges: %d\n' % (nv,nf,numberEles[2])

        # assume there is no space line in the file
        varray = []
        for k in range(0,nv):
            # vline = re.findall(retrieveNumberRule,allLines[k + 2].strip())
            vline = allLines[k + 2].strip().split()
            varray.append([float(x) for x in vline] )
        if not quiet:
            print '%d vertices read' % len(varray)

        farray = []
        for k in range(0,nf):
            # fline = re.findall(retrieveNumberRule,allLines[k + 2 + nv].strip())
            fline = allLines[k + 2 + nv].split()
            farray.append([int(x) for x in fline[1:4]])
        if not quiet:
            print '%d faces read' % len(farray)
        return varray,farray

def writeOff(filename,varray,farray,quiet=True):
    with open(filename,'w') as file:
        allLines = []
        allLines.append('OFF\n')
        allLines.append('%d %d %d\n' % (len(varray), len(farray), 0))

        for k in range(0,len(varray)):
            allLines.append('%.8f %.8f %.8f\n' % (varray[k][0],varray[k][1],varray[k][2]))
        for k in range(0,len(farray)):
            allLines.append('3 %d %d %d\n' % (farray[k][0],farray[k][1],farray[k][2]))

        file.writelines(allLines)
        if not quiet:
            print 'finished writting off mesh into %s' % filename



def operateMesh(fn):
    # eval is useful
    readMode = r'^read.*'
    writeMode = r'^write.*'
    name = fn.__name__
    mr = re.findall(readMode,name)
    mw = re.findall(writeMode,name)
    if not mr and not mw:
        raise ValueError('Mesh operation functions must start with \'read\' or \'write\': %s' % name)
    def operate(filename,quiet=True):
        return fn(filename,quiet)
    return operate

@operateMesh
def readMesh(filename,quiet=True):
    name, ext = os.path.splitext(filename)
    if '.ply' == ext:
        readPly(filename, quiet)
    elif '.off' == ext:
        readOff(filename, quiet)
    else:
        raise ValueError('Meshio does not suppot reading mesh in %s format' % ext)
@operateMesh
def writeMesh(filename,quiet=True):
    name, ext = os.path.splitext(filename)
    if '.ply' == ext:
        writePly(filename, quiet)
    elif '.off' == ext:
        writeOff(filename, quiet)
    else:
        raise ValueError('Meshio does not suppot writting mesh in %s format' % ext)
