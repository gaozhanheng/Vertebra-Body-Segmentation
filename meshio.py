#!/usr/bin/env python
#encoding: utf-8
'''
to read and write 3d mesh file with format: off

'''

import re # regular expression for processing string

retrieveNumberRule = r'-*\d+\.?\d*'
def readOff(filename,quiet=True):
    with open(filename,'r') as file:
        if not quiet:
            print 'reading file: %s' % filename
        allLines = file.readlines()
        if not allLines[0].strip().lower() == 'off':
            print 'current file is not an off file: %s' % filename
            file.close()
            return None
        numberStr = re.findall(retrieveNumberRule,allLines[1].strip())
        numberEles = []
        for s in numberStr:
            numberEles.append(int(s))
        nv = numberEles[0]
        nf = numberEles[1]
        if not quiet:
            print 'number of vertices: %d\nnumber of faces: %d\nnumber of edges: %d\n' % (nv,nf,numberEles[2])

        # assume there is no space line in the file
        varray = []
        for k in range(0,nv):
            vline = re.findall(retrieveNumberRule,allLines[k + 2].strip())
            varray.append([float(x) for x in vline])
        if not quiet:
            print '%d vertices read' % len(varray)

        farray = []
        for k in range(0,nf):
            fline = re.findall(retrieveNumberRule,allLines[k + 2 + nv].strip())
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




