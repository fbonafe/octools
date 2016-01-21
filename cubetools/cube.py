"""
These functions extract all the useful information from the Gaussian CUBE file
and analize it.
"""

import sys
import numpy as np
from matplotlib.mlab import griddata

def readCube(filename):
    class CubeFile:
        def __init__(self,coordsx,coordsy,coordsz,x0,y0,z0,pointsx,pointsy,pointsz,voldx,voldy,voldz,isovals):
            self.x0=x0
            self.y0=y0
            self.z0=z0
            self.coords_x=coordsx
            self.coords_y=coordsy
            self.coords_z=coordsz
            self.nx=pointsx
            self.ny=pointsy
            self.nz=pointsz
            self.dx=voldx
            self.dy=voldy
            self.dz=voldz
            self.dv=voldx*voldy*voldz
            self.isovals=isovals
         
        def __add__(self,other): #untested
            result = self
            result.isovals = [sum(x) for x in zip(self.isovals,other.isovals)]
            return result
        
        def __add__(self,other): #untested
            result = self
            result.isovals = [x-y for x,y in zip(self.isovals,other.isovals)]
            return result
        
    
    f = open(filename, "r")
    bohr_to_angst = 0.529177
    #skipping comment lines
    f.readline()
    f.readline()
    
    #reading the number of atoms and the origin of the cube
    l = f.readline().split()
    n_atoms = int(l[0])
    x_origin = float(l[1])
    y_origin = float(l[2])
    z_origin = float(l[3])
    
    #reading number of volume elements and their volume
    cubic_box = True
    l = f.readline().split()
    points_x = int(l[0])
    vol_dx = float(l[1])
    if float(l[2])!=0.0 or float(l[3])!= 0.0:
        cubic_box = False
    l = f.readline().split()
    points_y = int(l[0])
    vol_dy = float(l[2])
    if float(l[1])!=0.0 or float(l[3])!= 0.0:
        cubic_box = False
    l = f.readline().split()
    points_z = int(l[0])
    vol_dz = float(l[3])
    if float(l[1])!=0.0 or float(l[2])!= 0.0:
        cubic_box = False
    if cubic_box == False:
        print "Non-cubic box, cannot continue!"
        sys.exit()
#    volume_element = vol_dx * vol_dy * vol_dz
    
    #reading atomic coordinates
    coords_x = []
    coords_y = []
    coords_z = []
    for i in range(n_atoms):
        co = f.readline().split()
        coords_x.append(bohr_to_angst * float(co[2]))
        coords_y.append(bohr_to_angst * float(co[3]))
        coords_z.append(bohr_to_angst * float(co[4]))
    
    #memory consuming but easy
    isovalues = [] 
    for line in f:
        spl = line.split()
        for v in spl:
            isovalues.append(float(v))
            
    return CubeFile(coords_x,coords_y,coords_z,x_origin,y_origin,z_origin,points_x,points_y,points_z,vol_dx,vol_dy,vol_dz,isovalues)
    
def genGrid(filename):
    class Grid:
        def __init__(self,x,y,z,isoval,x0,y0,z0,pointsx,pointsy,pointsz,voldx,voldy,voldz):    
            self.x=x
            self.y=y
            self.z=z
            self.isovals=isoval
            self.x0=x0
            self.y0=y0
            self.z0=z0
            self.nx=pointsx
            self.ny=pointsy
            self.nz=pointsz
            self.dx=voldx
            self.dy=voldy
            self.dz=voldz
            self.dv=voldx*voldy*voldz
            
        def __add__(self,other): #untested
            result = self
            result.isovals = self.isovals + other.isovals
            return result
        
        def __sub__(self,other): #untested
            if self.isovals.size == other.isovals.size:
                result = self
                result.isovals = self.isovals - other.isovals
                return result
            else:
                print 'Error: not equal number of isovalues, keeping first variable as output' 
                return self
    
    CubeFile = readCube(filename)
    Coords = np.array([CubeFile.coords_x,CubeFile.coords_y,CubeFile.coords_z])
    
    bohr_to_angst = 0.529177
    xs = []
    ys = []
    zs = []
    for ix in range(CubeFile.nx):
        for iy in range(CubeFile.ny):
            for iz in range(CubeFile.nz):
                x = bohr_to_angst * (CubeFile.x0 + ix * CubeFile.dx)
                y = bohr_to_angst * (CubeFile.y0 + iy * CubeFile.dy)
                z = bohr_to_angst * (CubeFile.z0 + iz * CubeFile.dz)
                xs.append(x)
                ys.append(y)
                zs.append(z)
    return Grid(np.array(xs),np.array(ys),np.array(zs),np.array(CubeFile.isovals),CubeFile.x0,CubeFile.y0,CubeFile.z0, \
               CubeFile.nx,CubeFile.ny,CubeFile.nz,CubeFile.dx,CubeFile.dy,CubeFile.dz),Coords

def intZ(CubeFile):
    class Grid:
        def __init__(self,x,y,isoval):    
            self.x=x
            self.y=y
            self.isoval=isoval
            
    bohr_to_angst = 0.529177
    isovalues=CubeFile.isovals
#    isovalues.reverse()
    xs = []
    ys = []
    dens = []
    i = 0
    for ix in range(CubeFile.nx):
        for iy in range(CubeFile.ny):
            x = bohr_to_angst * (CubeFile.x0 + ix * CubeFile.dx)
            y = bohr_to_angst * (CubeFile.y0 + iy * CubeFile.dy)
            intdens=0.
            for iz in range(CubeFile.nz):
#                val = isovalues.pop()
#                intdens += val
                intdens += isovalues[i]
                i += 1
            intdens=intdens*CubeFile.dz
            xs.append(x)
            ys.append(y)
            dens.append(intdens)
#    CubeFile.isovals = aux
    return Grid(np.array(xs),np.array(ys),np.array(dens))

def transDipole(CubeFile1,CubeFile2):

    bohr_to_angst = 0.529177
    iso1 = CubeFile1.isovals
    iso1.reverse()
    iso2 = CubeFile2.isovals
    iso2.reverse()
    transDipArr = []
    for ix in range(CubeFile1.nx):
        for iy in range(CubeFile1.ny):
            for iz in range(CubeFile1.nz):
                x = bohr_to_angst * (CubeFile1.x0 + ix * CubeFile1.dx)
                y = bohr_to_angst * (CubeFile1.y0 + iy * CubeFile1.dy)
                z = bohr_to_angst * (CubeFile1.z0 + iz * CubeFile1.dz)
                prod = iso1.pop() * iso2.pop()
                transDip = prod**2 * (x**2 + y**2 + z**2)
                transDipArr.append(transDip)
    transDipArr.reverse()
    return transDipArr
    
def writeDens(a,fileout):
    fmt='15.8e'
    asTxtLst = map(lambda val: format(val, fmt), a)
    size = len(asTxtLst)
    ncol = 4
    f = open(fileout, 'w+')
    for i in range(size/ncol):
        f.write(' '.join(asTxtLst[i*ncol:(i+1)*ncol])+'\n')
    if ncol*(size/ncol) < size:
        f.write(' '.join(asTxtLst[ncol*(size/ncol):])+'\n')
    f.close()
    return 0
