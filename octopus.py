
import glob
import numpy as np
from scipy import interpolate
from cubetools.cube import *
from cubetools.plotters import *


hbar = 0.658
dens0 = np.empty(0)
outevery = 50 # In general

class dynamics:
    def __init__(self,time,field,dipole,forces,torque,alfa,omega,theta,charge,energy,mudot,muddot,torquepersite,muabs,muphi):
        self.time = time
        self.field = field
        self.mu = dipole
        self.forces = forces
        self.torque = torque
        self.alfa = alfa
        self.omega = omega
        self.theta = theta
        self.charge = charge
        self.energy = energy
        self.mudot = mudot
        self.muddot = muddot
        self.torquepersite = torquepersite
        self.muabs = muabs
        self.muphi = muphi

def Interpol(datax,datay,npoints):
    f1 = interpolate.interp1d(datax,datay,kind='linear')
    xnew = np.linspace(datax.min(),datax.max(),npoints)
    ynew = f1(xnew)
    return xnew,ynew


# This loads all the relevant data
def loadData(out_every,prepath='./'):
    
    time,field,mu,charge = dipoleFieldVec()
    dt = time[1]-time[0]
    npoints = time.shape[0]
    outevery = out_every 
    
    force,torque,torquePerSite,alfa = forces(prepath,time,dt*out_every)
    
    muabs = np.sqrt(mu[0]**2+mu[1]**2+mu[2]**2)
    
    muphi = np.angle(mu[0]+mu[1]*1j)
    #muphi = np.array([x+2.*np.pi for x in muphi if x < 0])
    for i in range(1,muphi.shape[0]):
        if muphi[i] < muphi[i-1]:
            muphi[i:] = 2*np.pi + muphi[i:]
        
    mudot = muDot(time,muabs)
    muddot = muDDot(time,muabs)
    omega,theta = alphaOmegaTheta(time,alfa)
    energ = energy()
    
    return dynamics(time,field,mu,force,torque,alfa,omega,theta,charge,energ,mudot,muddot,torquePerSite,muabs,muphi)
   

#################################################################
# This quantities are NOT automatically calculated by loadData,
# should be used separately.

def time():
    t = []
    path = './td.general/'
    with open(path+'multipoles','r') as multipoles:
        for line in multipoles:
            sl = line.strip()
            if not sl.startswith("#"):
                l = line.split()
                t.append(float(l[1]))
    t = np.array(t)*hbar
    return t # time in fs
                
                            
def ePot(efield,mu):
    epot = mu[0]*efield[0]+mu[1]*efield[1]+mu[2]*efield[2]
    return epot

def torqueNuc(efield,mu):
    torque = efield[0]*mu[1]-mu[0]*efield[1]
    return torque
#################################################################

def energy():
    path = './td.general/'
    e = []
    with open(path+'energy','r') as laser:
        for line in laser:
            sl = line.strip()
            if not sl.startswith("#"):
                l = line.split()
                e.append(float(l[2]))
    e = np.array(e)
    return e

#Dipole Moment
def dipoleFieldVec():
    path = './td.general/'
    ex, ey, ez, mux, muy, muz, time, charge = [], [], [], [], [], [], [], []
    
    with open(path+'multipoles','r') as multipoles:
        for line in multipoles:
            sl = line.strip()
            if not sl.startswith("#"):
                l = line.split()
                time.append(float(l[1]))
                charge.append(float(l[2]))
                mux.append(float(l[3]))
                muy.append(float(l[4]))
                muz.append(float(l[5]))
                    
    with open(path+'laser','r') as laser:
        for line in laser:
            sl = line.strip()
            if not sl.startswith("#"):
                l = line.split()
                ex.append(float(l[2]))
                ey.append(float(l[3]))
                ez.append(float(l[4]))
    
    time = np.array(time)*hbar
    charge = np.array(charge)
    efield = np.array([ex[:len(time)],ey[:len(time)],ez[:len(time)]])
    mu = np.array([mux,muy,muz])
    
    return time,efield,mu,charge


#DipoleMomentDervivative

def muDot(time,mu):
    tck,uout = interpolate.splprep([time,mu],s=0.,k=2,per=False)
    dx,dy = interpolate.splev(uout,tck,der=1)
    mudot = dy/dx
    return mudot
    
def muDDot(time,mu):
    tck,uout = interpolate.splprep([time,mu],s=0.,k=2,per=False)
    dx,dy = interpolate.splev(uout,tck,der=1)
    mudot = dy/dx
    tck,uout = interpolate.splprep([time,mudot],s=0.,k=2,per=False)
    ddx,ddy = interpolate.splev(uout,tck,der=1)
    muddot = ddy/ddx
    return muddot
    
#Force

def forces(path,realtime,dt): # dt in hbar/eV
    forcefiles = sorted(glob.glob(path+'td.0*/forces.xsf')) # VERY IMPORTANT TO INCLUDE sorted
    natoms = sum(1 for line in open(path+'td.0000000/forces.xsf')) - 1
    el, idn, x, y, z, fx0, fy0, fz0, fx, fy, fz, torque = [[None]*natoms for i in range(12)]
    nsteps = len(forcefiles)
    
    forcePerStep, torquePerStep, time = np.zeros((nsteps)), np.zeros((nsteps)), np.zeros((nsteps))
    torquePerSite = np.zeros((nsteps,natoms))
    

    with open(forcefiles[0],'r+') as f:
        f.readline()
        for j,line in enumerate(f):
            l = line.split()
            el[j] = l[0]
            if el[j] == 'C':
                idn[j] = 6
            elif el[j] == 'H':
                idn[j] = 1
            elif el[j] == 'Na':
                idn[j] = 9
            elif el[j] == 'Li':
                idn[j] = 3
            else:
                print 'Element mass not known'
            x[j] = float(l[1])
            y[j] = float(l[2])
            z[j] = float(l[3])
            fx0[j] = float(l[4]) 
            fy0[j] = float(l[5])
            fz0[j] = float(l[6])        
    
    # Create AXSF file to be read with XCRYSDEN
    g = open('netforces.axsf','w')
    g.write('ANIMSTEPS %i \n' % len(forcefiles))

    for i,ffile in enumerate(forcefiles):
        time[i] = i*dt
        g.write('ATOMS %i \n' % (i+1))
        with open(ffile,'r+') as f:
            f.readline()
            rForce = [0.,0.,0.]
            totalTorque = 0.
            for j,line in enumerate(f):
                l = line.split()
                fx[j] = float(l[4]) - fx0[j]
                fy[j] = float(l[5]) - fy0[j]
                fz[j] = float(l[6]) - fz0[j]   
                # Write to AXSF file
                g.write("%i %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n" % (idn[j],x[j],y[j],z[j],fx[j],fy[j],fz[j]))
                # Calculate net force
                rForce += [fx[j],fy[j],fz[j]]
                # Calculate torque
                torque[j] = x[j]*fy[j]-y[j]*fx[j]
                torquePerSite[i,j] = torque[j]
                totalTorque += torque[j]
            
            totalForce = np.sqrt(rForce[0]**2+rForce[1]**2+rForce[2]**2)
            forcePerStep[i] = totalForce
            torquePerStep[i] = totalTorque
            
    g.close()
    mInertia = momentInertia(natoms,el,x,y)
    
    npoints = realtime.shape[0]
    time2,force = Interpol(time,forcePerStep,npoints)
    time2,torque = Interpol(time,torquePerStep,npoints)
    alfa = torque/mInertia
    
    torquePerSite2 = np.zeros((npoints,natoms))
    for i in range(natoms):
        timenew,torquePerSite2[:,i] = Interpol(time,torquePerSite[:,i],npoints)
    
    return force,torque,torquePerSite,alfa
    
def alphaOmegaTheta(time,alfa): #$\alpha$ aceleracion angular
    def integ(x, tck, constant=0):
        x = np.atleast_1d(x)
        out = np.zeros(x.shape[0], dtype=x.dtype)
        for n in xrange(len(out)):
            out[n] = interpolate.splint(0, x[n], tck)
    #    out += constant
        return out

    tck = interpolate.splrep(time,alfa,s=0)
    omega = integ(time,tck)   # r'$\omega$ velocidad angular'
    
    tck2 = interpolate.splrep(time,omega,s=0)
    theta = integ(time,tck2) # $\theta$ angulo de rotacion
    return omega,theta

def momentInertia(natoms,name,x,y):
    mass = [1.00794, 12.0107, 22.99, 6.941]
    m_inertia = 0.0
    for i in range(natoms):
        if name[i] == 'C':
            m_inertia += mass[1]*(x[i]**2+y[i]**2)
        elif name[i] == 'H':
            m_inertia += mass[0]*(x[i]**2+y[i]**2)
        elif name[i] == 'Na':
            m_inertia += mass[2]*(x[i]**2+y[i]**2)  
        elif name[i] == 'Li':
            m_inertia += mass[3]*(x[i]**2+y[i]**2) 
        else:
            print 'Element mass not known'
            raise(SystemExit)
    moment_inertia = m_inertia*293.227 # In hbar^2/eV
    return moment_inertia
                  
def currentXY(path,Time,z0,npts=100):
    times = time()
    timestep = 50 * int(Time/(times[1]-times[0]) / 50)
    timename = str(timestep).zfill(7)
    realpath = path+'td.'+timename
                  
    cube1x,coords = genGrid(realpath+'/current-sp1-x.cube')
    cube1y,coords = genGrid(realpath+'/current-sp1-y.cube')
    cube2x,coords = genGrid(realpath+'/current-sp2-x.cube')
    cube2y,coords = genGrid(realpath+'/current-sp2-y.cube')
    currx = cube1x + cube2x
    curry = cube1y + cube2y
    tol = 0.05
    x = []
    y = []
    Ix = []
    Iy = []              
    for i in range(currx.z.size):
        if abs(currx.z[i] - z0) < tol:
            x.append(currx.x[i])
            y.append(currx.y[i])
            Ix.append(currx.isovals[i])
            Iy.append(curry.isovals[i])
   
    xi = np.linspace(min(x),max(x),npts)
    yi = np.linspace(min(y),max(y),npts)
    Ixi = griddata(x, y, Ix, xi, yi, interp='linear')
    Iyi = griddata(x, y, Iy, xi, yi, interp='linear')
    return xi,yi,Ixi,Iyi

                  
def densityXY(path,Time,z0,npts=100):
    global dens0
    times = time()
    timestep = 50 * int(Time/(times[1]-times[0]) / 50)
    timename = str(timestep).zfill(7)
    realpath = path+'td.'+timename
    
    if dens0.size == 0:
        cube10,coords = genGrid(path+'td.0000000/density-sp1.cube')
        cube20,coords = genGrid(path+'td.0000000/density-sp2.cube')
        cube0 = cube10 + cube20
        dens0 = cube0.isovals
        print 'Cube data at t=0 loaded'
        
    cube1,coords = genGrid(realpath+'/density-sp1.cube')
    cube2,coords = genGrid(realpath+'/density-sp2.cube')
    cube = cube1 + cube2
    cube.isovals = cube.isovals - dens0
    tol = 0.05
    x = []
    y = []
    dens = []
    for i in range(cube.z.size):
        if abs(cube.z[i] - z0) < tol:
            x.append(cube.x[i])
            y.append(cube.y[i])
            dens.append(cube.isovals[i])
   
    xi = np.linspace(min(x),max(x),npts)
    yi = np.linspace(min(y),max(y),npts)
    densi = griddata(x, y, dens, xi, yi, interp='linear')
    return xi,yi,densi


def coords(realpath):
    cube0,coords = genGrid(realpath+'td.0000000/density-sp1.cube')
    return coords

def toalCurrent(path,timestep,out_every,npts=100):
    """ current as 3D-vector as a function of 3D-coordinates. Not finished. 
    """
    timename = str(timestep*out_every).zfill(7)
    realpath = path+'td.'+timename
    cube1x,coords = genGrid(realpath+'/current-sp1-x.cube')
    cube1y,coords = genGrid(realpath+'/current-sp1-y.cube')
    cube1z,coords = genGrid(realpath+'/current-sp1-z.cube')
    cube2x,coords = genGrid(realpath+'/current-sp2-x.cube')
    cube2y,coords = genGrid(realpath+'/current-sp2-y.cube')
    cube2z,coords = genGrid(realpath+'/current-sp2-z.cube')
    
    currx = cube1x + cube2x
    curry = cube1y + cube2y
    curry = cube1z + cube2z 
    
    return 0 #unfinished
