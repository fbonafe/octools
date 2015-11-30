
import glob
import numpy as np
from scipy import interpolate

hbar = 0.658

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
def loadData(out_every):
    
    time,field,mu,charge = DipoleFieldVec()
    dt = time[1]-time[0]
    npoints = time.shape[0]
    
    force,torque,torquePerSite,alfa = Forces(time,dt*out_every)
    
    muabs = np.sqrt(mu[0]**2+mu[1]**2+mu[2]**2)
    
    muphi = np.angle(mu[0]+mu[1]*1j)
    #muphi = np.array([x+2.*np.pi for x in muphi if x < 0])
    for i in range(1,muphi.shape[0]):
        if muphi[i] < muphi[i-1]:
            muphi[i:] = 2*np.pi + muphi[i:]
        
    mudot = MuDot(time,muabs)
    muddot = MuDDot(time,muabs)
    omega,theta = AlphaOmegaTheta(time,alfa)
    energy = Energy()
    
    return dynamics(time,field,mu,force,torque,alfa,omega,theta,charge,energy,mudot,muddot,torquePerSite,muabs,muphi)
   

#################################################################
# This quantities are NOT automatically calculated by loadData,
# should be used separately.

def Epot(efield,mu):
    epot = mu[0]*efield[0]+mu[1]*efield[1]+mu[2]*efield[2]
    return epot

def TorqueNuc(efield,mu):
    torque = efield[0]*mu[1]-mu[0]*efield[1]
    return torque
#################################################################

def Energy():
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
def DipoleFieldVec():
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

def MuDot(time,mu):
    tck,uout = interpolate.splprep([time,mu],s=0.,k=2,per=False)
    dx,dy = interpolate.splev(uout,tck,der=1)
    mudot = dy/dx
    return mudot
    
def MuDDot(time,mu):
    tck,uout = interpolate.splprep([time,mu],s=0.,k=2,per=False)
    dx,dy = interpolate.splev(uout,tck,der=1)
    mudot = dy/dx
    tck,uout = interpolate.splprep([time,mudot],s=0.,k=2,per=False)
    ddx,ddy = interpolate.splev(uout,tck,der=1)
    muddot = ddy/ddx
    return muddot
    
#Force

def Forces(realtime,dt): # dt in hbar/eV
    forcefiles = sorted(glob.glob('td.0*/forces.xsf')) # VERY IMPORTANT TO INCLUDE sorted
    natoms = sum(1 for line in open('td.0000000/forces.xsf')) - 1
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
            
            totalForce = np.sqrt(rForce[1]**2+rForce[2]**2+rForce[3]**2)
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
    
def AlphaOmegaTheta(time,alfa): #$\alpha$ aceleracion angular
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
    mass = [1.00794, 12.0107, 22.99]
    m_inertia = 0.0
    for i in range(natoms):
        if name[i] == 'C':
            m_inertia += mass[1]*(x[i]**2+y[i]**2)
        elif name[i] == 'H':
            m_inertia += mass[0]*(x[i]**2+y[i]**2)
        elif name[i] == 'Na':
            m_inertia += mass[2]*(x[i]**2+y[i]**2)  
        else:
            print 'Element mass not known'
            raise(SystemExit)
    moment_inertia = m_inertia*293.227 # In hbar^2/eV
    return moment_inertia
