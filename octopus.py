
import glob
import numpy as np
from scipy import interpolate

#Dipole Moment
def DipoleFieldVec():
    path = './td.general/'
    ex = ey = ez = mux = muy = muz = time = []
    
    with open(path+'multipoles','r') as multipoles:
        for line in multipoles:
            sl = line.strip()
            if not sl.startswith("#"):
                l = line.split()
                mux.append(l[2])
                muy.append(l[3])
                muz.append(l[4])
                    
    with open(path+'laser','r') as laser:
        for line in laser:
            sl = line.strip()
            if not sl.startswith("#"):
                l = line.split()
                time.append(l[1])
                ex.append(l[2])
                ey.append(l[3])
                ez.append(l[4])
    
    time = np.array(time[:len(mux)])
    efield = np.array(ex[:len(mux)],ey[:len(mux)],ez[:len(mux)])
    mu = np.array(mux,muy,muz)
    
    return time,efield,mu


def Epot(efield,mu):
    epot = mu[1]*efield[1]+mu[2]*efield[2]+mu[3]*efield[3]
    return epot

def TorqueNuc(efield,mu):
    torque = efield[1]*mu[2]-mu[1]*efield[2]
    return torque
    
def Interpol(datax,datay,npoints):
    f1 = interpolate.interp1d(datax,datay,kind='linear')
    xnew = np.linspace(datax.min(),datax.max(),npoints)
    ynew = f1(xnew)
    return xnew,ynew

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

def Forces():
    dt = 0.001
    forcefiles = glob.glob('td.0*/forces.xsf')
    natoms = sum(1 for line in open('td.0000000/forces.xsf')) - 1
    el = idn = x = y = z = fx0 = fy0 = fz0 = fx = fy = fz = torque = []*natoms
    
    forcePerStep = torquePerStep = time = np.zeros((len(forcefiles)))
    torquePerSite = np.zeros((len(forcefiles),natoms))

    with open(forcefiles[0],'r+') as f:
        f.readline()
        for j,line in enumerate(f):
            l = f.readline().split()
            el[j] = l[0]
            if el[j] == 'C':
                idn[j] = 6
            elif el[j] == 'H':
                idn[j] == 1
            elif el[j] == 'Na':
                idn[j] == 9
            x[j] = float(l[1])
            y[j] = float(l[2])
            z[j] = float(l[3])
            fx0[j] = float(l[4]) 
            fy0[j] = float(l[5])
            fz0[j] = float(l[6])        
    
    # Create AXSF file to be read with XCRYSDEN
    g = open('netforces.axsf','w')
    g.write('ANIMSTEPS %i \n' & len(forcefiles))

    for i,ffile in enumerate(forcefiles):
        time[i] = i*dt
        g.write('ATOMS %i \n' & i)
        with open(ffile,'r+') as f:
            f.readline()
            rForce = [0.,0.,0.]
            totalTorque = 0.
            for j,line in enumerate(f):
                l = f.readline().split()
                fx[j] = float(l[4]) - fx0[j]
                fy[j] = float(l[5]) - fy0[j]
                fz[j] = float(l[6]) - fz0[j]   
                # Write to AXSF file
                g.write("%i %10.6 %10.6 %10.6    %10.6 %10.6 %10.6" % (idn[j],x[j],y[j],z[j],fx[j],fy[j],fz[j]))
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
    
    return time,forcePerStep,torquePerStep,torquePerSite


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



#Torque