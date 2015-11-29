# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 10:48:13 2015

@author: franco
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def integ(x, tck, constant=0):
    x = np.atleast_1d(x)
    out = np.zeros(x.shape[0], dtype=x.dtype)
    for n in xrange(len(out)):
        out[n] = interpolate.splint(0, x[n], tck)
#    out += constant
    return out

muxx = np.genfromtxt('mux.dat')
muyy = np.genfromtxt('muy.dat')
exx = np.genfromtxt('ex.dat')
eyy = np.genfromtxt('ey.dat')
torquee = np.genfromtxt('torque.dat')

time = muxx[:,0]
mux = muxx[:,1]
muy = muyy[:,1]
ex = exx[:,1]
ey = eyy[:,1]
timetor = torquee[:,0]
torque = torquee[:,1]
alfa = torquee[:,2]

total_time=time.max()
npoints=4000
#cut=0.5 #Where to slice the arrays
timenew = np.linspace(time.min(),time.max(),npoints)
f1 = interpolate.interp1d(timetor,torque,kind='linear')
f2 = interpolate.interp1d(time,mux,kind='linear')
f3 = interpolate.interp1d(time,muy,kind='linear')
f4 = interpolate.interp1d(time,ex[:mux.shape[0]],kind='linear')
f5 = interpolate.interp1d(time,ey[:mux.shape[0]],kind='linear')

torquenew = f1(timenew)
muxnew = f2(timenew)
muynew = f3(timenew)
exnew = f4(timenew)
eynew = f5(timenew)
#mutot = np.sqrt(muxnew**2+muynew**2)
#
#tck = interpolate.splrep(timetor,alfa,s=0)
#torqueint= integ(timetor,tck)
#
#tck2 = interpolate.splrep(timetor,torqueint,s=0)
#torqueint2= integ(timetor,tck2)

#plt.plot(timetor,torqueint2)
outputevery = 0.02
nplotsanim = 3010
plot_every=int(npoints*outputevery/total_time)

timenew=timenew*0.658
cut = nplotsanim*outputevery/total_time
#for i in range(int(npoints*cut/plot_every)):
for i in range(3000,nplotsanim):
    plt.figure(figsize=(8, 12), dpi=200)
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    ax1.plot(timenew[:i*plot_every],exnew[:i*plot_every], linewidth = 2, color='b')
    ax1.plot(timenew[:i*plot_every],eynew[:i*plot_every], linewidth = 2, color='r')
    #ax1.quiver(exarr, eyarr, eex, eey, color='b')
    #ax1.quiver(exarr, eyarr, exy[:,0], exy[:,1], color='b')
    #ax1.set_xlabel('time / ($\hbar/eV$)',fontsize=18)
    ax1.set_ylabel('electric field / $(V/\AA)$',fontsize=20)
    ax1.tick_params(axis = 'both', labelsize = 12)
    ax1.axis([timenew.min(),timenew.max()*cut,exnew.min(),exnew.max()])
#    ax1.legend(['$E_x$','$E_y$'],fontsize=30)
    #ax1.set_title('$-\cos(x)\hat{i}+\sin(x)\hat{j}$', fontsize=20)
    #ax1.ticklabel_format(axis='both', style = 'sci', scilimits=(0,0))
    #plt.savefig('eyvsex')
    #plt.subplot(121)
    #plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='r')
    ax2.plot(timenew[:i*plot_every],muxnew[:i*plot_every], linewidth = 2, color='b')
    ax2.plot(timenew[:i*plot_every],muynew[:i*plot_every], linewidth = 2, color='r')
    #ax2.quiver(muxarr, muyarr, muux, muuy, units='x', color='r')
    #ax2.set_xlabel('time / $fs$',fontsize=18)
    ax2.set_ylabel('dipole moment / $(e \AA)$',fontsize=20)
    ax2.tick_params(axis = 'both', labelsize = 12)
    ax2.legend(['$\mu_x$','$\mu_y$'],fontsize=20)
    #ax2.ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
    #ax2.set_autoscale_on(False)
    ax2.axis([timenew.min(),timenew.max()*cut,muxnew[:nplotsanim*plot_every].min(),muxnew[:nplotsanim*plot_every].max()])
    ax3.plot(timenew[:i*plot_every],torquenew[:i*plot_every], linewidth = 2, color='r')
    #ax2.quiver(muxarr, muyarr, muux, muuy, units='x', color='r')
    ax3.set_xlabel('time / $fs$',fontsize=20)
    ax3.set_ylabel('torque / $(eV/rad)$',fontsize=20)
    ax3.tick_params(axis = 'both', labelsize = 12)
    ax3.ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
    ax3.axis([timenew.min(),timenew.max()*cut,torquenew[:nplotsanim*plot_every].min(),torquenew[:nplotsanim*plot_every].max()])
    #ax2.set_autoscale_on(False)
#    ax3.axis([0,80,-8e-5,2e-5])
    #plt.suptitle('$(\cos(\omega t)\hat{i}-\sin(\omega t)\hat{j})\sin^2(\pi t/ \tau)$', fontsize=30)
    plt.savefig('dyplots'+str(i).zfill(4)+'.jpg')
    plt.cla()
    plt.clf()
    plt.close('all')


#mplt.plot_2_legends(time[:100000],ex[:100000],time[:100000],ey[:100000],'time / ($\hbar/eV)$','electric field / ($V/A$)','$E_x$','$E_y$','electric_field.pdf')
#mplt.plot_2_legends(time[:100000],mux[:100000],time[:100000],muy[:100000],'time / ($\hbar/eV)$','dipole moment / ($e \AA$)','$\mu_x$','$\mu_y$','dipole_moments.pdf')
#mplt.plot_default(timetor[:400],torque[:400],'time / ($\hbar/eV$)','torque / ($eV/rad$)','torque.pdf')
#plt.plot(mutot,torquenew)
#plt.plot(muynew,torquenew)

