#!/home/ar6116/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os.path

density=np.genfromtxt("density.txt",skip_header=2, dtype=None)
density1=np.genfromtxt("density1.txt",skip_header=2, dtype=None)
fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
#plt.plot(bin_vec,density,lw=1.5,ls='-',label='w=0.001')
plt.plot(density,lw=1.5,ls='--')
plt.plot(density1,lw=1.5,ls='-')
plt.xlabel(r'y', fontsize=12)
plt.ylabel(r'$\rho(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=10)
plt.legend()
fig.tight_layout()
fig.savefig('rho.pdf',dpi=fig.dpi)
plt.show()

velocity=np.genfromtxt("velocity.txt",skip_header=2, dtype=None)
velocity_total=np.genfromtxt("velocity1.txt",skip_header=2, dtype=None)
fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
plt.plot(velocity,lw=1.5,ls='-')
plt.plot(velocity_total,lw=1.5,ls='-')
plt.xlabel(r'y', fontsize=12)
plt.ylabel(r'$v_x(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=10)
plt.legend()
fig.tight_layout()
fig.savefig('v.pdf',dpi=fig.dpi)
plt.show()

#if os.path.isfile("velocity_y.txt"):
#    velocity_y=np.genfromtxt("velocity_y.txt",skip_header=2, dtype=None)
#    fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
#    plt.plot(velocity_y,lw=1.5,ls='-')
#    plt.xlabel(r'y', fontsize=12)
#    plt.ylabel(r'$v_y(y)$', fontsize=12)
#    plt.tick_params(axis='both', labelsize=10)
#    plt.legend()
#    fig.tight_layout()
#    fig.savefig('vy.pdf',dpi=fig.dpi)
#    plt.show()

#if os.path.isfile("velocity_div.txt"):
#    velocity_div=np.genfromtxt("velocity_div.txt",skip_header=2, dtype=None)
#    fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
#    plt.plot(velocity_div,lw=1.5,ls='-')
#    plt.xlabel(r'y', fontsize=12)
#    plt.ylabel(r'$div_y(v_y)(y)$', fontsize=12)
#    plt.tick_params(axis='both', labelsize=10)
#    plt.legend()
#    fig.tight_layout()
#    fig.savefig('div_vy.pdf',dpi=fig.dpi)
#    plt.show()


temperature=np.genfromtxt("temperature.txt",skip_header=2, dtype=None)
fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
plt.plot(temperature,lw=1.5,ls='-')
plt.xlabel(r'y', fontsize=12)
plt.ylabel(r'$T(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=10)
plt.legend()
fig.tight_layout()
fig.savefig('temp.pdf',dpi=fig.dpi)
plt.show()


stressK22=np.genfromtxt("stressK22_1.txt",skip_header=2, dtype=float)
stressV22=np.genfromtxt("stressV22_1.txt",skip_header=2, dtype=float)
stressT22=stressV22+stressK22

stressK22_total=np.genfromtxt("stressK22.txt",skip_header=2, dtype=float)
stressV22_total= np.genfromtxt("stressV22.txt",skip_header=2, dtype=float)
stressT22_total=stressV22_total+stressK22_total

fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
plt.plot(stressV22,lw=1.5,ls='-',label=r'$\sigma^v_{yy}(y)$')
plt.plot(stressK22,lw=1.5,ls='-',label=r'$\sigma^k_{yy}(y)$')
plt.plot(stressT22,lw=1.5,ls='-',label=r'$\sigma^t_{yy}(y)$')

plt.plot(stressV22_total,lw=1.5,ls='--',label=r'$\sigma^v_{yy}(y)$')
plt.plot(stressK22_total,lw=1.5,ls='--',label=r'$\sigma^k_{yy}(y)$')
plt.plot(stressT22_total,lw=1.5,ls='--',label=r'$\sigma^t_{yy}(y)$')

plt.ylim (-2.5,1)
plt.xlabel(r'y', fontsize=12)
plt.ylabel(r'$\sigma(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=10)
plt.legend(ncol=3)
fig.tight_layout()
fig.savefig('stress22.pdf',dpi=fig.dpi)
plt.show()


stressK33=np.genfromtxt("stressK33_1.txt",skip_header=2, dtype=float)
stressV33=np.genfromtxt("stressV33_1.txt",skip_header=2, dtype=float)
stressT33=stressV33+stressK33

stressK33_total=np.genfromtxt("stressK33.txt",skip_header=2, dtype=float)
stressV33_total= np.genfromtxt("stressV33.txt",skip_header=2, dtype=float)
stressT33_total=stressV33_total+stressK33_total

fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
plt.plot(stressV33,lw=1.5,ls='-',label=r'$\sigma^v_{zz}(y)$')
plt.plot(stressK33,lw=1.5,ls='-',label=r'$\sigma^k_{zz}(y)$')
plt.plot(stressT33,lw=1.5,ls='-',label=r'$\sigma^t_{zz}(y)$')

plt.plot(stressV33_total,lw=1.5,ls='--',label=r'$\sigma^v_{zz}(y)$')
plt.plot(stressK33_total,lw=1.5,ls='--',label=r'$\sigma^k_{zz}(y)$')
plt.plot(stressT33_total,lw=1.5,ls='--',label=r'$\sigma^t_{zz}(y)$')

plt.ylim (-2.5,2)
plt.xlabel(r'y', fontsize=12)
plt.ylabel(r'$\sigma(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=10)
plt.legend()
fig.tight_layout()
fig.savefig('stress33.pdf',dpi=fig.dpi)
plt.show()

#if os.path.isfile("stressK11.txt") and os.path.isfile("stressV11.txt") and os.path.isfile("stressK22.txt") and os.path.isfile("stressV22.txt") and os.path.isfile("stressK33.txt") and os.path.isfile("stressV33.txt") :
#    stressK11=np.genfromtxt("stressK11.txt",skip_header=2, dtype=float)
#    stressV11=np.genfromtxt("stressV11.txt",skip_header=2, dtype=float)
#    stressT11=stressV11+stressK11
#    stressK22=np.genfromtxt("stressK22.txt",skip_header=2, dtype=float)
#    stressV22=np.genfromtxt("stressV22.txt",skip_header=2, dtype=float)
#    stressT22=stressV22+stressK22
#    stressK33=np.genfromtxt("stressK33.txt",skip_header=2, dtype=float)
#    stressV33=np.genfromtxt("stressV33.txt",skip_header=2, dtype=float)
#    stressT33=stressV33+stressK33
#
#    fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
#    plt.plot( density,lw=1.5,ls='-',label=r'$\rho(y)$')
#    plt.plot(stressT11,lw=1.5,ls='-',label=r'$\sigma^t_{xx}(y)$')
#    plt.plot(stressT22,lw=1.5,ls='-',label=r'$\sigma^t_{yy}(y)$')
#    plt.plot(stressT33,lw=1.5,ls='-',label=r'$\sigma^t_{zz}(y)$')
#    plt.xlabel(r'y', fontsize=12)
#    plt.ylabel(r'$\sigma(y)$', fontsize=12)
#    plt.tick_params(axis='both', labelsize=10)
#    plt.legend()
##    plt.xlim(-5,5)
#    plt.ylim(-0.5,0.5)
#    fig.tight_layout()
#    fig.savefig('stressN.pdf',dpi=fig.dpi)
#    plt.show()


stressK12=np.genfromtxt("stressK12_1.txt",skip_header=2, dtype=float)
stressV12=np.genfromtxt("stressV12_1.txt",skip_header=2, dtype=float)
stressT12=stressV12+stressK12

stressK12_total=np.genfromtxt("stressK12.txt",skip_header=2, dtype=float)
stressV12_total= np.genfromtxt("stressV12.txt",skip_header=2, dtype=float)
stressT12_total=stressV12_total+stressK12_total

fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
plt.plot(stressV12,lw=1.5,ls='-',label=r'$\sigma^v_{yy}(y)$')
plt.plot(stressK12,lw=1.5,ls='-',label=r'$\sigma^k_{yy}(y)$')
plt.plot(stressT12,lw=1.5,ls='-',label=r'$\sigma^t_{yy}(y)$')

plt.plot(stressV12_total,lw=1.5,ls='--',label=r'$\sigma^v_{yy}(y)$')
plt.plot(stressK12_total,lw=1.5,ls='--',label=r'$\sigma^k_{yy}(y)$')
plt.plot(stressT12_total,lw=1.5,ls='--',label=r'$\sigma^t_{yy}(y)$')

plt.ylim (-0.25,0.05)
plt.xlabel(r'y', fontsize=12)
plt.ylabel(r'$\sigma(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=10)
plt.legend()
fig.tight_layout()
fig.savefig('stress12.pdf',dpi=fig.dpi)
plt.show()



velocity_grad=np.genfromtxt("velocity_grad.txt",skip_header=2, dtype=None)
viscosity=stressT12/velocity_grad
viscosity_total=stressT12_total/velocity_grad

fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
plt.plot(viscosity,lw=1.5,ls='-')
plt.plot(viscosity_total,lw=1.5,ls='--')
plt.xlabel(r'y', fontsize=12)
plt.ylabel(r'$eta_x(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=10)
plt.ylim(0,10)
plt.legend()
fig.tight_layout()
fig.savefig('viscosity.pdf',dpi=fig.dpi)
plt.show()

#    density=np.genfromtxt("density.txt",skip_header=2, dtype=None)
#    fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
#    plt.plot(density,viscosity,lw=1.5,ls='-')
#    plt.xlabel(r'y', fontsize=12)
#    plt.ylabel(r'$eta_x(y)$', fontsize=12)
#    plt.tick_params(axis='both', labelsize=10)
#    plt.xlim(0,2)
#    plt.ylim(0,15)
#    plt.legend()
#    fig.tight_layout()
#    fig.savefig('viscosity-density.pdf',dpi=fig.dpi)
#    plt.show()


