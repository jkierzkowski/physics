import numpy as np
import matplotlib.pyplot as plt

def max_min(U,plot=False):
    Uabs =np.abs(U)
    m = np.max(Uabs,axis=1)
    x_idx = np.argmax(Uabs,axis=1)
    mask = np.diff(m)>0
    mask = np.where(np.diff(mask.astype(int)) ==-1)[0]+1
    if plot:
        plt.scatter(np.linspace(0,len(m),len(m)),m,s=3)
        plt.scatter(np.linspace(0,len(m),len(m))[mask],m[mask],s=3,marker='^')
        plt.show()
    return x_idx[mask],mask
    
def freq(t,tx):
    return 1/(np.diff(t[tx]))*np.pi

def vel(x,t,x_idx,tx):
    return np.diff(x[x_idx])/np.diff(t[tx])

def measure_all(U,x,t,plot=False):
    x_idx, tx = max_min(U,plot)
    v = vel(x,t,x_idx,tx)
    gamma = 1/np.sqrt(1-v**2)
    f = freq(t,tx)*gamma
    fmean = np.average(f)
    fmask = f<fmean*1.5
    return (f[fmask].mean(),v[fmask].mean())

def measure_shift(U,x,t,before_col,after_col,x0,plot=False):
    mid = len(x)//2
    U1 = U[:,:mid]
    U2 = U[:,mid:]
    m1before, t1before = max_min(U1[:before_col,:])
    m2before, t2before = max_min(U2[:before_col,:])
    v01 =np.average(vel(x[:mid],t[:before_col],m1before, t1before))
    v02 = np.average(vel(x[mid:],t[:before_col],m2before, t2before))
    m2after, t2after = max_min(U1[after_col:,:])
    m1after, t1after = max_min(U2[after_col:,:])
    t1after += after_col
    t2after += after_col
    d1 = np.abs(-np.full_like(t1after,x0)+v01*t[t2after]-x[m1after+mid])
    d2 = np.abs(np.full_like(t2after,x0)+v02*t[t1after]-x[m2after])
    return d1,d2

def u_t(U,dt):
    Ut = np.empty_like(U)
    Ut[1:-1,:] = (U[2:,:] - U[:-2,:]) / (2*dt)
    Ut[0,:]    = (U[1,:] - U[0,:]) / dt 
    Ut[-1,:]   = (U[-1,:] - U[-2,:]) / dt  
    return Ut

def u_x(U,dx):
    Ux = np.empty_like(U)
    Ux[:,1:-1] = (U[:,2:] - U[:,:-2]) / (2*dx)
    Ux[:,0]    = (U[:,1] - U[:,0]) / dx
    Ux[:,-1]   = (U[:,-1] - U[:,-2]) / dx  
    #return Ux 
    return (np.roll(U, -1, axis=1) - np.roll(U, 1, axis=1)) / (2*dx)

def energy_density(U,V,dx,dt):
    return  0.5*(u_x(U,dx)**2+u_t(U,dt)**2)+V