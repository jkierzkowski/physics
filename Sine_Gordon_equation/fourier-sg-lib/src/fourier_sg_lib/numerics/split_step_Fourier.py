import numpy as np
#import mkl 
from scipy.interpolate import approximate_taylor_polynomial

def linear_half_step(u,v,c,ws,s_over_w):
    Uh, Vh = np.fft.fft(u), np.fft.fft(v)
    Uh, Vh = c*Uh + s_over_w*Vh, -ws*Uh + c*Vh
    return (np.fft.ifft(Uh).real, np.fft.ifft(Vh).real)

def sg_linear_propagator(k,dt):
    w = np.abs(k)
    c = np.cos(w*dt/2.0)
    s_over_w = np.zeros_like(w)
    m        = w>0
    s_over_w[m]  = np.sin(w[m]*dt/2.0)/w[m]
    s_over_w[~m] = dt/2.0
    ws = w*np.sin(w*dt/2.0)
    return c, s_over_w, ws

#removes higher modes removing noise
def filtered(arr):
    N = arr.size
    cutoff = N // 3
    Ah = np.fft.fft(arr)
    Ah[cutoff:-cutoff] = 0
    return np.fft.ifft(Ah).real

def sg_nonlinear_step(u,v,t,dt,kfilter):
    return filtered(v-dt*np.sin(u)) if kfilter else v-dt*np.sin(u)

def approx_nonlinear_step(u,v,t,dt,kfilter,order=15):
    sin_taylor = approximate_taylor_polynomial(np.sin, 0, order, 1,
                                               order=order)
    return filtered(v-dt*sin_taylor(u)) if kfilter else v-dt*sin_taylor(u)

def perturb_sg(u,v,t,dt,kfilter,eps1=0,eps2=0):
    r =v-dt*np.sin(u)-eps1*u-eps2*u**3
    return filtered(r) if kfilter else r

def strang_split_step(u, v, dt, k,kfilter,propagator,nonlinear_step=sg_nonlinear_step,
                      lin_split_step=linear_half_step):
    c, s_over_w, ws =propagator
    t=0
    u,v = lin_split_step(u,v,c,ws,s_over_w)
    v  = nonlinear_step(u,v,t,dt,kfilter)
    u,v = lin_split_step(u,v,c,ws,s_over_w)
    t+=dt
    return u, v



def step_4th_order(u,v,dt,k,kfilter,propagator,nonlinear_step=sg_nonlinear_step,
                      lin_split_step=linear_half_step):
    g1 = 1.0/(2-2**(1/3))
    g2 = -2**(1/3)/(2-2**(1/3))
    prop1 = sg_linear_propagator(k,g1*dt)
    prop2 = sg_linear_propagator(k,g2*dt)
    u,v = strang_split_step(u,v,g1*dt,k,kfilter,prop1,nonlinear_step,
                      lin_split_step)
    u,v = strang_split_step(u,v,g2*dt,k,kfilter,prop2,nonlinear_step,
                      lin_split_step)
    u,v = strang_split_step(u,v,g1*dt,k,kfilter,prop1,nonlinear_step,
                      lin_split_step)
    return u,v

def integration_loop(u0,v0,steps,save_every,dt,k,kfilter,strang_split=strang_split_step,nonlinear_step=sg_nonlinear_step,
                      lin_split_step=linear_half_step):
    timeline = []
    time = []
    u, v = u0.copy(),v0.copy()
    t = 0.0
    propagator = sg_linear_propagator(k,dt)
    for n in range(steps+1):
        if n % save_every == 0:
            timeline.append(u.copy())
            time.append(t)
        u, v = strang_split(u, v, dt, k,kfilter,propagator,nonlinear_step,
                      lin_split_step)
        t += dt
    return (np.array(timeline),np.array(time))

def spike_filter(timeline,spike = 0.9,cutoff=5):
    diff = np.array([np.abs(np.diff(t)) for t in timeline])
    indx = [np.argwhere(row>spike)[len(np.argwhere(row>spike))//2-1:len(np.argwhere(row>spike))//2+1] for row in diff ]
    indx[0] = np.array([[0],[len(timeline[0])]])
    indx = np.array(indx)
    tl_arr =[]
    for i, frame in enumerate(timeline):
        pad = np.full_like(timeline[0],np.nan)
        pad[indx[i][0][0]+cutoff+1:indx[i][1][0]-cutoff] = frame[indx[i][0][0]+cutoff+1:indx[i][1][0]-cutoff]
        tl_arr.append(pad)
    tl_arr= np.array(tl_arr)
    return tl_arr

def sum_error(U,U_a):
    return np.sum(np.abs(U_a-U),axis=1)/len(U[0])

