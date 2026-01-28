import numpy as np
import json
from fourier_sg_lib.numerics import split_step_Fourier as split_f
def init_sol(file_path="params.json"):
    with open(file_path) as params_file:
        data = json.load(params_file)
    #adding x, dx & k to space grid
    Lx = data['spaceGrid']['lx']
    Nx = data['spaceGrid']['nx']
    data['spaceGrid']['x'] =  np.linspace(-Lx/2, Lx/2, Nx, endpoint=False)
    data['spaceGrid']['dx']  = data['spaceGrid']['x'][1]-data['spaceGrid']['x'][0]
    data['spaceGrid']['k'] = 2*np.pi*np.fft.fftfreq(Nx, d=data['spaceGrid']['dx'])
    #adding steps to time grid
    data['timeGrid']['steps'] = int(data['timeGrid']['t']/data['timeGrid']['dt'])
    #method of integration
    m = data['method']
    data['method'] = split_f.strang_split_step if m=='strangSplitStep' else split_f.step_4th_order
    #adding gamma to soliton params
    print("Initial data for calculations:\n",data)
    return data
def marker(x, time, v,x0):
    m = []

    for t in time:
        # Position of the travelling point at time t
        pos =x0+ v * t
        
        # Index of x closest to pos
        idx = np.argwhere(x >= pos)[0][0]
        m.append(x[idx])
    return m
def kink(x, t, v, x0):
    g = 1.0/np.sqrt(1.0 - v*v)
    return 4.0*np.arctan(np.exp(g*(x - x0 - v*t)))

def antikink(x, t, v, x0):
    g = 1.0/np.sqrt(1.0 - v*v)
    return 4.0*np.arctan(np.exp(-g*(x - x0 - v*t)))

def kink_t(x, t, v, x0):
    g = 1.0/np.sqrt(1.0 - v*v)
    cosh_arg = g*(x-x0-v*t)
    return -2*g*v/np.cosh(cosh_arg)

def antikink_t(x, t, v, x0):
    return -kink_t(x, t, x0, v)

def kink_kink(x, t, v):
    gamma = 1/np.sqrt(1-v**2)
    frac = v*np.sinh(x*gamma)/np.cosh(v*t*gamma)
    return 4*np.arctan(frac)

def kink_antikink(x, t, v,t0=0):
    gamma = 1/np.sqrt(1-v**2)
    frac = np.sinh(gamma*v*(t-t0))/v/np.cosh(gamma*x)
    return -4*np.arctan(frac)

def kink_antikink_t0(x,v,t0=0):
    gamma = 1 / np.sqrt(1 - v**2)

    num = 4 * gamma * v**2 * np.cosh(gamma * v * t0) * np.cosh(gamma * x)
    den = v**2 * np.cosh(gamma * x)**2 + np.sinh(gamma * v * t0)**2

    return -num / den
def breather(x, t,v,w):
    gamma = 1/np.sqrt(1-v**2)
    sq = np.sqrt(1-w**2)
    frac = sq*np.cos(w*gamma*(t-v*x))/w/np.cosh(sq*gamma*(x-v*t))
    return 4*np.arctan(frac)

def breather_t(x,t,v,w):
    gamma = 1/np.sqrt(1-v**2)
    sq = np.sqrt(1.0 - w**2)
    alpha = w * gamma * v * x
    beta  = sq * gamma * x
    F0 = (sq * np.cos(-w * gamma * v * x)) / (w * np.cosh(beta))

    numer = 4.0 * sq * gamma * (sq * v * np.cos(alpha) * np.tanh(beta) + w * np.sin(alpha))
    denom = w * np.cosh(beta)

    return  numer / (denom * (1.0 + F0**2))

def two_breathers(x,t,v,w,x0,v0=False):
    gamma = 1/np.sqrt(1-v**2)
    al1 = np.sqrt((1-v)/(1+v))
    al2 = np.sqrt((1+v)/(1-v))
    cosq = np.sqrt(1-w**2/gamma**2)
    sinq = w/gamma
    x1,x2 = x+x0,x-x0
    a1, a2 = (x1-v*t)*gamma*cosq, (x2+v*t)*gamma*cosq
    b1, b2 = w*(t-v*x1), w*(t+v*x2)
    
    B, C ,D = 4*al1**2*al2**2*(1-2*w**2/gamma**2), -al1**2+al2**2, al1**2+al2**2
    
    u1,w1 = al1*cosq,al1*w/gamma
    u2,w2 = al2*cosq, al2*w/gamma

    fi = -u1*w2*C**2*(np.cosh(a1)*np.cos(b2) + np.cosh(a2)*np.cos(b1)) + 4*(u1*w2)**2*C*(
        np.sinh(a1)*np.sin(b2) - np.sinh(a2)*np.sin(b1))
    
    fr = -w1*w2*(D**2+B)*np.cosh(a1)*np.cosh(a2)+u1*u2*(D**2-B)*np.cos(b2)*np.cos(b1) + 4*(u1*w2)**2*D*(
        np.sinh(a1)*np.sinh(a2)+np.sin(b1)*np.sin(b2))
    if v0:
        a1_t = -v * gamma * cosq
        a2_t =  v * gamma * cosq
        b1_t = w 
        b2_t = w 
        fi_t = (
        -u1*w2*C**2 * (
            np.sinh(a1)*a1_t*np.cos(b2)
            + np.cosh(a1)*(-np.sin(b2))*b2_t
            + np.sinh(a2)*a2_t*np.cos(b1)
            + np.cosh(a2)*(-np.sin(b1))*b1_t
        )
        + 4*(u1*w2)**2*C * (
            np.cosh(a1)*a1_t*np.sin(b2)
            + np.sinh(a1)*np.cos(b2)*b2_t
            - np.cosh(a2)*a2_t*np.sin(b1)
            - np.sinh(a2)*np.cos(b1)*b1_t
        )
        )

        fr_t = (
            -w1*w2*(D**2+B) * (
            np.sinh(a1)*a1_t*np.cosh(a2)
            + np.cosh(a1)*np.sinh(a2)*a2_t
            )
            + u1*u2*(D**2-B) * (
                -np.sin(b2)*b2_t*np.cos(b1)
                + np.cos(b2)*(-np.sin(b1))*b1_t
            )
            + 4*(u1*w2)**2*D * (
                np.cosh(a1)*a1_t*np.sinh(a2)
            + np.sinh(a1)*a2_t*np.cosh(a2)
            + np.cos(b1)*b1_t*np.sin(b2)
            + np.sin(b1)*np.cos(b2)*b2_t
            )
        )
        return  4*np.arctan(fi/fr),4*(fr*fi_t - fi*fr_t) / (fr**2 + fi**2)
    return 4*np.arctan(fi/fr)

def pos_shift(v,w):
    gamma = 1/np.sqrt(1-v**2)
    al1 = np.sqrt((1-v)/(1+v))
    al2 = np.sqrt((1+v)/(1-v))
    cosq = np.sqrt(1-w**2/gamma**2)
    u1 = al1*cosq
    u2 = al2*cosq
    B,D = 4*al1**2*al2**2*(1-2*w**2/gamma**2), al1**2+al2**2
    return np.arctanh(4*u1*u2*D/(D**2+B))


def breather_plus_breather(x,t,v,w,x0=15):
    return breather(x-x0,t,-v,w)+breather(x+x0,t,v,w)

def kink_plus_antikink(x,t,v,x0):
   return kink(x, t, -x0,  v) + antikink(x, t,  x0, -v) - 2*np.pi

def kink_plus_kink(x,t,v,x0):
    return kink(x, t, -x0,  v) + kink(x, t,  x0, -v)-2*np.pi