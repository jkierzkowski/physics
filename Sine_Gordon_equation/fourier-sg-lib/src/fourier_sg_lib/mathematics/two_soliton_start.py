import numpy as np
from . import sg_soliton_functions as ssf
def kink_kink_start(x,v,x0):
    u0 = ssf.kink_kink(x,t=0,v=v)
    v0 = np.zeros_like(u0)  
    return(u0,v0,'kink_kink_filtered.gif','kink_kink',ssf.kink_kink,True)

def kink_start(x,v,x0):
    u0 = ssf.kink(x,t=0,v=v,x0=x0)
    v0 = ssf.kink_t(x, t=0, x0=x0, v=v)  
    return(u0,v0,'kink_filtered.gif','kink_kink',ssf.kink,True)

def kink_antikink_start(x,v,t0):
    u0 = ssf.kink_antikink(x,0,v,t0)
    v0 = ssf.kink_antikink_t0(x,v,t0)
    
    return(u0,v0,'kink_antikink.gif','kink_antikink',ssf.kink_antikink,True)

def stat_breather_start(x,w):
    u0 = ssf.breather(x,0,v=0,w=w)
    v0 = np.zeros_like(u0)
    return(u0,v0,'stationary_breather.gif','breather',ssf.breather,False)

def breather_start(x,v,w):
    u0 = ssf.breather(x,0,v,w)
    v0 = ssf.breather_t(x,0,v,w)
    return(u0,v0,'moving_breather.gif','moving breather',ssf.breather,True)

def breather_breather_start(x,v,w,x0):
    u0,v0 = ssf.two_breathers(x,0,v,w,x0,v0=True)
    return(u0,v0,'breather_collision.gif','breather_breather',ssf.two_breathers,False)