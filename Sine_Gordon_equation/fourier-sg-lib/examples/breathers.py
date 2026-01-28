import numpy as np
import fourier_sg_lib.drawing_tool as gg
import fourier_sg_lib.mathematics.two_soliton_start as two_s
import fourier_sg_lib.numerics.split_step_Fourier as split_f
import fourier_sg_lib.numerics.measure as msr
def starting_point():
    # Grid
    Nx, Lx = 1028, 100.0
    x = np.linspace(-Lx/2, Lx/2, Nx, endpoint=False)
    dx = x[1]-x[0]
    k = 2*np.pi*np.fft.fftfreq(Nx, d=dx)
    x0=20
    
    # Time
    dt, T = 0.01, 60.0
    steps = int(T/dt)
    save_every = max(1,20)
    speed = 0.5
    w_vec = [0.2,0.5,0.8]
    t0=10
    for w in w_vec:
    # Starting point
        u0,v0,name,title,f_a,kfilter = two_s.breather_start(x+x0,speed,w)
    # Integration
        method = split_f.strang_split_step
        Ut_approx,time = split_f.integration_loop(u0,v0,steps,save_every,dt,k,kfilter,method,nonlinear_step=split_f.sg_nonlinear_step)
        U_a = np.array([f_a(x+x0,t,speed,w) for t in time])
        #gg.plot_3d(x,[Ut_approx],time,f'breather_v{speed}_f{w}.png',azm=10,el=25)
        gg.plot_cm(x,[Ut_approx],time,f'brether_cm_v{speed}_f{w}.png')
        delta = split_f.sum_error(U_a,Ut_approx)
        print(f"{w} average error: {np.average(delta)}, sumed up: {np.sum(delta)} ")
starting_point()