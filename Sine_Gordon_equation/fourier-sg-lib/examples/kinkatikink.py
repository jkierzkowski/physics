import numpy as np
import fourier_sg_lib.drawing_tool as gg
import fourier_sg_lib.mathematics.two_soliton_start as two_s
import fourier_sg_lib.numerics.split_step_Fourier as split_f

def starting_point():
    # Grid
    Nx, Lx = 1024, 100.0
    x = np.linspace(-Lx/2, Lx/2, Nx, endpoint=False)
    dx = x[1]-x[0]
    k = 2*np.pi*np.fft.fftfreq(Nx, d=dx)
    x0=0
    
    # Time
    dt, T = 0.01, 60.0
    steps = int(T/dt)
    save_every = max(1,20)
    speed = 0.5
    w = 0.7
    t0=T/2
    # Starting point
    u0,v0,name,title,f_a,kfilter = two_s.kink_antikink_start(x,speed,t0)
    # Integration
    method = split_f.strang_split_step
    Ut_approx,time = split_f.integration_loop(u0,v0,steps,save_every,dt,k,kfilter,method,nonlinear_step=split_f.sg_nonlinear_step)
    #np.savetxt('two_breather.txt', Ut_approx, fmt="%.10e")
    #Ut_approx = np.loadtxt('two_breather.txt')
    #x_idx,tx = msr.max_min(Ut_approx)
    #V = 1-np.cos(Ut_approx)
    #plt.plot(time,msr.energy_density(Ut_approx,V,dx,dt))
    #plt.show()
    #time = np.linspace(0,T,steps//save_every+1)
    #print(Ut_approx.shape)
    #U_a = np.array([f_a(x,t,speed,x0) for t in time])
    gg.plot_3d(x,[Ut_approx],time,['numerical outcome'],f'Kink anti-kink collision at velocity {speed}',f"kink_antikink_v{speed}.png")
    U_a = np.array([f_a(x,t,speed,t0) for t in time])
    delta = split_f.sum_error(U_a,Ut_approx)
    print(np.average(delta), np.sum(delta))
    #gg.generate_gif(x,Ut_approx,time,name,U_a=U_a,title=title,y_min=-10,y_max=10,
                    #msr=False)
    #msr.measure_all(Ut_approx,x,time)
    #Ut,time = split_f.integration_loop(u0,v0,steps,save_every,dt,k,kfilter,method)
    #print(Ut.shape)
    #Clearing error from FFT periodic BC mainly for kink_kink interactions
    #timeline = split_f.spike_filter(timeline)
    # Analitic partx,t,w,v,gamma
    
    
   

    #np.savetxt('Ua_breather.txt', U_a, fmt="%.10e")
    #msr.energy_density(U_a,1-np.cos(U_a),dx,dt,time)
    #U_a = np.loadtxt('Ua_breather.txt')
    
    #e_a =msr.energy_density(U_a,1-np.cos(U_a),dx,dt)
    #e_t = msr.energy_density(Ut_approx,V,dx,dt)
    #plt.plot(time,np.trapz(e_a,x=x,axis=1)[:-1],label='Energy analytic')
    #plt.plot(time,np.trapz(e_t,x=x,dx=x,axis=1)[:-1],linestyle='--',label='Energy numerical')
    #plt.legend()
    #plt.show()
    # Making a gif
   
    
    #print(timeline[0:4])
starting_point()