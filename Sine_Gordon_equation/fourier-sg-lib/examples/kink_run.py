import numpy as np
import fourier_sg_lib.drawing_tool as gg
import fourier_sg_lib.mathematics.two_soliton_start as two_s
import fourier_sg_lib.numerics.split_step_Fourier as split_f
import fourier_sg_lib.numerics.measure as msr
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
    save_every = 1#max(1,20)
    speed = 0.0
    w = 0.7
    t0=10
    # Starting point
    u0,v0,name,title,f_a,kfilter = two_s.kink_start(x,speed,x0)
    # Integration
    method = split_f.strang_split_step
    Ut_approx,time = split_f.integration_loop(u0,v0,steps,save_every,dt,k,kfilter,method,nonlinear_step=split_f.sg_nonlinear_step)
    #np.savetxt('two_breather.txt', Ut_approx, fmt="%.10e")
    #Ut_approx = np.loadtxt('two_breather.txt')
    x_idx,tx = msr.max_min(Ut_approx)
    #V = 1-np.cos(Ut_approx)
    #plt.plot(time,msr.energy_density(Ut_approx,V,dx,dt))
    #plt.show()
    #time = np.linspace(0,T,steps//save_every+1)
    #print(Ut_approx.shape)
    U_a = np.array([f_a(x,t,speed,x0) for t in time])
    print(len(time)//1.5)
    gg.plot_a_frame(x,[Ut_approx,U_a],time,nframe=int(len(time)//6),labels=['numeric','analytic'],
               title='Comparision of analytic and numeric kink',filename='kink_mid.png',ls=['-',':'],ylim=[-7,10])
    gg.plot_a_frame(x,[Ut_approx,U_a],time,nframe=1,labels=['numeric','analytic'],
               title='Comparision of analytic and numeric kink',filename='kink_start.png',ls=['-',':'],ylim=[-7,10])    
    gg.plot_a_frame(x,[Ut_approx,U_a],time,nframe=-1,labels=['numeric','analytic'],
               title='Comparision of analytic and numeric kink',filename='kink_end.png',ls=['-',':'],ylim=[-7,10])
    #gg.generate_gif(x,Ut_approx,time,name,U_a=U_a,title=title,y_min=-10,y_max=10,
                    #msr=True, Umsr = Ut_approx[tx,x_idx],xmsr=x[x_idx])
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