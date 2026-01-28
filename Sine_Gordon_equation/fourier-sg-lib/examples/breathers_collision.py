import numpy as np
import fourier_sg_lib.drawing_tool as gg
import fourier_sg_lib.mathematics.sg_soliton_functions as ssf
import fourier_sg_lib.mathematics.two_soliton_start as two_s
import fourier_sg_lib.numerics.split_step_Fourier as split_f
from functools import partial
import fourier_sg_lib.numerics.measure as msr

def shift(Ut_approx,U_a,space_grid,time,soliton_params,speed,w,gamma):
    before_col = np.argwhere(time>(soliton_params['x0']/speed-5))[0][0]
    after_col =np.argwhere(time>(soliton_params['x0']/speed+6))[0][0]
    d1,d2 = msr.measure_shift(Ut_approx,space_grid['x'],time,before_col,after_col,soliton_params['x0'])
    d1a,d2a = msr.measure_shift(U_a,space_grid['x'],time,before_col,after_col,soliton_params['x0'])
    print(f"Shift measured on analytic data: {d1a.mean()},{d2a.mean()}")                          
    print(f"Shift perturb: {d1.mean()}, {d2.mean()}")
    print(f"Shift theoretical: {ssf.pos_shift(speed,w/gamma)}")

def starting_point():
    params = ssf.init_sol("examples/params.json")
    time_grid = params['timeGrid']
    # Space grid
    space_grid = params['spaceGrid']
    soliton_params = params['solitonParams']
    method = params['method']
    #method = split_f.strang_split_step
    for w,speed in zip(soliton_params['omegaVec'],soliton_params['speedVec']):
        gamma = 1/np.sqrt(1-speed**2)
        # Starting point
        u0,v0,name,title,f_a,kfilter = two_s.breather_breather_start(space_grid['x'],speed,w,soliton_params['x0'])
        # Integration
        Ut_approx,time = split_f.integration_loop(u0,v0,time_grid['steps'],time_grid['saveEvery'],
                                                time_grid['dt'],space_grid['k'],kfilter,method,
                                                nonlinear_step=partial(split_f.perturb_sg,eps1=0,eps2=1e-5))
        #time = np.linspace(0,T,steps//save_every)
        U_a = np.array([f_a(space_grid['x'],t,speed,w,soliton_params['x0']) for t in time])
        U_marker = ssf.marker(space_grid['x'],time,speed,-soliton_params['x0'])
        shift(Ut_approx,U_a,space_grid,time,soliton_params,speed,w,gamma)
        
        gg.plot_cm(space_grid['x'],[Ut_approx],time,U_marker=U_marker,two_breahters=True,xlim=[-20,20],filename=f'examples/breathers_colisions_cm_pertub_v{speed}_f{w}.png')
        gg.plot_cm(space_grid['x'],[U_a],time,U_marker=U_marker,two_breahters=True,xlim=[-20,20],filename=f'examples/breathers_colisions_cm_analytic_v{speed}_f{w}.png')
        #gg.plot_3d(x,[Ut_approx],time,'breathers_colisions_3d.png')
        gg.generate_gif(space_grid['x'],Ut_approx,time,U_a=U_a,filename=f'examples/breather_collision_perturb_v{speed}_f{w}.gif',title='',
                        mark=True,Umark=U_marker)
starting_point()