import numpy as np
import fourier_sg_lib.mathematics.sg_soliton_functions as ssf
import fourier_sg_lib.mathematics.two_soliton_start as two_s
import fourier_sg_lib.numerics.split_step_Fourier as split_f
import fourier_sg_lib.drawing_tool as gg
import pytest
import os
params = ssf.init_sol("tests/params.json")
soliton_params = params['solitonParams']
time_grid = params['timeGrid']
# Space grid
space_grid = params['spaceGrid']
method = params['method']
# Starting point
u0,v0,name,title,f_a,kfilter = two_s.breather_breather_start(space_grid['x'],soliton_params['speedVec'][0],soliton_params['omegaVec'][0],soliton_params['x0'])
# Integration
Ut,time = split_f.integration_loop(u0,v0,2,time_grid['saveEvery'],
                                            time_grid['dt'],space_grid['k'],kfilter,method)
def test_integrate():
    calculated = np.column_stack((time,Ut))
    assert Ut.shape[0] == time.shape[0]
    assert np.isfinite(Ut).all()
    assert np.isfinite(time).all()
    ref_path = "tests/numerics.txt"
    if os.path.exists(ref_path):
        loaded = np.loadtxt("tests/numerics.txt")
        assert np.allclose(
            loaded,
            calculated,
            rtol=1e-12,
            atol=1e-14
            )
def test_giff():
    gg.generate_gif(space_grid['x'],Ut,time,name)
    assert os.path.exists(name)
    os.remove(name)