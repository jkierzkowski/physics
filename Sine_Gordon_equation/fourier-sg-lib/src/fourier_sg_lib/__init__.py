# __init__.py

# --- Public imports from submodules ---
from .drawing_tool import (
    generate_gif,
    giff_comp,
    plot_a_frame,
    plot_3d,
    plot_cm,
    plt_eps,

)
from .mathematics.two_soliton_start import (
    kink_kink_start,
    kink_start,
    kink_antikink_start,
    stat_breather_start,
    breather_start,
    breather_breather_start,


)
from .mathematics.sg_soliton_functions import (
    kink,
    kink_t,
    antikink,
    antikink_t,
    breather,
    breather_t,
    two_breathers,
    pos_shift,
)

from .numerics.split_step_Fourier import (
    integration_loop,
    strang_split_step,
    step_4th_order,
    approx_nonlinear_step,
)

# --- Metadata ---
__version__ = "1.1"
__author__ = "Jan Kierzkowski"

# --- Public API ---
__all__ = [
    "__version__",
    "__author__",
    "generate_gif",
    "giff_comp",
    "plot_a_frame",
    "plot_3d",
    "plot_cm",
    "plt_eps",
    "kink_kink_start",
    "kink_start",
    "kink_antikink_start",
    "stat_breather_start",
    "breather_start",
    "breather_breather_start",
    "kink",
    "kink_t",
    "antikink",
    "antikink_t",
    "breather",
    "breather_t",
    "two_breathers",
    "pos_shift",
    "integration_loop",
    "strang_split_step",
    "step_4th_order",
    "approx_nonlinear_step",
]
