from . import Gauss_Legendre as GL
from .BaseMethod import (
    CPN_Converg_1d_Focus_Newt, 
    CPN_Converg_2d_Newt,
    CPNLS_SemiImpLin,
    CPNLS_SemiImpLinfu
)

class Lap_Type():
    peri_opt = 'peri'
    Diri_opt = 'Dirichlet'

__all__ = [
    'CPN_Converg_1d_Focus_Newt',
    'GL',
    'CPN_Converg_2d_Newt',
    'CPNLS_SemiImpLin',
    'CPNLS_SemiImpLinfu'
]