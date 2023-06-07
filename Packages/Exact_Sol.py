from ngsolve import *
import numpy as np

sech = lambda x_: 1/cosh(x_)
tanh = lambda x_: sinh(x_)/cosh(x_)
sol_M = 1.2
sol_N = 1
J = np.arctanh(2*sol_M*sol_N/(sol_M**2+sol_N**2))

chx = (exp(x)+exp(-x))/2
Exact_u_Dict = {
    'Order2_Soliton': lambda t: (4+3*(exp(8*1j*t)-1)/chx**2)/(
            2*chx+3/4*(cos(8*t)-1)/chx**3
        )*exp(1j*t),
    'Standard_Soliton': lambda t: 2/(exp(x+4*t)+exp(-x-4*t))*exp(-1j*(2*x+3*t)),
    'Bi_Soliton_T_Inv': lambda t: (
            exp(-1j*sol_M**2*(2-t))*sol_M*sech(sol_M*x) - exp(-1j*sol_N**2*(2-t))*sol_N*sech(sol_N*x))/(
            np.cosh(J) - np.sinh(J)*(tanh(sol_M*x)*tanh(sol_N*x) + cos((sol_M**2-sol_N**2)*(2-t))*sech(sol_M*x)*sech(sol_N*x))
        ),
    'Bi_Soliton': lambda t: (
            exp(1j*sol_M**2*t)*sol_M*sech(sol_M*x) - exp(1j*sol_N**2*t)*sol_N*sech(sol_N*x))/(
            np.cosh(J) - np.sinh(J)*(tanh(sol_M*x)*tanh(sol_N*x) + cos((sol_M**2-sol_N**2)*t)*sech(sol_M*x)*sech(sol_N*x))
        )
}