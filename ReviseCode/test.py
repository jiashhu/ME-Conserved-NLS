import numpy as np
sech = lambda x_: 1/np.cosh(x_)
tanh = lambda x_: np.sinh(x_)/np.cosh(x_)
sol_M = 1.2
sol_N = 1
J = np.arctanh(2*sol_M*sol_N/(sol_M**2+sol_N**2))
func = lambda t,x: (
            np.exp(1j*sol_M**2*t)*sol_M*sech(sol_M*x) - np.exp(1j*sol_N**2*t)*sol_N*sech(sol_N*x))/(
            np.cosh(J) - np.sinh(J)*(tanh(sol_M*x)*tanh(sol_N*x) 
                                     + np.cos((sol_M**2-sol_N**2)*t)*sech(sol_M*x)*sech(sol_N*x))
        )
print(sech(0))