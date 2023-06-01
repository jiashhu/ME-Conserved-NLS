import numpy as np
import sympy as sym

x = sym.Symbol('x')
y = sym.Symbol('y')
t = sym.Symbol('t')
exact_u_Sym = 4/(2*sym.exp(x+y-2*np.sqrt(2)*t)+sym.exp(-(x+y-2*np.sqrt(2)*t)))*sym.exp(sym.I*(t+(x+y)*np.sqrt(2)/2))
Norm_u = 4/(2*sym.exp(x+y-2*np.sqrt(2)*t)+sym.exp(-(x+y-2*np.sqrt(2)*t)))

Equ_Left = sym.I*exact_u_Sym.diff(t) + exact_u_Sym.diff(x).diff(x) + exact_u_Sym.diff(y).diff(y)
Nonlinear_term = 2*Norm_u**2*exact_u_Sym

residue = Equ_Left + Nonlinear_term

numerical_res = sym.lambdify((x,y,t),residue)

print('good')

# exact_u_2d = lambda t: 4/(2*exp(x+y-2*t)+exp(-(x+y-2*t)))*exp(1j*(t+(x+y)/2))