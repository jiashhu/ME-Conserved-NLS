import sympy as sym
import numpy as np
try: 
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except:
    pass

aim_code_set = [
    'draw amplitude',
    'check exact solution'
]
aim_code = 'check exact solution'

x = sym.Symbol('x',real=True)
t = sym.Symbol('t',real=True)
chx = (sym.exp(x)+sym.exp(-x))/2
exact_u_ng = (
    4+3*(sym.exp(8*sym.I*t)-1)/chx**2
)/(
    2*chx+3/4*(sym.cos(8*t)-1)/chx**3
)*sym.exp(sym.I*t)

diff_terms = sym.lambdify((x,t),exact_u_ng.diff(t)*sym.I + exact_u_ng.diff(x).diff(x))
u_numpy = sym.lambdify((x,t),exact_u_ng)
L = 5

## draw the amplitude (np.pi/8 periodic)
try:
    if aim_code == 'draw amplitude':
        x_arr = np.linspace(-L,L,2**14)
        t_arr = np.linspace(0,np.pi/4,400)
        X, Y = np.meshgrid(x_arr,t_arr)
        Z = abs(u_numpy(X,Y))
        figure = plt.figure()
        axes = Axes3D(figure)
        axes.plot_surface(X,Y,Z,cmap = 'rainbow')
        plt.show()
    elif aim_code == 'check exact solution':
        x_arr = np.linspace(-L,L,2**14)
        t_arr = np.linspace(0,np.pi/4,2**10)
        res = []
        for t_val in t_arr:
            residue = diff_terms(x_arr,t_val) + 2*abs(u_numpy(x_arr,t_val))**2*u_numpy(x_arr,t_val)
            res.append(max(abs(residue)))
        print(np.linalg.norm(res))
except:
    pass