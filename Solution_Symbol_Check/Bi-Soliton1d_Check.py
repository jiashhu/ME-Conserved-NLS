import numpy as np
import sympy as sym
try: 
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except:
    pass

sech = lambda x_: 1/np.cosh(x_)
aim_code_set = [
    'draw amplitude',
    'check exact solution'
]
aim_code = 'check exact solution'
M = 1.2
N = 1

x1 = sym.Symbol('x')
x2 = x1 + 0.4
t = sym.Symbol('t')
S = (M**2-N**2)*t
J = np.arctanh(2*M*N/(M**2+N**2))

exact_u_ng = (
    sym.exp(sym.I*M**2*t)*M*sym.sech(M*x1) - \
    sym.exp(sym.I*N**2*t)*N*sym.sech(N*x2)
    )/(
    sym.cosh(J) - sym.sinh(J)*(sym.tanh(M*x1)*sym.tanh(N*x2) + sym.cos(S)*sym.sech(M*x1)*sym.sech(N*x2))
)

diff_terms = sym.lambdify((x1,t),exact_u_ng.diff(t)*sym.I + exact_u_ng.diff(x1).diff(x1),modules=['numpy',{'sech':sech}])

u_numpy = sym.lambdify((x1,t),exact_u_ng,modules=['numpy',{'sech':sech}])
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
        t_arr = np.linspace(0,1,2**10)
        res = []
        for t_val in t_arr:
            residue = diff_terms(x_arr,t_val) + 2*abs(u_numpy(x_arr,t_val))**2*u_numpy(x_arr,t_val)
            res.append(max(abs(residue)))
        print(np.linalg.norm(res))
except:
    pass
