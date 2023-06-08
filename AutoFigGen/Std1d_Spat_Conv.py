import json
import matplotlib.pyplot as plt
from IO_Func import *

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle()

f_path = '../Numerical_Tests/1d_Standard_Soliton/Spat_Conv/Sum_Spat_Conv.json'
# Read Convergece Errors and Plot
with open(f_path, "r") as json_file:
    Res = json.load(json_file)

Res['Params'] = {
        "Legendlabel": "$\\tau$",
        "xlabel": "h",
        "T": 1,
        "$\\tau$": "1e-3",
        "L": 40,
        "N_iter_tol": 1e-10,
        "Name": "StandardSoliton_Spatial",
        "Err": "err",
        "Param": "param",
        "Param_Total": 80,
        "add_order": 0,
        "var": "p"
    }
Legendlabel = Res['Params']['Legendlabel']

_order_key = ['1','2','3']
Zoom_xy_scale = [
    {'x': [1/2,2], 'y': [1/2,2]},
    {'x': [1/3,3], 'y': [1/3,3]},
    {'x': [1/4,4], 'y': [1/4,4]}
]
Val_List = [
    [[0.02,0.05,0.2], [3e-2,6e-2,1e-1,2e-1,3e-1]],
    [[0.01,0.1,1],    [8e-4,4e-3,1e-2,4e-2,0.08]],
    [[0.01,0.1,1,10], [5e-5,1e-4,1e-3,1e-2,0.05]]
]
XYTick_scale = GenerateXYTick(Val_List)

Zoom=dict(zip(_order_key,Zoom_xy_scale))
XYTick=dict(zip(_order_key,XYTick_scale))
proposed_md = '{}={}'.format(Legendlabel,Res['Params'][Legendlabel])
Legend_scale = [
    [proposed_md,'$\mathcal{O}(h)$'],
    [proposed_md,'$\mathcal{O}(h^2)$'],
    [proposed_md,'$\mathcal{O}(h^3)$'],
]
Legend=dict(zip(_order_key,Legend_scale))

ERR_Plot = CanonicalErrPlot(Res,Legend=Legend,Zoom=Zoom,opt='centeral',adj_p=0.7)
ERR_Plot.ReadParam()
fig = plt.figure(figsize=(5*ERR_Plot.n_order,5))
ERR_Plot.DrawParam(fig)
plt.show()