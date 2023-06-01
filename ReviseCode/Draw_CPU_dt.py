import numpy as np
import matplotlib.pyplot as plt
from Package_MyCode import *

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=10,fontsize=16)

PP_CPU_time = np.array([
    0.1130828857421880, 0.26062631607055700, 0.4572722911834720, 0.722142219543457, 
    1.4711434841156000, 3.074068069458010, 6.351077079772950, 12.293888092041000, 
    24.52998375892640, 49.70548748970030
])

Time_Step = np.array([
    0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 
    0.0009765625, 0.00048828125, 0.000244140625	
])

CPU_ALL = np.array([
    6.0, 12.0, 19.0, 30.0, 51.0, 104.0, 167.0, 335.0, 671.0, 1035.0
])

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.semilogx(Time_Step,PP_CPU_time,'o-',label='The post-processing correction step')
ax.semilogx(Time_Step,CPU_ALL-PP_CPU_time,'^--',label='The standard collocation step')
plt.legend(loc='upper right')
plt.grid(True,ls='--')
ax.set_xlabel('$\\tau$')
ax.set_ylabel('CPU time(s)')
plt.show()

SavePath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Pic'
fig.savefig(os.path.join(SavePath,'{}.pdf'.format('Example1CPUvsdt')),dpi=600,format='pdf')