import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

data = pd.read_csv(r'eMotor\FEM Studies\EM\Efficiency map\Efficiency map_efficiency_map.txt', sep='\t', header=None)
print(data)

RF_pow = 1.0
P_Lox = [RF_pow*143, RF_pow*1760]
n_Lox = [20260, 38800]
T_Lox = [P_Lox[i] / (n_Lox[i] * 2 * math.pi / 60) for i in range(0, len(n_Lox))]
P_prop = [RF_pow*138.7, RF_pow*1403]
n_prop = [18930, 38000]
T_prop = [P_prop[i] / (n_prop[i] * 2 * math.pi / 60) for i in range(0, len(n_prop))]



T_min = 0
T_max = 0.6038
T_length = data.shape[0] 
Tdata = np.linspace(T_min, T_max, T_length)
n_min = 0
n_max = 75000
n_length = data.shape[1]
ndata = np.linspace(n_min, n_max, n_length)
contour_levels = np.arange(start=60,stop=96,step=0.5)
fig, ax = plt.subplots()
cf = ax.contourf(ndata, Tdata, 100 * data, cmap='cividis', levels = contour_levels)
cs = ax.contour(ndata, Tdata, 100 * data, [60,70,80,85,90,92,93,94,95], colors='k')
ax.clabel(cs, inline=True, fontsize=10)

lox = ax.plot(n_Lox, T_Lox, 'o', ls='-', ms=5, color='b', linewidth=3, label='LOx')
prop = ax.plot(n_prop, T_prop, 'o', ls='-', ms=5, color='r', linewidth=3, label='Propane')
lns = lox+prop
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=1)

cbar = fig.colorbar(cf)
cbar.ax.set_ylabel('Efficiency [%]')

ax.set_title('Efficiency Map of the Motor')
ax.set_xlabel('Speed [rpm]')
ax.set_ylabel('Torque [Nm]')


plt.savefig(r'D:\Python_Prejects\GitHub\e-Pump-RPII\eMotor\FEM Studies\EM\Eff_Map.pdf')
plt.show()