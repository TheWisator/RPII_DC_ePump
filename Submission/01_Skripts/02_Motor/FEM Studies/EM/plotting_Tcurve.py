import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data1 = pd.read_csv(r'eMotor\FEM Studies\EM\Efficiency map\Efficiency map_torque_curve.txt', sep='\t', header=None)
data2 = pd.read_csv(r'eMotor\FEM Studies\EM\Efficiency map\Efficiency map_power_curve.txt', sep='\t', header=None)

fig, ax = plt.subplots()
p1 = ax.plot(data1.iloc[0,:], data1.iloc[1,:], label='Torque', color='k')
ax.set_xlim(0,75000)
ax.set_ylim(0,0.65)
ax.set_xlabel('Speed [rpm]')
ax.set_ylabel('Torque [Nm]')
ax.set_title('Torque over Speed Characteristic')
ax2 = ax.twinx()
p2 = ax2.plot(data2.iloc[0,:], data2.iloc[1,:], label='Input Power')
p3 = ax2.plot(data2.iloc[0,:], data2.iloc[2,:], label='Output Power')
ax2.set_ylabel('Power [W]')
ax2.set_ylim(0,3500)
lns = p1+p2+p3
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=4)

plt.savefig(r'D:\Python_Prejects\GitHub\e-Pump-RPII\eMotor\FEM Studies\EM\T_over_n.pdf')
plt.show()