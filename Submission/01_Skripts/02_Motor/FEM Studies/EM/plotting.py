import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv(r'eMotor\FEM Studies\EM\Optim_Results_1.csv',  sep = '\t')

print(data)

fig, ax = plt.subplots()

scatter = ax.scatter(data.P, data.Mass, c=data.EFF, cmap="cividis")
ax.scatter(2445, 0.75, s=80, c='red', marker=(5, 0))

legend1 = ax.legend(*scatter.legend_elements(num=6),
                    loc="upper left", title="Efficnency [%]")
ax.add_artist(legend1)

ax.set_title('Results of Optimization')
ax.set_ylabel('Motor mass [kg]')
ax.set_xlabel('Motor power [W]')

ax.grid(True)

plt.savefig(r'D:\Python_Prejects\GitHub\e-Pump-RPII\eMotor\FEM Studies\EM\Optim_Scatter.pdf')
plt.show()