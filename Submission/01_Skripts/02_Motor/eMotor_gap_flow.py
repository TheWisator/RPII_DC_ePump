## Python package of functions to calculate losses due to gap flow in electrical machines
# source: Bilgen & Boulos: Functional Dependence of Torque Coefficient of Coaxial Cylinders on Gap Width and Reynolds Numbers
# DOI: https://doi.org/10.1115/1.3446944 


import CoolProp.CoolProp as CP
import math
import numpy as np
import matplotlib.pyplot as plt

def fluid_loss_CP(speed_rpm:float, r_rotor_m:float, l_gap_m:float, l_motor_m:float, fluid_name, T_fluid_K:float, p_fluid_bar:float) -> float:
    ## get fluid data with coolprop
    p_Pa = p_fluid_bar * 1e5
    kin_vis_Pas = CP.PropsSI('V', 'P', p_Pa, 'T', T_fluid_K, fluid_name)
    rho_kg_m3 = CP.PropsSI('D', 'P', p_Pa, 'T', T_fluid_K, fluid_name)
    dyn_vis_m2_s = kin_vis_Pas / rho_kg_m3

    # get angular velocity
    omega_rad_s = speed_rpm / (2 * math.pi * 60)
    # reynolds number of gap
    Re = omega_rad_s * r_rotor_m * l_gap_m / dyn_vis_m2_s
    # taylor number
    Ta = omega_rad_s * rho_kg_m3 * r_rotor_m * l_gap_m / kin_vis_Pas * math.sqrt(l_gap_m / r_rotor_m)
    # critical taylor number
    Ta_c = 41.3
    # ratio of gap length and rotor radius
    LR = l_gap_m / r_rotor_m 
    # "constant calculated from the theory"
    D = 1.4472

    ## estimate losses of the fluid flow with correlations for different renolds numbers
    if Re < 64:
        # linear regime
        cm = 10 * LR ** 0.3 * Re ** (-1)
    elif Re >= 64 and Re < 500 and LR > 0.07:
        # transition regime case 1
        cm = 10 * LR ** 0.3 * Re ** (-0.6)
        print('Transition regime case 1')
    elif Re >= 64 and Re < 500 and LR <= 0.07:
        # transition regime case 2
        psi = (2 + LR) * (1 + D * (1 - (Ta_c / Ta) ** (2)))
        cm = 2 * psi / Re
        print('Transition regime case 2')
    elif Re >= 500 and Re < 1e4:
        # turbulent regime 1
        cm = 1.03 * LR ** 0.3 * Re ** (-0.5)
    elif Re >= 1e4:
        # turbulent regime case 2
        cm = 0.065 * LR ** 0.3 * Re ** (-0.2)
    else:
        return("ERROR: reynolds number out of range!")

    P_loss_W = 0.5 * math.pi * omega_rad_s ** 3 * rho_kg_m3 * r_rotor_m ** 4 * l_motor_m * cm

    return [P_loss_W, Re, Ta]


# motor specification
l_gap = 0.00055
d_rot = 0.0192
l_mot = 0.03
omega = np.linspace(1,45000)

res_holder = [fluid_loss_CP(i, d_rot/2, l_gap, l_mot, 'oxygen', 90, 30) for i in omega]
P_holder_1 = [res_holder[i][0] for i in range(0,len(res_holder))]

res_holder = [fluid_loss_CP(i, d_rot/2, l_gap, l_mot, 'propane', 90, 32) for i in omega]
P_holder_2 = [res_holder[i][0] for i in range(0,len(res_holder))]
#Re_holder = [res_holder[i][1] for i in range(0,len(res_holder))]
#Ta_holder = [res_holder[i][2] for i in range(0,len(res_holder))]
#print(res_holder)
#
fig, ax = plt.subplots()
p1 = ax.plot(omega, P_holder_1, label='LOx')
p2 = ax.plot(omega, P_holder_2, label='Propane')

lns = p1+p2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=2)

ax.set(xlabel='Speed [rpm]', ylabel='Losses [W]',
       title='Viscous losses in machine air gap')
ax.grid()

plt.savefig(r'D:\Python_Prejects\GitHub\e-Pump-RPII\eMotor\Visc_Losses.pdf')
plt.show()