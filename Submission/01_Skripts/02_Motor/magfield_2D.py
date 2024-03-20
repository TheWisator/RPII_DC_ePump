import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

"""
Script to solve the 2D magnetic field distribution of permanent magnet (PM) excited synchronour inner rotor (IR) electrical machines with surface permanent magnets (SPM)

sources:
[1]: Bolte, E: Elektrische Maschinen (2018), pp. 643 ff.
[2]: Mezani, S. et. al: Sizing Electrical Machines Using OpenOffice Calc, 9th EUCASS (2022)
[3]: Zhu, Z. K. et al: Instantaneous Magnetic Field Distribution in Brushless Permanent Magnet dc Motors, Part I: Open-circuit Field  (1993)


(c) Thilo WITZEL, TU Munich, 2023

"""


MU_0 = 4 * math.pi * 1e-07      # permeability of vacuum in H/m
FMAX = 100                     # maximum amount of Fourier summation terms

def F2(mu:int, mu1:float, mu2:float, r1_m:float):
    # Bolte page 167, equation 2.123
    return (mu1 + mu2) / (mu1 - mu2) * r1_m ** (-2 * mu)

def F4(mu:int, mu4:float, mu5:float, mu6:float, r4_m:float, r5_m:float):
    # Bolte page 169, equation 2.136
    b = (mu6 - mu5) / (mu6 + mu5) * (r4_m / r5_m) ** (2*mu)
    a = (mu4 * (b - 1)) / (mu5 * (b + 1))

    return (1 + a) / (1 - a) * r4_m ** (-2 * mu)

def EFG3(mu:int, F2:float, F4:float, r2_m:float, r3_m:float, mu2:float, mu3:float, mu4:float):
    # Bolte page 168, equations 2.128 ff.

    # useful constants for calculation
    r2_mu = r2_m ** (-mu + 1) / (mu ** 2 - 1)
    r3_mu = r3_m ** (-mu + 1) / (mu ** 2 - 1)
    mu_32 = mu3 / mu2
    mu_34 = mu3 / mu4
    f2 = (F2 + r2_m ** (-2 * mu)) / (F2 - r2_m ** (-2 * mu))
    f4 = (F4 + r3_m ** (-2 * mu)) / (F4 - r3_m ** (-2 * mu))


    Det = -r2_m ** (-2 * mu) * (f2 / mu_32 + 1) * (f4 / mu_34 - 1) + r3_m ** (-2 * mu) * (f4 / mu_34 + 1) * (f2 / mu_32 - 1) 
    DetG3 = -r2_m ** (-2 * mu) * (f2 / mu_32 + 1) * r3_mu * (1 - f4 / (mu * mu_34)) + r3_m ** (-2 * mu) * (f4 / mu_34 + 1) * r2_mu * (1 - f2 / (mu * mu_32))
    DetE3 = r2_mu * (1 - f2 / (mu * mu_32)) * (f4 / mu_34 - 1) - r3_mu * (1 - f4 / (mu * mu_34)) * (f2 / mu_32 - 1)
    
    E3 = DetE3 / Det
    F3 = DetE3 / DetG3
    G3 = DetG3 / Det
    return [E3, F3, G3]

def E4(mu:int, pf:int, F4:float, r3_m:float, G3:float, E3:float, M:float, rMag_m:float, alpha:float, Br_T):
    # bolte page 168, equation 2.131
    r3_mu = r3_m ** (-mu + 1) / (mu ** 2 - 1)
    
    return ((F4 + r3_m ** (-2 * mu)) ** (-1) * (G3 + E3 * r3_m ** (-2 * mu) + r3_mu) * MU_0 * mu * mu_M(mu, pf, rMag_m, alpha, Br_T))

def mu_M(mu:int, pf:int, rMag_m:float, alpha:float, Br_T:float):
    # Bolte page 156, equations 2.216 ff.
    
    # assumption for magnetization (may be not appropriate if magnet is operated at large negative H, only applicable if PM has approximately constant J(H) over H curve) 
    Mref = Br_T / MU_0
    tau_pf = math.pi * rMag_m / pf
    yf = alpha * tau_pf                     # pole coverage alpha needs to be between 0 and 1
    bf = 0.25 * (1-alpha) * yf
    
    
    return Mref * pf * 4 / (mu * math.pi) * np.sinc(mu * bf / (pf * 2 * tau_pf)) * np.sin(math.pi * mu * yf / (pf * 2 * tau_pf))

def A3(r:float, phi:float, pf:int, mu_rotor:float, mu_stator:float, mu_mag:float, r1:float, r2:float, r3:float, h_gap:float, r5:float, alpha:float, Br_T:float):
    A3 = 0
    for a in range(0, FMAX):
        mu = pf * (2* a + 1)
        f2 = F2(mu, 1, mu_rotor, r1)
        f4 = F4(mu, 1, mu_stator, 1, r3 + h_gap, r5)
        ef3 = EFG3(mu, f2, f4, r2, r3, mu_rotor, mu_mag, 1)
        m = mu_M(mu, pf, r3, alpha, Br_T)
        
        A3 += (ef3[0] * (ef3[1] * r ** (mu) + r ** (-mu)) + r / (mu ** 2 - 1)) * MU_0 * mu * m * math.sin(mu * phi)
        
    return A3

def A4(r:float, phi:float, pf:int, mu_rotor:float, mu_stator:float, mu_mag:float, r1:float, r2:float, r3:float, h_gap:float, r5:float, alpha:float, Br_T:float):
    A4 = 0
    for a in range(0, FMAX):
        mu = pf * (2* a + 1)
        f2 = F2(mu, 1, mu_rotor, r1)
        f4 = F4(mu, 1, mu_stator, 1, r3 + h_gap, r5)
        efg3 = EFG3(mu, f2, f4, r2, r3, mu_rotor, mu_mag, 1)
        m = mu_M(mu, pf, r3, alpha, Br_T)
        e4 = E4(mu, pf, f4, r3, efg3[2], efg3[0], m, ((r2 + r3)/2), alpha, Br_T)
        
        
        A4 += e4 * (f4 * r ** mu + r ** (-mu)) * math.sin(mu*phi)
    return A4

def Br3Hat_T(r:float, pf:int, mu_rotor:float, mu_stator:float, mu_mag:float, r1:float, r2:float, r3:float, h_gap:float, r5:float, alpha:float, Br_T:float):
    Br3 = 0
    for a in range(0, FMAX):
        mu = pf * (2* a + 1)
        f2 = F2(mu, 1, mu_rotor, r1)
        f4 = F4(mu, 1, mu_stator, 1, r3 + h_gap, r5)
        ef3 = EFG3(mu, f2, f4, r2, r3, mu_rotor, mu_mag, 1)
        m = mu_M(mu, pf, r3, alpha, Br_T)
        
        Br3 += (ef3[0] * (ef3[1] * r ** (mu) + r ** (-mu)) + r / (mu ** 2 - 1)) * MU_0 * mu ** 2 * m / r
        
    return Br3

def Br3_T(phi:float, r:float, pf:int, mu_rotor:float, mu_stator:float, mu_mag:float, r1:float, r2:float, r3:float, h_gap:float, r5:float, alpha:float, Br_T:float):
    Br3 = 0
    for a in range(0, FMAX):
        mu = pf * (2* a + 1)
        f2 = F2(mu, 1, mu_rotor, r1)
        f4 = F4(mu, 1, mu_stator, 1, r3 + h_gap, r5)
        ef3 = EFG3(mu, f2, f4, r2, r3, mu_rotor, mu_mag, 1)
        m = mu_M(mu, pf, r3, alpha, Br_T)
        
        Br3 += (ef3[0] * (ef3[1] * r ** (mu) + r ** (-mu)) + r / (mu ** 2 - 1)) * MU_0 * mu ** 2 * m / r * math.cos(mu * phi)
        
    return Br3

def Br4Hat_T(r:float, pf:int, mu_rotor:float, mu_stator:float, mu_mag:float, r1:float, r2:float, r3:float, h_gap:float, r5:float, alpha:float, Br_T:float):
    Br4 = 0
    for a in range(0, FMAX):
        mu = pf * (2* a + 1)
        f2 = F2(mu, 1, mu_rotor, r1)
        f4 = F4(mu, 1, mu_stator, 1, r3 + h_gap, r5)
        efg3 = EFG3(mu, f2, f4, r2, r3, mu_rotor, mu_mag, 1)
        m = mu_M(mu, pf, r3, alpha, Br_T)
        e4 = E4(mu, pf, f4, r3, efg3[2], efg3[0], m, ((r2 + r3)/2), alpha, Br_T)
        
        Br4 += e4 * (f4 * r ** mu + r ** (-mu)) / r
        
    return Br4

def Br4_T(phi:float, r:float, pf:int, mu_rotor:float, mu_stator:float, mu_mag:float, r1:float, r2:float, r3:float, h_gap:float, r5:float, alpha:float, Br_T:float):
    Br4 = 0
    for a in range(0, FMAX):
        mu = pf * (2* a + 1)
        f2 = F2(mu, 1, mu_rotor, r1)
        f4 = F4(mu, 1, mu_stator, 1, r3 + h_gap, r5)
        efg3 = EFG3(mu, f2, f4, r2, r3, mu_rotor, mu_mag, 1)
        m = mu_M(mu, pf, r3, alpha, Br_T)
        e4 = E4(mu, pf, f4, r3, efg3[2], efg3[0], m, ((r2 + r3)/2), alpha, Br_T)
        
        Br4 += e4 * (f4 * r ** mu + r ** (-mu)) / r * math.cos(mu * phi)
        
    return Br4

def PhiPf(L_m:float, pf:int, mu_rotor:float, mu_stator:float, mu_mag:float, r1:float, r2:float, r3:float, h_gap:float, r5:float, alpha:float, Br_T:float):
    r_eval = r3+h_gap / 2 # evaluate integral in the middle of the airgap
    
    def integration_function(x:float, r:float, pf:int, mu_rotor:float, mu_stator:float, mu_mag:float, r1:float, r2:float, r3:float, h_gap:float, r5:float, alpha:float, Br_T:float):
        return abs(Br3_T(x, r, pf, mu_rotor, mu_stator, mu_mag, r1, r2, r3, h_gap, r5, alpha, Br_T) * r_eval)
    
    
    return sp.integrate.quad(integration_function, -math.pi / (2*pf), math.pi / (2*pf), args=(r_eval, pf, mu_rotor, mu_stator, mu_mag, r1, r2, r3, h_gap, r5, alpha, Br_T))[0] / (math.pi * r_eval / pf)

def Brms_T(pf:int, Br_T:float, r2_m:float, r3_m:float, h_gap:float):
    # source [2] and [3]
    r = r3_m + h_gap
    
    num = (4*Br_T * pf * (pf-1) * r3_m**(2*pf) - (pf + 1) * r2_m**(2*pf) + 2*r3_m**(pf-1)*r2_m**(pf+1))
    den = (math.pi * math.sqrt(2) * (pf**2 - 1) * (r**(2*pf) - r2_m**(2*pf)))
    return num / den * ((r) / r3_m)**(pf-1)

def rmsValue(arr, n):
    square = 0
    mean = 0.0
    root = 0.0
     
    #Calculate square
    for i in range(0,n):
        square += (arr[i]**2)
     
    #Calculate Mean 
    mean = (square / (float)(n))
     
    #Calculate Root
    root = math.sqrt(mean)
     
    return root

def torque2D(A_rms_Am:float, pf:int, roRotor_m:float, roMag_m:float, h_gap:float, Brem_T:float, L_m: float, k_fe:float):
    Brms = Brms_T(pf, Brem_T, roRotor_m, roMag_m, h_gap)
    #print(f'RMS air gap flux density: {Brms} [T]')
    
    # assumption: winding factor ~0.85...0.95 -> winding factor 0.9
    # assumption: power factor = 0.9

    return 2.0 * Brms * A_rms_Am * math.pi * (roMag_m + h_gap)**2 * L_m * k_fe * 0.9 * 0.9

def length_for_torque(Treq_Nm:float, A_rms_Am:float, pf:int, roRotor_m:float, roMag_m:float, h_gap:float, Brem_T:float, k_fe:float):
    Brms = Brms_T(pf, Brem_T, roRotor_m, roMag_m, h_gap)
    #print(f'RMS air gap flux density: {Brms} [T]')
    #print(f'sigma = {Brms * A_rms_Am} [N/m²]')

    # assumption: winding factor ~0.85...0.95 -> winding factor 0.9
    # assumption: power factor = 0.9
        
    return Treq_Nm / (2.0 * 0.9 * 0.9 * k_fe* Brms * A_rms_Am * math.pi * (roMag_m + h_gap)**2)

def Mn(Br:float, alpha_p:float, n:int) -> float:
    ## function to calculate normal magnet magnetization acc to [3] eqn 7
    # @ inputs: 
    #   Br: remanence flux density in T
    #   alpha_p: magnet pole arc to pole pitch ratio
    #   n: integer (1, 2, ...)
    # @output: Magnetization in T
    return 4 * Br * alpha_p / MU_0 * math.sin(n * math.pi * alpha_p / 2) / (n * math.pi * alpha_p)

def B_norm_gap(r, phi, rs_m, h_gap_m, h_mag_m, Br_T, alpha_p, ppairs, mu_rmag):
    # function to calculate normal air gap field distribution acc. to [3] eqn 43
    Rm = rs_m - h_gap_m
    Rr = rs_m - h_gap_m - h_mag_m
    rrm = Rr / Rm
    rrs = Rr / rs_m
    rms = Rm / rs_m
    rs = r / rs_m
    mr = Rm / r
    mplus = (mu_rmag + 1) / mu_rmag
    mminus = (mu_rmag - 1) / mu_rmag

    B = 0
    for i in range(1, FMAX):
        n = (2 * i - 1)
        np = (2 * i - 1) * ppairs
        np2 = 2 * np        
        B += MU_0 * Mn(Br_T, alpha_p, n) * np / (mu_rmag * (np ** 2 - 1)) * ((np - 1) + 2*rrm**(np + 1) - (np + 1) * rrm**np2) / (mplus * (1 - rrs**np2)-mminus * (rms ** np2 - rrm ** np2)) * (rs ** (np - 1) * rms**(np + 1) + mr**(np + 1)) * math.cos(np * phi)

    return B

def B_norm_gap_s(phi, rs_m, h_gap_m, h_mag_m, Br_T, alpha_p, ppairs, mu_rmag):
    # function to calculate normal air gap field distribution acc. to [3] eqn 43
    Rm = rs_m - h_gap_m
    Rr = rs_m - h_gap_m - h_mag_m
    rrm = Rr / Rm
    rsm = rs_m / Rm
    mplus = (mu_rmag + 1) / mu_rmag
    mminus = (mu_rmag - 1) / mu_rmag

    B = 0
    for i in range(1, FMAX):
        n = (2 * i - 1)
        np = (2 * i - 1) * ppairs
        np2 = 2 * np        
        B += 2 * MU_0 * Mn(Br_T, alpha_p, n) * np / (mu_rmag * (np ** 2 - 1)) * rsm ** (np - 1) * ((np - 1) * Rm ** np2 + 2* Rr ** (np+1) * Rm ** (np-1) - (np+1) * Rr**np2) / (mplus * (rs_m ** np2 - Rm ** np2) - mminus * (Rm ** np2 - rs_m ** np2 * rrm ** np2)) * math.cos(np * phi)
    return B

def plot_theor_gap_density(r1, r2, r3, r4, r5, L, ppairs, Brem, mu_iron):
    h_gap = r4 - r3
    mu_mag = 1.07
    Phi = PhiPf(L, ppairs, mu_iron, mu_iron, mu_mag, r1, r2, r3, h_gap, r5, 0.9, Brem)
    B2 = Brms_T(ppairs, Brem, r2, r3, h_gap)
    B3 = Br4Hat_T(r4, ppairs, mu_iron, mu_iron, mu_mag, r1, r2, r3, h_gap, r5, 0.9, Brem)

    phi_range = np.linspace(0,math.pi*2,1000)
    Br3 = [Br3_T(i, r3, ppairs, mu_iron, mu_iron, mu_mag, r1, r2, r3, h_gap, r5, 0.9, Brem) for i in phi_range]
    Br4 = [Br4_T(i, r4, ppairs, mu_iron, mu_iron, mu_mag, r1, r2, r3, h_gap, r5, 0.9, Brem ) for i in phi_range]
    Bavg = [(Br3[i] + Br4[i]) / 2 for i in range(0, len(phi_range))]

    plt.plot(phi_range, Br4)
    plt.plot(phi_range, Br3)
    plt.plot(phi_range, Bavg)
    plt.show()

def rotor_flux(rs_m, h_gap_m, h_mag_m, Br_T, alpha_p, ppairs, mu_rmag) -> float:
    def integrand(phi, r, rs_m, h_gap_m, h_mag_m, Br_T, alpha_p, ppairs, mu_rmag):
        return r * abs(B_norm_gap(r, phi, rs_m, h_gap_m, h_mag_m, Br_T, alpha_p, ppairs, mu_rmag))
    
    I = sp.integrate.quad(integrand, 0, 2*math.pi, args=(rs_m, rs_m, h_gap_m, h_mag_m, Br_T, alpha_p, ppairs, mu_rmag))
    return I[0]

def plot_radial_field(r_stator, h_gap, h_mag, Br_T, alpha_p, ppairs, mu_rmag):
    phi_range = np.linspace(0,math.pi*2,1000)
    r = r_stator + h_gap/2

    B = [B_norm_gap(r, ph, r_stator, h_gap, h_mag, Br_T, alpha_p, ppairs, mu_rmag) for ph in phi_range]
    #B2 = [B_norm_gap_s(ph, r_stator, h_gap, h_mag, Br_T, alpha_p, ppairs, mu_rmag) for ph in phi_range]
    plt.plot(phi_range, B)
    #plt.plot(phi_range, B2)
    print(rmsValue(B, len(B)))
    #print(rmsValue(B2, len(B2)))
    plt.show()

def torque2DFlux(A_rms_Am:float, pf:int, roRotor_m:float, roMag_m:float, h_gap:float, Brem_T:float, L_m: float, k_fe:float):
  
    flux = rotor_flux(roMag_m + h_gap, h_gap, roMag_m - roRotor_m, Brem_T, 0.95, pf, 1)

    print('Total flux: ', flux)
    return 0.9 * flux * 2 * math.pi * (roMag_m + h_gap) * A_rms_Am


'''
r1 = 0.005
r2 = 0.008
r3 = 0.0126
h_gap = 0.004
r4 = r3 + h_gap
r5 = 0.05
L = 0.1
ppairs = 4
mu_iron = 1000
mu_mag = 1.07
Brem = 1.15

Phi = PhiPf(L, ppairs, mu_iron, mu_iron, mu_mag, r1, r2, r3, h_gap, r5, 0.9, Brem)
B2 = Brms_T(ppairs, Brem, r2, r3, h_gap)
B3 = Br4Hat_T(r4, ppairs, mu_iron, mu_iron, mu_mag, r1, r2, r3, h_gap, r5, 0.9, Brem)
#B3 = rmsValue(Phi, len(Phi))

print(Phi)
print(B2)
print(B3)

phi_range = np.linspace(0,math.pi*2,1000)
Br3 = [Br3_T(i, r3, ppairs, mu_iron, mu_iron, mu_mag, r1, r2, r3, h_gap, r5, 0.9, Brem) for i in phi_range]
Br4 = [Br4_T(i, r4, ppairs, mu_iron, mu_iron, mu_mag, r1, r2, r3, h_gap, r5, 0.9, Brem ) for i in phi_range]
Bavg = [(Br3[i] + Br4[i]) / 2 for i in range(0, len(phi_range))]

plt.plot(phi_range, Br4)
plt.plot(phi_range, Br3)
plt.plot(phi_range, Bavg)
plt.show()
'''