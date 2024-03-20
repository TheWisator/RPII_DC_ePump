"""
Script for analytical presizing of permanent magnet (PM) excited synchronour inner rotor (IR) electrical machines

source: Bolte, E: Elektrische Maschinen (2018), pp. 643 ff.

(c) Thilo WITZEL, TU Munich, 2023

"""

# imports
import math
from eMotor import Motor
import matplotlib.pyplot as plt
import numpy as np

## Define input parameters for dimensioning of the motor
r1_m = 0.005            # inner radius (shaft) of the rotor in m
nMax_rpm = 40000        # maximum rotational speed in rpm
Preq_W = 2000           # required shaft power @ max speed in W
Treq_Nm = Preq_W * 60 / (math.pi * nMax_rpm * 2)    # requred torque to achieve power at shaft speed in Nm
Udc_V = 55             # available DC voltage of the power supply
m = 3                   # number of phases
p_init = 2              # initial value for pole pairs

# build a motor class with the respective initial parameters
# material parameters have to be set in the __init__ function of the motor class in eMotor.py
motor1 = Motor(nMax_rpm, Preq_W, Udc_V, r1_m*2, p_init)
# calculate motor for 30 bar operational pressure with drowned rotor
motor1.calculate_drowned_motor(30)
print(f'Motor mass = {motor1.calculate_mass():0.2f} [kg]')
print(f'Motor power density = {(motor1.Preq_W / motor1.m_kg):0.2f} [W/kg]')
