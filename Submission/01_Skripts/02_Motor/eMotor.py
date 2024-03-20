"""
Class for analytical presizing of permanent magnet (PM) excited synchronour inner rotor (IR) electrical machines

sources:
[1]: Bolte, E: Elektrische Maschinen (2018), pp. 643 ff.
[2]: Mezani, S. et. al: Sizing Electrical Machines Using OpenOffice Calc, 9th EUCASS (2022)
[3]: Zhu, Z. K. et al: Instantaneous Magnetic Field Distribution in Brushless Permanent Magnet dc Motors, Part I: Open-circuit Field  (1993)
[4]: Skaar, S. E. et al: Distribution, coil-span and winding factors for PM machines with concentrated windings (2006)
[5]: Binder, A. et al: Fixation of buried and surface mounted magnets in high-speed permanent magnet synchronous motors 

(c) Thilo WITZEL, TU Munich, 2023

"""
import math
import magfield_2D as m2d

class Motor:
    def __init__(self, nMax_rpm:float, Preq_W:float, UDC_V:float, dShaft:float, p:int) -> None:
        self.Treq_Nm = Preq_W * 60 / (math.pi * nMax_rpm * 2)       # required torque in Nm
        self.nMax_rpm = nMax_rpm                                    # required rpm at load point
        self.Preq_W = Preq_W                                        # required power at load point
        self.r1_m = dShaft / 2                                      # shaft radius in m
        self.p = p                                                  # pole pair number
        self.Udc_V = UDC_V                                          # DC voltage
        
        # definitions materal / thermal limits
        self.delta_m = 0.0005                                        # fluid (non magnetic) gap thickness in m
        self.Arms_Am = 50000                                        # rms current loading (thermal considerations) in A/m
        self.Jn_Am2 = 10e6                                           # conductor current (thermal considerations) density in A/m²
        self.rho_copper = 8600                                      # copper density in kg/m³
        # permanent magnet
        self.BrMag_T = 1.41                                         # remanent flux density of magnet in T
        self.mu_Rmag = 1.07                                         # relative permeability 
        self.alpha = 0.95                                           # relative pole pitch of the magnets
        self.rho_mag = 7500                                         # magnet density in kg/m³
        # magnetic steel
        self.Bmax_T = 2.20                                          # maximum allowable flux density of the iron cores in T
        self.k_fe = 0.95                                            # iron stacking factor
        self.Rp02_MPa = 90                                          # sleeve material yield strength in MPa (GF reinforced PEEK currently)
        self.Rm_Mpa = 90                                           # sleeve material tensile strength in MPa (GF reinforced PEEK currently)
        self.rho_steel = 7800                                       # steel density in kg/m³                                
        ## numerical parameters                             
        self.REL_TOL = 1e-14                                        # relative tolerance for iteration stop criterion
                                                                      
        pass
    
    def iterate_magnet_thickness(self, b_req:float, By_T:float, delta_g_m:float) -> float:
        ## relation to calculate magnet thickness iterative solution for given b-ratio
        # @inputs: 
        #   b_req: required ratio of Bm / Br
        #   By_T: desired flux density in the yoke in T
        #   delta_g_m: gemoetric air gap length in m
        #   dMag_init_m: magnet thickness of previous iteration in m
        # @output:
        #   updated magnet thickness in m

        # initial guess:
        if b_req > 0.9:
            dM_old = self.mu_Rmag * delta_g_m / (1 / b_req - 1) 
        else:
            dM_old = 9 * self.mu_Rmag * delta_g_m

        delta_dM = 1
        iters = 1
        while abs(delta_dM) > self.REL_TOL:
            dM_new = self.magnet_thickness_relation(By_T, b_req * self.BrMag_T, delta_g_m, dM_old)
            delta_dM = (dM_new - dM_old)/ dM_new
            dM_old = dM_new
            iters = iters + 1
        print(f' magnet thickness iterations: {iters}')
        return (dM_new)
    
    def magnet_thickness_relation(self, By_T:float, Bm_T:float, delta_g_m:float, dMag_init_m:float) -> float:
        ## relation to calculate magnet thickness (for iterative solution)
        # @inputs: 
        #   Bm_T: desired flux density at the PM/yoke iron interface in T
        #   By_T: desired flux density in the yoke in T
        #   delta_g_m: gemoetric air gap length in m
        #   dMag_init_m: magnet thickness of previous iteration in m
        # @output:
        #   updated magnet thickness in m
        # !! current (wrong) assumption: carter factor neglected for effective air gap !!
        b = Bm_T / self.BrMag_T    
        pfk2 = math.pi * self.alpha * self.BrMag_T * b / (2 * self.k_fe * By_T)
        k1 = pfk2 / (self.p - pfk2)

        dMag = b * self.r1_m * (1 + k1) * (math.log(1 + dMag_init_m / (self.r1_m * (1 + k1))) + self.mu_Rmag * math.log(1 + delta_g_m / (self.r1_m * (1 + k1) + dMag_init_m)))
        return dMag
    
    def determine_rotor_yoke(self, By_T:float, Bm_T:float) -> float:
        ## function to calculate rotor yoke height of a surface PM machine
        # @inputs: 
        #   Bm_T: desired flux density at the PM/yoke iron interface in T
        #   By_T: desired flux density in the yoke in T
        # @output:
        #   rotor yoke height in m
        b = Bm_T / self.BrMag_T    
        pfk2 = math.pi * self.alpha * self.BrMag_T * b / (2 * self.k_fe * By_T)
    
        return (self.r1_m * pfk2 / (self.p - pfk2))
    
    def calculate_rotor_geometry(self, delta_g_m:float) -> None:
        self.h_yr = self.determine_rotor_yoke(self.Bmax_T, 0.5 * self.Bmax_T)     # rotor yoke height
        self.r2_m = self.r1_m + self.h_yr
        self.h_mag = self.iterate_magnet_thickness(0.8, self.Bmax_T, delta_g_m) # permanent magnet height
        self.r3_m = self.r2_m + self.h_mag
        self.r4_m = self.r3_m + delta_g_m
        self.length = m2d.length_for_torque(self.Treq_Nm, self.Arms_Am, self.p, self.r2_m, self.r3_m, delta_g_m, self.BrMag_T, self.k_fe)  
        self.Brms_T = m2d.Brms_T(self.p, self.BrMag_T, self.r2_m, self.r3_m, delta_g_m)
        pass

    def calculate_concentrated_winding_factor(self, slots:int, poles, phases, UPh_V:float):
        ## function to calculate fundamental winding pactor of a 3 phase concentrated winding electrical machine
        # @inputs: 
        #   number of slots, poles and phases
        # @output:
        #   fundamental winding factor
        # calculstions acc to [4]
        
        q = slots / phases / poles                                                          # slots per pole per phase
        z = slots / math.gcd(slots, poles * phases)
        sigma = 60 * math.pi / 180                                                          # phase spread angle
        gamma_s = math.pi / q / phases                                                      # slot pitch angle
        eps = math.pi - gamma_s                                                             # coil span angle

        # some variables to check for design rules acc to. [4]
        F = math.gcd(round(slots), round(poles/2))                                                        # section F
        Pp = poles / (2 * F)                                                                # number of pole pairs in one section
        Na = slots / F                                                                      # number of slots per section

        k_mn = lambda n : math.sin(0.5 * sigma * n) / z / math.sin(0.5 * sigma * n / z)     # distribution factor
        k_en = lambda n : math.cos(0.5 * eps * n)                                           # coil span factor

        ## implement design rules for concentrated windings acc. to [4]:
        if q < 0.25:
            print('Too few slots per pole per phase for concentrated winding!')
            return [0, 0]
        if q > 0.5:
            print('Too many slots per pole per phase for concentrated winding!')
            return [0, 0]
        if (Pp % phases == 0):
            print('Unbalanced Winding!')
            return [-1, 0]
        if (poles == slots):
            print('Poles = number of slots -> unfeasible design')
            return [0, 0]  
        if (slots % phases != 0):
            print('Number of slots incompatible with number of phases!')
            return [-1, 0]  

        kwf = k_mn(1) * k_en(1)                                                             # fundamental winding factor

        # calculate the number of turns / phase N with back EMF
        # assumption 1: EMF must not exceed 0.9 * phase voltage [2]
        # assumption 2: fundamental winding factor is about 0.9 [2]
        N = UPh_V / (2 * kwf * self.r4_m * self.length * self.Brms_T * self.nMax_rpm * 2 * math.pi / 60)
        # number of coils per phase assuming a double layer winding
        Nc_Ph = slots / phases
        # number of turns per coil
        Nc = round(N / Nc_Ph)

        return [kwf, Nc]

    def calculate_drowned_rotor_geometry(self, delta_g_m:float, p_MPa:float) -> None:
        self.h_yr = self.determine_rotor_yoke(self.Bmax_T, 0.5 * self.Bmax_T)     # rotor yoke height
        self.r2_m = self.r1_m + self.h_yr

        # iterate over sleeve thickness to determine magnet height and air gap length
        dSleeve = 1
        gap_length = delta_g_m
        tSleeve_old = 0
        tJacket_m = 0
        iters = 0
        print('## iterating...')
        while abs(dSleeve) > self.REL_TOL:
            tSleeve_m = self.calculate_sleeve_thickness(p_MPa, self.r2_m + delta_g_m + tSleeve_old / 2)
            gap_length = delta_g_m + tSleeve_m + tJacket_m
            self.h_mag = self.iterate_magnet_thickness(0.8, self.Bmax_T, delta_g_m) # permanent magnet height
            self.r3_m = self.r2_m + self.h_mag + tJacket_m
            self.r4_m = self.r3_m + gap_length
            # check for magnet glue normal stress
            if self.calculate_magnet_centrifugal_force()[1] < 20e6:
                print("Magnet glue bond stress exceeds 20 [PMa] -> Polymer sleeve for mechanical integrity necessary!")
                tJacket_m = self.calculate_rotor_sleeve()
            else:
                tJacket_m = 0
            
            self.length = m2d.length_for_torque(self.Treq_Nm, self.Arms_Am, self.p, self.r2_m, self.r3_m, delta_g_m, self.BrMag_T, self.k_fe) 
            self.Brms_T = m2d.Brms_T(self.p, self.BrMag_T, self.r2_m, self.r3_m, delta_g_m)
            # relative sleeve thckness change per iteration
            dSleeve = (tSleeve_m - tSleeve_old) / tSleeve_m
            tSleeve_old = tSleeve_m
            iters +=1
            
        print(f' sleeve thickness calculations: {iters} \n')
        print('## Rotor geometry calculation finished ##')
        self.tSleeve_m = tSleeve_m
        pass
    
    def export_rotor_geometry(self) -> None:
        print(f'Rotor dimensions: outer D = {((self.r3_m) * 2000):0.1f} [mm], inner D = {((self.r1_m) * 2000):0.1f} [mm]')
        print(f'Magnet thickness = {(self.h_mag * 1000):0.1f} [mm]')
        print(f'Iron length = {(self.length * 1000):0.1f} [mm]')
        print(f'RMS air gap flux density = {self.Brms_T:0.2f} [T]')
        pass

    def calculate_distributed_winding(self, UPh_V:float, q:int):
        ## calculate 3 phase distributed winding for the machine according to [2]
        # @inputs:
        #   UPh_V: Peak phase voltage in V
        #   q: number of slots per pole per phase, should be even acc. to [2] to reduce 5th and 7th harmonic
        # @output: [Nc, Ns]
        #   Nc: number of turns per coil 
        #   Ns: number of slots

        # calculate the number of turns / phase N with back EMF
        # assumption 1: EMF must not exceed 0.9 * phase voltage [2]
        # assumption 2: fundamental winding factor is about 0.9 [2]
        N = UPh_V / (2 * 0.9 * self.r4_m * self.length * self.Brms_T * self.nMax_rpm * 2 * math.pi / 60)
        # calculate number of slots
        Ns = 6 * self.p * q
        # number of coils per phase assuming a double layer winding
        Nc_Ph = Ns / 3
        # number of turns per coil
        print(N / Nc_Ph)
        Nc = round(N / Nc_Ph)

        return [Nc, Ns] 
        
    def calculate_stator_geometry(self) -> None:
        ## function to estimate stator geometry based on assumptions acc. to [2]

        # calculate phase voltage in star connected 3 phase system
        Uph_V = self.Udc_V / math.sqrt(3) / math.sqrt(2)

        # check if the winding should be distributed or concentrated with pole and slot number based on slots per pole per phase
        s_min = round(0.25 * self.p * 2 * 3)
        s_max = round(0.5 * self.p * 2 * 3) + 1
        s_range = range(s_min, s_max)
        kw_holder = 0
        for slot in s_range:
            print('Slots: ', slot)
            kw_i = self.calculate_concentrated_winding_factor(slot, 2*self.p, 3, Uph_V)[0]
            if kw_i < 0.1 or slot < 4:
                continue
            elif not(kw_i < kw_holder):
                self.slots = slot
                kw_holder = kw_i
        print('Fundamental winding factor: ', kw_holder)
        print(self.slots)

        if kw_holder < 0.8 or self.slots < 5:
            print('Using a distributed winding')
            self.wTpe = 'distributed'
            # concentrated winding unfesible / undesirable
            winding = self.calculate_distributed_winding(Uph_V, 2)
            self.slots = winding[1]
            self.Nc = winding[0]
        else:
            print('Using a concentrated winding')
            self.wTpe = 'concentrated'
            winding = self.calculate_concentrated_winding_factor(self.slots, 2*self.p, 3, Uph_V)
            self.Nc = winding[1]

        # calculate required rms current in one conductor in A
        Irms = self.Arms_Am * 2 * math.pi * self.r4_m / (self.Nc * self.slots)
        # calculate required conductor area in m2
        Aw_m2 = Irms / self.Jn_Am2
        print(f'Phase current {Irms:0.1f} [Arms]')
        # calculate required slot area assuming a 40% copper fill factor and a double layer winding
        self.Aslot_m2 = Aw_m2 * 2 * self.Nc / 0.4
        print(f'Slot area = {self.Aslot_m2 * 1e6:0.1f} [mm2]')
        ## next section: calculate lamination dimensions (tooth, slot, yoke) acc. to [2] with 10% margin
        Ltooth_m = 2 * math.sqrt(2) * math.pi * self.r4_m * self.Brms_T / (self.slots * self.Bmax_T) * 1.025
        print(f'Tooth = {Ltooth_m * 1000:0.2f} [mm]')
        # angle per slot
        tau_u = 2 * math.pi * self.r4_m / self.slots

        # calculate slot depth acc. to [2] by finding the roots of the following function (OBSOLETE: TW 9.12.23)
        '''def f(x, n, tau, l, S):
            return x ** 2 + n/math.pi * (tau - l) * x - S * n / math.pi
        hSlot_m = max(fsolve(f, [0,1], args=(self.slots, tau_u, self.length, Aslot_m2))) * 1.1
        print('Slot Solve  [m]:', fsolve(f, [0,1], args=(self.slots, tau_u, self.length, Aslot_m2)))'''
        # calculte slot depth alternatively
        hSlot_m = max(self.calculate_tooth_height(self.slots, self.r4_m, Ltooth_m, self.Aslot_m2))
        print(f'Slot Solve = {1000*hSlot_m:0.2f} [mm]')
        # calculate yoke height acc. to [2]
        hYoke_m = math.pi * self.r4_m * self.Brms_T / (math.sqrt(2) * self.p * self.Bmax_T) * 1.025
        print(f'Yoke {1000*hYoke_m:0.2f} [mm]')
        # assemble to estimate outer radius of stator
        self.r5_m = self.r4_m + hSlot_m + hYoke_m
                
        pass
        
    def export_stator_geometry(self) -> None:
        print(f'Stator dimensions: outer D = {((self.r5_m) * 2000):0.1f} [mm], inner D = {((self.r4_m) * 2000):0.1f} [mm]')
        print(f'Iron length = {(self.length * 1000):0.1f} [mm]')
        print(f'Winding data: 3 Phases, {self.slots} slots, {self.Nc} turns per coil! Type: {self.wTpe} | Fundamental frequency= {round(self.nMax_rpm / 60 * self.p)} [Hz] !')
        pass

    def calculate_tooth_height(self, n_slots, ri_m, w_tooth_m, Areq_m2):
        phi = 2*math.pi / n_slots
        h1 = (w_tooth_m / phi - ri_m) + math.sqrt((ri_m - w_tooth_m)**2 + 2 * Areq_m2 / phi)
        h2 = (w_tooth_m / phi - ri_m) - math.sqrt((ri_m - w_tooth_m)**2 + 2 * Areq_m2 / phi)
        return [h1, h2]

    def calculate_sleeve_thickness(self, p_MPa:float, rm_m:float):
        # calculate minimum sleeve thicknesa acc. to barlows formula
        sLimit_mm = 0.5         # manufacturing limit for sleeve in mm
        # calculation for yield (proof pressure = MoS x MEOP)
        MoSY = 1.5
        sMinY_m = MoSY * p_MPa * rm_m / (self.Rp02_MPa)
    
        # calculation for burst (burst pressure = MoS x MEOP)
        MoSB = 2.5
        sMinB_m = MoSB * p_MPa * rm_m / (self.Rm_Mpa)
    
        sMin_m = max([sMinB_m, sMinY_m])
        if sMin_m < sLimit_mm/1000:
            print(f'Sleeve thickness: {sMin_m*1000:0.3f} [mm] is below the manufacturing limit of {sLimit_mm} [mm] ! ')
            return sLimit_mm/1000
        else:
            return sMin_m
    
    def calculate_motor(self) -> None:
        print('### Electric PMSM v1, regular rotor, Thilo Witzel, TUM ###')
        print("## Inputs for the simulation: ##")
        print(f'Load point: Power = {self.Preq_W} [W] | speed = {self.nMax_rpm} [rpm] | torque = {self.Treq_Nm:0.3f} [N*m]')
        print(f'Power supply voltage = {self.Udc_V} [VDC]')
        print(f'Shaft diameter = {self.r1_m * 2 * 1000} [mm]')
        print(f'Fundamental frequency= {self.nMax_rpm / 60 * self.p} [Hz] ! \n')
        print('# Calculating Rotor...')
        self.calculate_rotor_geometry(self.delta_m)
        print('# Finished! \n')
        print('# Calculating Stator...')
        self.calculate_stator_geometry()
        print('# Finished! \n')
        print('## Displaying results...')
        self.export_rotor_geometry()
        self.export_stator_geometry()
        print('### Programme finished sucessfully! ### \n \n')
        pass

    def calculate_drowned_motor(self, p_bar:float) -> None:
        print('### Electric PMSM, drowned rotor, v1 Thilo Witzel, TUM ###')
        print("## Inputs for the simulation: ##")
        print(f'Load point: Power = {self.Preq_W} [W] | speed = {self.nMax_rpm} [rpm] | torque = {self.Treq_Nm:0.3f} [N*m]')
        print(f'Power supply voltage = {self.Udc_V} [VDC]')
        print(f'Shaft diameter = {self.r1_m * 2 * 1000} [mm]')
        print(f'Operational pressure for sleeve dimensioning = {p_bar} [bar] \n')
        print('# Calculating Rotor...')
        self.calculate_drowned_rotor_geometry(self.delta_m, p_bar / 10)
        print('# Finished! \n')
        print('# Calculating Stator...')
        self.calculate_stator_geometry()
        print('# Finished! \n')
        print('## Displaying results...')
        self.export_drowned_rotor_geometry()
        self.export_stator_geometry()
        print('### Programme finished sucessfully! ### \n \n')
        pass

    def export_drowned_rotor_geometry(self) -> None:
        print(f'Rotor dimensions: outer D = {((self.r3_m) * 2000):0.1f} [mm], inner D = {((self.r1_m) * 2000):0.1f} [mm]')
        print(f'Magnet thickness = {(self.h_mag * 1000):0.1f} [mm]')
        print(f'RMS air gap flux density = {self.Brms_T:0.2f} [T]')
        print(f'Iron length = {(self.length * 1000):0.1f} [mm]')
        print(f'Sleeve thickness = {(self.tSleeve_m * 1000):0.1f} [mm]')
        print(f'\t resulting in total gap = {((self.tSleeve_m + self.delta_m) * 1000):0.1f} [mm]')
        pass
    
    def calculate_mass(self) -> float:
        self.m_kg = 0
        
        # rotor mass
        self.m_kg += math.pi * (self.r2_m**2 - self.r1_m**2) * self.length * self.rho_steel
        # magnet mass
        self.m_kg += math.pi * (self.r3_m**2 - self.r2_m**2) * self.length * self.rho_mag * self.alpha
        # calculate stator winding volume
        V_wind = self.Aslot_m2 * self.slots * self.length
        # winding mass assuming 40% copper fill factor
        self.m_kg += V_wind * 0.4 * self.rho_copper
        # calculate stator iron mass
        self.m_kg += (math.pi * (self.r5_m**2 - self.r4_m**2) * self.length - V_wind) * self.rho_steel
        
        return self.m_kg

    def plot_theor_gap_density(self) -> None:
        #m2d.plot_theor_gap_density(self.r1_m, self.r2_m, self.r3_m, self.r4_m, self.r5_m, self.length, self.p, self.BrMag_T, 1000)
        alt_torque = m2d.torque2DFlux(self.Arms_Am, self.p, self.r2_m, self.r3_m, self.delta_m, self.BrMag_T, self.length, self.k_fe)
        print('Alternative Torque: ', alt_torque, ' [Nm]')
        m2d.plot_radial_field(self.r4_m, self.delta_m, self.h_mag, self.BrMag_T, self.alpha, self.p, self.mu_Rmag)

    def __del__(self):
        print('Destructor of class Motor has been called. Object deleted!')
        pass

    def calculate_magnet_centrifugal_force(self):
        mMag_kg = math.pi * (self.r3_m ** 2 - self.r2_m ** 2) * self.rho_mag / (2 * self.p) * self.alpha    # mass of single magnet in kg/m
        Fc_N =  mMag_kg * (self.r2_m + self.h_mag / 2) * (2 * math.pi * self.nMax_rpm / 60) ** 2                          # centripedal force in N/m calculated with single magnet mass at magnet middle diameter
        sigmaC_Pa = Fc_N / (self.r2_m * self.alpha * 2 * math.pi / (2 * math.pi))

        print(f'Rotor surface speed = {(self.nMax_rpm / 60 * 2 * math.pi * self.r3_m):0.2f} [m/s]')

        return [Fc_N, sigmaC_Pa]
    
    def calculate_rotor_sleeve(self):
        # function to calculate retaining sleeve thickness for high speed rotors acc. to [5]
        w_max = 1.2 * 2 * math.pi * self.nMax_rpm / 60          # maximum overspeed with 20% margin
        pCont_Pa = 10e6                                            # safe remaining contact pressure at max speed in MPa

        sigmaSleeve_Pa = 1100 * 1e6                             # maximum allowable sleeve stress, CF 55% acc. to [5]
        rhoSleeve_kg_m3 = 1850                                  # sleeve material density in kg /m3
        YMSleeve_Pa = 132000 *1e6                               # sleeve young's modulus in MPa, CF 55% acc. to [5]

        hS_m_cache = 0
        hS_m = 0.25 * self.h_mag                                # initial sleeve thickness

        while (abs((hS_m - hS_m_cache) / hS_m) > self.REL_TOL):
            hS_m_cache = hS_m
            rmag = (self.r2_m + self.r3_m) / 2  # middle magnet radius
            rsl = self.r3_m + hS_m / 2          # middle sleeve radius
            sig_tw = w_max ** 2 * rsl ** 2 * rhoSleeve_kg_m3 # tangential sleeve stress due to sleeve rotation, thin shell, [5] eqn. 5a
            sig_tpre = sigmaSleeve_Pa - sig_tw # prestress chosen so that allowable stress is not exceeded
            pwm = rmag * self.rho_mag * w_max ** 2 * self.h_mag
            pws = rsl * rhoSleeve_kg_m3 * w_max ** 2 * hS_m
            p_pre = pCont_Pa + pwm + pws # prestress contact pressure chosen so that defined safe pressure is reached at overspeed
            hS_m = self.r3_m * (p_pre)/(sig_tpre)# - p_pre / 2)
            print(f'Sleeve iter: h = {(hS_m*1000):0.3f} [mm]')

        # check results
        rmag = (self.r2_m + self.r3_m) / 2  # middle magnet radius
        rsl = self.r3_m + hS_m / 2          # middle sleeve radius
        sig_tw = w_max ** 2 * rsl ** 2 * rhoSleeve_kg_m3 # tangential sleeve stress due to sleeve rotation, thin shell, [5] eqn. 5a
        sig_max = sig_tpre + sig_tw
        dD = sig_tpre / YMSleeve_Pa * 2 * self.r3_m
        print(f'Sleeve stress [MPa] = {sig_max*1e-6} | diameter stretch [mm]: {dD * 1000} / [%]:{dD / 2 / self.r3_m * 100}')
        
        return hS_m