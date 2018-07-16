# -*- coding: utf-8 -*-
"""
This file contains all the required functions to run the ugrid economical model.

Created on Tue Jun 21 10:13:16 2016

@author: queralt
"""

from __future__ import division
import numpy as np
import pandas as pd


def EconomicInputs(m_propane,m_propane_genset,peakload,E_load,E_load_th,BattkWh,TESkWh,PVkW,gensetkW,ORCkW,CSP_A,E_bat):
    inputs = {}
    inputs['lifetime'] = 15
    inputs['term'] = 15
    inputs['interest_rate'] = 0.05
    inputs['LEC'] = 0.1
    
    Batt_lifecycle = 1750       # The capacity of the battery is approximately 1750 times its C10 rating OPG2-2000.pdf
    inputs['Batt_life_yrs'] = np.floor((BattkWh*Batt_lifecycle)/(E_bat+0.01))     # Years of battery life before replacement is necessary, rounded down to an integer
    
    
    # UPDATE THESE PARAMETERS FOR A COMMUNITY OF INTEREST
    num_nodes = 90              # Number of connections / meters on the network
    Dist_km = 5                 # km from ViPor
    num_step_up_trans = 1       # from ViPor
    num_pole_trans = 3          # Number of transformers on network from the ViPor topology
    num_poles = Dist_km/0.05    # Number of poles, based on 50 meter pole-to-pole distance
    
    
    # COMPONENT COST FUNCTIONS
    Cost_dist_wire = 0.5*1000           # 0.5 USD/m
    Cost_batt = 130                     # USD/kWh
    Cost_panels = PVkW*1000             # PV price via Alibaba 2014
    Cost_control = 5000                 # STG Build or PLC
    Cost_Pole = 40                      # Transmission Pole Prices 2016 from treatedpoles.co.za
    Cost_Pole_Trans = 150               # USD, Alibaba 20kVA single phase 11kV/.22kV
    Cost_Step_up_Trans = 1000           # USD, Alibaba 63kVA single phase 11kV/.22kV
    
    Cost_Smartmeter = 50*num_nodes      # Iometer
    Cost_MPesa = 70*peakload            # Estimate for merchant services with vodacom
    Cost_inv = peakload*800             # USD/kW peak
    
    
    inputs['Cost_bank'] = BattkWh*Cost_batt
    inputs['Cost_Propane_yr'] = m_propane*1.24       # 2.037USD/gallon  0.00378541178m³/gallon  2.01kg/m³    1.24 USD/kg, RSA prices 2014 (1USD/kg)
    inputs['Cost_Propane_genset_yr'] = m_propane_genset*1.24
    inputs['Cost_Propane_burner_yr'] = (m_propane-m_propane_genset)*1.24
    
    Cost_dist = Cost_dist_wire*Dist_km+Cost_Step_up_Trans*num_step_up_trans+Cost_Pole_Trans*num_pole_trans+Cost_Pole*num_poles
    Cost_BOS = Cost_inv + Cost_control + Cost_dist + Cost_Smartmeter + Cost_MPesa # USD, Balance Of System
    
    inputs['C1_PV'] = Cost_panels + Cost_BOS
    inputs['C1_LPG'] = -10354.1143 + 6192.606*np.log(peakload)    # USD, propane genset costs based on general lineup
    if ORCkW>0:
        inputs['C1_CHP'] = ORCkW*2000 + CSP_A*150 + 25*TESkWh     #  ORC: $2/W 
        inputs['C1_CSP'] = 0
    else:
        inputs['C1_CHP'] = 0    
        inputs['C1_CSP'] = CSP_A*150 + 25*TESkWh        # Collector: $150/m², TES: $25/kWh
    
    inputs['dailyload'] = E_load/365 # kWh
    inputs['dailyload_tot'] = (E_load_th+E_load)/365
    inputs['E_load'] = E_load
    inputs['E_load_th'] = E_load_th
    
    return inputs


def EconomicInputs2(m_propane,m_propane_genset,peakload,E_load,E_load_th,BattkWh,TESkWh,PVkW,gensetkW,ORCkW,CSP_A,E_bat,A_FPC):
    '''
    Economic inputs for a system with Flat Plate Collectors instead of ORC
    '''    
    inputs = {}
    inputs['lifetime'] = 15
    inputs['term'] = 15
    inputs['interest_rate'] = 0.05
    inputs['LEC'] = 0.1
    
    Batt_lifecycle = 1750   # The capacity of the battery is approximately 1750 times its C10 rating OPG2-2000.pdf
    inputs['Batt_life_yrs'] = np.floor((BattkWh*Batt_lifecycle)/(E_bat+0.01))     # Years of battery life before replacement is necessary, rounded down to an integer
    
    
    # UPDATE THESE PARAMETERS FOR A COMMUNITY OF INTEREST
    num_nodes = 90              # Number of connections / meters on the network
    Dist_km = 5                 # km from ViPor
    num_step_up_trans = 1       # from ViPor
    num_pole_trans = 3          # Number of transformers on network from the ViPor topology
    num_poles = Dist_km/0.05    # Number of poles, based on 50 meter pole-to-pole distance
    
    
    # COMPONENT COST FUNCTIONS
    Cost_dist_wire = 0.5*1000           # 0.5 USD/m
    Cost_batt = 130                     # USD/kWh
    Cost_panels = PVkW*1000             # PV price via Alibaba 2014
    Cost_control = 5000                 # STG Build or PLC
    Cost_Pole = 40                      # Transmission Pole Prices 2016 from treatedpoles.co.za
    Cost_Pole_Trans = 150               # USD, Alibaba 20kVA single phase 11kV/.22kV
    Cost_Step_up_Trans = 1000           # USD, Alibaba 63kVA single phase 11kV/.22kV
    
    Cost_Smartmeter = 50*num_nodes      # Iometer
    Cost_MPesa = 70*peakload            # Estimate for merchant services with vodacom
    Cost_inv = peakload*800             # USD/kW peak
    
    
    inputs['Cost_bank'] = BattkWh*Cost_batt
    inputs['Cost_Propane_yr'] = m_propane*1.24       # USD/kg, RSA prices 2014 (1USD/kg)
    inputs['Cost_Propane_genset_yr'] = m_propane_genset*1.24
    inputs['Cost_Propane_burner_yr'] = (m_propane-m_propane_genset)*1.24
    
    Cost_dist = Cost_dist_wire*Dist_km+Cost_Step_up_Trans*num_step_up_trans+Cost_Pole_Trans*num_pole_trans+Cost_Pole*num_poles
    Cost_BOS = Cost_inv + Cost_control + Cost_dist + Cost_Smartmeter + Cost_MPesa # USD, Balance Of System
    
    inputs['C1_PV'] = Cost_panels + Cost_BOS
    inputs['C1_LPG'] = -10354.1143 + 6192.606*np.log(peakload)    # USD, propane genset costs based on general lineup
    # Case with flat plate solar thermal collectors
    inputs['C1_CHP'] = 0                            # No ORC
    inputs['C1_CSP'] = 25*TESkWh + A_FPC*100        # FPC: $100/m², TES: $25/kWh
    
    inputs['dailyload'] = E_load/365 # kWh
    inputs['dailyload_tot'] = (E_load_th+E_load)/365
    inputs['E_load'] = E_load
    inputs['E_load_th'] = E_load_th
    
    return inputs


  
def LCOE(lifetime, Batt_life_yrs, term, loadkWh, loadkWh_th, interest_rate, LEC, C1_PV, C1_CHP, C1_CSP, C1_LPG, Cost_bank, Cost_Propane_yr, Cost_Propane_genset_yr, Cost_Propane_burner_yr):
    '''
    Function that calculates the levelized cost of energy.
    :param lifetime:          Project lifetime, years
    :param term:              not used (assumed equal to lifetime)
    :param loadkWh:           kWh, yearly load, kWh
    :param interest_rate:     - 
    :param C1_PV:             PV Investment costs, USD
    :param C1_CSP:            CSP Investment costs, USD
    :param C1_LPG:            genset Investment costs, USD
    :param Cost_bank:         battery investmen costs, USD
    :param Cost_Propane_yr:   yearly fuel costs, USD
    :return out:              Detailed (yearly financial information)
    :return LCOE_ThEl:        Global tariff for electricity and heat, USD/kWh
    :return LCOE:             Tariff for electricity, USD/kWh
    '''
    
    years = range(1,lifetime+1)
    
    # Maintenance costs over each year of the project:
    # Factors for distributing maintenance costs in time as a function of capex:
    # (function developed in Orosz IMechE paper in 2012) 
    f_PV = 0.25
    a_PV = 0.25
    f = 1.25
    a = 0.25
    M = [f_PV*a_PV*C1_PV/lifetime + (f_PV*(1-a_PV)*C1_PV/lifetime**2)*(2*y-1) + f*a*C1_LPG/lifetime + (f*(1-a)*C1_LPG/lifetime**2)*(2*y-1) + f*a*C1_CHP/lifetime + (f*(1-a)*C1_CHP/lifetime**2)*(2*y-1) for y in range(1,lifetime+1)]
    M_ThEl = [f_PV*a_PV*C1_PV/lifetime + (f_PV*(1-a_PV)*C1_PV/lifetime**2)*(2*y-1) + f*a*C1_LPG/lifetime + (f*(1-a)*C1_LPG/lifetime**2)*(2*y-1) + f*a*C1_CSP/lifetime + (f*(1-a)*C1_CSP/lifetime**2)*(2*y-1) + f*a*C1_CHP/lifetime + (f*(1-a)*C1_CHP/lifetime**2)*(2*y-1) for y in range(1,lifetime+1)]
    
    # Operating costs:
    O_ThEl = [Cost_Propane_yr] * lifetime
    O = [Cost_Propane_genset_yr] * lifetime    
    
    # Investment costs:
    I = np.zeros(lifetime)
    I_ThEl = np.zeros(lifetime)
    I[0] = C1_PV + C1_LPG + C1_CHP
    I_ThEl[0] = C1_PV + C1_CSP + C1_CHP + C1_LPG
    for y in range(lifetime):
        if y % Batt_life_yrs==0:
            I[y] = I[y] + Cost_bank
            I_ThEl[y] = I_ThEl[y] + Cost_bank

    # Levelized cost of energy
    LCOE_ThEl = np.sum([(M_ThEl[y-1] + O_ThEl[y-1] + I_ThEl[y-1])/(1+interest_rate)**y for y in years])/np.sum([(loadkWh+loadkWh_th)/(1+interest_rate)**y for y in years])
    # Levelized cost of electricity
    LCOE = np.sum([(M[y-1] + O[y-1] + I[y-1])/(1+interest_rate)**y for y in years])/np.sum([loadkWh/(1+interest_rate)**y for y in years])    
    
    out = pd.DataFrame()
    out['Year'] = years
    out['M'] = M
    out['O'] = O
    out['I'] = I
    
    return out,LCOE_ThEl,LCOE
    
    
#%% The MIT License
'''
SPDX short identifier: MIT

Copyright  2017  Queralt Altés Buch, Matthew Orosz

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''