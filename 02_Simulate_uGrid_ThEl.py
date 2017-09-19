# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 16:49:02 2016

@author: queralt

Test function operation from uGrid Technical Model
"""
from __future__ import division
import sys
from ugrid_tools_ThEl import generate_results_table
sys.path.append('Functions')
import time
start_time = time.clock()


# Storing results in an Excel requires some seconds. As default, results are not stored to make the simulation faster.
# Want to store results? Change this line to True:
store_results = False
plot_results = False

################################## PARAMETERS ###################################

### COMMUNITY ### - Need to be adapted for the studied community. Ha Nkau:
timestep=1
NOCT=45         # [C] Nominal Open Circuit Temperature of the PV panel
Pmax_Tco=-0.4   # [%/C]
lat=-29
lon=28
smart=1         # use a charging strategy that attempts to minimize the usage of the genset
Pmax_load=30    # maximum power output of the load curve [kW]


### SYSTEM ###

#Optimal sizing
BattkWh=262     # kWh
Pnom_PV=66      # kW
CSP_A=66        # m²
Pnom_ORC=0      # kW
TESkWh=86       # kWh
heating_TES = True
heating_ORC = False



############################ RUN THE TECHNICAL MODEL ############################
from technical_model_ThEl import operation
(Batt_kWh_tot,Propane, DNI,hour,P_genset,P_genset_ToBat,P_PV_ToBat,P_bat_inv,SOC_bat,SOC_TES,dayHour,P_load,E_loadleft,State,beam,T_ambient,P_PV,P_CSP,P_ORC_ToBat,ORC_Eta,P_ORC,ChargeTES,TES_loss,P_bat,P_dumped,Batt_frac,Genset_fuel,Fuel_kW,P_load_th,P_heating,P_heating_backup,P_dumped_th,eta_genset,m_burner,m_genset) = operation(Pmax_Tco, NOCT, smart, lat, lon, CSP_A, Pnom_ORC, Pnom_PV, TESkWh, BattkWh, Pmax_load, heating_TES, heating_ORC)

up = generate_results_table(P_load,P_genset,P_genset_ToBat,P_PV,P_PV_ToBat,P_ORC,P_ORC_ToBat,
                           P_bat_inv,P_bat,P_dumped,P_CSP,SOC_bat,SOC_TES,dayHour,E_loadleft,
                           State,beam,T_ambient,ORC_Eta,ChargeTES,TES_loss,Batt_frac,Genset_fuel,
                           Fuel_kW,P_load_th,P_heating,P_heating_backup,P_dumped_th,m_burner,m_genset)

### STORE RESULTS IN EXCEL ###
if store_results:
    up.to_excel('Results/Py_ThEl_results_1.xlsx')


### PLOT RESULTS ###
if plot_results: 
    from ugrid_tools_ThEl import dispatch_plot, dispatch_ThEl_plot
    dispatch_plot(up,rng=range(5400,5700))
    dispatch_ThEl_plot(up,rng=range(5400,5700))


################################ ECONOMIC STUDY #################################
from ugrid_tools_ThEl import results_analysis
from economic_tools_ThEl import EconomicInputs,LCOE
results = results_analysis(up)

# Economic study:
econ = EconomicInputs(results['m_propane'],results['m_propane_genset'],results['Pmax_load'],results['E_load'],results['E_load_th'],BattkWh,TESkWh,Pnom_PV,results['Pmax_load'],Pnom_ORC,CSP_A,results['E_bat'])

# Calculation of the levelized cost of electricity:
cash,LCOE_thel,LCOE_el = LCOE(econ['lifetime'], econ['Batt_life_yrs'], econ['term'], econ['E_load'], econ['E_load_th'], econ['interest_rate'], econ['LEC'], econ['C1_PV'], econ['C1_CHP'], econ['C1_CSP'], econ['C1_LPG'], econ['Cost_bank'], econ['Cost_Propane_yr'],econ['Cost_Propane_genset_yr'], econ['Cost_Propane_burner_yr'])

print 'LCOE_thel: ' + str(LCOE_thel) + ', LCOE: ' + str(LCOE_el)
print time.clock() - start_time, "seconds"





################################ COMMENTS #################################

#heating_TES = True
#heating_ORC = False

# CASE 1: EVERYTHING
#BattkWh=90      # kWh
#Pnom_PV=35      # kW
#CSP_A=230       # m²
#Pnom_ORC=10     # kW
#TESkWh=700      # kWh

## CASE 2: JUST PV and BATTERIES
#BattkWh=220     # kWh
#Pnom_PV=70      # kW
#CSP_A=0         # m²
#Pnom_ORC=0      # kW
#TESkWh=0        # kWh
#
## CASE 3: JUST ORC and TES
#BattkWh=0       # kWh
#Pnom_PV=0       # kW
#CSP_A=10        # m²
#Pnom_ORC=4      # kW
#TESkWh=2        # kWh
#
## CASE 4: JUST ORC-TES and BATTERIES
#BattkWh=50      # kWh
#Pnom_PV=0       # kW
#CSP_A=730       # m²
#Pnom_ORC=10     # kW
#TESkWh=2600     # kWh

## CASE 5: JUST ORC bottoming on GENSET
#BattkWh=0       # kWh
#Pnom_PV=0       # kW
#CSP_A=0         # m²
#Pnom_ORC=10     # kW
#TESkWh=0        # kWh
#
## CASE 6: JUST ORC bottoming on GENSET and BATTERIES
#BattkWh=80      # kWh
#Pnom_PV=0       # kW
#CSP_A=0         # m²
#Pnom_ORC=10     # kW
#TESkWh=0        # kWh
#
## CASE 7: JUST GENSET and BATTERIES
#BattkWh=70      # kWh
#Pnom_PV=0       # kW
#CSP_A=0         # m²
#Pnom_ORC=0      # kW
#TESkWh=0        # kWh
#
## CASE 8: JUST PV and ORC bottoming on GENSET
#BattkWh=30       # kWh
#Pnom_PV=40       # kW
#CSP_A=0          # m²
#Pnom_ORC=10      # kW
#TESkWh=0         # kWh