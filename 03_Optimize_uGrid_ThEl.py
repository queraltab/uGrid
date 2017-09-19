# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 00:31:24 2016
@author: queralt

Main optimization function

"""


from __future__ import division
import sys
sys.path.append('Functions')
from technical_model_ThEl import operation
from ugrid_tools_ThEl import results_analysis,generate_results_table
from economic_tools_ThEl import EconomicInputs,LCOE
import numpy as np


############################## CONSTANT PARAMETERS #############################
timestep=1
NOCT=45                                 # [C] Nominal Open Circuit Temperature of the PV panel
Pmax_Tco=-0.4                           # [%/C]
lat=-29
lon=28
smart=1                                 # use a charging strategy that attempts to minimize the usage of the genset
Pmax_load=30                            # maximum power output of the load curve [kW]

heating_TES = True
heating_ORC = False

############################## OBJECTIVE FUNCTION #############################
def objective_function(params):
    print params
    # Imposing a minimum of zero for each parameter in case the optimizer is not bounded
    [BattkWh,Pnom_PV,CSP_A,Pnom_ORC,TESkWh] = [max(0,param) for param in params]    
    # Penalty for negative values    
    penalty = np.maximum(0,-params[0])+np.maximum(0,-params[1])+np.maximum(0,-params[2])+10*np.maximum(0,-params[3])+np.maximum(0,-params[4])
    # Run the technical model
    (Batt_kWh_tot,Propane,DNI,hour,P_genset,P_genset_ToBat,P_PV_ToBat,P_bat_inv,SOC_bat,SOC_TES,dayHour,P_load,E_loadleft,State,beam,T_ambient,P_PV,P_CSP,P_ORC_ToBat,ORC_Eta,P_ORC,ChargeTES,TES_loss,P_bat,P_dumped,Batt_frac,Genset_fuel,Fuel_kW,P_load_th,P_heating,P_heating_backup,P_dumped_th,eta_genset,m_burner,m_genset) = operation(Pmax_Tco, NOCT, smart, lat, lon, CSP_A, Pnom_ORC, Pnom_PV, TESkWh, BattkWh, Pmax_load, heating_TES, heating_ORC)
    up = generate_results_table(P_load,P_genset,P_genset_ToBat,P_PV,P_PV_ToBat,P_ORC,P_ORC_ToBat,
                           P_bat_inv,P_bat,P_dumped,P_CSP,SOC_bat,SOC_TES,dayHour,E_loadleft,
                           State,beam,T_ambient,ORC_Eta,ChargeTES,TES_loss,Batt_frac,Genset_fuel,
                           Fuel_kW,P_load_th,P_heating,P_heating_backup,P_dumped_th,m_burner,m_genset)

    # Economic study
    results = results_analysis(up,verbose=False)
    econ = EconomicInputs(results['m_propane'],results['m_propane_genset'],results['Pmax_load'],results['E_load'],results['E_load_th'],BattkWh,TESkWh,Pnom_PV,results['Pmax_load'],Pnom_ORC,CSP_A,results['E_bat'])
    cash,LCOE_thel,LCOE_el = LCOE(econ['lifetime'], econ['Batt_life_yrs'], econ['term'], econ['E_load'], econ['E_load_th'], econ['interest_rate'], econ['LEC'], econ['C1_PV'], econ['C1_CHP'], econ['C1_CSP'], econ['C1_LPG'], econ['Cost_bank'], econ['Cost_Propane_yr'],econ['Cost_Propane_genset_yr'], econ['Cost_Propane_burner_yr'])
    print "LCOE_thel: " + str(LCOE_thel)
    print "Penalty: " + str(penalty)
    return LCOE_thel + penalty
    

################################# OPTIMIZATION ################################
# Testing the function with the basic parameters:    
    
# Guess values for optimization parameters for Nelder-Mead:
BattkWh_guess=190     # kWh
Pnom_PV_guess=70      # kW
CSP_A_guess=230       # m²
Pnom_ORC_guess=5      # kW
TESkWh_guess=700      # kWh

tariff_test = objective_function([BattkWh_guess,Pnom_PV_guess,CSP_A_guess,Pnom_ORC_guess,TESkWh_guess])


### CHOOSE THE ALGORITHM ####
use_PSO = True

if not use_PSO:
    # optimization using scipy and nelder-mead:
    import scipy.optimize as optimize
    
    # Constraints:
    cons = ({'type': 'ineq', 'fun': lambda x: x[0]},
            {'type': 'ineq', 'fun': lambda x: x[1]})
    bnds = ((0, 1), (0, 10))
        
    guesses = [BattkWh_guess,Pnom_PV_guess,CSP_A_guess,Pnom_ORC_guess,TESkWh_guess]
    
    result = optimize.minimize(objective_function, guesses, method='Nelder-Mead', tol=1e-3).values()


else:
    # to install pyswarm: pip install pyswarm (as root)
    from pyswarm import pso

    lb = [0, 0, 0, 0, 0]            # lower bound
    ub = [300, 120, 600, 70, 2000]  # upper bound
        
    xopt, fopt = pso(objective_function, lb, ub, swarmsize=20)







