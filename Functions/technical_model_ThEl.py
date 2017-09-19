# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 11:04:42 2016

@author: queralt
"""
from __future__ import division
import numpy as np
from scipy.interpolate import interp1d
from ugrid_tools_ThEl import load_xls
import time

start_time = time.clock()



########################  Load DATA  ########################

P_loadT = load_xls('Data/DATA.xlsx','LoadKWT')
FullYearEnergy = load_xls('Data/DATA.xlsx','FullYearEnergy')
Heating_P_load = load_xls('Data/DATA.xlsx','Heating')

tmy = load_xls('Data/IncidentRadiation.xlsx','Sheet1')
timedeltas = tmy.index - tmy.index[0]
x=timedeltas.total_seconds()/3600
f_Tamb = interp1d(x,tmy['T'].values)
f_Ib=interp1d(x,tmy['Direct with tracking'].values)
f_I_PV=interp1d(x,tmy['poa_global'].values)


########################  FUNCTIONS  ########################

def P_bat_minmax(T_amb, Pmax_bat, SOC_bat_in, BattkWh, timestep):
    '''
    Calculation of maximum power that can be charged/discharged in a time step. Takes into account the 
    state of charge and the charging/discharge limit characteristic of the batteries.
    :param T_amb:           Ambient temperature, °C
    :param Pmax_bat:        Limit value for battery charging/discharging, kW
    :param SOC_bat_in:      State of Charge out from the previous time step, kWh
    :param BattkWh:         Installed batteries full capacity, kWh
    :param timestep:        Time step, hours
    :return P_bat_min:      Maximum discharge power in curent timestep, kW
    :return P_bat_max:      Maximum charge power in current time step, kW
    :return high_trip:      Limit above which batteries are considered full and cannot be charged beyond it (95 % of batteries full capacity + effect of temperature considered), kWh
    :return SOC_bat_out:    Batteries' state of charge after accounting for self discharge, kWh
    '''
    Self_discharge=timestep*3600*BattkWh*3.41004573*1e-9*np.exp(0.0693146968*T_amb) # Effect of temperature on self discharge of the battery bank from 5% per month at 25C and doubling every 10C
    SOC_frac=0.711009126 + 0.0138667895*T_amb - 0.0000933332591*T_amb**2            # Effect of temperature on the capacity of the battery bank from IEEE 485 Table 1
    SOC_bat_inout=SOC_bat_in-Self_discharge

    high_trip=0.95*BattkWh*SOC_frac    
    
    P_bat_max = np.minimum((high_trip - SOC_bat_inout)/timestep,Pmax_bat)   
    P_bat_min = np.maximum(-SOC_bat_inout/timestep,-Pmax_bat)                      
    
    return P_bat_min,P_bat_max,high_trip,SOC_bat_inout

def batt_bal(P_bat_max_discharge, P_bat_max_charge, P_PV_ToBat, P_bat_dis, P_genset_ToBat, P_ORC_ToBat, SOC_bat_inout, BattkWh, timestep):
    '''
    Calculation of the energy balance of the batteries bank, accounting for temperature and 
    depth of discharge effects on the bank.
    :param P_bat_max_discharge: Maximum discharge power in curent timestep, kW
    :param P_bat_max_charge:    Maximum charge power in current time step, kW
    :param P_PV_ToBat:          Load powered from PV, kW
    :param P_bat_dis:           Load powered from the inverter, kW
    :param P_genset_ToBat:      Load powered from the genset, kW
    :param P_ORC_ToBat:         Load powered from the ORC, kW
    :param SOC_bat_inout:       Batteries' state of charge after accounting for self discharge, kWh
    :param BattkWh:             Installed batteries full capacity, kWh
    :param timestep:            Time step, hours
    :return SOC_bat_out:        Batteries' state of charge in the current time step after accounting for self discharge and charging/discharging, kWh
    :return Charge:             Load charged to the batteries, kWh
    :return P_bat_charging: 
    :return P_dumped:           Dumped load, kWh
    :return high_trip:          Limit above which batteries are considered full and cannot be charged beyond it (95 % of batteries full capacity - effect of temperature considered), kWh
    '''    
    # Possible charge/discharge mechanisms accounted for:
    P_bat = P_PV_ToBat + P_bat_dis + P_genset_ToBat + P_ORC_ToBat
    
    charging = True if P_bat>0 else False
    discharging = False if P_bat>0 else True
    
    if BattkWh==0:
        # If there are no batteries, all energy gets dumped:
        P_bat = 0
        P_dumped = P_PV_ToBat + P_genset_ToBat + P_ORC_ToBat
        SOC_bat_out = 0
        P_bat_charging = 0
        if P_bat_dis<0:
            print 'WARNING: System has no Battery but Inverter is operational -> ERROR!!!'
    else:
        # Batteries cannot be charged beyond available capacity, extra goes to a dump
        if charging:
            if P_bat > P_bat_max_charge:
                P_dumped = P_bat - P_bat_max_charge
                P_bat = P_bat_max_charge
            else:
                P_dumped = 0
            SOC_bat_out = SOC_bat_inout + P_bat*timestep*0.9
            P_bat_charging = P_bat
        
        if discharging:
            P_dumped = 0
            if P_bat < P_bat_max_discharge:
                P_bat = P_bat_max_discharge
            SOC_bat_out = SOC_bat_inout + P_bat*timestep/0.9
            P_bat_charging = 0
    
    return P_bat_charging,P_bat,SOC_bat_out,P_dumped


def storage_levels(SOC_bat,low_trip,high_trip,SOC_TES,lowT_trip,mediumT_trip,highT_trip):
    '''
    Definition of storage levels in a certain time step by classifying them into different cases
    :param SOC_bat:         Batteries' State of Charge, kWh
    :param low_trip:        Limit below which batteries are considered too low (5 % of batteries full capacity), kWh
    :param high_trip:       Limit above which batteries are considered full and cannot be charged beyond it (95 % of batteries full capacity + effect of temperature considered), kWh
    :param SOC_TES:         TES State of Charge, kWh
    :param lowT_trip:       Limit below which TES level is considered too low to serve heating demand, kWh
    :param mediumT_trip:    Limit below which TES level is considered too low to run the ORC (approximative amount of energy withdrawn in one time step if the ORC is ON - set for minimum possible ORC efficiency of 10%), kWh
    :param highT_trip:      Limit above which TES is considered full and cannot be charget beyond it (95 % of TES full capacity), kWh
    :return batt_state:     Batteries storage level, classified in: 'empty' - too low level, 'part charge'- batteries are charged but not full, 'full' - batteries are full
    :return TES_state:      TES storage level, classified in: 'too low' - TES temperature is too low and cannot be used, 'low' - TES is charged but not full and cannot run the ORC, 'high' - TES is charged but not full and can run the ORC, 'too high' - TES is full (temperature reached its maximum)
    '''
    from ugrid_tools_ThEl import batt_level, TES_level
    batt_state = batt_level(SOC_bat,low_trip,high_trip)
    TES_state = TES_level(SOC_TES,lowT_trip,mediumT_trip,highT_trip)
    
    return batt_state, TES_state



def set_state_elec(P_bat_max_discharge,P_bat_max_charge,there_is_ORC,there_is_CSP,there_is_batt,heating_TES,heating_ORC,timestep,T_amb,P_PV,P_CSP,P_load,TES_state,batt_state,Pnom_ORC,Pmax_load,SOC_bat,low_trip,E_loadleft,P_load_th,SOC_TES,highT_trip,TESkWh,lowT_trip,T_w_return,there_is_heating_demand,TES_tooSmallforORC):
    ''' 
    Power dispatching for an electrical driven control, i.e. the following cases:
    - when heat demand is provided by TES,
    - when heat demand is served by the ORC and there is no CSP-TES, i.e. ORC running with genset exhaust gases,
    - when heat demand is served by the ORC and there is CSP-TES but there is no heating demand at this time step.
    :param P_bat_max_discharge:     Batteries' maximum discharge power in curent timestep, kW
    :param P_bat_max_charge:        Batteries' maximum charge power in current time step, kW
    :param there_is_ORC:            Boolean variable that is True: if there is ORC installed in the system, False: on the contrary
    :param there_is_CSP:            Boolean variable that is True: if there is CSP installed in the system, False: on the contrary
    :param there_is_batt:           Boolean variable that is True: if there is a battery bank installed in the system, False: on the contrary
    :param heating_TES:             Boolean variable that is True: if the heating system is the TES, False: on the contrary
    :param heating_ORC:             Boolean variable that is True: if the heating system is the ORC, False: on the contrary
    :param timestep:                Time step, hours
    :param T_amb:                   Ambient temperature, °C
    :param P_PV:                    PV power generation, kW
    :param P_CSP:                   CSP power generation, kW
    :param P_load:                  Electrical demand that must be supplied, kW
    :param TES_state:               TES storage level, classified in: 'too low' - TES temperature is too low and cannot be used, 'low' - TES is charged but not full and cannot run the ORC, 'high' - TES is charged but not full and can run the ORC, 'too high' - TES is full (temperature reached its maximum)
    :param batt_state:              Batteries storage level, classified in: 'empty' - too low level, 'part charge'- batteries are charged but not full, 'full' - batteries are full
    :param Pnom_ORC:                Installed ORC capacity, kW
    :param Pmax_load:               Maximum power output of the load curve, kW
    :param SOC_bat:                 Batteries' state of charge after accounting for self discharge, kWh
    :param low_trip:                Limit below which batteries are considered too low (5 % of batteries full capacity), kWh
    :param E_loadleft:              Amount of load defined by smart charging (during daylight: a value higher than batteries capacity, otherwise: load remaining from premidnight until dawn), kWh
    :param P_load_th:               Thermal demand that must be supplied, kW
    :param SOC_TES:                 TES State of Charge, kWh
    :param highT_trip:              Limit above which TES is considered full and cannot be charget beyond it (95 % of TES full capacity), kWh
    :param TESkWh:                  Installed TES capacity, kWh
    :param lowT_trip:               Limit below which TES level is considered too low to serve heating demand, kWh
    :param T_w_return:              Water return temperature from heating network, ºC
    :param there_is_heating_demand: Boolean variable that is True: if there heating demand in the current time step, False: on the contrary
    :param TES_tooSmallforORC:      Boolean variable that is True: if the TES is too small to run the ORC (i.e. its full capacity is smaller than the energy needed to run the ORC for a timestep), False: on the contrary    
    :return residual_load:          Residual load from substracting the PV power generation from the total load, kW
    :return P_PV:                   PV power generation, kW
    :return P_ORC:                  ORC power generation, kW
    :return P_genset_min:           Genset power generation for supplying the demand, kW
    :return nextSOC_TES:            TES SOC at the end of the time-step, kWh
    :return Eta_ORC:                Efficiency of the ORC, -
    :return CSP_Pow:                CSP power generation, kW
    :return ChargeTES:              Power charged in TES, kW
    :return SOC_TES_out:            TES SOC in current time step after accounting for charging and self discharge, kWh
    :return TES_loss:               TES losses as a function of ambient temperature and TES State of Charge, kW
    :return PV_on:                  Boolean variable computed at each time step that is True: if PV is generating power, False: otherwise
    :return ORC_on:                 Boolean variable computed at each time step that is True: if the ORC is generating power, False: otherwise
    :return Genset_on:              Boolean variable computed at each time step that is True: if the Genset is generating power, False: otherwise
    :return heating_on:             Boolean variable computed at each time step that is True: if the heating system is generating power, False: otherwise
    :return Pnom_ORC_set:           Actual nominal power produced by ORC depending on the heat source, kW
    :return P_heating:              Heating system power generation, kW
    :return heating_backup_on:      Boolean variable computed at each time step that is True: if the heating backup system is generating power, False: otherwise
    :return P_heating_backup:       Heating backup system power generation, kW
    :return P_dumped_th:            Curtailed power from the heating system, kW
    '''
    
    # Initializations
    ORC_on = False
    Genset_on = False
    Pnom_ORC_set = Pnom_ORC
    
    # Define if PV output is ON/OFF (Recall low-light cutoff and PV_kW=0 checks in solar_ambient)
    PV_on = True if P_PV>0 else False    
    
    # Check if there is residual load, i.e. PV power generation is < than demand
    residual_load = P_load - P_PV
    
    # Electrical power dispatching
    if residual_load>0:
        if not there_is_batt:
            if there_is_ORC and there_is_CSP and (TES_state=='too high' or TES_state=='high'):
                # TES has adequate charge to run ORC. Compute potential power generation to see if residual load is already covered by the ORC
                ORC_on = True
                CSP_Pow,Charge_TES,SOC_TES_out,TES_loss,T_htf=TES_energy_balance(lowT_trip,highT_trip,SOC_TES,TESkWh,P_CSP,timestep,T_amb,Genset_on,Pmax_load)
                P_ORC, Eta_ORC, P_ORC_th = ORC_on_TES_performances(Pnom_ORC, T_amb, T_htf, there_is_heating_demand, heating_ORC, T_w_return)
                if residual_load > P_ORC:
                    Genset_on = True
                else:
                    Genset_on = False
            elif there_is_ORC and not there_is_CSP:
                # ORC bottoming on genset, no TES in system             
                Genset_on = True
                ORC_on = True
                Pnom_ORC_set=0.14*Pmax_load
            else:
                # Not adequate charge on TES / No ORC
                Genset_on = True
                ORC_on = False
                
        elif there_is_batt and (batt_state == 'full' or batt_state == 'part charge' and SOC_bat-low_trip > E_loadleft) and residual_load <= -P_bat_max_discharge:
            # Enough discharge power in batteries to supply residual load
            ORC_on = False
            Genset_on = False
                
        else:
            # Catches cases where there_is_batt and (batt_state == 'empty' or batt_state == 'part charge' and SOC_bat-low_trip < E_loadleft) and/or residual_load > -P_bat_max_discharge
            # Battery is empty / too low to supply load left
            if not there_is_ORC or there_is_ORC and there_is_CSP and (TES_state=='low' or TES_state=='too low'):
                # There is no ORC, or TES is too low to run the ORC
                Genset_on = True
                ORC_on = False
            elif there_is_ORC and not there_is_CSP:
                # ORC bottoming on genset, no TES in system
                Genset_on = True
                ORC_on = True
                Pnom_ORC_set=0.14*Pmax_load
            elif P_CSP > 0 and (TES_state == 'high' or TES_state == 'too high'):
                # When sunny, run ORC not only if TES is 'too high' but also 'high' because sun continues filling it
                ORC_on = True
                # Compute potential power generation to see if residual load is already covered by the ORC or ORC + batteries
                Genset_on = False
                CSP_Pow,Charge_TES,SOC_TES_out,TES_loss,T_htf=TES_energy_balance(lowT_trip,highT_trip,SOC_TES,TESkWh,P_CSP,timestep,T_amb,Genset_on,Pmax_load)
                P_ORC, Eta_ORC, P_ORC_th = ORC_on_TES_performances(Pnom_ORC, T_amb, T_htf, there_is_heating_demand, heating_ORC, T_w_return)
                if -P_bat_max_discharge >= (residual_load - P_ORC):
                    Genset_on = False
                else:
                    Genset_on = True
            elif there_is_ORC and there_is_CSP and (TES_state=='high' or TES_state=='too high'):
                Genset_on = False
                CSP_Pow,Charge_TES,SOC_TES_out,TES_loss,T_htf=TES_energy_balance(lowT_trip,highT_trip,SOC_TES,TESkWh,P_CSP,timestep,T_amb,Genset_on,Pmax_load)
                P_ORC, Eta_ORC, P_ORC_th = ORC_on_TES_performances(Pnom_ORC, T_amb, T_htf, there_is_heating_demand, heating_ORC, T_w_return)
                if -P_bat_max_discharge >= (residual_load - P_ORC):
                    # TES has adequate charge to run ORC and batteries are not too low
                    ORC_on = True
                else:
                    # Save TES for later
                    Genset_on = True
                    ORC_on = False
            else:
                # Adequate charge in TES but battery is too low to power load from ORC
                ORC_on = False
                Genset_on = True


    else:
        # If residual_load < 0, PV is already supplying all load and more
        Genset_on = False
        if not there_is_batt:
            # Do not run ORC/Genset because generated power will be dumped
            ORC_on = False
            Genset_on = False
        elif there_is_batt and batt_state!='full':
            # Check if there is available charge capacity. If so, run ORC if adequate charge in TES. If not, do not run ORC because all generated power will be dumped, defocusing instead
            P_PV_ToBat = -residual_load
            if P_PV_ToBat < P_bat_max_charge and there_is_ORC and there_is_CSP and (TES_state == 'too high' or TES_state == 'high'):
                ORC_on = True
        else:
            ORC_on = False
    
    
    # Initialize
    Eta_ORC = -9e7
    
    # ORC and TES usage
    CSP_Pow,ChargeTES,SOC_TES_out,TES_loss,T_htf=TES_energy_balance(lowT_trip,highT_trip,SOC_TES,TESkWh,P_CSP,timestep,T_amb,Genset_on,Pmax_load)    

    if ORC_on and there_is_CSP:
        # The ORC is running with TES
        P_ORC, Eta_ORC, P_ORC_th = ORC_on_TES_performances(Pnom_ORC, T_amb, T_htf, there_is_heating_demand, heating_ORC, T_w_return)
        nextSOC_TES = SOC_TES_out - P_ORC*timestep/Eta_ORC
    else:
        # The ORC is not running (see that the case when it runs with genset exhaust is considered later)
        P_ORC = 0
        P_ORC_th = 0
        nextSOC_TES = SOC_TES_out
    
    # Genset
    if Genset_on and (there_is_CSP or not there_is_ORC):
        P_genset_min = residual_load - P_ORC
    elif Genset_on and there_is_ORC and not there_is_CSP: 
        # Means ORC is ON, P_ORC depend on genset
        P_genset_min = residual_load / (1 + (0.146775411598832 + 2.66715750826574e-7*T_amb**3 - 2.25242051687131e-5*T_amb**2))
        (P_ORC, P_ORC_th) = ORC_on_Genset_performances(T_amb, P_genset_min, there_is_heating_demand, heating_ORC, T_w_return)
        nextSOC_TES = SOC_TES_out
    else:
        P_genset_min = 0
        

    # Initializations
    heating_on = False
    P_heating = 0
    heating_backup_on = False
    P_heating_backup = 0    
    P_dumped_th = 0
    
    # Maximum available power from TES
    P_TES_max = (nextSOC_TES-lowT_trip)/timestep 
    
    # Heating power dispatching
    if P_load_th>0 and heating_TES:
        residual_load_th = P_load_th - P_TES_max
        if TES_state!='too low':
            heating_on = True
            if residual_load_th>0:
                P_heating = P_TES_max
                heating_backup_on = True
                P_heating_backup = residual_load_th
            else:
                P_heating = P_load_th
                heating_backup_on = False
        else:
            heating_on = False
            heating_backup_on = True
            P_heating_backup = P_load_th
            
    elif P_load_th>0 and heating_ORC:
        if ORC_on:
            residual_load_th = P_load_th - P_ORC_th
            heating_on = True
            if residual_load_th>0:
                P_heating = P_ORC_th
                heating_backup_on = True
                P_heating_backup = residual_load_th
            else:
                P_heating = P_load_th
                heating_backup_on = False
        else:
            heating_on = False
            heating_backup_on = True
            P_heating_backup = P_load_th
        
    else:
        heating_on = False
        heating_backup_on = False
    
    if heating_on and heating_TES:
        nextSOC_TES -= P_heating*timestep

    return residual_load,P_PV,P_ORC,P_genset_min,nextSOC_TES,Eta_ORC,CSP_Pow,ChargeTES,SOC_TES_out,TES_loss,PV_on,ORC_on,Genset_on,heating_on,Pnom_ORC_set,P_heating,heating_backup_on,P_heating_backup,P_dumped_th


def set_state_heat(P_bat_max_discharge,P_bat_max_charge,there_is_ORC,there_is_CSP,there_is_batt,heating_TES,heating_ORC,timestep,T_amb,P_PV,P_CSP,P_load,TES_state,batt_state,Pnom_ORC,Pmax_load,SOC_bat,low_trip,E_loadleft,P_load_th,SOC_TES,highT_trip,TESkWh,lowT_trip,T_w_return,there_is_heating_demand):
    '''  
    Power dispatching for a heat driven control, i.e. the following case:
    - when heat demand is served by the ORC and there is CSP-TES and there is heating demand at this time step
    :param P_bat_max_discharge:     Maximum discharge power in curent timestep, kW
    :param P_bat_max_charge:        Maximum charge power in current time step, kW
    :param there_is_ORC:            Boolean variable that is True: if there is ORC installed in the system, False: on the contrary
    :param there_is_CSP:            Boolean variable that is True: if there is CSP installed in the system, False: on the contrary
    :param there_is_batt:           Boolean variable that is True: if there is a battery bank installed in the system, False: on the contrary
    :param heating_TES:             Boolean variable that is True: if the heating system is the TES, False: on the contrary
    :param heating_ORC:             Boolean variable that is True: if the heating system is the ORC, False: on the contrary
    :param timestep:                Time step, hours
    :param T_amb:                   Ambient temperature, °C
    :param P_PV:                    PV power generation, kW
    :param P_CSP:                   CSP power generation, kW
    :param P_load:                  Electrical demand that must be supplied, kW
    :param TES_state:               TES storage level, classified in: 'too low' - TES temperature is too low and cannot be used, 'low' - TES is charged but not full and cannot run the ORC, 'high' - TES is charged but not full and can run the ORC, 'too high' - TES is full (temperature reached its maximum)
    :param batt_state:              Batteries storage level, classified in: 'empty' - too low level, 'part charge'- batteries are charged but not full, 'full' - batteries are full
    :param Pnom_ORC:                Installed ORC capacity, kW
    :param Pmax_load:               Maximum power output of the load curve, kW
    :param SOC_bat:                 Batteries' state of charge after accounting for self discharge, kWh
    :param low_trip:                Limit below which batteries are considered too low (5 % of batteries full capacity), kWh
    :param E_loadleft:              Amount of load defined by smart charging (during daylight: a value higher than batteries capacity, otherwise: load remaining from premidnight until dawn), kWh
    :param P_load_th:               Thermal demand that must be supplied, kW
    :param SOC_TES:                 TES State of Charge, kWh
    :param highT_trip:              Limit above which TES is considered full and cannot be charget beyond it (95 % of TES full capacity), kWh
    :param TESkWh:                  Installed TES capacity, kWh
    :param lowT_trip:               Limit below which TES level is considered too low to serve heating demand, kWh
    :param T_w_return:              Water return temperature from heating network, ºC
    :param there_is_heating_demand: Boolean variable that is True: if there heating demand in the current time step, False: on the contrary
    :return residual_load:          Residual load from substracting the PV power generation from the total load, kW
    :return P_PV:                   PV power generation, kW
    :return P_ORC:                  ORC power generation, kW
    :return P_genset_min:           Genset power generation for supplying the demand, kW
    :return nextSOC_TES:            TES SOC at the end of the time-step, kWh
    :return Eta_ORC:                Efficiency of the ORC, -
    :return CSP_Pow:                CSP power generation, kW
    :return ChargeTES:              Power charged in TES, kW
    :return SOC_TES_out:            TES SOC in current time step after accounting for charging and self discharge, kWh
    :return TES_loss:               TES losses as a function of ambient temperature and TES State of Charge, kW
    :return PV_on:                  Boolean variable computed at each time step that is True: if PV is generating power, False: otherwise
    :return ORC_on:                 Boolean variable computed at each time step that is True: if the ORC is generating power, False: otherwise
    :return Genset_on:              Boolean variable computed at each time step that is True: if the Genset is generating power, False: otherwise
    :return heating_on:             Boolean variable computed at each time step that is True: if the heating system is generating power, False: otherwise
    :return Pnom_ORC_set:           Actual nominal power produced by ORC depending on the heat source, kW
    :return P_heating:              Heating system power generation, kW
    :return heating_backup_on:      Boolean variable computed at each time step that is True: if the heating backup system is generating power, False: otherwise
    :return P_heating_backup:       Heating backup system power generation, kW
    :return P_dumped_th:            Curtailed power from the heating system, kW
    '''
       
    # Initializations
    ORC_on = False
    Pnom_ORC_set = Pnom_ORC
    heating_on = False
    P_heating = 0
    heating_backup_on = False
    P_heating_backup = 0
    P_dumped_th = 0
    
    # Check ORC potential power
    Genset_on = False
    CSP_Pow,Charge_TES,SOC_TES_out,TES_loss,T_htf=TES_energy_balance(lowT_trip,highT_trip,SOC_TES,TESkWh,P_CSP,timestep,T_amb,Genset_on,Pmax_load)
    P_ORC_potential, Eta_ORC, P_ORC_th_potential = ORC_on_TES_performances(Pnom_ORC, T_amb, T_htf, there_is_heating_demand, heating_ORC, T_w_return)
    
    # Heating power dispatching: set ORC_on
    if P_load_th>0:
        if TES_state == 'high' or TES_state == 'too high':
            ORC_on = True
            heating_on = True
            P_heating = P_ORC_th_potential
            residual_load_th = P_load_th - P_heating
            if residual_load_th>0:
                heating_backup_on = True
                P_heating_backup = residual_load_th
            else:
                # even if residual_load<0, ORC produces the same. Some thermal energy will be dumped
                P_dumped_th = -residual_load_th
        else:
            heating_backup_on = True
            P_heating_backup = P_load_th
    
    # Define if PV output is ON/OFF (Recall low-light cutoff and PV_kW=0 checks in Solar Ambient)
    PV_on = True if P_PV>0 else False
    
    # Check if there is residual load, i.e. PV power generation is < than demand
    residual_load = P_load - P_PV
    
    # Electrical power dispatching
    if residual_load>0:
        if ORC_on:
            # Heat driven control decides to run the ORC, if there is residual load must be supplied by genset
            if residual_load > P_ORC_potential:
                Genset_on = True
            else:
                Genset_on = False
        else:
            # If the heat driven control does not decide to run the ORC, run an electrical driven control
            if not there_is_batt:
                if there_is_ORC and there_is_CSP and (TES_state=='too high' or TES_state=='high'):
                    # TES has adequate charge to run ORC. Check if residual load is already covered by the ORC
                    ORC_on = True
                    if residual_load > P_ORC_potential:
                        Genset_on = True
                    else:
                        Genset_on = False
                elif there_is_ORC and not there_is_CSP:
                    # ORC bottoming on genset, no TES in system
                    Genset_on = True
                    ORC_on = True
                    Pnom_ORC_set=0.14*Pmax_load
                else:
                    # Not adequate charge on TES / No ORC
                    Genset_on = True
                    ORC_on = False
                
            elif there_is_batt and (batt_state == 'full' or batt_state == 'part charge' and SOC_bat-low_trip > E_loadleft) and residual_load <= -P_bat_max_discharge:
                # Enough discharge power in batteries to supply residual load
                ORC_on = False
                Genset_on = False
                
            else:
                # Catches cases where there_is_batt and (batt_state == 'empty' or batt_state == 'part charge' and SOC_bat-low_trip < E_loadleft) and/or residual_load > -P_bat_max_discharge
                # Battery is empty / too low to supply load left
                if not there_is_ORC or there_is_ORC and there_is_CSP and (TES_state=='low' or TES_state=='too low'):
                    # There is no ORC, or TES is too low to run the ORC
                    Genset_on = True
                    ORC_on = False
                elif there_is_ORC and not there_is_CSP:
                    # ORC bottoming on genset, no TES in system
                    Genset_on = True
                    ORC_on = True
                    Pnom_ORC_set=0.14*Pmax_load
                elif P_CSP > 0 and (TES_state == 'high' or TES_state == 'too high'):
                    # When sunny, run ORC even if TES is not 'too high' but just 'high' because sun continues filling it
                    ORC_on = True
                    # Check if residual load is already covered by the ORC or ORC + batteries
                    if -P_bat_max_discharge >= (residual_load - P_ORC_potential):
                        Genset_on = False
                    else:
                        Genset_on = True
                elif there_is_ORC and there_is_CSP and (TES_state=='high' or TES_state=='too high') and -P_bat_max_discharge >= (residual_load - P_ORC_potential):
                    # TES has adequate charge to run ORC and batteries are not too low
                    ORC_on = True
                    Genset_on = False
                else:
                    # Adequate charge in TES but battery is too low to power load from ORC
                    ORC_on = False
                    Genset_on = True


    else:
        # If residual_load < 0, PV is already supplying all load and more
        if ORC_on:
            # Heat driven control decides to run the ORC, if there is residual load must be supplied by genset
            Genset_on = False
        else:
            # If the heat driven control does not decide to run the ORC, run an electrical driven control
            Genset_on = False
            if not there_is_batt:
                # Do not run ORC/Genset because generated power will be dumped
                ORC_on = False
                Genset_on = False    
            elif there_is_batt and batt_state!='full':
                # Check if there is available charge capacity. If so, run ORC if adequate charge in TES. If not, do not run ORC because all generated power will be dumped, defocusing instead
                P_PV_ToBat = -residual_load
                if P_PV_ToBat < P_bat_max_charge and there_is_ORC and there_is_CSP and (TES_state == 'too high' or TES_state == 'high'):
                    ORC_on = True
                else:
                    ORC_on = False
    
    
    # Initialize
    Eta_ORC = -9e7
    
    # ORC and TES usage
    CSP_Pow,ChargeTES,SOC_TES_out,TES_loss,T_htf=TES_energy_balance(lowT_trip,highT_trip,SOC_TES,TESkWh,P_CSP,timestep,T_amb,Genset_on,Pmax_load)
    
    if ORC_on and there_is_CSP:
            P_ORC, Eta_ORC, P_ORC_th = ORC_on_TES_performances(Pnom_ORC, T_amb, T_htf, there_is_heating_demand, heating_ORC, T_w_return)
            nextSOC_TES = SOC_TES_out - P_ORC*timestep/Eta_ORC
    else:
        P_ORC = 0
        nextSOC_TES = SOC_TES_out
    
    if Genset_on:
        P_genset_min = residual_load - P_ORC
    else:
        P_genset_min = 0    
    
    return residual_load,P_PV,P_ORC,P_genset_min,nextSOC_TES,Eta_ORC,CSP_Pow,ChargeTES,SOC_TES_out,TES_loss,PV_on,ORC_on,Genset_on,heating_on,Pnom_ORC_set,P_heating,heating_backup_on,P_heating_backup,P_dumped_th


def set_state(ORC_on,Genset_on,heating_on):
    '''
    Define current state based on which power generators are used. They indicate the direction of the power flows
    under various configurations as follows (for generators in order of ORC, Genset, Heating system):
        STATE B    = 0 0 0 
        STATE H    = 0 0 1 
        STATE G    = 0 1 0
        STATE GH   = 0 1 1
        STATE O    = 1 0 0
        STATE OH   = 1 0 1
        STATE OG   = 1 1 0
        STATE OGH  = 1 1 1
    :param ORC_on:          Boolean variable that defines if ORC is True: ON, False: OFF, -
    :param Genset_on:       Boolean variable that defines if genset is True: ON, False: OFF, -
    :param heating_on:      Boolean variable that defines if the heating system is True: generating, False: otherwise, -
    :return current_state:  System state in the current timestep, -
    '''
    if (ORC_on,Genset_on,heating_on)==(0,0,0):
        current_state = 'B'
    elif (ORC_on,Genset_on,heating_on)==(0,0,1):
        current_state = 'H'
    elif (ORC_on,Genset_on,heating_on)==(0,1,0):
        current_state = 'G'
    elif (ORC_on,Genset_on,heating_on)==(0,1,1):
        current_state = 'GH'     
    elif (ORC_on,Genset_on,heating_on)==(1,0,0):
        current_state = 'O'
    elif (ORC_on,Genset_on,heating_on)==(1,0,1):
        current_state = 'OH'
    elif (ORC_on,Genset_on,heating_on)==(1,1,0):
        current_state = 'OG'
    elif (ORC_on,Genset_on,heating_on)==(1,1,1):
        current_state = 'OGH'
    return current_state

def TES_energy_balance(lowT_trip,highT_trip,SOC_TES_in,TESkWh,P_CSP_potential,timestep,T_amb,Genset_on,Pmax_genset):
    '''
    Energy balance on TES tank. Calculate addition of energy from Genset WHR and from solar collectors.
    Calculate TES losses as a function of T_amb and TES SOC (as a proxy for average Tank T) extracted from a
    dataset exercising the Schuman model (10-node) for a TES of up to 40 MWh.
    Calculation of the TES exhaust temperature.
    :param lowT_trip:       Limit below which TES level is considered too low to serve heating demand, kWh
    :param highT_trip:      Limit above which TES is considered full and cannot be charget beyond it (95 % of TES full capacity), kWh
    :param SOC_TES_in:      TES State of Charge out from the previous time step, kWh
    :param TESkWh:          Installed TES capacity, kWh
    :param P_CSP_potential: CSP power generation (function of irradiation and CSP efficiency), kW
    :param timestep:        Time step, hours
    :param T_amb:           Ambient temperature, °C
    :param Genset_on:       Boolean variable that defines if genset is ON/OFF, -
    :param Pmax_genset:     Genset maximum power, defined as maximum power output of the load curve, kW
    :return CSP_Pow:        CSP power thermal addition to TES, kW
    :return ChargeTES:      Total thermal energy additions to TES from CSP and/or Genset, kW
    :return SOC_TES_out:    TES State of charge in current time step after accounting for charging and self discharge, kWh
    :return TES_loss:       TES losses as a function of ambient temperature and TES State of Charge, kW
    :return T_htf_ex_TES:   Heat transfer fluid temperature at the TES outlet, ºC
    '''    
    if TESkWh>0:
        if Genset_on and SOC_TES_in < highT_trip:
            # Thermal potential for WHR from the genset exhaust derived from LP GENSET.EES model created for the Shell project, fit using Eureqa R2=0.997
            Genset_Pow = 0.932071139000018 + 0.00521289941525483*T_amb - 0.0801630128318715*SOC_TES_in/TESkWh - 3.87675754014623e-5*T_amb**2
        else:
            Genset_Pow=0
        
        if SOC_TES_in >= highT_trip:
            CSP_Pow=0                                       # TES cannot be charged beyond full capacity, collectors are defocused if TES is full (above high cutoff)
        else:
            CSP_Pow=P_CSP_potential
        
        ChargeTES = CSP_Pow + Genset_Pow                    # Total thermal energy additions to TES
        
        TES_loss = 0.000831064446466728*TESkWh + 3.49875188913509*np.arccosh(2.0551369082496 + 0.00148277899794591*TESkWh) - 4.37187301750623 - 7.37880204435673e-6*TESkWh*T_amb - 2.14136361650572e-9*TESkWh**2
        
        P_charge_max = (highT_trip - SOC_TES_in)/timestep   # Limit charge capacity of TES
        if ChargeTES > P_charge_max:
            ChargeTES = P_charge_max                        # Charge to the limit
        
        SOC_TES_out = SOC_TES_in + (ChargeTES-TES_loss)*timestep
        
        T_htf_max_TES = 180
        T_htf_min_TES = 130
        T_htf_ex_TES = T_htf_min_TES + (T_htf_max_TES - T_htf_min_TES)*(SOC_TES_out - lowT_trip)/(highT_trip - lowT_trip)
    
    else:
        ChargeTES = 0
        SOC_TES_out = 0
        TES_loss = 0
        CSP_Pow = 0
        T_htf_ex_TES = 0
    
    return CSP_Pow,ChargeTES,SOC_TES_out,TES_loss,T_htf_ex_TES



    
def ORC_on_TES_performances(Pnom_ORC, T_amb, T_htf, there_is_heating_demand, heating_ORC, T_w_return):
    '''
    Get characteristics of the ORC based on an efficiency law and a power law (fitted from a detailed 
    model in EES) for when the ORC heat source is the TES.
    :param Pnom_ORC:                Installed ORC capacity, kW
    :param T_amb:                   Ambient temperature, ºC
    :param T_htf:                   Heat Transfer Fluid Temperature at the inlet of ORC evaporator, ºC 
    :param there_is_heating_demand: Boolean variable that is True: if there heating demand in the current time step, False: on the contrary
    :param heating_ORC:             Boolean variable that is True: if the heating system is the ORC, False: on the contrary
    :param T_w_return:              Water return temperature from heating network, ºC
    :return P_ORC:                  ORC electrical output power when the heat source is the TES, kW
    :return eta_cycle:              Efficiency of ORC, -
    :return P_ORC_th:               ORC thermal output power when the heat source is the TES, kW
    '''
    if there_is_heating_demand and heating_ORC:
        # Water condenser
        T_cf_su_cd = T_w_return
        eta_carnot = 1 - (T_cf_su_cd +273.15)/(T_htf + 273.15)
        eta_II = -5.88069776E-01+4.28358532E-03*T_cf_su_cd-4.76586742E-05*T_cf_su_cd**2+1.41450496E-02*T_htf-5.67157252E-05*T_htf**2
        X_Q_dot_heating = (-5.49687048E+04-9.20619550E+01*T_cf_su_cd+5.88049918E+02*T_htf)/5000
        X_W_dot_net = (-5.17739828E+03-2.87192524E+01*T_cf_su_cd+6.19595654E+01*T_htf)/5000
    else:
        # No heating demand / Heating by TES: Air condenser
        T_cf_su_cd = T_amb
        eta_carnot = 1 - (T_cf_su_cd +273.15)/(T_htf + 273.15)
        eta_II = -1.83049346E-01+1.95586114E-03*T_amb-1.12194277E-05*T_amb**2+8.90005113E-03*T_htf-3.84177089E-05*T_htf**2
        X_W_dot_net = (-6.37092564E+03-1.71805264E+01*T_amb+6.78053231E+01*T_htf)/5000
        X_Q_dot_heating = 0
 
    eta_cycle = eta_carnot*eta_II
    P_ORC = Pnom_ORC*X_W_dot_net
    P_ORC_th = Pnom_ORC*X_Q_dot_heating

    return P_ORC, eta_cycle, P_ORC_th


    
def ORC_on_Genset_performances(T_amb, P_genset, there_is_heating_demand, heating_ORC, T_w_return):
    '''
    Get characteristics of the ORC based on a power law (fitted from a detailed model in EES)
    for when the ORC is bottoming on the genset.
    :param T_amb:                   Ambient temperature, ºC
    :param P_genset:                Output power from genset, kW
    :param there_is_heating_demand: Boolean variable that is True: if there heating demand in the current time step, False: on the contrary
    :param heating_ORC:             Boolean variable that is True: if the heating system is the ORC, False: on the contrary
    :param T_w_return:              Water return temperature from heating network, ºC
    :return P_ORC:                  ORC electrical output power when the heat source is the genset, kW
    :return P_ORC_th:               ORC thermal output power when the heat source is the genset, kW   
    '''
    # Thermal potential for WHR from the genset exhaust derived from LP GENSET.EES model created for the Shell project, fit using Eureqa R2=0.997:
    P_ORC = P_genset * (0.146775411598832 + 2.66715750826574e-7*T_amb**3 - 2.25242051687131e-5*T_amb**2)
    
    if there_is_heating_demand and heating_ORC:
        T_cf_su_cd = T_w_return
        T_hf = 442
        eta_carnot = 1 - (T_cf_su_cd +273.15)/(T_hf + 273.15)
        eta_II = -5.88069776E-01+4.28358532E-03*T_cf_su_cd-4.76586742E-05*T_cf_su_cd**2+1.41450496E-02*T_hf-5.67157252E-05*T_hf**2
        eta_cycle = eta_carnot*eta_II
        P_ORC_th = P_ORC / eta_cycle - P_ORC
    else:
        P_ORC_th = 0
    
    return P_ORC,P_ORC_th



def P_bat_in_out(P_genset,residual_load,P_ORC,Genset_on,heating_on,there_is_batt,Pmax_genset,P_bat_max_charge):
    '''    
    Calculation of the power flows to and/or from the batteries and the total genset output power.
    :param P_genset:            Output power from genset, kW
    :param residual_load:       Residual load from substracting the PV power generation from the total load, kW
    :param P_ORC:               ORC electrical output power when the heat source is the genset, kW
    :param Genset_on:           Boolean variable that defines if genset is ON/OFF, -
    :param heating_on:          Boolean variable that defines if the heating system is True: generating, False: otherwise, -
    :param there_is_batt:       Boolean variable that is True: if there is a battery bank installed in the system, False: on the contrary
    :param Pmax_genset:         Genset maximum power, defined as maximum power output of the load curve, kW
    :param P_bat_max_charge:    Batteries maximum charge power in current time step, kW
    :return P_PV_ToBat:         PV power to charging the batteries, kW
    :return P_bat_dis:          Load powered from the inverter, kW
    :return P_genset_ToBat:     Genset power to charging the batteries, kW
    :return P_ORC_ToBat:        ORC power to charging the batteries, kW
    :return P_genset:           Genset power output, kW
    '''
    # Batteries and Inverter
    if residual_load<0:
        # If residual_load is negative there is power available
        P_PV_ToBat = -residual_load     # PV charging batteries
        P_bat_dis = 0                   # PV is powering entire load
        P_ORC_ToBat = P_ORC
    else:
        # The difference is positive, the load is more than the PV can supply
        P_PV_ToBat = 0                  # There is no power surplus from the PV
        if P_ORC > residual_load:
            # ORC has excess available beyond supplying load
            P_ORC_ToBat = P_ORC - residual_load
            P_bat_dis = 0
        else:
            # PV + ORC not enough to supply load, also need inverter
            P_ORC_ToBat = 0
            if not Genset_on and there_is_batt:
                P_bat_dis = P_ORC - residual_load
            else:
                P_bat_dis = 0           # The inverter is not supplying any loads if genset is on
            if not Genset_on and not there_is_batt:
                print 'WARNING: ORC is not producing at nominal power so small part of load is not beign supplied? P_ORC =', P_ORC
    
    # Determine whether and how much the genset is charging the battery bank based on presence of a battery bank, residual genset capacity, and battery charge rate limit
    if there_is_batt and Genset_on:
        if P_genset < Pmax_genset:
            # If loads fall below its peak output the generator will charge the battery bank up to the full charge rate
            Difference = Pmax_genset - P_genset                     # Residual genset capacity
            if Difference > P_bat_max_charge - P_ORC_ToBat:
                P_genset_ToBat = P_bat_max_charge - P_ORC_ToBat     # If capacity is higher than the charge rate limit, the genset charges at the limit
            else:
                P_genset_ToBat = Difference                         # The genset charges up to its full residual capacity after supplying the load
            P_genset = P_genset + P_genset_ToBat
        else:
            # Supplying peak load
            P_genset_ToBat = 0
    else:
        P_genset_ToBat=0
            
    return P_PV_ToBat,P_bat_dis,P_genset_ToBat,P_ORC_ToBat,P_genset



def genset(Pmax,P):
    '''
    Efficiency curve for the LPG genset 
    "RV generator set Quiet Gasoline TM Series RV QG 4000"
    Datasheet:
    http://coloradostandby.com/media/custom/upload/document_technical_a-1399.pdf
    The fuel is asusmed to be LPG
    According to the datasheet, the maximum efficiency is quite low: 19.87%
    Refer to the attached EES file to see the details of the regression
    :param Pmax:        Nominal genset power, kW
    :param P:           Current requested power, kW
    :return P_el:       Generated power, kW
    :return Mdot_lpg:   Current fuel consumption, kg/s 
    :return eta:        Current genset efficiency
    '''
    eta_max = 0.1987
    partload = np.minimum(1,np.maximum(0,P/Pmax))
    P_el = partload * Pmax
    Qdot_lpg_max = Pmax/eta_max
    fuel_ratio = 0.384615 + 0.923077*partload - 0.307692*partload**2
    Qdot_lpg = fuel_ratio * Qdot_lpg_max
    eta = P_el/Qdot_lpg
    Mdot_lpg_kgh = Qdot_lpg/46e6*3600
    
    return P_el,eta,Qdot_lpg,Mdot_lpg_kgh



def fuel_calcs(P_genset,Pmax_load,timestep,P_heating_backup):
    '''
    Determine fuel consumption based on genset power output and taking into account
    part load effect on its efficiency
    :param P_genset:    Genset power output, kW
    :param Pmax_load:   Genset maximum power, defined as maximum power output of the load curve, kW
    :param timestep:    Current time step, hours
    :return Fuel_kW:    Fuel consumption, kW
    :return Fuel_kg:    Fuel massic consumption, kg propane
    '''  
    if P_genset>0:
        P_el,Eta_genset,Qdot_genset,Mdot_lpg_kgh = genset(Pmax_load,P_genset)
    else:
        Qdot_genset=0
        Eta_genset=-1
    
    # Burner
    eta_burner = 0.82                              # This type of burners have very high efficiency even when working in partload
    Qdot_burner = P_heating_backup/eta_burner      # kW

    Qdot_fuel = Qdot_burner + Qdot_genset
    m_burner = Qdot_burner*timestep*3600/46E3
    m_genset = Qdot_genset*timestep*3600/46E3
    # Fuel consumption
    Fuel_kJ=Qdot_fuel*timestep*3600                 # kJ
    Fuel_kg=Fuel_kJ/46E3                            # Converts kJ to kg propane
    
    return Qdot_fuel,Fuel_kg,Eta_genset,m_burner,m_genset


def solar_ambient(th_Hour,NOCT,Pmax_Tco,Pnom_PV,SOC_TES,TESkWh,CSP_A,lowlightcutoff):
    '''
    Determine ambient conditions and perform the calculations related to irradiance and temperature.
    Determine CSP and PV power outputs for those conditions, i.e. the available power.
    :param th_Hour:             Hour corresponding to current time step, hours
    :param NOCT:                Nominal Operating Cell Temperature on the PV panel, ºC
    :param Pmax_Tco:            Outdoor temperature coefficient (TCo) of the maximum power (Pmax), %/C
    :param Pnom_PV:             Installed PV capacity, kW
    :param SOC_TES:             TES STate of Charge, kWh
    :param TESkWh:              Installed TES capacity, kWh
    :param CSP_A:               Installed CSP surface, m²
    :param lowlightcutoff:      DNI level under which solar technologies not responsive, W/m²
    :return I_b:                Beam irradiation, W/m²
    :return T_amb:              Ambient temperature, ºC
    :return P_PV_potential:     PV available power, kW
    :return P_CSP_potential:    CSP available power, kW
    '''
    from ugrid_tools_ThEl import CSP_power, PV_power
    
    T_amb = f_Tamb(th_Hour)
    
    if CSP_A>0 or Pnom_PV>0:
        I_b = f_Ib(th_Hour)
    else:
        I_b = -999
    
    T_HTF_col_su = 150+30*(SOC_TES/TESkWh)
    
    P_CSP_potential = CSP_power(CSP_A,I_b,lowlightcutoff,T_HTF_col_su,T_amb)
    P_PV_potential = PV_power(Pnom_PV,NOCT,Pmax_Tco,f_I_PV(th_Hour),lowlightcutoff,T_amb)
    
    return I_b,T_amb,P_PV_potential,P_CSP_potential
    

def operation(Pmax_Tco, NOCT, smart, lat, lon, CSP_A, Pnom_ORC, Pnom_PV, TESkWh, BattkWh, Pmax_load, heating_TES, heating_ORC):
    '''
    :param Pmax_Tco:        Outdoor temperature coefficient (TCo) of the maximum power (Pmax), %/C
    :param NOCT:            Nominal Operating Cell Temperature on the PV panel, ºC
    :param smart:           Smart charging indicator (1: ON, 0: OFF), -
    :param lat:             Latitude, º
    :param lon:             Longitude, º
    :param CSP_A:           Installed CSP surface, m²
    :param Pnom_ORC:        Installed ORC capacity, kW
    :param Pnom_PV:         Installed PV capacity, kW
    :param TESkWh:          Installed TES capacity, kWh
    :param BattkWh:         Installed batteries full capacity, kWh
    :param Pmax_load:       Maximum power output of the load curve, kW
    :return Batt_kWh_tot:   Total energy charged to the batteries, kWh
    :return Propane:        Total propane massic consumption, kg
    :return DNI:            
    :return hour:           Array of hour of the year at each timestep, h
    :return P_genset:       Array of genset output power at each timestep, kW
    :return P_genset_ToBat: Array of genset power to charging the batteries, kW
    :return P_PV_ToBat:     Array of PV power to charging the batteries, kW
    :return P_bat_inv:      Array of load powered from the inverter, kW     
    :return SOC_bat:        Array of batteries SOC, kWh
    :return SOC_TES:        Array of TES SOC, kWh
    :return dayHour:        Array of hour of the day, h
    :return P_load:         Array of electrical demand that is being supplied, kW
    :return E_loadleft:     Array of amount of load defined by smart charging (during daylight: a value higher than batteries capacity, otherwise: load remaining from premidnight until dawn), kWh
    :return State:          Array of current system state, -
    :return beam:           Array of beam irradiation at each timestep, W/m^2
    :return T_ambient:      Array of ambient temperature at each timestep, ºC
    :return P_PV:           Array of PV power generation, kW
    :return P_CSP:          Array of CSP power generation, kW
    :return P_ORC_ToBat:    Array of ORC power to charging the batteries, kW
    :return ORC_Eta:        Array of ORC efficiency, -
    :return P_ORC:          Array of ORC electrical output power, kW
    :return ChargeTES:      Array of TES charge power, kW      
    :return TES_loss:       Array of TES power losses, kW
    :return Charge:         Array of batteries Charged power at each timestep, kW
    :return P_dumped:       Array of curtailed power because of the battery limit at each timestep, kW
    :return Batt_frac:      Array of the fraction of batteries SOC on their full capacity, -
    :return Genset_fuel:    Array of propane massic consumption at each timestep, kg
    :return Fuel_kW:        Array of propane consumption at each timestep, kW
    '''
    # Define system characteristics 
    there_is_PV = True if Pnom_PV>0 else False
    there_is_ORC = True if Pnom_ORC>0 else False
    there_is_CSP = True if CSP_A>0 else False
    there_is_batt = True if BattkWh>0 else False
        
    
    if there_is_CSP and not there_is_ORC and not heating_TES:
        Propane = -9999
        print 'WARNING: Pnom_ORC CANNOT BE ZERO WHEN CSP_A>0 and heating is not provided by TES --> ERROR!!!'
    
    if not there_is_ORC and heating_ORC:
        print 'WARNING: Pnom_ORC is ZERO and HEATING SYSTEM IS ORC --> USING BURNER ALL THE TIME!!!'
    if not TESkWh>0 and heating_TES:
        print 'WARNING: THERE IS NO TES and HEATING SYSTEM IS TES --> USING BURNER ALL THE TIME!!!'
    
    # Initialize timestep:
    timestep=1
    
    # Days of the year:
    year_day_number = 365
    hmax=int(round(year_day_number*24/timestep))
        
    # Parametrizations of year
    hour=np.linspace(0,hmax-1,num=hmax,endpoint=True)        
          
    T_w_return = 45
    # Initialize variables
    Batt_kWh_tot = 0
    DNI=0
    Propane=0
    fuel_kg=0
    Emax_bat=BattkWh/3          # Battery can only be charged at rate up to arbitrary 1/5 of its full capacity rating
    Pmax_bat = Emax_bat/timestep
    lowlightcutoff=100          # DNI level under which solar technologies not responsive
    
    # Battery:    
    high_trip = 0.95*BattkWh    # kWh
    low_trip = 0.05*BattkWh     # kWh
    I_b=0
    P_genset=np.zeros(hmax)
    m_burner=np.zeros(hmax)
    m_genset=np.zeros(hmax)
    Eta_genset=np.zeros(hmax)
    P_genset_ToBat=np.zeros(hmax)
    P_PV_ToBat=np.zeros(hmax)
    P_bat_dis=np.zeros(hmax)
    SOC_bat=np.zeros(hmax)
    SOC_TES=np.zeros(hmax)
    dayHour=np.zeros(hmax)
    P_load=np.zeros(hmax)
    E_loadleft=np.zeros(hmax)
    P_load_th=np.zeros(hmax) 
    P_heating=np.zeros(hmax)
    P_heating_backup=np.zeros(hmax)
    P_dumped_th=np.zeros(hmax)
    
    # Create cumulative variables array for later update
    State=[]
    beam=np.zeros(hmax)
    T_ambient=np.zeros(hmax)
    P_PV=np.zeros(hmax)
    P_CSP=np.zeros(hmax)
    P_ORC_ToBat=np.zeros(hmax)
    ORC_Eta=np.zeros(hmax)
    P_ORC=np.zeros(hmax)
    ChargeTES=np.zeros(hmax)
    TES_loss=np.zeros(hmax)
    Charge=np.zeros(hmax)
    P_dumped=np.zeros(hmax)
    Batt_frac=np.zeros(hmax)
    Genset_fuel=np.zeros(hmax)
    Fuel_kW=np.zeros(hmax)
    
    # TES utility range
    if TESkWh==0:
        there_is_CSP = False
        
    mediumT_trip = Pnom_ORC*timestep/0.1                    # kWh, approx amount of energy withdrawn in one timestep if the ORC is ON (set for minimum possible ORC efficiency of 10%)
    highT_trip = TESkWh*0.95                                # kWh
    #lowT_trip = mediumT_trip - (highT_trip-mediumT_trip)   # kWh
    lowT_trip = - highT_trip                                # kWh
    if there_is_CSP and mediumT_trip > highT_trip:
        print 'WARNING: TES is too small to run ORC!'
        TES_tooSmallforORC = True        # ORC will bottom on genset
    else:
        TES_tooSmallforORC = False

        
    # If battery bank is included in generation equipment, initialize it to half full / else no battery bank:
    SOC_bat[0] = BattkWh/2 if there_is_batt else 0
    
    # If thermal storage is included in the generation equipment, initialize it to cold start / else there is no TES tank:
    SOC_TES[0] = TESkWh/2 if TESkWh>0 else 0
    
    
    for h in range(hmax):
        ##### Analyze external conditions #####       
        
        # Establish the time of day in hours:       
        dayHour[h] = np.mod(hour[h],24)
#        day = np.ceil(hour[h]/24.)
        
        if there_is_CSP or there_is_PV or there_is_ORC:
            # Assess solar avaiability
            (I_b,T_amb,P_PV_potential,P_CSP_potential)=solar_ambient(hour[h],NOCT,Pmax_Tco,Pnom_PV,SOC_TES[h],TESkWh,CSP_A,lowlightcutoff)
        else:
            T_amb = f_Tamb(h)
            P_PV_potential=0
            P_CSP_potential=0
            I_b=-9999
        
        # LOAD (kW) for this timestep
        # calculate the row value in the Nkautest simulated load dataset (based on Hobhouse NRS) corresponding to the simulation timestep
        Load_Clinic=0.5*(0.690641672 - 0.11362214*dayHour[h] - 0.0000869340922*dayHour[h]**2 + 0.00596194972*dayHour[h]**3 - 0.000750723835*dayHour[h]**4 + 0.0000343476381*dayHour[h]**5 - 5.45160266E-07*dayHour[h]**6)
        loadcounter=int(hour[0]/timestep) + h
        
        P_load[h] = P_loadT.Load[loadcounter] + Load_Clinic
        
        # Smart charging strategy - the generator will only continue charging the battery bank if the SOC of the bank
        # is less than the SOC which would permit inverter operation until dawn (when the PV can ideally be prioritized
        # to recharge the bank at a lower operating cost than genset charging):
        if smart == 1:
            # estimate load remaining until dawn for nighttimes (set zero for daytimes)
            if dayHour[h]>17 or dayHour[h]<7:
                # determine the amount of load in kWh remaining from premidnight until dawn
                E_loadleft[h] = FullYearEnergy.FYE[loadcounter]
            else:
                E_loadleft[h] = 0             # Daytime load prediction not part of algorithm
        else:
            E_loadleft[h] = 10000 + 2*BattkWh # Number much higher than can possibly be stored in battery bank
            # This forces genset to stay on (not smart)
        
        # Thermal load
        if heating_TES or heating_ORC:
            heating = load_xls('Data/Heating_1h.xlsx','Sheet1')
            P_load_th[h] = heating.Q_dot_W[h] / 1000
            
        if P_load_th[h] > 0:
            there_is_heating_demand = True
        else:
            there_is_heating_demand = False
        
        # Determine dis/charging capacity of batteries
        P_bat_max_discharge,P_bat_max_charge,high_trip,SOC_bat_inout = P_bat_minmax(T_amb, Pmax_bat, SOC_bat[h], BattkWh, timestep)        
        
        # Determine how full are thermal - chemical storage:
        batt_state,TES_state=storage_levels(SOC_bat_inout,low_trip,high_trip,SOC_TES[h],lowT_trip,mediumT_trip,highT_trip)
        
        # Determine state:
        if TES_tooSmallforORC:
            there_is_CSP = False
        if heating_TES or heating_ORC and not there_is_CSP or heating_ORC and there_is_CSP and not there_is_heating_demand:
            control = 'elec driven'
            (residual_load,P_PV_h,P_ORC_h,P_genset_min,nextSOC_TES,eta,CSP_Pow,ChTES,SOC_TES_out,TESloss,PV_on,ORC_on,Genset_on,heating_on,Pnom_ORC_set,P_heating_h,heating_backup_on,P_heating_backup_h,P_dumped_th_h) = set_state_elec(P_bat_max_discharge,P_bat_max_charge,there_is_ORC,there_is_CSP,there_is_batt,heating_TES,heating_ORC,timestep,T_amb,P_PV_potential,P_CSP_potential,P_load[h],TES_state,batt_state,Pnom_ORC,Pmax_load,SOC_bat_inout,low_trip,E_loadleft[h],P_load_th[h],SOC_TES[h],highT_trip,TESkWh,lowT_trip,T_w_return,there_is_heating_demand,TES_tooSmallforORC)
        if heating_ORC and there_is_CSP and there_is_heating_demand:
            control = 'heat driven'
            (residual_load,P_PV_h,P_ORC_h,P_genset_min,nextSOC_TES,eta,CSP_Pow,ChTES,SOC_TES_out,TESloss,PV_on,ORC_on,Genset_on,heating_on,Pnom_ORC_set,P_heating_h,heating_backup_on,P_heating_backup_h,P_dumped_th_h) = set_state_heat(P_bat_max_discharge,P_bat_max_charge,there_is_ORC,there_is_CSP,there_is_batt,heating_TES,heating_ORC,timestep,T_amb,P_PV_potential,P_CSP_potential,P_load[h],TES_state,batt_state,Pnom_ORC,Pmax_load,SOC_bat_inout,low_trip,E_loadleft[h],P_load_th[h],SOC_TES[h],highT_trip,TESkWh,lowT_trip,T_w_return,there_is_heating_demand)
        if TES_tooSmallforORC:
            there_is_CSP = True
        
        #(current_state,PV_on,ORC_on,Genset_on,heating_on,Pnom_ORC_set,P_ORC_th,P_heating,heating_backup_on,P_heating_backup) = set_state(there_is_ORC,there_is_CSP,there_is_batt,heating_TES,heating_ORC,timestep,T_amb,Pmax_bat,P_PV_potential,CSP_Avail_Pow,P_load[h],TES_state,batt_state,Pnom_ORC,BattkWh, Pmax_load,SOC_bat[h],low_trip,E_loadleft[h],P_load_th[h],SOC_TES[h],highT_trip,TESkWh)
        current_state = set_state(ORC_on,Genset_on,heating_on)
        # Set system variables:
        PV_BCP,Inv_BDP,Gen_BCP,ORC_BCP,P_genset_h = P_bat_in_out(P_genset_min,residual_load,P_ORC_h,Genset_on,heating_on,there_is_batt,Pmax_load,P_bat_max_charge)
        
        # Balance energy storage in batteries:
        P_bat_charging,Ch,SOC_bat_out,dump = batt_bal(P_bat_max_discharge, P_bat_max_charge, PV_BCP, Inv_BDP, Gen_BCP, ORC_BCP, SOC_bat_inout, BattkWh, timestep)
        
        # Calculate fuel usage:
        if P_genset_h>0 and P_heating_backup_h>0:
            FkW,fuel_kg,Eta_genset_h,m_burner_h,m_genset_h = fuel_calcs(P_genset_h,Pmax_load,timestep,P_heating_backup_h)
        elif P_genset_h>0 and P_heating_backup_h<=0:
            FkW,fuel_kg,Eta_genset_h,m_burner_h,m_genset_h = fuel_calcs(P_genset_h,Pmax_load,timestep,P_heating_backup_h)
            m_burner_h=0
        elif P_genset_h<=0 and P_heating_backup_h>0:
            FkW,fuel_kg,Eta_genset_h,m_burner_h,m_genset_h = fuel_calcs(P_genset_h,Pmax_load,timestep,P_heating_backup_h)
            FkW=0
            m_genset_h = 0
            Eta_genset_h = -1
        else:
            FkW=0
            fuel_kg=0
            Eta_genset_h = -1
            m_burner_h=0
            m_genset_h=0
        
        
        ##### Update cumulative variables #####
        
        # State variables:
        #State[h]=current_state
        State.append(current_state)
        
        # Env parameters:
        beam[h]=I_b                             # to track irradiance
        DNI = DNI + beam[h]*timestep/1000
        T_ambient[h]=T_amb                      # to track ambient temperature
        
        # Solar:
        P_PV[h] = P_PV_h                   # kW electrical
        P_CSP[h] = CSP_Pow                 # kW thermal
        P_PV_ToBat[h] = PV_BCP
        P_ORC_ToBat[h] = ORC_BCP
        ORC_Eta[h] = eta
        P_bat_dis[h] = Inv_BDP
        P_ORC[h] = P_ORC_h
        
        # TES:
        ChargeTES[h] = ChTES
#        SOC_TES[h] = SOC_TES_out       #commented because what I really wanna see is available power and see if model works OK, instead of level of TES after being discharged
        TES_loss[h] = TESloss
        if h!=hmax-1:
            SOC_TES[h+1] = nextSOC_TES
        
        # Batteries:
        Charge[h] = Ch
        Batt_kWh_tot += P_bat_charging*timestep
        P_dumped[h] = dump
        SOC_bat[h]=SOC_bat_out
        if BattkWh>0:
            Batt_frac[h]=SOC_bat[h]/BattkWh
        else:
            Batt_frac[h]=0
        if h!=hmax-1:
            SOC_bat[h+1]=SOC_bat[h]
        
        # Genset:
        m_burner[h] = m_burner_h
        m_genset[h] = m_genset_h
        Eta_genset[h] = Eta_genset_h
        P_genset[h] = P_genset_h
        P_genset_ToBat[h] = Gen_BCP
        Genset_fuel[h]=fuel_kg
        Fuel_kW[h] = FkW
        Propane = Propane + Genset_fuel[h]      # Cumulative genset fuel consumption in kg Propane
        
        # Heating:
        P_heating[h] = P_heating_h
        P_heating_backup[h] = P_heating_backup_h
        P_dumped_th[h] = P_dumped_th_h
    

    if max(SOC_TES)<highT_trip:
        print 'TES is never full'
    
    return Batt_kWh_tot,Propane,DNI,hour,P_genset,P_genset_ToBat,P_PV_ToBat,P_bat_dis,SOC_bat,SOC_TES,dayHour,P_load,E_loadleft,State,beam,T_ambient,P_PV,P_CSP,P_ORC_ToBat,ORC_Eta,P_ORC,ChargeTES,TES_loss,Charge,P_dumped,Batt_frac,Genset_fuel,Fuel_kW,P_load_th,P_heating,P_heating_backup,P_dumped_th,Eta_genset,m_burner,m_genset



print time.clock() - start_time, "seconds"



