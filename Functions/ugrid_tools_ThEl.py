# -*- coding: utf-8 -*-
"""
This file contains all the required functions to run the ugrid model.
These functions have been adapted for the .lib file

Created on Sun May  8 18:05:45 2016

@author: Queralt Altés Buch
"""

from __future__ import division
import numpy as np

def results_analysis(up,verbose=True):
    results = {}
    results['E_load'] = up.P_load.sum()
    results['E_PV'] = up.P_PV.sum()
    results['E_ORC'] = up.P_ORC.sum()
    results['E_genset'] = up.P_genset.sum()
    results['E_dumped'] = up.P_dumped.sum()
    results['E_bat'] = np.abs(up.P_bat).sum()/2
    results['E_TES_losses'] = up.TES_loss.sum()
    results['m_propane'] = up.Genset_fuel.sum()
    results['eta_genset'] = results['E_genset']/(up.Fuel_kW.sum()+0.0000001)
    results['Pmax_load'] = up.P_load.max()
    results['E_load_th'] = up.P_load_th.sum()
    results['m_propane_burner'] = up.m_burner.sum()
    results['m_propane_genset'] = up.m_genset.sum()
    results['E_heating'] = up.P_heating.sum()
    results['E_burner'] = up.P_heating_backup.sum()
    if verbose:
        print 'Yearly values: '
        print 'Maximum load: ' + str(results['Pmax_load']) + ' kW'
        print 'Total load: ' + str(results['E_load']) + ' kWh'
        print 'Total PV generation: ' + str(results['E_PV']) + ' kWh'
        print 'Total ORC generation: ' + str(results['E_ORC']) + ' kWh'
        print 'Total Genset generation: ' + str(results['E_genset']) + ' kWh'
        print 'Total energy curtailed (dumped): ' + str(results['E_dumped']) + ' kWh'
        print 'Total TES losses: ' + str(results['E_TES_losses']) + ' kWh'
        print 'Total thermal load: ' + str(results['E_load_th']) + ' kWh'
        print 'Total heating system generation: ' + str(results['E_heating']) + 'kWh'
        print 'Total burner generation: ' + str(results['E_burner']) + 'kWh'
        print 'Propane consumption from genset: '  + str(results['m_propane_genset']) + ' kg'        
        print 'Propane consumption from burner: '  + str(results['m_propane_burner']) + ' kg'
        print 'Total propane consumption: ' + str(results['m_propane']) + ' kg'
        print 'Average genset efficiency: ' + str(results['eta_genset']*100) + '%'
    return results
    



def dispatch_plot(dispatch,rng=[]):
    '''
    Plotting the results of the dispatch
    Parameters:
        dispatch (dataframe): Values of the main dispatch vectors
        rng (pd datetimeindex): selected index for plotting
    '''
    import matplotlib.pyplot as plt
    import pandas as pd
    
    alpha = '0.3'
    index = dispatch.index
    
    
    demand = dispatch['P_load']
    
    pdrng = rng

    pv = dispatch['P_PV']
    BatteryConsumption = pd.Series(np.maximum(0,dispatch['P_bat']),index=index)
    BatteryGeneration = pd.Series(np.maximum(0,-dispatch['P_bat']),index=index)
    LevelOfCharge = pd.Series(dispatch['SOC_bat'],index=index)

    vec1 = pd.Series(- dispatch['P_dumped'] - BatteryConsumption,index=index)
    vec2 = pd.Series(- BatteryConsumption,index=index)
    vec2b = pd.Series(BatteryGeneration,index=index)
    vec3 = pd.Series(BatteryGeneration+dispatch['P_genset'],index=index)
    vec4 = pd.Series(BatteryGeneration+pv + dispatch['P_genset'],index=index)
    vec5 = pd.Series(BatteryGeneration+pv + dispatch['P_genset'] + dispatch['P_ORC'],index=index)


    fig = plt.figure(figsize=(13,7))
    
    # Create left axis:
    ax = fig.add_subplot(111)
    ax.plot(pdrng,demand[pdrng],color='k')
    
    plt.fill_between(pdrng,vec1[pdrng],vec2[pdrng],color='r',alpha=alpha,hatch="x")
    plt.fill_between(pdrng,vec2[pdrng],0,color='g',alpha=alpha,hatch="x")
    plt.fill_between(pdrng,0,vec2b[pdrng],color='g',alpha=alpha)
    plt.fill_between(pdrng,vec2b[pdrng],vec3[pdrng],color='y',alpha=alpha)
    plt.fill_between(pdrng,vec3[pdrng],vec4[pdrng],color='b',alpha=alpha,hatch="//")
    plt.fill_between(pdrng,vec4[pdrng],vec5[pdrng],color='r',alpha=alpha,hatch="//")

    ax.set_ylabel('Power [kW]')
    ax.yaxis.label.set_fontsize(16)

    # Create right axis:
    ax2 = fig.add_subplot(111, sharex=ax, frameon=False,label='aa')
    ax2.plot(pdrng,LevelOfCharge[pdrng],color='k',alpha=0.3,linestyle='--')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('Battery SOC [kWh]')
    ax2.yaxis.label.set_fontsize(16)

    # Legend:
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines
    a = mpatches.Patch(color='red',alpha=0.3,hatch='x',label='Dumped')
    b = mpatches.Patch(color='green',alpha=0.3,hatch='x',label='Battery charging')
    c = mpatches.Patch(color='green',alpha=0.3,label='Battery generation')
    d = mpatches.Patch(color='yellow',alpha=0.3,label='Genset')
    e = mpatches.Patch(color='blue',alpha=0.3,hatch='//',label='PV')
    f = mpatches.Patch(color='red',alpha=0.3,hatch='//',label='ORC')
    line_demand = mlines.Line2D([], [], color='black',label='Load')
    line_SOC = mlines.Line2D([], [], color='black',alpha=0.3,label='SOC',linestyle='--')

    plt.legend(handles=[a,b,c,d,e,f,line_demand,line_SOC],loc=4)
    
    return True
    
    
def dispatch_ThEl_plot(dispatch,rng=[]):
    '''
    Plotting the results of the thermal dispatch
    Parameters:
        dispatch (dataframe): Values of the main dispatch vectors
        rng (pd datetimeindex): selected index for plotting
    '''
    import matplotlib.pyplot as plt
    import pandas as pd
    
    alpha = '0.3'
    index = dispatch.index
    
    
    demand = dispatch['P_load_th']
    
    pdrng = rng

    dumped = dispatch['P_dumped_th']
    heating = dispatch['P_heating']
    backup = dispatch['P_heating_backup']
    LevelOfCharge = pd.Series(dispatch['SOC_TES'],index=index)

    vec0 = pd.Series(-dumped,index=index)
    vec1 = pd.Series(backup,index=index)
    vec2 = pd.Series(backup + heating,index=index)

    fig = plt.figure(figsize=(13,7))
    
    # Create left axis:
    ax = fig.add_subplot(111)
    ax.plot(pdrng,demand[pdrng],color='k')
    
    plt.fill_between(pdrng,vec0[pdrng],0,color='r',alpha=alpha,hatch="x")
    plt.fill_between(pdrng,0,vec1[pdrng],color='y',alpha=alpha)
    plt.fill_between(pdrng,vec1[pdrng],vec2[pdrng],color='g',alpha=alpha,hatch="x")

    ax.set_ylabel('Power [kW]')
    ax.yaxis.label.set_fontsize(16)

    # Create right axis:
    ax2 = fig.add_subplot(111, sharex=ax, frameon=False,label='aa')
    ax2.plot(pdrng,LevelOfCharge[pdrng],color='k',alpha=0.3,linestyle='--')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('TES SOC [kWh]')
    ax2.yaxis.label.set_fontsize(16)

    # Legend:
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines
    aa = mpatches.Patch(color='red',alpha=0.3,hatch='x',label='Dumped')
    a = mpatches.Patch(color='yellow',alpha=0.3,label='Backup')
    b = mpatches.Patch(color='green',alpha=0.3,hatch='x',label='Heating')
    line_demand = mlines.Line2D([], [], color='black',label='Load')
    line_SOC = mlines.Line2D([], [], color='black',alpha=0.3,label='SOC',linestyle='--')

    plt.legend(handles=[aa,a,b,line_demand,line_SOC],loc=4)
    
    return True


def scale_series(data,index):
    ''' 
    Function that takes a pandas Dataframe as input and interpolates it to the proper index
    '''
    import pandas as pd
    aa = pd.Series(index=index)
    bb = data.join(pd.DataFrame(aa),how='outer')    
    cc = bb.interpolate()
    dd = cc.loc[index,:]
    return dd


def load_xls(filename,sheet,TempPath='.pickle',header=0,skiprows=[],skip_footer=0,index_col=None,parse_col=None,parse_dates=False):
    '''
    Function that loads an xls sheet into and save a temporary pickle version of it. If the pickle is newer
    than the sheet, do no load the sheet again.
    :param file_excel: path to the excel file
    :param sheet: excel sheet to load
    :param TempPath: path to store the temporary data files
    '''
    import os
    import pandas as pd
    
    filepath_pandas = TempPath + '/' + filename.replace('/','-') + '-' + sheet + '.p'
    if not os.path.isdir(TempPath):
        os.mkdir(TempPath)
    if not os.path.isfile(filepath_pandas):
       time_pd = 0
    else:
        time_pd = os.path.getmtime(filepath_pandas)
    if os.path.getmtime(filename) > time_pd:
        data = pd.read_excel(filename,sheet,header=header,skiprows=skiprows,skip_footer=skip_footer,index_col=index_col,parse_col=parse_col,parse_dates=parse_dates)
        data.to_pickle(filepath_pandas)
    else:
        data= pd.read_pickle(filepath_pandas)
    return data


  
    
def I_tilt_year(Lat,Lon):
    '''
    DNI Modifier to account for latitutude tilt
    Derive factor for latitude tilted PV panel based on I_tilt to DNI ratio from NASA SSE dataset
    Handbook of Photovoltaic Science and Engineering, Hegedus/Luque 2011  pg 801
    This function loads the data_I_tilt.xlsx file and interpolates the data for the given month, latitude and longitude
    :param Lat: Latitude
    :param Lon: Longitude
    :param NumMonth: Month number (from 1 to 12)
    :return: Latitude factor (DNI multiplier)
    '''
    import pandas as pd
    from scipy import interpolate
    
    data = pd.read_excel('Data/data_I_tilt.xlsx','YEAR')
    
    f = interpolate.interp2d(data.Lat, data.Lon, data.Val, kind='linear')
    
    return f(Lat, Lon)[0]
   
   
   
def CSP_power(CSP_A,I_b,lowlightcutoff,T_HTF_col_su,T_amb):
    '''
    Calculation of the CSP generation based on an efficiency law (fitted from a detailed EES model)
    :param CSP_A: Area of the CSP system, m^2
    :param I_b: Beam Irradiation, W/m^2
    :param lowlightcutoff: Minimum irradiation for power generation, W/m^2
    :param T_HTF_col_su: Inlet temperature of the heat transfer fluid in the collectors, °C
    :param T_amb: Ambient temperature, °C
    :return: Power generation, kW
    '''
    if I_b>lowlightcutoff:
        #If CSP in equipment
        if CSP_A > 0:
            # CSP Efficiency based on SORCE2.1_collector_eckerd.EES (Table 3) exercised across a range of T_amb, T_HTF_su and I_b parameters (constant vwind), functionalized in Eureqa R2=0.995
            PTC_ETA = 0.734327868107202 + (6.93950871301196 + 0.252496799296039*T_amb - 0.320896097990017* T_HTF_col_su )/I_b
            CSP_Avail_Pow=CSP_A*I_b*PTC_ETA/1000    # CSP kW output, simplified CSP efficiency based on TES_SOC assuming range 150-200C
        else:
            CSP_Avail_Pow=0    
    else:
        CSP_Avail_Pow=0    
    return CSP_Avail_Pow
    
    
    
def PV_power(PVkW,NOCT,Pmax_Tco,I_b,lowlightcutoff,T_amb):
    '''
    Calculation of the PV generation based on the effect of temperature on the PV cell efficiency.
    From Handbook of Photovoltaic Science and Engineering, Hegedus/Luque 2011 pg 801 with modified DNI to account for latitutude tilt
    :param PVkW: Peak power , kW
    :param NOCT: Nominal Operating Cell Temperature of the PV panel, °C
    :param Pmax_Tco: Temperature coefficient of rated power, %/C
    :param I_b: Beam Irradiation, W/m^2
    :param lowlightcutoff: Minimum irradiation for power generation, W/m^2
    :param T_amb: Ambient temperature, °C
    :return: Power generation, kW
    '''
    if I_b>lowlightcutoff:
        if PVkW > 0:
            # Determine the effect of temperature on the PV cell efficiency
            T_cell=T_amb+(NOCT-20)*I_b/800       # T_Cell in relation to NOCT, irradiance and T_amb
            P_nom=(100+(T_cell-25)*Pmax_Tco)/100            # Nominal power
            PV_Avail_Pow=PVkW*P_nom*I_b/1000     # cale with I_b and adjust for lat tilt to DNI ratio, and include P_norm to account for temperature coefficient
        else:
            PV_Avail_Pow=0
    else:
        PV_Avail_Pow=0    
    return PV_Avail_Pow
    
    
    
def batt_level(Batt_SOC,low_trip,high_trip):
    '''
    Calculation of batteries storage level. This is classified as
    'empty': battery is too low to supply load, 
    'part charge': enough energy in batteries to supply load,
    'full': batteries are full.
    :param BattSOC: Battery state of charge, kWh
    :param low_trip: 5% of battery full capacity, kWh
    :param high_trip: 95% of battery full capacity, kWh
    :return: Batteries storage level ('empty','part charge','full' - low to high level)
    '''
    if Batt_SOC > high_trip:
        batt_state = 'full'
    elif Batt_SOC < low_trip:
        batt_state = 'empty'
    else:
        batt_state = 'part charge'
    return batt_state



def TES_level(TES_SOC,lowT_trip,mediumT_trip,highT_trip):
    '''
    Calculation of TES storage level. This is classified as 
    'too low' - TES temperature is too low and cannot be used, 
    'low' - TES is charged but not full and cannot run the ORC, 
    'high' - TES is charged but not full and can run the ORC, 
    'too high' - TES is full (temperature reached its maximum)
    :param TES_SOC: TES state of charge
    :param lowT_trip: Limit below which TES level is considered too low to serve heating demand, kWh
    :param mediumT_trip: Limit below which TES level is considered too low to run the ORC (approximative amount of energy withdrawn in one time step if the ORC is ON - set for minimum possible ORC efficiency of 10%), kWh
    :param highT_trip: Limit above which TES is considered full and cannot be charget beyond it (95 % of TES full capacity), kWh
    :return: TES storage level (from 1 to 3 - low to high)
    '''
    if TES_SOC >= highT_trip:
        TES_state = 'too high'
    elif TES_SOC <= lowT_trip:
        TES_state = 'too low'
    else:
        if TES_SOC > mediumT_trip:
            TES_state = 'high'
        else:
            TES_state = 'low'
    return TES_state


def generate_results_table(P_load,P_genset,P_genset_ToBat,P_PV,P_PV_ToBat,P_ORC,P_ORC_ToBat,
                           P_bat_inv,P_bat,P_dumped,P_CSP,SOC_bat,SOC_TES,dayHour,E_loadleft,
                           State,beam,T_ambient,ORC_Eta,ChargeTES,TES_loss,Batt_frac,Genset_fuel,
                           Fuel_kW,P_load_th,P_heating,P_heating_backup,P_dumped_th,m_burner,m_genset):
    import pandas as pd    
    up = pd.DataFrame()
    up['P_load'] = P_load
    up['P_genset'] = P_genset
    up['P_genset_ToBat'] = P_genset_ToBat
    up['P_PV'] = P_PV
    up['P_PV_ToBat'] = P_PV_ToBat
    up['P_ORC'] = P_ORC
    up['P_ORC_ToBat'] = P_ORC_ToBat
    up['P_bat_inv'] = P_bat_inv
    up['P_bat'] = P_bat
    up['P_dumped'] = P_dumped
    up['P_CSP'] = P_CSP
    up['SOC_bat'] = SOC_bat
    up['SOC_TES'] = SOC_TES
    up['dayHour'] = dayHour
    up['E_loadleft'] = E_loadleft
    up['State'] = State
    up['beam'] = beam
    up['T_ambient'] = T_ambient
    up['ORC_Eta'] = ORC_Eta
    up['ChargeTES'] = ChargeTES
    up['TES_loss'] = TES_loss
    up['Batt_frac'] = Batt_frac
    up['Genset_fuel'] = Genset_fuel
    up['Fuel_kW'] = Fuel_kW
    up['P_load_th'] = P_load_th
    up['P_heating'] = P_heating
    up['P_heating_backup'] = P_heating_backup
    up['P_dumped_th'] = P_dumped_th
    up['m_burner'] = m_burner
    up['m_genset'] = m_genset
    return up    


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