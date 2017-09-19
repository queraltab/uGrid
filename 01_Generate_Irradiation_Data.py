# -*- coding: utf-8 -*-
"""
This scripts computes the direct normal radiation (DNI) incident on a north-south tracking collector and the incident 
total radiation on a titled surface. 

The input data is the typical meteorological year (TMY) from EnergyPlus, read as an epw file.

The script makes use of the following libraries: 
- pyepw for to ready the typical meteorological year data
- pvlib to compute the incident irradiation

"""

# Loading libraries:
import numpy as np
from pyepw.epw import EPW
import datetime
import pandas as pd
from pvlib import irradiance
from pvlib import solarposition
from pvlib.location import Location
from pvlib import tracking
from ugrid_tools_ThEl import scale_series


# Loading the Typical Meteorological year" for Johannesburg 
datafile="Data/ZAF_Johannesburg.683680_IWEC.epw"
# Specify coordinates:
latitude = -29
longitude = 28
# Define the orientation of the PV panel (28° is the optimum for south africa):
surface_tilt = 28
surface_azimuth = 0

# Read the TMY file and write the relevant results in numpy arrays:
epw = EPW()
epw.read(datafile)

N = len(epw.weatherdata)
I_hor = np.zeros(N)
I_d_hor = np.zeros(N)
I_0 = np.zeros(N)
DNI = np.zeros(N)
albedo = np.zeros(N)
T = np.zeros(N)

for i in range(N):
    wd = epw.weatherdata[i]
    print wd.year, wd.month, wd.day, wd.hour, wd.minute, wd.dry_bulb_temperature
    I_hor[i] = wd.global_horizontal_radiation
    I_d_hor[i] = wd.diffuse_horizontal_radiation
    DNI[i] = wd.direct_normal_radiation
    albedo[i] = wd.albedo
    T[i] = wd.dry_bulb_temperature
    I_0[i] = wd.extraterrestrial_direct_normal_radiation

# use the pvlib library to define the characteristics of the location:
tus = Location(latitude, longitude, 'Africa/Johannesburg', 700, 'Johannesburg')

# Define the pandas index (arbitrarily set to 2015):
first = epw.weatherdata[0]
last = epw.weatherdata[-1]
times = pd.date_range(start=datetime.datetime(2015,first.month,first.day),
                      end=datetime.datetime(2015,last.month,last.day,23,59,59), freq='h')

# Compute the sun position at each hour of the year:
ephem_data = solarposition.get_solarposition(times,
                                             latitude=latitude,
                                             longitude=longitude,
                                             method='nrel_numpy')

# Compute the diffuse irradiance on the panel, reflected from the ground:
S_d_reflect = irradiance.grounddiffuse(surface_tilt, I_hor, albedo=0.25, surface_type=None)

# Compute the diffuse irradiance on the panel, from the sky:
S_d_sky = irradiance.klucher(surface_tilt, surface_azimuth, I_d_hor, I_hor, ephem_data['zenith'], ephem_data['azimuth'])

# Compute the angles between the panel and the sun:
aoi = irradiance.aoi(surface_tilt, surface_azimuth, ephem_data['zenith'], ephem_data['azimuth'])

# Compute the global irradiance on the panel:
S = irradiance.globalinplane(aoi, DNI, S_d_sky, S_d_reflect)

# Second case: with tracking (axis is supposed to be north-south):

S_track = tracking.singleaxis(ephem_data['apparent_zenith'], ephem_data['azimuth'], axis_tilt=0, axis_azimuth=0, max_angle=360, backtrack=True)

S['Direct with tracking'] = DNI * np.cos(np.radians(S_track.aoi))

S = S.fillna(0)
S.plot()
#S.to_excel('data/incidentPV.xlsx')

# Adding temperature to the dataframe:
S['T'] = T

# print S with w 15min time resolution:
index15 = pd.DatetimeIndex(start='2015-01-01 00:00',end='2015-12-31 23:59:00',freq='15min')
S15 = scale_series(S,index15)
del S15[0]
S15.to_excel('Data/IncidentRadiation.xlsx')

print 'Total incident radiation without tracking: ' + str(S['poa_global'].sum()/1E3) + ' kWh/m²'
print 'Total incident direct radiation with tracking: ' + str(S['Direct with tracking'].sum()/1E3) + ' kWh/m²'




