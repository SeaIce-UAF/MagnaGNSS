# -*- coding: utf-8 -*-
# pylint: skip-file
# -1: [file-ignored]
"""
Created on Mon Jan 13 15:37:42 2025

@author: apinzner
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:48:25 2025

@author: apinzner
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:27:52 2024

@author: acapelli
"""


path=r'C:\Users\AcCap\Google Drive\myISOSP2\OTI_ISOPS2\Magnaprobe\data\250112_Parkinglot\mg2_0130.25' # Achilles path. Anika just comment this if you use the script



filename_cont="Mg2_0130.pos"
filename_events="mg2_0130.25__SBF_ExtEvent1.txt"
file_magna="GEODEL_0.dat"


file_save='MagnaMergedA_mg2_0130.25.csv'
file_saveB='MagnaMergedB_mg2_0130.25.csv'

DeltaS=16.605
tolerance=0.5
window=0.8


events,eventsB=merge(filename_events,filename_cont, file_magna, file_save=file_save, file_saveB=file_saveB, 
      path=path, plot=True, 
      dt_min=0.5, window=window, tolerance=tolerance, DeltaH=8, DeltaS=DeltaS)




#%% Plots


# latlon
pl.figure()
pl.plot(events.lon,events.lat,'-x')
pl.ylabel('lat')
pl.xlabel('lon')

pl.figure()
pl.plot(events.time,events.lon,'-x')
pl.xlabel('time')
pl.ylabel('lon')



# plot heights
leg1=events.query('time> "2025-01-13 03:03:23" and time < "2025-01-13 03:15:46" ')  # the times in the strimg are used to cut out the data of the single transects
plot_heigts(leg1 )  # function in mergeMagna.py

leg2=events.query('time> "2025-01-13 03:03:23" and time < "2025-01-13 03:15:46" ')
plot_heigts(leg2)