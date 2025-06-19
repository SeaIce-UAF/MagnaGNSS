# -*- coding: utf-8 -*-
# pylint: skip-file
# -1: [file-ignored]
"""
Created on Tue Dec 17 11:27:52 2024

@author: acapelli
"""




path=r'C:\Users\AcCap\MagnaGNSS\examples\withSyncEvents' # Achille's path

filename_cont="mg1_merged_339_340.pos" 
filename_events="mg1_3390_3400_ExtEvent1.txt"
file_magna="DataMagnaprobe.dat"


file_save='MagnaMergedA.csv'
file_saveB='MagnaMergedB.csv'

DeltaS=16.62
tolerance=0.5
window=1


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


plot_heigts_time(events )


# plot heights
leg1=events.query('time> "2024-12-04 22:04" and time < "2024-12-04 23:03" ') # the times in the strimg are used to cut out the data of the single transects
plot_heigts(leg1 ) # function in mergeMagna.py # the times in the strimg are used to cut out the data of the single transects
plot_heigts_quality(leg1, title='SW')



# plot heights
leg2=events.query('time> "2024-12-04 23:15"  and time < "2024-12-05 00:34"  ') # the times in the strimg are used to cut out the data of the single transects
plot_heigts(leg2 ) # function in mergeMagna.py # the times in the strimg are used to cut out the data of the single transects
plot_heigts_quality(leg2, title='SE')


# plot heights
leg3=events.query('time> "2024-12-05 00:40"  and time < "2024-12-05 00:55" ') # the times in the strimg are used to cut out the data of the single transects
plot_heigts(leg3 ) # function in mergeMagna.py # the times in the strimg are used to cut out the data of the single transects
plot_heigts_quality(leg3, title='S')

