# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:54:33 2024

Modified: Dec 17 2024


@author: AcCap
"""
import sys
import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from datetime import timedelta
import math 

try:
    from geopy import distance
except ModuleNotFoundError as e:
    print('You need to install geopy with "pip install geopy" ')   
    raise e

# %% functions
def read_events(filename_events,path):
    
    eventsB=pd.DataFrame([])
    if filename_events.split('.')[-1]=='txt':
        events=pd.read_csv(path+'\\'+filename_events,header=5, sep=',',
                           names=['TOW','WNc','Source','Polarity','Offset' ,'Offset_TOW','RxClkBias' ,'PVTAge'])
        

        events['time']=[gps_to_datetime(events.WNc[i], events.TOW[i]) for i in range(len(events))]
        
        eventsB=events.query("Source=='EventB' ")
        events=events.query("Source=='EventA' ")
        
        events.reset_index(drop=True, inplace=True)
        eventsB.reset_index(drop=True, inplace=True)
        
    elif filename_events.split('.')[-1]=='pos':
        events=pd.read_csv(path+'\\'+filename_events,header=9, sep='\s+', 
                           names=['date', 'GPST', 'lat', 'lon', 'height', 'Q', 'ns',
                                  'sdn', 'sde', 'sdu', 'sdne', 'sdeu', 'sdun', 'age', 'ratio',])
        events['time']=pd.to_datetime(events['date']+' '+events['GPST'])
    else:
        raise ValueError("Can't open events file!")
    
    return events,eventsB

    

def gps_to_datetime(gps_week, time_of_week):
    # GPS epoch start date (January 6, 1980)
    gps_epoch = pd.Timestamp("1980-01-06")
    
    # Calculate the total time in weeks and seconds
    gps_time = gps_epoch + timedelta(weeks=float(gps_week), seconds=time_of_week)
    
    return gps_time   

# get deviation from reference data
def RMSE_Ref(ref, events,data=[],d_radious=0.25):
     # radious to look for nearby data points in meters 
    dist=[]
    d_h=[]
    index2=[]
    dist_ref=[]
    d_h_ref=[]
    
    for i in ref.index:
        index=[]
        dist2=[]
        d_h2=[] 
        for j in events.index:
            d=distance.distance([ref.loc[i,'lat'],ref.loc[i,'lon']],[events.loc[j,'lat'],events.loc[j,'lon']]).m
    
            if  d< d_radious:
                
                if len(data)>0:
                    # go through data and find same timestamp
                    for k in data.index:
                        if events.loc[j,'time']==data.loc[k,'time']:
                            index.append(k)
                            index2.append(k)
                            d=distance.distance([ref.loc[i,'lat'],ref.loc[i,'lon']],[data.loc[k,'lat'],data.loc[k,'lon']]).m
                            dist.append(d)
                            dist2.append(d)
                            d_h.append(ref.loc[i,'height']-data.loc[k,'height'])
                            d_h2.append(ref.loc[i,'height']-data.loc[k,'height'])
                    
                else:
                    index.append(j)
                    index2.append(j)
                    dist.append(d)
                    dist2.append(d)
                    d_h.append(ref.loc[i,'height']-events.loc[j,'height'])
                    d_h2.append(ref.loc[i,'height']-events.loc[j,'height'])
            
        dist_ref.append(np.mean(dist2))
        d_h_ref.append(np.mean(d_h2))
    
    
    d_h=np.array(d_h)
    dist=np.array(dist)
    
    Delta_hor=np.mean(np.abs(dist))
    Delta_h=np.mean(np.abs(d_h))
    Delta_hor_std=np.std(dist)
    Delta_h_std=np.std(np.abs(d_h))
    RMSE_h=np.sqrt(np.sum(d_h**2)/len(d_h))
    RMSE_x=np.sqrt(np.sum(dist**2)/len(dist))
    print('Delta_h= {:.3f} +-{:.4f} m'.format(Delta_h,Delta_h_std))
    print('Delta_x_hor={:.3f} +-{:.4f} m'.format(Delta_hor,Delta_hor_std))
    print('RMSE(h)= {:.3f} m'.format(RMSE_h))
    print('RMSE(x)= {:.3f} m'.format(RMSE_x))

    return  RMSE_h,RMSE_x,index2,d_h,dist,dist_ref,d_h_ref



# filter events too near
def filterDoubleClick(events,dt_min=0.3):
    index2=[events.index[0]]
    n=0
    dts=[]
    
    print('\nTrowing out double click events\n------------------------------------')
    for i in events.index[1:]:
        
        if (events.loc[i,'time']- events.loc[i-1,'time'])<pd.Timedelta(seconds=dt_min):
            index2.append(i)
            print(events.loc[i-1,'time'], ': ',(events.loc[i,'time']- events.loc[i-1,'time']))
            n+=1
            dts.append((events.loc[i,'time']- events.loc[i-1,'time']).total_seconds())
            
    index2=np.array(index2)
    events=events.drop(index=index2)
    events.set_index(pd.RangeIndex(0,len(events)),inplace=True)
    
    print('Mean dt double points: {:.3f} s'.format(np.mean(dts)))
    print('Std dt double points: {:.3f} s'.format(np.std(dts)))
    
    return events

# get magnaprobe data
def add_magna(events,magna,dt_magnaprobe=0,tolerance=0.1):
    ind=magna.copy()
    if isinstance(dt_magnaprobe,pd.Timedelta):
        ind['time']+dt_magnaprobe
    else:
        ind['time']=ind['time']+pd.Timedelta(seconds=dt_magnaprobe)
    ind.set_index('time',inplace=True)
    ii=ind.index.get_indexer(events.time,method='nearest',tolerance=pd.Timedelta(seconds=tolerance))
    
    ii2=ii.copy()
    ii2[ii==-1]=0  # set events not found as first value. They will be set to nan later
    
    events.loc[:,['Counter', 'DepthCm', 'BattVolts']]=magna.loc[ii2,['Counter', 'DepthCm', 'BattVolts']].values
    
    # set magnaprobe values to nan for evebts not found in Magnaprobe files
    events.loc[ii==-1,['Counter', 'DepthCm', 'BattVolts']]=np.nan
    
    return events


# get mean of continous data
def cont_to_events(cont,events,dt_mean=1):
    """
    

    Parameters
    ----------
    cont : TYPE
        DESCRIPTION.
    events : TYPE
        DESCRIPTION.
    dt_mean : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    lines=[]
    ind=[]
    for i in events.index:
        t1=events.loc[i,'time']-pd.Timedelta(seconds=dt_mean)
        t2=events.loc[i,'time']+pd.Timedelta(seconds=dt_mean)
        a=cont.loc[(cont['time'] >= t1) & (cont['time'] <= t2),['lat', 'lon', 'height', 'Q', 'ns', 'sdn', 'sde', 'sdu',
               'sdne', 'sdeu', 'sdun', 'age', 'ratio', 'time']]
        line=a.mean()
        line['lat_std']=a.loc[:,'lat'].std()
        line['lon_std']=a.loc[:,'lon'].std()
        line['height_std']=a.loc[:,'height'].std()
        
        
        line['time']=events.loc[i,'time']
        
        if len(a)>0:   # trow away events if there are no corresponding countinous events!
            lines.append(line)
            ind.append(i)
    
    return pd.DataFrame(lines,index=ind)





def plot_heigts_time(events2 ):
    """
    Plot snow depth and ice and snow surface elevations  as a function of time

    Parameters
    ----------
    events2 : pandas
        Merged data.


    """
    
    fig,[ax1,ax2]=pl.subplots(2,1,sharex=True,figsize=(8,6))
    ax1.plot(events2.time,events2.height-events2.height.mean()+events2.DepthCm/100,'x:',label='snow surface')
    ax1.plot(events2.time,events2.height-events2.height.mean(),'x:',label='ice surface')
    ax1.set_ylabel('h (m)')
    # ax1.set_xlabel('time')
    ax1.legend()
    
    ax2.plot(events2.time,events2.DepthCm/100,'x:',label='snow depth')
    ax2.set_ylabel('snow depth (m)')
    ax2.set_xlabel('time')
    ax2.legend()



def plot_heigts(events2 ):
    """
    Plot snow depth and ice and snow surface elevations as a function of distance form the first datapoint

    Parameters
    ----------
    events2 : pandas
        Merged data.


    """


    lat0=events2['lat'].iloc[0]
    lon0=events2['lon'].iloc[0]
    d=[]
    for i in events2.index:
        d.append(distance.distance((lat0,lon0),(events2.loc[i,'lat'],events2.loc[i,'lon'])).m)
        
    d=np.array(d)
    
    fig,[ax1,ax2]=pl.subplots(2,1,sharex=True,figsize=(8,6))
    ax1.plot(d,events2.height-events2.height.mean()+events2.DepthCm/100,'x:',label='snow surface')
    ax1.plot(d,events2.height-events2.height.mean(),'x:',label='ice surface')
    ax1.set_ylabel('h (m)')
    ax1.set_xlabel('distance (m)')
    ax1.legend()
    
    ax2.plot(d,events2.DepthCm/100,'x:',label='snow depth')
    ax2.set_ylabel('snow depth (m)')
    ax2.set_xlabel('distance (m)')
    ax2.legend()




def plot_heigts_quality(events, title='',xaxis='distance' ):
    """
    Plot snow depth and ice and snow surface elevations as a function of distance from the first datapoint or time.
    
    

    Parameters
    ----------
    events: pandas
        Merged data.
    xaxis: string if ='distance' (default) use distance as xaxis. otherwise time.

    """


    lat0=events['lat'].iloc[0]
    lon0=events['lon'].iloc[0]
    
    if xaxis='distance':
        d=[]
        for i in events.index:
            d.append(distance.distance((lat0,lon0),(events.loc[i,'lat'],events.loc[i,'lon'])).m)
        d=np.array(d)
        x_label='distance (m)'
    else:
        d=events.time
        x_label='time'
    
    fig,[ax1,ax2,ax3,ax4]=pl.subplots(4,1,sharex=True,figsize=(8,12))
    ax1.plot(d,events.height-events.height.mean()+events.DepthCm/100,'x:',label='snow surface')
    ax1.plot(d,events.height-events.height.mean(),'x:',label='ice surface')
    ax1.set_ylabel('h (m)')
    # ax1.set_xlabel('distance (m)')
    ax1.legend()
    
    ax2.plot(d,events.DepthCm/100,'x:',label='snow depth')
    ax2.set_ylabel('snow depth (m)')
    # ax2.set_xlabel('distance (m)')
    ax2.legend()

    ax3.plot(d,events.sdn,'x:',label='Std North')
    ax3.plot(d,events.sde,'x:',label='Std East')
    ax3.plot(d,events.sdu,'x:',label='Std Up')
    ax3.set_ylabel('mean std GPS solution (m)')
    # ax3.set_xlabel('distance (m)')
    ax3.legend()
    
    ax4.plot(d,events.Q,'x:',label='Fix quality')
    ax4.plot([],[],'x:g',label='AR ratio')
    ax4b=ax4.twinx()
    ax4b.plot(d,events.ratio,'x:g')
    ax4b.set_ylabel('AR ratio ()')
    ax4.set_ylabel('Fix type ()')
    ax4.set_xlabel(x_label)
    ax4.legend(loc=1)

    # Set custom y-axis ticks and labels
    y_ticks = np.linspace(-1, 1, 5)  # Set 5 ticks between -1 and 1
    y_labels = ['Low', 'Below Zero', 'Zero', 'Above Zero', 'High']  # Custom labels for the ticks

    # Apply the custom ticks and labels
    ax4.set_yticks(ticks=(1,2,3), labels=['fix','float','3d'])


    if len(title)>0:
        ax1.set_title(title)




def get_XY(events,origin=0):  
    """
    Transform lat,lon to local coordinates around origin point and add x,y data to events dataframe

    Parameters
    ----------
    events : pandas dataframe with Magnaprobe data
        events data with lat lon colums.
    origin : TYPE
        if int: index of coordinate to se as origin of local coordinate system
        if: [lat,lon] coordinates of origin of local reference system.

    Returns
    -------
    x : local x-coordinate in meters (west-east)
    y : local y-coordinate in meters (south- north)

    """
     
    
    origin=np.array(origin)
    if len(origin.shape)==0:  
        p0=[events.lat[origin],events.lon[origin]]
    elif len(origin.shape)==1:  
            p0=origin    
    
    d=[]
    bearing=[]
    for lat,lon in zip(events.lat,events.lon):
        d.append(distance.distance(p0,[lat,lon]).m)
        bearing.append(get_bearing(p0,[lat,lon]))
    
    d=np.array(d)
    bearing=np.array(bearing)
    
    y=d*np.cos(bearing/360*2*np.pi)
    x=d*np.sin(bearing/360*2*np.pi)
    
    events['x']=x
    events['y']=y
    
    return x,y



def get_bearing(start, end):
    # if (type(start) != tuple) or (type(end) != tuple):
    #     raise TypeError("Only tuples are supported as arguments")

    lat1 = math.radians(start[0])
    lat2 = math.radians(end[0])

    diffLong = math.radians(end[1] - start[1])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
                                           * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(x, y)

    # Now we have the initial bearing but math.atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


def get_d(events,loc0):
    
    lat0=loc0[0]
    lon0=loc0[1]   
    d=[]
    for i in events.index:
        d.append(distance.distance((lat0,lon0),(events.loc[i,'lat'],events.loc[i,'lon'])).m)
    events['d']=np.array(d)
    
    
    

def merge(filename_events, filename_cont, file_magna,
         file_save='PPK_Magna.csv', file_saveB='PPK_Magna_B.csv',
         path='\\', plot=True,
         dt_min=0.1, window=1, tolerance=1, DeltaH=8, DeltaS=16.6,
         useEvents=True):
    """ 
    Parameters
    ----------
    filename_events : TYPE
        Event data GNSS.   .txt from RXTools or .pos form Emlid studio
    filename_cont : TYPE
        Continus data GNSS.  .pos form Emlid studio
    file_magna : TYPE
        File data magnaprobe.
    file_save : TYPE, optional
        Filename for saving data. The default is 'PPK_Magna.csv'.
    file_saveB : TYPE, optional
        Filename for saving events from button B. The default is 'PPK_Magna_B.csv'.
    path: string, optional
        path of files. The default is '\\'.
    plot : TYPE, optional
        Produce plot to check sync. The default is True.
    dt_min : TYPE, optional
        filter events nearer then dt_min. The default is 0.1.
    window : TYPE, optional
        Window size for averaging continous data around event in 's'. The default is 0.5.
    tolerance : TYPE, optional
        tolerance in 's'for syncing Magnaprobe and GNSS. The default is 0.1.
    DeltaH : TYPE, optional
        Time shift Magnaprobe in hours. GPS time to local time. The default is 8.
    DeltaS : TYPE, optional
        Time shift between Magnaprobe and GNSS in seconds (leap seconds+ error). The default is 16.62
    useEvents: bool, optional
            Default true. If false the GNSS Timestamps are ignored and the data are synced based on the Magnaprobe time. 
            To be used in case the conenction between the GNSS and Magnaprobe didn't work. 
                
    Returns
    -------
    None.
    
    """
    
    dt_min=float(dt_min) 
    dt_mean=float(window)/2
    tolerance=float(tolerance) 
    DeltaH=float(DeltaH)
    DeltaS=float(DeltaS)
    eventsA_empty=False
    
    
    # read files
    events,eventsB=read_events(filename_events,path)
            
    cont=pd.read_csv(path+'\\'+filename_cont,header=9, sep='\s+', 
                       names=['date', 'GPST', 'lat', 'lon', 'height', 'Q', 'ns',
                              'sdn', 'sde', 'sdu', 'sdne', 'sdeu', 'sdun', 'age', 'ratio'])

    cont['time']=pd.to_datetime(cont['date']+' '+cont['GPST'])
    
    magna=pd.read_csv(path+'\\'+file_magna,skiprows=4,
                     names=["TIMESTAMP","RECORD","Counter","DepthCm","BattVolts",
                            "latitude_a","latitude_b","Longitude_a","Longitude_b","fix_quality",
                            "nmbr_satellites","HDOP","altitudeB","DepthVolts","LatitudeDDDDD","LongitudeDDDDD",
                            "month","dayofmonth","hourofday","minutes","seconds","microseconds"],
                     )
    try:
        magna['time']=pd.to_datetime(magna['TIMESTAMP'],format='ISO8601')
    except ValueError:
        magna['time']=pd.to_datetime(magna['TIMESTAMP'])
    
    
    # time difference
    magna['time']=magna['time']+pd.Timedelta(hours=DeltaH, minutes=00, seconds=DeltaS)
    
    
    if len(events)==0 or not useEvents:
        events=pd.DataFrame(magna['time'].values,columns=['time'])
        eventsA_empty=True
        print('eventsA is empty!! Using Magnaprobe timestamp instead!!')
    
    if plot:
        pl.figure()
        pl.plot(cont.time,cont.height,'x',label='GPSRTK GNSS continous')
        
        pl.scatter(magna.time,magna.altitudeB,marker='o',edgecolors='g',facecolor='None',label='Magna GPS')
        ylim=pl.gca().get_ylim()
        
        if eventsA_empty:
            for t in events.time:
                pl.plot([t,t],ylim,'--g')
                ax=pl.gca()
                pl.text(0.05,0.95,'eventsA is empty!! Using MAgnaprobe timestamp instead!!',
                        transform=ax.transAxes,        
                        ha='left', va='top', fontsize=9, color='red')
        else: 
            for t in events.time:
                pl.plot([t,t],ylim,'--k')
        for t in eventsB.time:
               pl.plot([t,t],ylim,'--r')
               
        pl.plot([t,t],ylim,'--k',label='EventsA')  
        pl.plot([t,t],ylim,'--r',label='EventsB')   
        if eventsA_empty:
            pl.plot([],[],'--g',label='GPS Magnaprobe timestamps')        
        pl.legend()
        pl.xlabel('time')
        pl.ylabel('h (m ASL)')   
        pl.title('Merge quality check')


    # filter events too near
    events2=filterDoubleClick(events,dt_min=dt_min)
    
    # get mean of continous data
    events2=cont_to_events(cont,events2,dt_mean=dt_mean)
    
    # get mean of continous data
    eventsB2=cont_to_events(cont,eventsB,dt_mean=dt_mean)
    
    
    # merge with magnaprobe data
    events2=add_magna(events2,magna,dt_magnaprobe=0,tolerance=tolerance)


    # round numbers
    keys=['lat', 'lon', 'height', 'Q',  'sdn', 'sde', 'sdu', 'sdne', 'sdeu',
           'sdun',  'lat_std', 'lon_std', 'height_std','age','ratio','ns']
    decimals = pd.Series([9,9,3,0,4,4,4,4,4,4,11,11,4,2,2,0], 
                         index=keys)
    
    #eventsA
    events2=events2.round(decimals)
    events2.Counter=np.int32(events2.Counter)
    events2.Q=np.int32(events2.Q)
    events2.ns=np.int32(events2.ns)
    #save files
    events2.to_csv(path+'\\'+file_save,index=False,
                   columns=['time','Counter','lat', 'lon', 'height','DepthCm',  
                            'sdn', 'sde', 'sdu','sdne', 'sdeu', 'sdun', 'age', 'ratio', 
                            'lat_std', 'lon_std', 'height_std','Q','ns',
                             'BattVolts'])
    
    # eventsB
    if len(eventsB2)>0:
        eventsB2=eventsB2.round(decimals)
        eventsB2.Q=np.int32(eventsB2.Q)
        eventsB2.ns=np.int32(eventsB2.ns)
        eventsB2.to_csv(path+'\\'+file_saveB,index=False,
                       columns=['time','lat', 'lon', 'height',  
                                'sdn', 'sde', 'sdu','sdne', 'sdeu', 'sdun', 'age', 'ratio', 
                                'lat_std', 'lon_std', 'height_std','Q','ns'])
    

    
    return events2,eventsB2
    
#%% Main


# if __name__ == "__main__":
       
        
#     if len(sys.argv) > 3:
#         # Print the arguments passed in from the command line
#         # print("Arguments passed in from the command line:")
#         kwargs = dict(arg.split('=') for arg in sys.argv[1:] if '=' in arg)
#         # Print the keyword arguments
#         if kwargs:
           
#             # print("Keyword arguments passed in from the command line:")
#             for key, value in kwargs.items():
#                 print(f"{key}: {value}")
#             merge(sys.argv[1],sys.argv[2],sys.argv[3], **kwargs)
#         else:
#             merge(sys.argv[1],sys.argv[2],sys.argv[3],)
#     elif len(sys.argv) == 3:
#         merge(sys.argv[1],sys.argv[2],sys.argv[3],)
    
        
#     else:
#         print("""Missing arguments!! Syntax: mergeMagna.py file_events file_cont file_Magna  plot=True, dt_min=0.1, dt_mean=0.5, tolerance=0.1, DeltaH=8, DeltaS=17 
#               Parameters
#               ----------
#               filename_events : TYPE
#                   Event data GNSS.   .txt from RXTools or .pos form Emlid studio
#               filename_cont : TYPE
#                   Continus data GNSS.  .pos form Emlid studio
#               file_magna : TYPE
#                   File data magnaprobe.
#               file_save : TYPE, optional
#                   Filename for saving data. The default is 'PPK_Magna.csv'.
#               file_saveB : TYPE, optional
#                   Filename for saving events from button B. The default is 'PPK_Magna_B.csv'.
#               path: string, optional
#                   path of files. The default is '\\'.
#               plot : TYPE, optional
#                   Produce plot to check sync. The default is True.
#               dt_min : TYPE, optional
#                   filter events nearer then dt_min. The default is 0.1.
#               window : TYPE, optional
#                   Window size for averaging continous data around event in 's'. The default is 0.5.
#               tolerance : TYPE, optional
#                   tolerance in 's'for syncing Magnaprobe and GNSS. The default is 0.1.
#               DeltaH : TYPE, optional
#                   Time shift Magnaprobe in hours. GPS time to local time. The default is 8.
#               DeltaS : TYPE, optional
#                   Time shift Magnaprobe in seconds (leap seconds+ error). The default is 16.62
#               """ )
