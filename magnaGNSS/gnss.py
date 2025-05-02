# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:54:33 2024

Modified: Dec 17 2024
Amended: May 1, 2025, by MegaVolts
@author: AcCap, MegaVolts

"""

from datetime import timedelta
import numpy as np
import matplotlib.pyplot as pl
import pandas as pd

try:
    from geopy import distance
except ModuleNotFoundError as e:
    print('You need to install geopy with "pip install geopy" ')
    raise e

# %% Default variable
F_PLOT = True
DT_MIN = 0.1
WINDOW = 1
TOL = 1
DELTA_H = 8
DELTA_S = 16.6
F_USEEVENT = True
F_PRESERVE_ORPHAN = True

# %% functions
def read_events(fp):
    """
    Parameters
    ----------
    fp : STRING
        Filepath to event or position data file(.dat)

    Returns
    -------
    event_df: pd.DataFrame()
        Dataframe containing the event list
    """
    if fp is None:
        event_df = pd.DataFrame()
    elif fp.split('.')[-1]=='txt':
        event_df = pd.read_csv(fp, header=5, sep=',', skip_blank_lines=True,
                               names=['TOW', 'WNc', 'Source', 'Polarity', 'Offset' , 'Offset_TOW', 'RxClkBias' ,
                                      'PVTAge'])
        event_df['time']=[gps_to_datetime(event_df.WNc[i], event_df.TOW[i]) for i in range(len(event_df))]

    elif fp.split('.')[-1]=='pos':
        event_df = pd.read_csv(fp, header=9, sep=r'\s+',
                           names=['date', 'GPST', 'lat', 'lon', 'height', 'Q', 'ns',
                                  'sdn', 'sde', 'sdu', 'sdne', 'sdeu', 'sdun', 'age', 'ratio',])
        event_df['time'] = pd.to_datetime(event_df['date'] + ' ' + event_df['GPST'])
        event_df['time'] = pd.to_datetime(event_df['date'] + ' ' + event_df['GPST'])
    else:
        raise ValueError("No event file to open")

    return event_df


def split_events(event_df):
    """
    :param event_df: pd.DataFrame()
        Dataframe containing events A and B from
    :return: eventA_df: pd.DataFrame()
        Dataframe containing the event A, corresponding to snow depth measurement with magnaprobe
    :return: eventB_df: pd.DataFrame()
        Dataframe containing the event B
    """
    if event_df.empty:
        event_Adf = pd.DataFrame()
        event_Bdf = pd.DataFrame()
    else:
        event_Adf=event_df.query("Source=='EventB' ")
        event_Bdf=event_df.query("Source=='EventA' ")

        event_Adf.reset_index(drop=True, inplace=True)
        event_Bdf.reset_index(drop=True, inplace=True)
    return event_Adf, event_Bdf


def gps_to_datetime(gps_week, time_of_week):
    # GPS epoch start date (January 6, 1980)
    gps_epoch = pd.Timestamp("1980-01-06")

    # Calculate the total time in weeks and seconds
    gps_time = gps_epoch + timedelta(weeks=float(gps_week), seconds=time_of_week)

    return gps_time


# get deviation from reference data
def RMSE_Ref(ref, events, data=None, d_radious=0.25):
    # radius to look for nearby data points in meters
    if data is None:
        data = []
    dist=[]
    d_h=[]
    index2=[]
    dist_ref=[]
    d_h_ref=[]

    for i in ref.index:
        index = []
        dist2 = []
        d_h2 = []
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

    delta_hor = np.mean(np.abs(dist))
    delta_h = np.mean(np.abs(d_h))
    delta_hor_std = np.std(dist)
    delta_h_std = np.std(np.abs(d_h))
    RMSE_h = np.sqrt(np.sum(d_h**2) / len(d_h))
    RMSE_x = np.sqrt(np.sum(dist**2) / len(dist))
    print(f'Delta_h= {delta_h:.3f} +-{delta_h_std:.4f} m')
    print(f'Delta_x_hor={delta_hor:.3f} +-{delta_hor_std:.4f} m')
    print(f'RMSE(h)= {RMSE_h:.3f} m')
    print(f'RMSE(x)= {RMSE_x:.3f} m')

    return  RMSE_h,RMSE_x,index2,d_h,dist,dist_ref,d_h_ref



# filter events too near
def filterDoubleClick(events, dt_min=DT_MIN):
    index2=[events.index[0]]
    n = 0
    dts = []

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

    if len(dts) == 0:
        print('No double click event found')
    else:
        print(f"Mean dt double points: {np.mean(dts):.3f} s")
        print(f"Std dt double points: {np.std(dts):.3f} s")

    return events


# get magnaprobe data
def add_magna(events, magna, dt_magnaprobe=0, tolerance=0.1, preserve_mg_data=False):
    ind = magna.copy()
    if isinstance(dt_magnaprobe,pd.Timedelta):
        ind['time']+dt_magnaprobe
    else:
        ind['time']=ind['time']+pd.Timedelta(seconds=dt_magnaprobe)
    ind.set_index('time',inplace=True)
    ii=ind.index.get_indexer(events.time,method='nearest',tolerance=pd.Timedelta(seconds=tolerance))

    ii2=ii.copy()
    ii2[ii==-1]=0  # set events not found as first value. They will be set to nan later

    # # OLD VERSION
    # events.loc[:,['Counter', 'DepthCm', 'BattVolts']]=magna.loc[ii2,['Counter', 'DepthCm', 'BattVolts']].values
    #
    # if not preserve_orphan:
    #     # set magnaprobe values to nan for evebts not found in Magnaprobe files
    #     events.loc[ii==-1,['Counter', 'DepthCm', 'BattVolts']]=np.nan
    #

    # MO: rather than adding Counter, DepthCm, BattVolts from magna to events, I suggest the other way around,
    # this allow to preserve all the magnaprobe data even if there is no GNSS data (or need to be added later
    # due to daily file rotation at midnight UTC.
    cols = ['lat', 'lon', 'height', 'Q', 'ns', 'sdn', 'sde', 'sdu', 'sdne', 'sdeu',
       'sdun', 'age', 'ratio', 'time', 'lat_std', 'lon_std', 'height_std']

    if not preserve_mg_data:
        magna = magna[['Counter', 'DepthCm', 'BattVolts', 'TIMESTAMP']]

    events = pd.concat([magna, events], axis=1)

    # This happens automatically with the concat
    #if not preserve_orphan:
    #    # set magnaprobe values to nan for events not found in Magnaprobe files
    #    events.loc[ii==-1,['Counter', 'DepthCm', 'BattVolts']]=np.nan

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



def merge(mg_fp, pos_fp, event_fp=None, f_plot=F_PLOT, fig_fp=None,
         dt_min=DT_MIN, window=WINDOW, tolerance=TOL, deltaH=DELTA_H, deltaS=DELTA_S, f_useEvents=F_USEEVENT):
    """ 
    Parameters
    ----------
    mg_fp : STRING
        Filepath to raw data file from magnaprobe (.dat)
    pos_fp : STRING
        Filepath to the continuous GNSS position data file form Emlid studio (.pos)
    event_fp : STRING, optional
        Filepath to GNSS event data file from RXTools (.txt) or form Emlid studio (.pos)
    f_plot : BOOL, optional
        Produce plot to check sync. The default is True.
    fig_fp : STRING, optional
        Filepath to figure output for the plot. Default is None.
    dt_min : FLOAT, optional
        Filter events nearer than dt_min. The default is 0.1 min.
    window : FLOAT, optional
        Duration in seconds of the period centered on an event used for averaging continuous data. The default is 0.5 s.
    tolerance : FLOAT, optional
        Tolerance in seconds for syncing GNSS time to maganprobe timestamp. The default is 0.1.
    deltaH : FLOAT, optional
        Time shift Magnaprobe in hours. GPS time to local time. The default is 8.
    deltaS : FLOAT, optional
        Time shift between Magnaprobe and GNSS in seconds (leap seconds+ error). The default is 16.62
    useEvents: bool, optional
            Default true. If false the GNSS Timestamps are ignored and the data are synced based on the Magnaprobe time. 
            To be used in case the conenction between the GNSS and Magnaprobe didn't work. 
                
    Returns
    -------
        None
    """
    
    dt_mean = float(window)/2
    
    # read files
    event_df = read_events(event_fp)
    event_Adf, event_Bdf = split_events(event_df)
    cont = pd.read_csv(pos_fp, header=9, sep=r'\s+',
                       names=['date', 'GPST', 'lat', 'lon', 'height', 'Q', 'ns', 'sdn', 'sde', 'sdu', 'sdne', 'sdeu', 'sdun', 'age', 'ratio'])

    cont['time'] = pd.to_datetime(cont['date']+' '+cont['GPST'])
    
    magna=pd.read_csv(mg_fp, skiprows=4,
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
    magna['time']=magna['time']+pd.Timedelta(hours=deltaH, minutes=00, seconds=deltaS)

    if len(events)==0 or not f_useEvents:
        events = pd.DataFrame(magna['time'].values,columns=['time'])
        eventsA_empty=True
        print('eventsA is empty!! Using Magnaprobe timestamp instead!!')

    if f_plot or fig_fp is not None:
        # Plot only for the overlapping time period between magnaprobe and gnss
        t_max = min(magna.time.max(), cont.time.max())
        t_min = max(magna.time.min(), cont.time.min())
        # Convert t_min in t_index
        t_start = magna.loc[magna.time <= t_min, 'time'].iloc[-1].value/1e9/60/60/24
        t_end = magna.loc[t_max <= magna.time, 'time'].iloc[0].value/1e9/60/60/24

        fig, axs = pl.subplots(2, 1, figsize=[8.5,11])
        for ii, ax in enumerate(axs):
            ax.plot(cont.time,cont.height,'x',label='GPSRTK GNSS continous')

            ax.scatter(magna.time,magna.altitudeB,marker='o',edgecolors='g',facecolor='None',label='Magna GPS')
            ylim=ax.get_ylim()

            if eventsA_empty:
                for t in events.time:
                    ax.plot([t,t],ylim,'--g')
                    if ii == 0:
                        ax.text(0.05,0.95,'eventsA is empty!! Using MAgnaprobe timestamp instead!!',
                                transform=ax.transAxes,
                                ha='left', va='top', fontsize=9, color='red')
            else:
                for t in events.time:
                    ax.plot([t,t],ylim,'--k')
            for t in eventsB.time:
                   ax.plot([t,t],ylim,'--r')

            ax.plot([t, t],ylim,'--k',label='EventsA')
            ax.plot([t, t],ylim,'--r',label='EventsB')
            if eventsA_empty:
                ax.plot([],[],'--g',label='GPS Magnaprobe timestamps')

            h_max = max([magna.altitudeB.max(), cont.height.max()])
            h_min = min([magna.altitudeB.min(), cont.height.min()])
            ax.set_ylim([h_min-0.2, h_max+0.2])

            if ii==1:
                index_min = magna.loc[magna.time <= t_min, 'time'].index[-1]
                index_max = magna.loc[t_max <= magna.time, 'time'].index[0]
                index_mid = np.round((index_min + index_max)/2).astype(int)
                t_start = magna['time'].iloc[index_mid-5].value/1e9/60/60/24
                t_end = magna['time'].iloc[index_mid+5].value/1e9/60/60/24
                ax.legend(ncol=3)

            ax.set_xlim([t_start, t_end])
            ax.set_ylabel('h (m ASL)')
            ax.set_xlabel('time')
        fig.suptitle('Merge quality check')
        fig.show()
        if fig_fp is not None:
            fig.savefig(fig_fp, dpi='150')

    # filter events too near
    events2 = filterDoubleClick(events,dt_min=dt_min)
    
    # get mean of continous data
    events2 = cont_to_events(cont,events2,dt_mean=dt_mean)
    
    # get mean of continous data
    eventsB2 = cont_to_events(cont,eventsB,dt_mean=dt_mean)

    # merge with magnaprobe data
    # Rather than add magnaprobe data, I prefer to add events2 data so that not to destroy any orphan magnaprobe data
    events2 = add_magna(events2,magna,dt_magnaprobe=0,tolerance=tolerance)

    # round numbers
    keys=['lat', 'lon', 'height', 'Q',  'sdn', 'sde', 'sdu', 'sdne', 'sdeu',
           'sdun',  'lat_std', 'lon_std', 'height_std','age','ratio','ns']
    decimals = pd.Series([9,9,3,0,4,4,4,4,4,4,11,11,4,2,2,0], 
                         index=keys)

    #eventsA
    events2=events2.round(decimals)

    events2 = pd.DataFrame(events2)
    eventsB2 = pd.DataFrame(eventsB2)

    if F_PLOT:
        return events2, eventsB2, fig
    else:
        return events2, eventsB2

def save(event_df, output_fp='PPK_Magna.csv',columns=None):
    """ 
    Parameters
    ----------
    event_df : pd.DataFrame()
        Dataframe containing data
    output_fp: STRING
        Filepath of the output file use to save the dataframe
    Returns
    -------
    Nothing

    """
    if columns is None:
        columns = ['time', 'Counter', 'lat', 'lon', 'height', 'DepthCm',
                            'sdn', 'sde', 'sdu', 'sdne', 'sdeu', 'sdun', 'age', 'ratio',
                            'lat_std', 'lon_std', 'height_std', 'Q', 'ns',
                            'BattVolts']
        for col in columns:
            if col not in event_df:
                columns.remove(col)

    # Number rounding
    round_dict = dict(zip(['lat', 'lon', 'height', 'Q',  'sdn', 'sde', 'sdu', 'sdne', 'sdeu', 'sdun',  'lat_std',
                           'lon_std', 'height_std','age','ratio','ns'],
                          [9, 9, 3, 0, 4, 4, 4, 4, 4, 4, 11, 11, 4, 2, 2, 0]))
    decimals = pd.Series(round_dict)
    for key in round_dict.keys():
        if key not in columns:
            decimals = decimals.drop(key)

    event_df = event_df.round(decimals)
    event_df.Q = np.int32(event_df.Q)
    event_df.ns = np.int32(event_df.ns)
    event_df.Counter = np.int32(event_df.Counter)
    # save files
    event_df.to_csv(output_fp, index=False,
                   columns=columns)