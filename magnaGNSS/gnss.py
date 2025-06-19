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
import magnaGNSS

# %% Default variable
F_PLOT = False
DT_MIN = 0.1
WINDOW = 1
TOLERANCE = 1
DELTA_H = 8
DELTA_S = 16.6
F_USE_EVENT = True

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
                               names=['TOW', 'WNc', 'Source', 'Polarity', 'Offset' , 'Offset_TOW',
                                      'RxClkBias' , 'PVTAge'])
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
        event_a_df = pd.DataFrame(columns=['datetime'])
        event_b_df = pd.DataFrame(columns=['datetime'])
    else:
        event_a_df=event_df.query("Source=='EventB' ")
        event_b_df=event_df.query("Source=='EventA' ")

        event_a_df.reset_index(drop=True, inplace=True)
        event_b_df.reset_index(drop=True, inplace=True)
    return event_a_df, event_b_df


def gps_to_datetime(gps_week, time_of_week):
    """
    :param gps_week:
        GPS week number
    :return: time_of_week: int
        Time of the week in second
    :return: gps_time: dt.datetime
        GPS time
    """
    # GPS epoch start date (January 6, 1980)
    gps_epoch = pd.Timestamp("1980-01-06")

    # Calculate the total time in weeks and seconds
    gps_time = gps_epoch + timedelta(weeks=float(gps_week), seconds=time_of_week)

    return gps_time


# get deviation from reference data
# def RMSE_Ref(ref, events, data=None, d_radious=0.25):
#     # radius to look for nearby data points in meters
#     if data is None:
#         data = []
#     dist=[]
#     d_h=[]
#     index2=[]
#     dist_ref=[]
#     d_h_ref=[]
#
#     for i in ref.index:
#         index = []
#         dist2 = []
#         d_h2 = []
#         for j in events.index:
#             d=distance.distance([ref.loc[i,'lat'],ref.loc[i,'lon']],[events.loc[j,'lat'],events.loc[j,'lon']]).m
#             if  d< d_radious:
#                 if len(data)>0:
#                     # go through data and find same timestamp
#                     for k in data.index:
#                         if events.loc[j,'time']==data.loc[k,'time']:
#                             index.append(k)
#                             index2.append(k)
#                             d=distance.distance([ref.loc[i,'lat'],ref.loc[i,'lon']],
#                                                 [data.loc[k,'lat'],data.loc[k,'lon']]).m
#                             dist.append(d)
#                             dist2.append(d)
#                             d_h.append(ref.loc[i,'height']-data.loc[k,'height'])
#                             d_h2.append(ref.loc[i,'height']-data.loc[k,'height'])
#                 else:
#                     index.append(j)
#                     index2.append(j)
#                     dist.append(d)
#                     dist2.append(d)
#                     d_h.append(ref.loc[i,'height']-events.loc[j,'height'])
#                     d_h2.append(ref.loc[i,'height']-events.loc[j,'height'])
#
#         dist_ref.append(np.mean(dist2))
#         d_h_ref.append(np.mean(d_h2))
#
#     d_h=np.array(d_h)
#     dist=np.array(dist)
#
#     delta_hor = np.mean(np.abs(dist))
#     delta_h = np.mean(np.abs(d_h))
#     delta_hor_std = np.std(dist)
#     delta_h_std = np.std(np.abs(d_h))
#     RMSE_h = np.sqrt(np.sum(d_h**2) / len(d_h))
#     RMSE_x = np.sqrt(np.sum(dist**2) / len(dist))
#     print(f'Delta_h= {delta_h:.3f} +-{delta_h_std:.4f} m')
#     print(f'Delta_x_hor={delta_hor:.3f} +-{delta_hor_std:.4f} m')
#     print(f'RMSE(h)= {RMSE_h:.3f} m')
#     print(f'RMSE(x)= {RMSE_x:.3f} m')
#
#     return  RMSE_h,RMSE_x,index2,d_h,dist,dist_ref,d_h_ref



# filter out events happening within dt_min=0.1 s of each other
def filter_double_click(event_df, dt_min=DT_MIN):
    """
    Parameters
    ----------
    event_df : pd.DataFrame()
        Dataframe containing events to filters
    dt_min : FLOAT, optional
        Filter events nearer than dt_min. The default is 0.1 min.

    Returns
    -------
    eventf_df : pd.DataFrame()
        Dataframe containing the filtered events
    """
    print('\nTrowing out double click events\n------------------------------------')
    # Find index of the second of two consecutive events, when both events happened within 0.1 s
    index_dc = event_df.loc[event_df.datetime.diff() < pd.Timedelta(seconds=dt_min)].index
    dts = event_df.datetime.diff().dt.total_seconds()[index_dc]

    event_df.drop(index=index_dc, inplace=True)
    event_df.reset_index(drop=True, inplace=True)

    if len(dts) == 0:
        print('No double click event found')
    else:
        print(f"Mean dt double points: {np.mean(dts):.3f} s")
        print(f"Std dt double points: {np.std(dts):.3f} s")

    return event_df


# get magnaprobe data
def add_magna(event_df, mg_df, dt_magnaprobe=0, tolerance=TOLERANCE):
    """

    Parameters
    ----------
    event_df : TYPE
        DESCRIPTION.
    mg_df : TYPE
        DESCRIPTION.
    dt_magnaprobe : TYPE, optional
        DESCRIPTION. The default is 1.
    tolerance: TYPE, optional
        DESCRIPTION. The default is 0.1.
    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    mg_df_ = mg_df.copy()
    if isinstance(dt_magnaprobe, pd.Timedelta):
        mg_df_['datetime'] = mg_df_['datetime'] + dt_magnaprobe
    else:
        mg_df_['datetime'] = mg_df_['datetime'] + pd.Timedelta(seconds=dt_magnaprobe)
    mg_df_.set_index('datetime',inplace=True)
    new_index = mg_df_.index.get_indexer(event_df.datetime,method='nearest',tolerance=pd.Timedelta(seconds=tolerance))
    del mg_df_
    # Reindex event_df with index from mg_df
    event_df.set_index(new_index)
    # Drop all event of event_df non-existing in mg_df
    event_df = event_df[event_df.index >= 0]

    # rename latitude and longitude in the magnaprobe dataframe to avoid double header
    mg_df.rename(columns={'latitude': 'latitude_mg', 'longitude': 'longitude_mg'}, inplace=True)

    # concatenate data frame together
    event_df = pd.concat([mg_df, event_df], axis=1)
    return event_df


# compute the mean of the GNSS position continous data around each event
def mean_position_of_events(pos_df, event_df, dt_mean=0.5):
    """

    Parameters
    ----------
    pos_df : pd.DataFrame()
        Dataframe containing the continuous GNSS position data.
    event_df : pd.DataFrame()
        Dataframe containing the magnaprobe measurement event data.
    dt_mean : float, optional
        Duration in seconds of the period centered on an event used for averaging continuous data.
        Default is 1, but

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    lines=[]
    indices=[]
    for ii in event_df.index:
        columns = ['lat', 'lon', 'height', 'Q', 'ns', 'sdn','sde', 'sdu', 'sdne', 'sdeu', 'sdun',
                   'age', 'ratio', 'datetime']
        t1 = event_df.loc[ii,'datetime'] - pd.Timedelta(seconds=dt_mean)
        t2 = event_df.loc[ii,'datetime'] + pd.Timedelta(seconds=dt_mean)
        a = pos_df.loc[(pos_df['datetime'] >= t1) & (pos_df['datetime'] <= t2), columns]

        line=a.mean()
        line['lat_std']=a.loc[:,'lat'].std()
        line['lon_std']=a.loc[:,'lon'].std()
        line['height_std']=a.loc[:,'height'].std()

        line['datetime']=event_df.loc[ii,'datetime']

        if len(a) > 0:   # trow away events if there are no corresponding countinous events!
            lines.append(line)
            indices.append(ii)

    return pd.DataFrame(lines, index=indices)


def merge(mg_fp, pos_fp, event_fp=None, f_plot=F_PLOT, fig_fp=None, dt_min=DT_MIN, window=WINDOW,
          tolerance=TOLERANCE, delta_h=DELTA_H, delta_s=DELTA_S, f_use_event=F_USE_EVENT):
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
        Duration in seconds of the period centered on an event used for averaging continuous data.
        The default is 0.5 s.
    tolerance : FLOAT, optional
        Tolerance in seconds for syncing GNSS time to maganprobe timestamp. The default is 0.1.
    delta_h : FLOAT, optional
        Time shift Magnaprobe in hours. GPS time to local time. The default is 8.
    delta_s : FLOAT, optional
        Time shift between Magnaprobe and GNSS in seconds (leap seconds+ error). The default is 16.62
    f_use_event: bool, optional
            Default true. If false the GNSS Timestamps are ignored and the data are synced based on
             the Magnaprobe time. To be used in case the conenction between the GNSS and Magnaprobe
             didn't work.
                
    Returns
    -------
        None
    """

    dt_mean = float(window)/2

    # read files
    event_df = read_events(event_fp)
    event_a_df, event_b_df = split_events(event_df)

    pos_df = pd.read_csv(pos_fp, header=9, sep=r'\s+', names=['date', 'GPST', 'lat', 'lon',
                                                              'height', 'Q', 'ns', 'sdn', 'sde',
                                                              'sdu', 'sdne', 'sdeu', 'sdun',
                                                              'age', 'ratio'])
    pos_df['datetime'] = pd.to_datetime(pos_df['date']+' '+pos_df['GPST'])

    mg_df = magnaGNSS.load.raw_data(mg_fp)
    # time difference
    mg_df['datetime']=mg_df['datetime']+pd.Timedelta(hours=delta_h, minutes=00, seconds=delta_s)

    if len(event_a_df)==0 or not f_use_event:
        event_a_df = pd.DataFrame(mg_df['datetime'].values,columns=['datetime'])
        event_a_empty = True
        print('eventsA is empty!! Using Magnaprobe timestamp instead!!')
    else:
        event_a_empty = False

    if f_plot or fig_fp is not None:
        # Plot only for the overlapping time period between magnaprobe and gnss
        t_max = min(mg_df.datetime.max(), pos_df.datetime.max())
        t_min = max(mg_df.datetime.min(), pos_df.datetime.min())
        # Convert t_min in t_index
        t_start = mg_df.loc[mg_df.datetime <= t_min, 'datetime'].iloc[-1].value/1e9/60/60/24
        t_end = mg_df.loc[t_max <= mg_df.datetime, 'datetime'].iloc[0].value/1e9/60/60/24

        fig, axs = pl.subplots(2, 1, figsize=[8.5,11])
        for ii, ax in enumerate(axs):
            ax.plot(pos_df.datetime,pos_df.height,'x',label='GPSRTK GNSS continous')

            ax.scatter(mg_df.datetime,mg_df.altitude,marker='o',edgecolors='g',facecolor='None',label='Magna GPS')
            ylim = ax.get_ylim()

            if event_a_empty:
                for t in event_a_df.datetime:
                    ax.plot([t,t],ylim,'--g')
                    if ii == 0:
                        ax.text(0.05,0.95,'eventsA is empty!! Using MAgnaprobe timestamp instead!!',
                                transform=ax.transAxes,
                                ha='left', va='top', fontsize=9, color='red')
            else:
                for t in event_df.datetime:
                    ax.plot([t,t],ylim,'--k')

            for t in event_b_df.datetime:
                ax.plot([t,t],ylim,'--r')

            ax.plot([t, t],ylim,'--k',label='EventsA')
            ax.plot([t, t],ylim,'--r',label='EventsB')
            if event_a_empty:
                ax.plot([],[],'--g',label='GPS Magnaprobe timestamps')

            h_max = max([mg_df.altitude.max(), pos_df.height.max()])
            h_min = min([mg_df.altitude.min(), pos_df.height.min()])
            ax.set_ylim([h_min-0.2, h_max+0.2])

            if ii==1:
                index_min = mg_df.loc[mg_df.datetime <= t_min, 'datetime'].index[-1]
                index_max = mg_df.loc[t_max <= mg_df.datetime, 'datetime'].index[0]
                index_mid = np.round((index_min + index_max)/2).astype(int)
                t_start = mg_df['datetime'].iloc[index_mid-5].value/1e9/60/60/24
                t_end = mg_df['datetime'].iloc[index_mid+5].value/1e9/60/60/24
                ax.legend(ncol=3)

            ax.set_xlim([t_start, t_end])
            ax.set_ylabel('h (m ASL)')
            ax.set_xlabel('time')
        fig.suptitle('Merge quality check')
        fig.show()
        if fig_fp is not None:
            fig.savefig(fig_fp, dpi='150')

    # filter events too near
    event_a_df = filter_double_click(event_a_df,dt_min=dt_min)

    # get mean of continous data for events A and B
    event_a_df = mean_position_of_events(pos_df, event_a_df, dt_mean=dt_mean)
    event_b_df = mean_position_of_events(pos_df,event_b_df,dt_mean=dt_mean)

    # merge with magnaprobe data
    events_df = add_magna(event_a_df, mg_df, dt_magnaprobe=0, tolerance=tolerance)

    if f_plot or fig_fp is not None:
        return events_df, event_b_df, fig
    return events_df, event_b_df

def save(event_df, output_fp='PPK_mg_df.csv', columns=None):
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
    event_df.to_csv(output_fp, index=False, columns=columns)
