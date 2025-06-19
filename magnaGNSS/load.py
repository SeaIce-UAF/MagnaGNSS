"""
Function to load data from the raw magnaprobe data file
"""
__all__ = ['raw_data']

import logging
import pandas as pd
# import pint
# import pint_pandas as pint

logger = logging.getLogger(__name__)

# Define
# ureg = pint.UnitRegistry(auto_reduce_dimensions=True)
# pint.set_application_registry(ureg)
#ureg = pint.get_application_registry()

def read_raw(fp, header_row=1):
    """
    Read raw data file from filepath fp.
    :param fp: string
        Filepath to the magnaprobe raw data file.
    :param header_row: int
    :return:
    """
    if not fp.endswith(".dat"):
        logger.warning("This may not be an original raw datafile, extension is not '.dat'")
    else:
        pass

    df = pd.read_csv(fp, header=header_row)
    return df


def format_col_headers(df):
    """
    Format column headers. We prioritise for sanitation:
    - lower case
    - singular form. E.g. second rather than seconds
    - use underscore (_), instead of spacing. E.g. fix_quality rather than Fix quality
    - short abbreviation. E.g. n_satellite rather than NbrSatellite

    :param df: pd.DataFrame()
        Input Dataframe containing the raw data, imported with read_row.
    :return:
        Dataframe with formatted and sanitized column headers
    """
    # sanitize headers
    header_rename_d = {'TIMESTAMP': 'datetime', 'RECORD': 'record', 'Counter': 'counter',
                       'DepthCm': 'depth_cm', 'Battvolts': 'batt_v', 'Longitude_a': 'longitude_a',
                       'Longitude_b': 'longitude_b', 'nmbr_satellites':'n_satellite',
                       'HDOP': 'hdop', 'altitudeB': 'altitude_b', 'Depthvolts': 'depth_v',
                       'LatitudeDDDDD': 'latitude_d', 'LongitudeDDDDD': 'longitude_d'}
    header_rename_d = {k: header_rename_d[k] for k in iter(header_rename_d.keys()) if k in df.columns}
    df = df.rename(columns=header_rename_d)

    return df


def remove_junk_row(df):
    """
    Remove junk rows, line 3 and 4 in the original magnaprobe raw datafile

    :param df: pd.DataFrame()
        Dataframe containing the raw data.
    :return: pd.DataFrame()
        Dataframe stripped from the 2 first junk rows.
    """

    df = df.drop(df.index[:2])
    return df


def check_datetime(df):
    """
    Format timestamp to ISO. If timestamps are invalid, artificial timestamps are generated
    starting at midnight local time, with 1-minute increments.

    :param df: pd.DataFrame()
        Dataframe containing the raw data.
    :param date: dt.datetime()
        Date at which the dataset was acquired
    :return:
        Dataframe with formatted timestamp.
    """

    try:
        df['datetime'] = pd.to_datetime(df['datetime'], format='ISO8601')
    except ValueError:
        df['datetime'] = pd.to_datetime(df['datetime'])
    else:
        pass
    return df


def compute_coordinate(df):
    """
    Compute latitude and longitude coordinate in degree from the integer (_a) and minute (_b)
     fields. Unnecessary fields are then dropped (lat/long -_a, lat/long -_b,'lat/long -ddddd')
     or renamed (Altitudeb > Altitude)

    :param df: pd.DataFrame()
        Dataframe containing the raw data.
    :return: pd.DataFrame()
        Dataframe containing the latitude, respectively longitude, in the Latitude, respectively
         Longitude columns, given in Decimal Degree.
    """
    df['latitude'] = df['latitude_a'].astype('float') + df['latitude_b'].astype('float')/60
    df['longitude'] = df['longitude_a'].astype('float') + df['longitude_b'].astype('float')/60

    # drop unnecessary coordinate columns
    df.drop(['latitude_a', 'latitude_b', 'longitude_a', 'longitude_b', 'latitude_d',
             'longitude_d'], axis=1, inplace=True)
    df.rename({'altitude_b': 'altitude'}, axis=1, inplace=True)
    df['altitude']= df['altitude'].astype(float)
    # TODO: set units
    # df['Latitude'] = df['Latitude'].astype('pint[degree]')
    # df['Longitude'] = df['Longitude'].astype('pint[degree]')
    return df


def convert_depth_to_m(df):
    """
    Convert snow depth from cm to m. Remove unnecessary `Snowcm` column field.

    :param df: pd.DataFrame()
        Dataframe containing the raw data.
    :return: pd.DataFrame()
        Dataframe containing the snow depth given only in meter (m)
    """
    df['depth'] = df['depth_cm'].astype(float) / 100.0
    df.drop(columns=['depth_cm'], inplace=True)
    # TODO: set units
    # df['Depth'] = df['Depth'].astype('pint[m]')
    return df


def drop_unnecessary_columns(df, columns=None):
    """
    :param df: pd.DataFrame()
        Dataframe containing the raw data.
    :return: pd.DataFrame()
        Dataframe containing the snow depth given only in meter (m)
    """
    if columns is None:
        columns = ['depth_v', 'latitude_d', 'longitude_d', 'month', 'dayofmonth',
                   'hourofday', 'minutes', 'seconds', 'microseconds']
    # check if column exist in dataframe
    columns = [c for c in columns if c in df.columns]
    # remove column
    df.drop(columns=columns, inplace=True)
    return df




def raw_data(fp):
    """
    Load and clean raw data from filepath by computing local UTM coordinate, removing unnecessary headers and sanitizing
    the other headers

    :param fp: string
        Input data filepath
    :return df: pd.DataFrame()
        Dataframe containing the clean raw data
    """
    df = read_raw(fp)
    df = format_col_headers(df)
    df = remove_junk_row(df)
    df = check_datetime(df)
    df = convert_depth_to_m(df)
    df = compute_coordinate(df)
    df = drop_unnecessary_columns(df)
    return df.reset_index(drop=True)
