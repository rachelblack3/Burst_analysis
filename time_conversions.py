
""" This module contains utilties for converting different datetime formats to Python datetime objects """
import datetime
import cdflib

def h5_time_conversion(bytes_time) -> datetime.datetime:
    """
    Converts ISO 8601 datetime, which contains leap seconds, into Python
    datetime object

    :param bytes_time: time in bytes to convert
    :return: datetime in Python datetime object
    """
    date_str, time_str = bytes_time.decode("UTF-8").replace("Z", "").split("T")
    year, month, day = date_str.split("-")
    date = datetime.datetime(int(year), int(month), int(day))
    no_ms_time, dot, ms = time_str.partition(".")
    hours, minutes, seconds = no_ms_time.split(":")
    date = date + datetime.timedelta(
        hours=int(hours), minutes=int(minutes), seconds=int(seconds)
    )
    if dot:
        date = date + datetime.timedelta(milliseconds=int(ms))
    return date

def cdflib_time_conversion(time_str: str) -> datetime.datetime:
    """
    Converts cdflib CDF_EPOCH16 (picoseconds since Year 0) to Python 
    datetime object

    :param time_str: time in CDF_EPOCH16 as string to convert
    :return: datetime in Python datetime object
    """
    breakdown = cdflib.cdfepoch.breakdown(time_str)

    modified_breakdown = breakdown[:6] + [1000*breakdown[6]+breakdown[7]]
    return datetime.datetime(*modified_breakdown)