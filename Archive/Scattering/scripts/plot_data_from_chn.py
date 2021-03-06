# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 16:02:02 2021

@author: David
"""

import spinmob   as _s
import time      as _t
import struct    as _struct
import numpy     as _n

def load_chn(path=None, **kwargs):
    """
    Loads a Maestro CHN file at the specified path. Will return a databox with all
    the information.
    Parameters
    ----------
    path=None
        If None, a dialog will pop up. Otherwise specify a valid path (string).
    Optional keyword arguments, e.g., delimeter=',', are sent to spinmob.data.databox()
    """
    if path==None: path = _s.dialogs.load(filters="*.Chn")
    if path==None: return

    # Create a databox
    d = _s.data.databox(**kwargs)

    # Open the file
    fh = open(path,mode="rb")

    # Read the buffer
    buffer = fh.read()

    # Unpack
    [type, unitnumber, segment_number] = _struct.unpack("hhh",buffer[0:6])
    ascii_seconds = buffer[6:8].decode()
    [real_time_20ms, live_time_20ms] = _struct.unpack("ii",buffer[8:16])
    start_date = buffer[16:24].decode()
    start_time = buffer[24:28].decode()
    [channel_offset, channel_count] = _struct.unpack("hh",buffer[28:32])

    # Get the counts data
    spectrum = _n.zeros(channel_count,dtype=int)
    for i in range(channel_count): [spectrum[i]] = _struct.unpack("I",buffer[32+4*i:36+4*i])

    # Get the byte offset of the data
    offset = 36+4*(channel_count-1)

    # Get the size of the "Sample description" string
    info_size = int.from_bytes( buffer[offset+320:offset+321], "big" )

    # Extract "Sample description" string
    info = buffer[offset+321:offset+321+info_size].decode()

    d['Channel'] = range(channel_count)
    d['Counts']  = spectrum

    # Get the date in a nice format
    if "1" == start_date[7]: century = 20
    else:                    century = 19
    start_RFC2822 = "%s %s %02d%s %s:%s:%s" % (start_date[0:2], start_date[2:5], century, start_date[5:7], start_time[0:2], start_time[2:4], ascii_seconds)
    s = _t.strptime(start_RFC2822,"%d %b %Y %H:%M:%S")
    start = dict(
        year     = s.tm_year,
        month    = s.tm_mon,
        day      = s.tm_mday,
        hour     = s.tm_hour,
        minute   = s.tm_min,
        second   = s.tm_sec,
        year_day = s.tm_yday)

    # Header info
    d.path=path
    d.h(description = info,
        start_time  = start,
        real_time   = 0.02*real_time_20ms,
        live_time   = 0.02*live_time_20ms,
        path        = path,)

    return d

d = load_chn('Ba-133 Calibration 000.Chn')

import matplotlib.pyplot as plt

plt.plot(d[0],d[1])




