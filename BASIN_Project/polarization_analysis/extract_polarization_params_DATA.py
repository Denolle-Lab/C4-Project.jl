#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# from google.colab import drive
# drive.mount('/content/drive', force_remount=True)
# if "Trace" not in locals():
#     get_ipython().system('pip install obspy')
import numpy as np
# from scipy.io import loadmat
from scipy.signal import sosfilt, hann, iirfilter, zpk2sos
from obspy import Trace, Stream
import re
import polarization_analysis
import pandas as pd
import h5py

# =======================================
# input  =========================
# =======================================
# PREPARING TO READ IN JULIAN s DATA -- Copy from Marine
sourcestations = ["BRE", "IPT", "RUS", "TA2", "CHN", "LPC", "SNO", "CJM", "PASC", "SRN", "SVD"]
# station at which source is placed
metadata_file = 'full_socal.csv'
# bandpass filter
freqmin = 0.1
freqmax = 0.35
stacktype = "linear"
dir0 = '/mnt/DATA0/BASIN/nodestack/'
winlen = 5.0  # seconds
step = 2.5   # seconds
method = "frequency"  # set to frequency for John Vidale's method or timedomain for Sin-Mei Wu's method
ma = 0  # moving average: Smoothing for the polarization parameters after polarization analysis.
# set to 0 for no smoothing
# =======================================
# end input  =========================
# =======================================
print("Please note: It is assumed that data are in folders called\
    {}/2017 and {}/2018 and {}/2019, and are named as CI.{}_2017.h5 etc.".format(
        dir0, dir0, dir0, sourcestations[0]))

for sourcestation in sourcestations:
    print(sourcestation)
    locs = pd.read_csv(metadata_file)
    latsrc = locs[locs.station == sourcestation].latitude.values[0]
    lonsrc = locs[locs.station == sourcestation].longitude.values[0]

    # reassample the data over the years
    
    lfile = dir0 + "2017/CI.{}_2017.h5".format(sourcestation)
    lfile2 = dir0 + "2018/CI.{}_2018.h5".format(sourcestation)
    lfile3 = dir0 +"2019/CI.{}_2019.h5".format(sourcestation)

    if 'f' in locals(): f.close()
    if 'f2' in locals(): f2.close()
    if 'f3' in locals(): f3.close()

    f=h5py.File(lfile,'r')
    f2=h5py.File(lfile2,'r')
    f3=h5py.File(lfile3,'r')
    print(f, f2, f3)
    print('list of receivers')
    cc_len, maxlag = f['meta']['cc_len'][...], f['meta']['maxlag'][...]

    L = list(f.keys())
    L2 = list(f2.keys())
    L3 = list(f3.keys())
    Lc = L.copy()
    Lc.extend(list(f2.keys()))
    Lc.extend(list(f3.keys()))
    L_all=sorted(list(set(Lc)))
    LL = list(f[L[0]]) ;LLL=list(f[L[0]][LL[0]])
    data=f[L[0]][LL[0]][LLL[0]][:]
    d_len = max(data.shape)
    fs = (d_len - 1) / 2 / maxlag
    print(fs)
    print(f["meta"].keys())
    print(L_all)
    print(LL)
    print(d_len)


    # In[ ]:


    # ------------------------------------------------------------------
    # with modifications from obspy, LGPL-3 license
    # Modifications: Do not apply the filter, but only return it.
    #
    #
    # Filename: filter.py
    #  Purpose: Various Seismogram Filtering Functions
    #   Author: Tobias Megies, Moritz Beyreuther, Yannik Behr
    #    Email: tobias.megies@geophysik.uni-muenchen.de
    #
    # Copyright (C) 2009 Tobias Megies, Moritz Beyreuther, Yannik Behr
    # --------------------------------------------------------------------
    def bandpass(freqmin, freqmax, df, corners=4):
        """
        Butterworth-Bandpass Filter.

        Filter data from ``freqmin`` to ``freqmax`` using ``corners``
        corners.
        The filter uses :func:`scipy.signal.iirfilter` (for design)
        and :func:`scipy.signal.sosfilt` (for applying the filter).

        :type data: numpy.ndarray
        :param data: Data to filter.
        :param freqmin: Pass band low corner frequency.
        :param freqmax: Pass band high corner frequency.
        :param df: Sampling rate in Hz.
        :param corners: Filter corners / order.
        :param zerophase: If True, apply filter once forwards and once backwards.
            This results in twice the filter order but zero phase shift in
            the resulting filtered trace.
        :return: Filtered data.
        """
        fe = 0.5 * df
        low = freqmin / fe
        high = freqmax / fe
        # raise for some bad scenarios
        if high > 1:
            high = 1.0
            msg = "Selected high corner frequency is above Nyquist. " +               "Setting Nyquist as high corner."
            print(msg)
        if low > 1:
            msg = "Selected low corner frequency is above Nyquist."
            raise ValueError(msg)
        z, p, k = iirfilter(corners, [low, high], btype='band',
                            ftype='butter', output='zpk')
        sos = zpk2sos(z, p, k)
        return sos


    # In[ ]:


    def get_distance(lat1, lon1, lat2, lon2):
        """
        (c) Jonas Igel
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians
        lon1 = np.deg2rad(lon1)
        lat1 = np.deg2rad(lat1)
        lon2 = np.deg2rad(lon2)
        lat2 = np.deg2rad(lat2)

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
        c = 2 * np.arcsin(np.sqrt(a))
        km = 6371 * c
        return km*1000

    # taper function
    def def_taper(npts, frac_taper=0.1):
        n_samples = int(npts * frac_taper)
        taper = np.ones(npts)
        taper[0: n_samples] = hann(2 * n_samples)[0: n_samples]
        taper[-n_samples: ] = hann(2 * n_samples)[-n_samples: ]
        return(taper)

    # prepare:
    # set up filter
    sos = bandpass(freqmin=freqmin, freqmax=freqmax, df=fs)

    # set up dataframe for output
    df = pd.DataFrame()

    # set up taper
    taper = def_taper((d_len - 1)//2)


    # In[ ]:


    # loop over receiver locations
    ix_mid = (d_len - 1) // 2
    noise_ampl = 0.1
    for component_source in ["E", "N", "Z"]:
        times = []
        azimuths = []
        ellipticities = []
        amplitudes = []
        pgvs = []
        latitude = []
        longitude = []
        dist_m = []
        stas_out = []
        for i, netsta in enumerate(L_all):
            
            #print()
            try:
                sta = netsta.split(".")[1]
                print(sta)
                if sta[0:2] == "B3":
                    sta = re.sub("B3", "B30", sta)
                elif sta[0:2] == "B2":
                    sta = re.sub("B2", "B20", sta)
                lat = locs[locs.station==sta].latitude.values[0]
                lon = locs[locs.station==sta].longitude.values[0]
            except IndexError:
                print("Problems with {}".format(netsta), ",")
                continue
            # print(sta, lat, lon)
            distance = get_distance(latsrc, lonsrc, lat, lon)
            # find the x, y, and z correlation
            # cut out the anticausal part
            if i % 20 == 0:
                print(i, end=",")
            vx = np.zeros(ix_mid)
            vy = np.zeros(ix_mid)
            vz = np.zeros(ix_mid)
            cnt_x = 0
            cnt_y = 0
            cnt_z = 0
            for fin in [f, f2, f3]:
                try:
                    comp = "{}E".format(component_source)
                    #vx += np.gradient(fin[netsta][comp]['robust'][:ix_mid][::-1], 1./fs)
                    vx += fin[netsta][comp][stacktype][:ix_mid][::-1]
                    cnt_x += 1
                except KeyError:
                    print("did not find E")
                    pass
                try:
                    comp = "{}N".format(component_source)
                    #vy += np.gradient(fin[netsta][comp]['robust'][:ix_mid][::-1], 1./fs)
                    vy += fin[netsta][comp][stacktype][:ix_mid][::-1]
                    cnt_y += 1
                except KeyError:
                    print("did not find N")
                    pass
                try:
                    comp = "{}Z".format(component_source)
                    # vz += np.gradient(fin[netsta][comp]['robust'][:ix_mid][::-1], 1./fs)
                    vz += fin[netsta][comp][stacktype][:ix_mid][::-1]
                    cnt_z += 1
                except KeyError:
                    print("did not find Z")
                    pass
            
            if 0.0 in [vdat.sum() for vdat in [vx, vy, vz]]:
                print("zero data")
                continue
           # dist_m.append(distance)
            latitude.append(lat)
            longitude.append(lon)
            stas_out.append(sta)

            # filter data
            vx = sosfilt(sos, taper * vx / cnt_z)
            vy = sosfilt(sos, taper * vy / cnt_y)
            vz = sosfilt(sos, taper * vz / cnt_z)
            
            # get azimuth, ellipticity
            tr = Stream()
            tr += Trace(data=vx)
            tr += Trace(data=vy)
            tr += Trace(data=vz)
            for t in tr:
                t.stats.sampling_rate = fs
            #  print(tr)
            ts, az, ps, ls, amps, pgs = polarization_analysis.pol_win(tr, winlen, step=step, offset=0, method=method, moving_average=ma)
            azimuths.append(az)
            times.append(ts)
            ellipticities.append(ps)
            amplitudes.append(amps)
            pgvs.append(pgs)
        azimuths = np.array(azimuths)
        ellipticities = np.array(ellipticities)
        amplitudes = np.array(amplitudes)
        pgvs = np.array(pgvs)

      
        # add to dataframe
        for i, t in enumerate(times[0]):
            print(i, end=",")
            key = "azimuth_S{}_{}s".format(component_source, t)
            df[key] = azimuths[:, i]
            key = "ellipticity_S{}_{}s".format(component_source, t)
            df[key] = ellipticities[:, i]
            key = "rms_S{}_{}s".format(component_source, t)
            df[key] = amplitudes[:, i]
            key = "pgv_S{}_{}s".format(component_source, t)
            df[key] = pgvs[:, i]

      #  df["dist_m"] = dist_m
        df["latitude_{}".format(component_source)] = latitude
        df["longitude_{}".format(component_source)] = longitude
        df["station_{}".format(component_source)] = stas_out
    df.to_csv("polarization_parameters_{}-{}_{}-{}Hz_{}s_{}s_nodestack_{}.csv".format(sourcestation, "DATAV", freqmin, freqmax, winlen, step, stacktype))



