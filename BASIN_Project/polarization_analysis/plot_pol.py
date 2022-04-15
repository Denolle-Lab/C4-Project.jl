#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

# In[2]:
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

# =======================================
# input  =========================
# =======================================
# station at which source is placed in simulation / virtual source station
sourcestations = ["DEV"]#"BRE", "SNO", "SVD", "IPT", "BRE", "CHN", "TA2", "SRN", "LPC", "CJM", "PASC", "RUS"]
linecodes = ["B4", "G1", "G2", "B1", "B2", "B3", "B4", "B5"]
#simulation_file = "/content/drive/Shared drives/BASIN/SIMULATIONS/3D+topo/7_DEV_{}.mat".format(sim_output)
metadata_file = "basin/full_socal.csv"
indir = "basin/"
# important: The following string needs to be edited to point us to the files produced
# by extract_polarization_params_DATA.py
descriptivestring = "0.1-0.35Hz_10.0s_2.5s"

source_comps = ["E", "N", "Z"] #["x", "y", "z"]#   # "Q"
# below important: Need to edit the amplitude scaling depending on the stacking type!
# 1 for simul, 100 for robust stack, 5000 for linear stack, 1.e4 for quakes LOW FREQ, 1 or 10 for quakes HIGH freq
ampfac = 5.e3
# use opacity or not: if True, low RMS windows are made more transparent and "disappear"
useopac = True
# color scale:
cm = plt.cm.PRGn
cm2 = plt.cm.PRGn
# lag limit: (seconds)
xlims = [0., 240.]

# if it does not look nice, we can try to restrict the opacity range below and edit lines 205 ff
amin = None
amax = 1.0

# this below just is a way to know which station is the "line end" for each line
refstadict = {"B1": "B1260", "G2": "G2050", "B4": "B4096", "G1": "G1060", "G3": "G3034", "B3": "B3084", "B5": "B5062", "B2": "B2060",
"G4": "G4010"}
# =======================================
# end input  =========================
# =======================================


subplstr = "1" + str(len(source_comps))

for sourcestation in sourcestations:
    polarization_file = "{}/polarization_parameters_{}-DATAV_{}_nodestack_linear.csv".format(indir, sourcestation, descriptivestring) # "msr_drive/polarization_parameters_DEV-V_0.1-0.35Hz_10.0s_3.0s.csv"# "msr_drive/polarization_parameters_SVD-DATA_0.1-0.35Hz_10.0s_3.0s.csv"
    for linecode in linecodes:    
        
        refstation_for_lin = refstadict[linecode] #"B1260"  #"B6033" #"B1260" #"G2050"#"B4096G2
        
        
        try:
            freqmin= float(os.path.basename(polarization_file).split("_")[3].split("-")[0])
            freqmax= float(os.path.basename(polarization_file).split("_")[3].split("-")[1][:-3])
        except ValueError:
            freqmin= float(os.path.basename(polarization_file).split("_")[4].split("-")[0])
            freqmax= float(os.path.basename(polarization_file).split("_")[4].split("-")[1][:-3])

        #times_pol = np.linspace(0.0, 90., 19)
        # plt.title("Azimuth (degree)"9
        if linecode == "B1":
            orientation  = "E-W"
        else:
            orientation = "S-N"

        pols = pd.read_csv(polarization_file)
        locs = pd.read_csv(metadata_file)
        print(pols.keys()
            )


        # get time array:
        times_pol = []
        for k in list(pols.keys()):
            if "rms" in k:
                times_pol.append(float(k.split("_")[-1][:-1]))
        #print(times_pol)
        times_pol = np.array(times_pol)
        nt = len(times_pol)
        if len(source_comps) > 1:
            times_pol = times_pol[:nt//3]
        print(times_pol)

        azimuths = []
        ellipticities = []
        lons_plot = []
        lats_plot = []
        sta_plot = []
        amps_plot = []

        for component_source in source_comps:

            az_part = []
            ell_part = []
            lon_part = []
            lat_part = []
            sta_part = []
            amp_part = []
            for i in range(len(pols)):

                #stations = locs[locs.latitude == pols.latitude.values[i]]
                #station = stations.station[stations.longitude == pols.longitude.values[i]].values[0]
                try:
                    station = pols["station_{}".format(component_source)].values[i]
                    print(station)
                except:
                    continue
                if station[0:2] != linecode:
                    continue

                if len(station) == 4:
                    station = station[0:2] + "0" + station[2:]

                #if station not in goodstations:
                #    continue
                print(station)
                lon_part.append(pols["longitude_{}".format(component_source)].values[i])
                lat_part.append(pols["latitude_{}".format(component_source)].values[i])
                sta_part.append(station)
                az = []
                ell = []
                amp = []
                for t in times_pol:
                    key = "azimuth_S{}_{}s".format(component_source, t)
                    az.append(pols[key].values[i])
                    key = "ellipticity_S{}_{}s".format(component_source, t)
                    ell.append(pols[key].values[i])
                    key ="pgv_S{}_{}s".format(component_source, t)
                    amp.append(pols[key].values[i]) 
                
                ell_part.append(ell)
                az_part.append(az)
                amp_part.append(amp)
            amps_plot.append(amp_part)
            lats_plot.append(lat_part)
            lons_plot.append(lon_part)
            azimuths.append(az_part)
            ellipticities.append(ell_part)
            sta_plot.append(sta_part)

        azimuths = np.array(azimuths)
        amps_plot = np.array(amps_plot)
        ellipticities = np.array(ellipticities)
        lons_plot = np.array(lons_plot)
        lats_plot = np.array(lats_plot)
        sta_plot = np.array(sta_plot)

        plot_dist = []
        for i, s_comp in enumerate(source_comps):
            ix_sta0 = np.where(sta_plot[i] == refstation_for_lin)[0][0]
            print(sta_plot[i, ix_sta0], lats_plot[i, ix_sta0], lons_plot[i, ix_sta0])
            dist_part = []
            for j in range(len(sta_plot[0])):
                
                dist_part.append(get_distance(lats_plot[i, j], lons_plot[i, j],
                                              lats_plot[i, ix_sta0],
                                              lons_plot[i, ix_sta0]) / 1000.)
            plot_dist.append(dist_part)
        plot_dist = np.array(plot_dist)
        print(plot_dist.shape)

        fig = plt.figure(figsize=(11, 4))
        axes = []
        for ix_subplot in range(len(source_comps)):


            a_plot = amps_plot[ix_subplot].copy()
            ells_plot = ellipticities[ix_subplot].copy()
            plot_dist_plot = plot_dist[ix_subplot].copy()

            ax = fig.add_subplot(subplstr+str(ix_subplot + 1))
            dt = times_pol[1]
            # print(ells_plot_sorted.shape)
            # if amin is None:
            #     amin = a_sorted.min()
            # else:
            #     amin = amin * a_sorted.min()
            # if amax is None:
            #     amax = a_sorted.max()
            # else:
            #     amax = amax * a_sorted.max()
            # opac = (a_sorted - amin) / (amax - amin)
            opac = (a_plot - a_plot.min()) / (a_plot.max() - a_plot.min())
            # opac[opac < 0.02] = 0.0
            xv, yv = np.meshgrid(times_pol, plot_dist_plot)
            p = ax.pcolormesh(xv, yv, ells_plot, cmap=cm, shading="gouraud", edgecolors="face", vmin=0., vmax=1.)
            if useopac:
                plt.savefig("dummy.png")
                print("MAX/MIN OPACITY: {}, {}".format(opac.max(), opac.min()))
                print("% OPACITY>0.95: {}".format(len(np.where(opac.flatten() > 0.95)[0]) / len(opac.flatten()) * 100))
                print("N OPACITY>0.95: {}".format(len(np.where(opac.flatten() > 0.95)[0])))

                for i, j in zip(p.get_facecolors(), opac.flatten()):
                    i[3] = j  # Set the alpha value of the RGBA tuple using m2
            axes.append(p)
            if ix_subplot == 0:
                plt.ylabel("Dist. {} (km)".format(orientation))

            for ix_amps in [len(plot_dist_plot)//4, len(plot_dist_plot)//2, len(plot_dist_plot)//4*3]:
                plt.plot(times_pol, ampfac * (a_plot[ix_amps]-a_plot[ix_amps].mean()) + plot_dist_plot[ix_amps], color="0.7")
            
            plt.xlim(xlims)
            plt.ylim([0., plot_dist_plot.max()])
            plt.xlabel("Time (s)")
            if source_comps[0] in ["E", "x"]:
                plt.title("Line: {}, Source: {}".format(linecode, ["E", "N", "Z"][ix_subplot]), fontweight="bold")
            else:
                plt.title("Source: {}".format(sourcestation))

            del p

        if source_comps[0] == "E":  
            plt.savefig("polar_{}_{}_{}-{}Hz_10s_DATAV.png".format(sourcestation, linecode, freqmin, freqmax), dpi=100)
        elif source_comps[0] == "x":
            plt.savefig("polar_{}_{}_{}-{}Hz_10s_SIML.png".format(sourcestation, linecode, freqmin, freqmax), dpi=100)
        elif source_comps[0] == "Q":
            plt.savefig("polar_{}_{}_{}-{}Hz_10s_QUAKE.png".format(sourcestation, linecode, freqmin, freqmax), dpi=100)

        plt.close()
        #plt.show()



