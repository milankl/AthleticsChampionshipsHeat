#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 12:59:47 2023

@author: navarro
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import datetime

from matplotlib.dates import DayLocator, HourLocator, DateFormatter
from matplotlib.ticker import FormatStrFormatter

import pandas as pd
import glob

import os, sys, inspect
HERE_PATH = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(HERE_PATH)
#from droplet import droplet
plt.rcParams['font.family'] = 'Arial'

path = ""

###############################################################################

## LOAD DATA
events = ["2007_osaka","2009_berlin","2011_daegu","2012_helsinki","2012_london",
          "2013_luzhniki","2014_zurich","2016_amsterdam","2018_berlin","2022_munich"]
events_long = ["WOC 2007 Osaka","WOC 2009 Berlin","WOC 2011 Daegu","EOC 2012 Helsinki","OG 2012 London",
          "WOC 2013 Moscow","EOC 2014 Zurich","EOC 2016 Amsterdam","EOC 2018 Berlin","EOC 2022 Munich"]
alphabet = ["a","b","c","d","e","f","g","h","i","j"]

t0 = datetime.datetime(1900,1,1)               # netCDF is in hours since

tstarts = [datetime.datetime(2007,8,24,0),     # start day of the event
           datetime.datetime(2009,8,15,0),
           datetime.datetime(2011,8,27,0),
           datetime.datetime(2012,6,27,0),
           datetime.datetime(2012,8,1,0),
           datetime.datetime(2013,8,10,0),
           datetime.datetime(2014,8,12,0),
           datetime.datetime(2016,7,6,0),
           datetime.datetime(2018,8,6,0),
           datetime.datetime(2022,8,11,0)]

tends = [datetime.datetime(2007,9,2,23),       # end day of the event
         datetime.datetime(2009,8,23,23),
         datetime.datetime(2011,9,4,23),
         datetime.datetime(2012,7,1,23),
         datetime.datetime(2012,8,12,23),
         datetime.datetime(2013,8,18,23),
         datetime.datetime(2014,8,17,23),
         datetime.datetime(2016,7,10,23),
         datetime.datetime(2018,8,12,23),
         datetime.datetime(2022,8,21,23)]

tends = [tend + datetime.timedelta(hours=1) for tend in tends]
#tstarts = [tstart -  datetime.timedelta(hours=10) for tstart in tstarts]

utcs = [9,2,9,3,1,3,2,2,2,2]


#UTCIs = []
#times = []
tabs=[]
for i,event in enumerate(events):
    files=glob.glob(event+'/*.nc')
    files.sort(key=os.path.getmtime)
    tab=pd.DataFrame()
    hour=0#+utcs[i]
    days=[]
    
    for j in range((tends[i]-tstarts[i]).days):
        days.append((tstarts[i]+datetime.timedelta(days=j)))
        
    for j in range((tends[i]-tstarts[i]).days):
        vn = nc.Dataset(files[j])
        utcis = vn.variables["utci"][:]
        utcis=np.array(utcis)
        for h in range(24):
            h=h-utcs[i]
            ligne=[]
            lignetime=tstarts[i]+datetime.timedelta(hours=hour)
            ligne.append(lignetime)
            ligne.append(utcis[h,1,1]-273.15)
            ligne.append(utcis[h,0,0]-273.15)
            ligne.append(utcis[h,1,0]-273.15)
            ligne.append(utcis[h,0,1]-273.15)
            ligne.append(utcis[h,1,2]-273.15)
            ligne.append(utcis[h,2,1]-273.15)
            ligne.append(utcis[h,0,2]-273.15)
            ligne.append(utcis[h,2,0]-273.15)
            ligne.append(utcis[h,2,2]-273.15)
            ligne.append(np.mean(utcis[h,:,:])-273.15)
            ligne.append(np.std(utcis[h,:,:]))
            ligne=pd.DataFrame(ligne)
            ligne=ligne.transpose()
            tab=pd.concat([tab,ligne])
            hour=hour+1

    tabs.append(tab.rename(columns={0:'Date',1:'point0',2:'point1',3:'point2',4:'point3',
                            5:'point4',6:'point5',7:'point6',8:'point7',9:'point8',
                            10:'Mean',11:'Std'}))
    # plt.figure()
    # plt.plot(np.array(tabs[i]['Mean'].astype(float)))
                                    
 
###############################################################################
test=[]
  
def temp_plotter(ax, dates, tlow, thigh, tstart, tend, utc):
    # these temperatures will be associated with the lower and upper end of the colormap
    clev = [10,32]
    tmin,tmax = (np.min(tlow),np.max(thigh))   # absurdly low and high temperatures  
    cmap = "rainbow"
    cmat = [[clev[0],clev[0]],[clev[1],clev[1]]]
    cmat_high = clev[1]*np.ones((2,2))
    cmat_low = clev[0]*np.ones((2,2))
    
    ax.contourf([tstart,tend-dt],clev,cmat,128,cmap=cmap)
    ax.contourf([tstart,tend-dt],[clev[1],tmax],cmat_high,128,vmin=clev[0],vmax=clev[1],cmap=cmap)
    ax.contourf([tstart,tend-dt],[tmin,clev[0]],cmat_low,128,vmin=clev[0],vmax=clev[1],cmap=cmap)
    ylim = ax.get_ylim()

    mask = [date >= tstart and date <= tend for date in dates]
    dates_masked = [date for date,m in zip(dates,mask) if m]
    tlow = tlow[np.array(mask)]
    thigh = thigh[np.array(mask)]
    test.append(dates_masked)
    test.append(np.array(mask))
    test.append(dates)
    n_tsteps = len(dates_masked)
    #dates_masked_utc = [dates-datetime.timedelta(hours=utc) for dates in dates_masked]
   
    # plot first half from the bottom, the first without transparency
    ax.fill_between(dates_masked,ylim[0]*np.ones(n_tsteps)-30,tlow,facecolor="w",alpha=1.,edgecolor="grey")
    ax.fill_between(dates_masked,thigh,ylim[1]*np.ones(n_tsteps)+30,facecolor="w",alpha=1.,edgecolor="grey")
    ax.xaxis.set_minor_locator(HourLocator(np.arange(6, 25, 6)))
    ax.set_xlim(tstart,tstart+datetime.timedelta(days=12))
    ax.set_xticklabels([])

###############################################################################
  
t0 = datetime.datetime(2000,1,1)
dt = datetime.timedelta(hours=0.1)

## PLOT
fig,axs = plt.subplots(10,1,figsize=(9,12),sharex=False)

# LABEL FORMATTING
for i,event in enumerate(events):
    axs[i].yaxis.set_major_formatter(FormatStrFormatter('%d'+u'\N{DEGREE SIGN}'+'C'))
    axs[i].grid(alpha=0.3)
    if (event != '2022_munich') and (event !='2016_amsterdam'):
        axs[i].text(0.84,0.85,events_long[i],transform=axs[i].transAxes,fontweight="bold",fontsize=10)
#    axs[i].text(0.98,0.80,alphabet[i],transform=axs[i].transAxes,fontweight="bold",fontsize=14)

# DATA PLOTTING
for i,event in enumerate(events):
    tstart=tstarts[i]
    tend=tends[i]-datetime.timedelta(hours=1)
    length_of_event = tends[i]-tstarts[i]
    ndays = length_of_event.days
    tab = tabs[i]
    tidays = tab["Date"]
    tab=tab.values
    utcim=np.mean(tab[:,1:9],1)
    utcimin=np.min(tab[:,1:9],1)
    utcimax=np.max(tab[:,1:9],1)
    utcim=utcim.astype(float)
    utcimin=utcimin.astype(float)
    utcimax=utcimax.astype(float)
    
    # temperature shaded curve
    temp_plotter(axs[i],tidays,utcimin,utcimax,tstart,tend,utcs[i])
    
    # Adapt y limits
   # axs[i].set_ylim(np.min(utcimin),np.max(utcimax)+(np.max(utcimax)-np.min(utcimin))*0.6)
    axs[i].set_ylim(np.mean(utcim)-17,np.mean(utcim)+27)
    
    # daily max as text
    for n in range(ndays):
        # find maximum
        utci_of_day_n = np.array([[u,t] for u,t in zip(utcim,tidays) if (t-tstarts[i]).days == n])
        utcimax_of_day_n = np.array([[u,t] for u,t in zip(utcimax,tidays) if (t-tstarts[i]).days == n])
        daily_max = np.max(utci_of_day_n[:,0])
        daily_max_max = np.max(utcimax_of_day_n[:,0])
        t_of_dailymax = utci_of_day_n[np.argmax(utci_of_day_n[:,0]),1]
        axs[i].text(t_of_dailymax,daily_max_max+1,"{:d}˚C".format(int(np.round(daily_max))),ha="center",fontsize=8)
        axs[i].xaxis.set_major_locator(DayLocator(interval=1))

axs[-1].xaxis.set_minor_locator(HourLocator(np.arange(6, 25, 6)))    # minor
axs[-1].xaxis.set_minor_formatter(DateFormatter("%Hh"))
axs[-1].xaxis.set_tick_params(which='minor', direction='out',pad=2,labelsize=6)


plt.setp(axs[-1].get_xticklabels(), ha="left")
axs[-1].xaxis.set_major_locator(DayLocator(interval=1))
axs[-1].xaxis.set_tick_params(which='major', direction='out',pad=10,labelsize=10)
axs[-1].set_xticklabels(['Day 1', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6', 'Day 7', 'Day 8', 'Day 9', 'Day 10', 'Day 11', 'Day 12'])

#axs[-1].xaxis.set_major_formatter(DateFormatter(" Day %d"))
axs[9].text(0.75,0.85,events_long[9],transform=axs[9].transAxes,fontweight="bold",fontsize=10)
axs[7].text(0.80,0.85,events_long[7],transform=axs[7].transAxes,fontweight="bold",fontsize=10)

#axs[3].set_yticks([10,20,30,40])

plt.tight_layout(pad=0)
plt.savefig("meteograms.pdf")