import xarray as xr
import numpy as np
from scipy.interpolate import interp1d
import glob
import os

import thermofeel
import datetime

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import gridspec
from matplotlib.dates import WeekdayLocator, DayLocator, HourLocator, DateFormatter, drange, date2num, num2date
from matplotlib.ticker import FormatStrFormatter

from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.image as mpimg

import os, sys, inspect
HERE_PATH = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(HERE_PATH)

from meteogram_plotting import clouds_plotter, utci_plotter, utci_contr_plotter, hourly_interpolated_data, hourly_dates, storm, wind_plotter, rain_plotter

from tweet_meteogram import tweet_meteogram

tweet_intro_tokyo = "Heat forecast for #TokyoOlympics: "
tweet_intro_sapporo = "Heat forecast for Sapporo (Marathon & race walking at #TokyoOlympics): "

tweet_outro = ". Stay hydrated! #Olympics #OlympicGames #Tokyo2020"

def create_meteogram(location,timestr):
    
    print("Creating meteogram...")
    
    init_date = timestr   # MMDDHHMM
    folder = "2021"+timestr[:4]
    target_date_format = "????????1"
    stream = "A1E"           # A1E for ensemble forecast ENS
    product = "pf"           # perturbed forecast, i.e. 50 members, "cf" for control
    which_location = location
    
    tokyo_npoints = 961              # distinguish tokyo vs sapporo by size, 31x31 grid points around tokyo
    sapporo_npoints = 1024           # vs 32x32 grid points around sapporo

    # dimensions
    n_time_steps = 145       # number of time steps of full forecast
    n_members = 50           # number of ensemble members

    # lat lon coordinates of Tokyo Stadium and Central Sapporo
    tokyo_latlon = (35.678,139.715)
    sapporo_latlon = (43.07,141.35)

    if which_location == "tokyo":
        latlon = tokyo_latlon
        npoints = tokyo_npoints
        title_label = "Tokyo, Olympic Stadium"
        tweet_intro = tweet_intro_tokyo
    elif which_location == "sapporo":
        latlon = sapporo_latlon
        npoints = sapporo_npoints
        title_label = "Sapporo, Central Station"
        tweet_intro = tweet_intro_sapporo
    
    # GET ENSEMBLE FILES
    ensemble_files = np.sort(glob.glob(os.path.join(folder,stream+init_date+target_date_format)))
    ds0 = xr.open_dataset(ensemble_files[0],engine="cfgrib",filter_by_keys={'numberOfPoints': npoints, 'dataType': product})


    # preallocate
    t2m = np.zeros((n_members,n_time_steps))
    u10 = np.zeros((n_members,n_time_steps))
    v10 = np.zeros((n_members,n_time_steps))
    sp = np.zeros((n_members,n_time_steps))
    tp = np.zeros((n_members,n_time_steps))
    d2m = np.zeros((n_members,n_time_steps))
    ssrd = np.zeros((n_members,n_time_steps))
    fdir = np.zeros((n_members,n_time_steps))
    lcc = np.zeros((n_members,n_time_steps))
    mcc = np.zeros((n_members,n_time_steps))
    hcc = np.zeros((n_members,n_time_steps))
    strd = np.zeros((n_members,n_time_steps))
    ssr = np.zeros((n_members,n_time_steps))
    strr = np.zeros((n_members,n_time_steps))

    timevec = [datetime.datetime.utcfromtimestamp(ds0.valid_time.data.tolist()/1e9),]*n_time_steps

    # read in file by file
    for istep,ensemble_file in enumerate(ensemble_files):
        ds = xr.open_dataset(ensemble_files[istep],engine="cfgrib",filter_by_keys={'numberOfPoints': npoints, 'dataType': product})
        latidx = np.argmin(np.abs(np.array(ds.latitude) - latlon[0]))
        lonidx = np.argmin(np.abs(np.array(ds.longitude) - latlon[1]))
        t2m[:,istep] = ds.t2m[:,latidx,lonidx]
        u10[:,istep] = ds.u10[:,latidx,lonidx]
        v10[:,istep] = ds.v10[:,latidx,lonidx]
        sp[:,istep] = ds.sp[:,latidx,lonidx]
        tp[:,istep] = ds.tp[:,latidx,lonidx]
        d2m[:,istep] = ds.d2m[:,latidx,lonidx]
        ssrd[:,istep] = ds.ssrd[:,latidx,lonidx]
        strd[:,istep] = ds.strd[:,latidx,lonidx]
        ssr[:,istep] = ds.ssr[:,latidx,lonidx]
        strr[:,istep] = ds.str[:,latidx,lonidx]
        fdir[:,istep] = ds.fdir[:,latidx,lonidx]
        lcc[:,istep] = ds.lcc[:,latidx,lonidx]
        mcc[:,istep] = ds.mcc[:,latidx,lonidx]
        hcc[:,istep] = ds.hcc[:,latidx,lonidx]
    
        timevec[istep] = datetime.datetime.utcfromtimestamp(ds.valid_time.data.tolist()/1e9)

    
    print("GRIB data loaded.")
    
    # wind speed
    wspd = np.sqrt(u10**2 + v10**2).clip(0,16.99)
    
    sec_in_step = [(timevec[i+1]-timevec[i]).seconds for i in range(len(timevec)-1)]
    hrs_in_step = [(timevec[i+1]-timevec[i]).seconds/3600 for i in range(len(timevec)-1)]
    
    # Solar zenith angle, integrated
    cosszai = np.zeros(n_time_steps-1)

    for i in range(1,n_time_steps):
        cosszai[i-1] = thermofeel.calculate_cos_solar_zenith_angle_integrated(lat=latlon[0],
                                            lon=latlon[1],
                                            y=timevec[i].year,
                                            m=timevec[i].month,
                                            d=timevec[i].day,
                                            h=timevec[i].hour,
                                            base=0,
                                            step=hrs_in_step[i-1])
        

    # mean radiant temperature
    mrt = np.zeros_like(t2m[:,1:])

    for j in range(n_members):
        mrt[j,:] = thermofeel.calculate_mean_radiant_temperature(
                                         np.diff(ssrd[j,:])/sec_in_step,
                                         np.diff(ssr[j,:])/sec_in_step,
                                         np.diff(fdir[j,:])/sec_in_step,
                                         np.diff(strd[j,:])/sec_in_step,
                                         np.diff(strr[j,:])/sec_in_step,
                                         cosszai)

    # relative humidity, saturation water vapour
    rh_pc = thermofeel.calculate_relative_humidity_percent(t2m, d2m)
    ehPa = thermofeel.calculate_saturation_vapour_pressure(t2m) * rh_pc / 100.0
    ehPa_50perc = thermofeel.calculate_saturation_vapour_pressure(t2m) / 2
    
    
    # UTCI
    utci = thermofeel.calculate_utci(t2m[:,1:], wspd[:,1:], mrt, ehPa[:,1:])

    # contributions by setting to defaults:
    # wind speed = 0.5m/s
    # relative humidity = 50%
    # mean radiant temperature = air temperature
    utci_wind = thermofeel.calculate_utci(t2m[:,1:], wspd[:,1:], t2m[:,1:], ehPa_50perc[:,1:])
    utci_rh = thermofeel.calculate_utci(t2m[:,1:], 0.5, t2m[:,1:], ehPa[:,1:])
    utci_mrt = thermofeel.calculate_utci(t2m[:,1:], 0.5, mrt, ehPa_50perc[:,1:])
    
    # TIME
    utc = 9                                     # japan time zone

    # start day of forecast
    tstart = datetime.datetime(2021,
                           int(init_date[:2]),
                           int(init_date[2:4]),
                           0)+datetime.timedelta(hours=24)

    # end day of the forecast 10days later
    tend = tstart+datetime.timedelta(hours=10*24-1)  
    utc_td = datetime.timedelta(hours=utc)      # convert utc from int to timedelta
    time = [t+utc_td for t in timevec]          # local time

    time_sub = [True,]
    i_sub = 0
    for i in range(1,len(time)):
        if (time[i]-time[i_sub]) >= datetime.timedelta(hours=6):
            use_i = True
            i_sub = i
        else:
            use_i = False
    
        time_sub.append(use_i)

    time_sub = np.array(time_sub)
    timevec_sub = np.array(time)[time_sub]
    
    ## plot
    fig = plt.figure(figsize=(10,6),facecolor="white")
    all_ax = gridspec.GridSpec(5, 1, height_ratios=[0.8,.8,.8,6,3],hspace=0)
    cloud_ax = plt.subplot(all_ax[0])
    rain_ax1 = plt.subplot(all_ax[1])
    rain_ax2 = plt.subplot(all_ax[2])
    ax = plt.subplot(all_ax[3])
    cont_ax = plt.subplot(all_ax[4])

    tend1h = tend-datetime.timedelta(hours=1)
    dt = datetime.timedelta(hours=0.9) #used to shift symbols left/right

    # RAIN&CLOUDS
    rain_plotter(rain_ax1,rain_ax2,time,tp,hrs_in_step)
    clouds_plotter(cloud_ax,time,hcc,mcc,lcc)

    # LABEL FORMATTING
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'+u'\N{DEGREE SIGN}'+'C'))
    cont_ax.yaxis.set_major_formatter(FormatStrFormatter('%d'+u'\N{DEGREE SIGN}'+'C'))
    cont_ax.xaxis.set_minor_locator(HourLocator(np.arange(6, 25, 6)))    # minor
    rain_ax1.xaxis.set_minor_locator(HourLocator(np.arange(6, 25, 6)))    # minor
    rain_ax2.xaxis.set_minor_locator(HourLocator(np.arange(6, 25, 6)))    # minor
    # cont_ax.xaxis.set_major_locator(HourLocator([0]))    # minor
    cont_ax.xaxis.set_minor_formatter(DateFormatter("%Hh"))
    cont_ax.get_xaxis().set_tick_params(which='minor', direction='out',pad=2,labelsize=6)
    cont_ax.xaxis.set_major_formatter(DateFormatter(" %a\n %d %b"))
    plt.setp(cont_ax.get_xticklabels(), ha="left")
    cont_ax.get_xaxis().set_tick_params(which='major', direction='out',pad=10,labelsize=10)
    ax.set_xticklabels([])

    for axi in [rain_ax1,rain_ax2,ax,cont_ax]:
        axi.grid(alpha=.2)

    rain_ax1.set_xticklabels([])
    rain_ax1.set_yticks([])
    rain_ax2.set_xticklabels([])
    rain_ax2.set_yticks([])
    rain_ax1.tick_params(axis="x",direction="in")
    rain_ax2.tick_params(axis="x",direction="in")
    rain_ax2.xaxis.tick_top()

    cloud_ax.set_xticks([])
    cloud_ax.set_yticks([])
    cont_ax.set_yticks([-10,-5,0,5,10])
    cont_ax.yaxis.set_label_position("right")
    cont_ax.yaxis.tick_right()

    # LIMITS
    ax.set_xlim(tstart,tend)
    cont_ax.set_xlim(tstart,tend)
    rain_ax1.set_xlim(tstart,tend)
    rain_ax2.set_xlim(tstart,tend)
    cloud_ax.set_xlim(tstart,tend)
    ax.set_ylim(np.percentile(utci,[0.1,99.9]))
    cont_ax.set_ylim(-9,9)
    rain_ax1.set_ylim(-0.1,1.1)
    rain_ax2.set_ylim(-0.1,1.1)
    cloud_ax.set_ylim(0, 1)

    # TITLES
    cloud_ax.set_title(title_label, loc="left",fontweight="bold")
    ax.text(tend1h,ax.get_ylim()[1]-.7,"Feels-like temperature"+" "*4,va="top",ha="right",fontweight="bold")
    
    # cont_ax.set_ylabel("Contribution")
    cont_ax.text(tend1h,cont_ax.get_ylim()[1]-.5,"perceived\nwarmer",
                 va="top",fontsize=8,fontweight="bold",ha="right",zorder=10)
    cont_ax.text(tend1h,cont_ax.get_ylim()[0]+.5,"perceived\n     colder",
                 fontsize=8,va="bottom",ha="right",fontweight="bold",zorder=10)
    
    # RAIN TITLES
    rain_ax1.text(tend1h,rain_ax1.get_ylim()[1]-.15,"worst-case rain",
                 va="top",fontsize=8,fontweight="bold",ha="right",zorder=10)
    rain_ax2.text(tend1h,rain_ax2.get_ylim()[1]-.15,"best-case rain",
                 fontsize=8,va="top",ha="right",fontweight="bold",zorder=10)
    
    # DATA PLOTTING
    utci_plotter(ax, time[1:], utci, tend)
    wind_plotter(ax, timevec_sub, wspd[:,time_sub])
    
    # # CONTRIBUTIONS TO UTCI
    cont_ax.plot(time,np.zeros_like(t2m[1,:]),"k",zorder=11)
    utci_contr_plotter(cont_ax,time[1:],utci_wind,t2m[:,1:]-273.15,color="steelblue",label="Wind")
    utci_contr_plotter(cont_ax,time[1:],utci_mrt,t2m[:,1:]-273.15,color="C1",alpha=0.007,label="Sunshine")
    utci_contr_plotter(cont_ax,time[1:],utci_rh,t2m[:,1:]-273.15,color="lightseagreen",alpha=0.01,label="Humidity")
    cont_ax.legend(loc=(0,0),ncol=3,framealpha=0,title="Contribution to feels-like temperature"+" "*17)
    
    # ECMWF logo
    logoax = fig.add_axes([0.783, 0.815, 0.16, 0.18], anchor='NE')
    logoax.imshow(plt.imread('ECMWF-logo.png'))
    logoax.axis('off')
    
    plt.tight_layout()
    outputpngpath = "meteograms/"+which_location+"_"+init_date+".png"
    plt.savefig(outputpngpath,dpi=200,transparent=False)
    print("Meteogram created for "+which_location+" at "+init_date)
    
    # next three days
    next3days = [tstart + datetime.timedelta(hours=i*24) for i in [0,1,2]]
    next3days_wkday = [day.strftime("%a") for day in next3days]
    
    # median of daily max
    dailymax = np.zeros(len(next3days),dtype=np.int64)
    for iday,day in enumerate(next3days):
        # all time steps of given day
        idx = [t > next3days[iday] and t < (next3days[iday]+datetime.timedelta(hours=24)) for t in time[1:]]
        dailymax[iday] = np.round(np.median(np.max(utci[:,idx],axis=1)))
        
    # temperature string, e.g. "Mon 31˚C, Tue 31˚C and Wed 35˚C"
    tempstring = "{day1} {max1}˚C, {day2} {max2}˚C and {day3} {max3}˚C".format(day1=next3days_wkday[0],
                                                              day2=next3days_wkday[1],
                                                              day3=next3days_wkday[2],
                                                              max1=dailymax[0],
                                                              max2=dailymax[1],
                                                              max3=dailymax[2])

    tweet_text = tweet_intro+tempstring+tweet_outro
    twitter_feedback = tweet_meteogram(outputpngpath,tweet_text)
    print("Meteogram tweeted, status code "+str(twitter_feedback))
    

    