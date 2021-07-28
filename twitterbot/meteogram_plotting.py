import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime

from matplotlib.dates import date2num, num2date
from matplotlib.patches import Ellipse

from droplet import droplet
from wind_sock import wind_sock

def clouds_plotter(axis,dates,highcloud,midcloud,lowcloud):
    """ Adds the different types of clouds to a given axis."""
    # add sun (and moon?)
    idate = datetime.datetime(dates[0].year,dates[0].month,dates[0].day,12)
    while idate < dates[-1]:
        #sun = Circle((dates[t], 0.5), 0.2, color='yellow', zorder=0)
        sun = Ellipse((idate, 0.5), 0.4/1.7, 0.5, angle=0.0, color='yellow', zorder=0)
        axis.add_artist(sun)
        idate = idate + datetime.timedelta(1)
        
    highcloud = np.median(highcloud,axis=0)
    midcloud = np.median(midcloud,axis=0)
    lowcloud = np.median(lowcloud,axis=0)
        
    totalcloud=(highcloud+midcloud+lowcloud)/3.
    totalcloudhalf=totalcloud/2.
    lowerbound=-totalcloudhalf+0.5
    upperbound=totalcloudhalf+0.5

    # don't plot clouds where totalcloud <= e.g. 0.05
    threshold=0.03

    # highcloud light grey, lowcloud dark grey
    axis.fill_between(dates, y1=lowerbound, y2=upperbound, color='0.95',zorder=1, alpha=0.8, edgecolor='none',where=totalcloud>=threshold)
    axis.fill_between(dates, y1=lowerbound, y2=upperbound-highcloud/3., color='0.7',zorder=2, alpha=0.6, edgecolor='none',where=totalcloud>=threshold)
    axis.fill_between(dates, y1=lowerbound, y2=lowerbound+lowcloud/3.,  color='0.4',zorder=3, alpha=0.3, edgecolor='none',where=totalcloud>=threshold)
    axis.set_facecolor('lightskyblue')
    
    
def utci_plotter(ax, dates, temperature, tend, color='white',alpha=0.07):
    
    temp = hourly_interpolated_data(temperature,dates)
    
    # sort ensemble members for plotting
    temp.sort(axis=0)
    temp = temp[3:-3,:]   # ~5-95% confidence interval
    
    # these temperatures will be associated with the lower and upper end of the colormap
    clev = [10,35]
    tmin,tmax = (0,60)   # absurdly low and high temperatures  
    cmap = "rainbow"
    cmat = [[clev[0],clev[0]],[clev[1],clev[1]]]
    cmat_high = clev[1]*np.ones((2,2))
    cmat_low = clev[0]*np.ones((2,2))
    
    ax.contourf([dates[0],dates[-1]],clev,cmat,128,cmap=cmap)
    ax.contourf([dates[0],dates[-1]],[clev[1],tmax],cmat_high,128,vmin=clev[0],vmax=clev[1],cmap=cmap)
    ax.contourf([dates[0],dates[-1]],[tmin,clev[0]],cmat_low,128,vmin=clev[0],vmax=clev[1],cmap=cmap)
    
    numtime = date2num(hourly_dates(dates))
    
    n_tsteps = len(numtime)
    n_ens_members = temp.shape[0]
    ylim = ax.get_ylim()
    
    # plot first half of ensemble members from the bottom, the first without transparency
    ax.fill_between(numtime,ylim[0]*np.ones(n_tsteps),temp[0,:],facecolor=color,alpha=1.)
    for i in range(1,n_ens_members//2):
        ax.fill_between(numtime,ylim[0]*np.ones(n_tsteps),temp[i,:],facecolor=color,alpha=alpha)
    
    # and the second half from the top, the last without transparency
    for i in range(n_ens_members//2,n_ens_members):
        ax.fill_between(numtime,temp[i,:],ylim[1]*np.ones(n_tsteps),facecolor=color,alpha=alpha)
 
    ax.fill_between(numtime,temp[-1,:],ylim[1]*np.ones(n_tsteps),facecolor=color,alpha=1.)
    
    # uncertainty example
    tend2h = tend-datetime.timedelta(hours=2)
    day10 = np.argmin(np.abs(numtime-date2num(tend2h)))
    ax.plot([tend2h,tend2h],[temp[0,day10],temp[-1,day10]],"k",marker="_",ms=6,alpha=.7)
    ax.text(tend2h,np.mean(temp[[0,-1],day10]),"uncertainty",rotation=90,va="center",ha="right",alpha=.7)
    
def utci_contr_plotter(ax, dates, utci_contr, temperature, color="C0", alpha=0.02,label="Contribution X"):
    
    contribution = utci_contr-temperature
    contr = hourly_interpolated_data(contribution,dates)
    
    # sort ensemble members for plotting
    contr.sort(axis=0)
    contr = contr[3:-3,:]    # ~5-95% confidence interval
        
    numtime = date2num(hourly_dates(dates))  
    n_tsteps = len(numtime)
    n_ens_members = contr.shape[0]

    for i in range(n_ens_members):
        ax.fill_between(numtime,contr[i,:],facecolor=color,alpha=alpha)
        
    # for legend only
    ax.fill_between([0,0],[0,0],facecolor=color,alpha=alpha**(1/6),label=label)
    
# interpolation of data
def hourly_interpolated_data(data, dates):
    
    numdates = date2num(dates)
    numdates_hourly = np.arange(numdates[0], numdates[-1], numdates[1]-numdates[0])
    data_hourly = np.empty((data.shape[0],len(numdates_hourly)))
    
    # loop over ensemble members
    for e in range(0, data.shape[0]):
        spline = interp1d(numdates, data[e,:], kind='cubic')
        data_hourly[e,:] = spline(numdates_hourly)
    return data_hourly

def hourly_dates(dates):
    numdates = date2num(dates)
    numdates_hourly = np.arange(numdates[0], numdates[-1], numdates[1]-numdates[0])
    return num2date(numdates_hourly)

def storm(spd):
#     threshold1 = 8.     # m/s   5 Beaufort 
#     threshold2 = 10.7   # 6 Beaufort
#     threshold3 = 13.8   # 7 Beaufort
    
    threshold1 = 6.     # m/s   4 Beaufort 
    threshold2 = 8.     # 5 Beaufort
    threshold3 = 10.7   # 6 Beaufort
    
    p_storm = np.mean(spd >= threshold1,axis=0)
    
    s_storm = np.zeros_like(p_storm)
    for i,s in enumerate(spd.T):
        st = s >= threshold1
        if np.sum(st):
            s_storm[i] = np.median(s[st])
    
    storm_strength = np.zeros((3,s_storm.shape[0]))
    storm_strength[0,:] = 1.*(np.logical_and(s_storm >= threshold1,s_storm < threshold2))
    storm_strength[1,:] = 1.*(np.logical_and(s_storm >= threshold2,s_storm < threshold3))
    storm_strength[2,:] = 1.*(s_storm >= threshold3)
     
    color=[0.71,0.12,0.12]
    
    # turn into alpha values
    colormat = np.zeros((p_storm.shape[0],4))
    colormat[:,0] = color[0]
    colormat[:,1] = color[1]
    colormat[:,2] = color[2]
    
    colormat[:,3] = p_storm
    
    return colormat,storm_strength.astype(bool)

def wind_plotter(ax,dates,wind):
    
    p_storm,storm_strength = storm(wind)
    
    windsock_weak = wind_sock(rot=50)
    windsock_medi = wind_sock(rot=70)
    windsock_stro = wind_sock(rot=90)
    
    q0,q1,q2 = storm_strength       # for readability
    
    y0,y1 = ax.get_ylim()
    
    ax.scatter([d for q,d in zip(q0,dates) if q],np.ones_like(dates)[q0]*y0+(y1-y0)*0.04,
               300,color=p_storm[q0,:],marker=windsock_weak)
    ax.scatter([d for q,d in zip(q1,dates) if q],np.ones_like(dates)[q1]*y0+(y1-y0)*0.04,
               350,color=p_storm[q1,:],marker=windsock_medi)
    ax.scatter([d for q,d in zip(q2,dates) if q],np.ones_like(dates)[q2]*y0+(y1-y0)*0.04,
               400,color=p_storm[q2,:],marker=windsock_stro)

def rain_plotter(ax1,ax2,time,tp,hrs_in_step):
    dropletpath = droplet(rot=-30)
    nlen = 140
    timesc = time[:nlen]
    
    # convert to mm/6h
    tp_diff = np.diff(tp,axis=1)/hrs_in_step*6*1000
    tp_diff = tp_diff.clip(0,np.max(tp_diff))
    tp_best = np.percentile(tp_diff,33.3,axis=0)
    tp_worst = np.percentile(tp_diff,66.6,axis=0)
    
    # used for zero transparency
    max_precip = 5   # mm/6h
    
    rgba_best = [(0.12,0.46,0.7,np.clip(tp_best[i]/max_precip,0,1)) for i in range(nlen)]
    rgba_worst = [(0.12,0.46,0.7,np.clip(tp_worst[i]/max_precip,0,1)) for i in range(nlen)]

    ax1.scatter(timesc[::4],[.2,]*(nlen//4),80,rgba_worst[::4],marker=dropletpath)
    ax1.scatter(timesc[1::4],[.6,]*(nlen//4),80,rgba_worst[1::4],marker=dropletpath)
    ax1.scatter(timesc[2::4],[.4,]*(nlen//4),80,rgba_worst[2::4],marker=dropletpath)
    ax1.scatter(timesc[3::4],[.8,]*(nlen//4),80,rgba_worst[3::4],marker=dropletpath)
    
    ax2.scatter(timesc[::4],[.2,]*(nlen//4),80,rgba_best[::4],marker=dropletpath)
    ax2.scatter(timesc[1::4],[.6,]*(nlen//4),80,rgba_best[1::4],marker=dropletpath)
    ax2.scatter(timesc[2::4],[.4,]*(nlen//4),80,rgba_best[2::4],marker=dropletpath)
    ax2.scatter(timesc[3::4],[.8,]*(nlen//4),80,rgba_best[3::4],marker=dropletpath)
