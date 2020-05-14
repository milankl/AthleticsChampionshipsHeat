import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib.dates as mdates
from matplotlib import gridspec
from matplotlib.dates import WeekdayLocator, DayLocator, HourLocator, DateFormatter, drange, date2num, num2date
from matplotlib.ticker import FormatStrFormatter
from pythermalcomfort.models import utci

path = "git/AthleticsChampionshipsHeat/data/"

ncfile = nc.Dataset(path+"berlin2009.nc")
time_hours = ncfile.variables["time"][:]
T = ncfile.variables["t2m"][:]-273.15       # 2m temperature [°C]
Td = ncfile.variables["d2m"][:]-273.15      # 2m dew point temperature [°C]
u = ncfile.variables["u10"][:]              # 10m wind u-component [m/s]
v = ncfile.variables["v10"][:]              # 10m wind v-component [m/s]

# daily accumulated radiation fluxes
ssr_acc = ncfile.variables["ssr"][:]        # net shortwave radiation [J/m^2]
str_acc = ncfile.variables["str"][:]        # net longwave radiation [J/m^2]
ssrd_acc = ncfile.variables["ssrd"][:]      # downward shortwave radiation [J/m^2]
strd_acc = ncfile.variables["strd"][:]      # downward longwave radiation [J/m^2]
lats = ncfile.variables["latitude"][:]
lons = ncfile.variables["longitude"][:]
ncfile.close()

# read the solar direct downward fdir from file

ncfile = nc.Dataset(path+"berlin2009_fdir_i.nc")
fdir = ncfile.variables["fdir"][:]          # downward direct shortwave radiation [J/m^2]
ncfile.close()

## convert time
t0 = datetime.datetime(1900,1,1)    # netCDF is in hours since
utc2 = datetime.timedelta(hours=2)  # Berlin summer time is UTC+2
time_utc = [t0+datetime.timedelta(hours=np.float64(i)) for i in time_hours]
time = [t+utc2 for t in time_utc]

tstart = datetime.datetime(2009,8,15,0)    # start at 8am berlin time on Aug 6
tend = datetime.datetime(2009,8,23,18)    # end at 6pm berlin time on Aug 12

julianday0 = datetime.datetime(2009,1,1)
juliandays = np.array([(t-julianday0).days for t in time_utc])
hours = np.array([(t-julianday0).seconds/3600 for t in time_utc])

## relative humidity from dew point
# Magnus formula
A1 = 17.625
B1 = 243.04
rh = 100*(np.exp(A1*Td / (B1+Td)) / np.exp(A1*T / (B1+T)))

## wind speed
wspd = np.sqrt(u**2 + v**2).clip(0.5,17)

## convert radiation from accumulated J/m^2 to W/m^2
one_hour = 3600 # seconds

def gradient(acc):
    out = np.zeros_like(acc[:-1,:,:])

    for i,hrs in enumerate(hours[:-1]):
        if hrs == 0:    # no diff across midnight
            out[i,:,:] = acc[i+1,:,:]/one_hour
        else:           # normal np.diff case
            out[i,:,:] = (acc[i+1,:,:] - acc[i,:,:])/one_hour

    return out

ssr = gradient(ssr_acc)
str = gradient(str_acc)
ssrd = gradient(ssrd_acc)
strd = gradient(strd_acc)

## radiation decomposition
# ssrd is diffusive+direct radiation, subtract the direct fdir to obtain diffusive only
ssrd_diffusive = ssrd - fdir[1:,:,:]/3600

# ssrd is downward, ssr is net (down+up), subtract for upward only
ssru = ssrd - ssr

# strd is downward, str is net (down+up), subtract for updward only
stru = strd - str

## solar coordinates (see Di Napoli, 2020 for more information)
def cosd(x):
    return np.cos(np.pi*x/180)
    
def sind(x):
    return np.sin(np.pi*x/180)

def tand(x):
    return np.tan(np.pi*x/180)
    
def acosd(x):
    return 180/np.pi*np.arccos(x)

def g(JD,hr):
    # angular fraction of the year in degrees
    # JD: Julian day number 0 for January 1st, 0 UTC
    # hr: hour of day (UTC)
    return 360/365.25*(JD+hr/24)

def h(hr,lambdaa,TC):
    # solar hour angle (local), NOAA 1997
    # hr: hour of day (UTC)
    # lambdaa: longitude (degrees)
    # TC: time correction
    return (hr-12)*15 + lambdaa + TC

def TC(g):
    # time correction, NOAA 1997
    return 0.004297 + 0.107029*cosd(g) + \
        -1.837877*sind(g) - 0.837378*cosd(2*g) + \
        -2.340475*sind(2*g)

def delta(g):
    # solar declination angle, Spencer 1971
    return 180/np.pi*(0.006918-0.399912*cosd(g) +
        0.070257*sind(g) - 0.006758*cosd(2*g) +
        0.000907*sind(2*g) - 0.002697*cosd(3*g) +
        0.001480*sind(3*g))
        
        
time_g = g(juliandays,hours)
declination = delta(time_g)
time_correction = TC(time_g)
hour_angle = np.array([h(hours,lon,time_correction) for lon in lons])

# hour_angle_sunrise
cos_h0 = np.array([-tand(declination)*tand(lat) for lat in lats])

# clip hour angles to sunrise-sunset
hminmax = np.zeros_like(fdir)
nosun = np.zeros_like(fdir)

for ilat,_ in enumerate(lats):
    for ilon,_ in enumerate(lons):
        for i,cos_h0i in enumerate(cos_h0[ilat,:]): 
            if hour_angle[ilon,i] < -acosd(cos_h0i):    # before sunrise
                hminmax[i,ilat,ilon] = -acosd(cos_h0i)
                nosun[i,ilat,ilon] = 1
            elif hour_angle[ilon,i] > acosd(cos_h0i):   # after sunset
                hminmax[i,ilat,ilon] = acosd(cos_h0i)
                nosun[i,ilat,ilon] = 1
            else:                                       # in between
                hminmax[i,ilat,ilon] = hour_angle[ilon,i]

# solar zenith angle
cos_phi0 = np.zeros_like(fdir[1:,:,:])     # integration over timestep
cos_phi = np.zeros_like(fdir)               

for ilat,lat in enumerate(lats):
    for ilon,_ in enumerate(lons):
        for i,decl in enumerate(declination[:-1]):
            # no sun in the time interval
            if nosun[i,ilat,ilon] == 1 and nosun[i+1,ilat,ilon] == 1:
                cos_phi0[i,ilat,ilon] = 0
            else:
                cos_phi0[i,ilat,ilon] = sind(decl)*sind(lat) + \
            1/(hminmax[i+1,ilat,ilon]-hminmax[i,ilat,ilon])* \
            cosd(decl)*cosd(lat)*(sind(hminmax[i+1,ilat,ilon])-
                                    sind(hminmax[i,ilat,ilon]))
            
            cos_phi[i,ilat,ilon] = sind(decl)*sind(lat) + \
                                cosd(decl)*cosd(lat)*cosd(hour_angle[ilon,i])

## elevation angle and projection factor
elevation_angle = 90-acosd(cos_phi[1:,:,:])
fp = 0.308*cosd(elevation_angle*(0.998-elevation_angle**2/50000))

## COMPUTE MRT
Istar = np.zeros_like(cos_phi0)

for i in range(len(time_utc)-1):
    for j in range(len(lons)):
        for k in range(len(lats)):
            if cos_phi0[i,j,k] > 0:
                # convert to W/m^2 from J/m^2
                Istar[i,j,k] = fdir[i+1,j,k]/cos_phi0[i,j,k]/one_hour


alpha_ir = 0.7      # absorption coefficient of a human body
eps_p = 0.97        # emissivity coefficient of a human body
sigma = 5.67e-8     # stefan-boltzmann constant
fa = 0.5            # angle factors

MRT = ((fa*strd + fa*stru + alpha_ir/eps_p*(fa*(ssrd_diffusive+ssru)+fp*Istar))/sigma)**0.25-273.15

## COMPUTE UTCI
UTCI = np.zeros_like(MRT)
UTCIwind = np.zeros_like(MRT)
UTCIrh = np.zeros_like(MRT)
UTCImrt = np.zeros_like(MRT)

for i in range(UTCI.shape[0]):
    for j in range(UTCI.shape[1]):
        for k in range(UTCI.shape[2]):
            if MRT[i,j,k] > T[i+1,j,k]+70:
                MRT[i,j,k] = T[i+1,j,k]+70
            
            UTCI[i,j,k] = utci(T[i+1,j,k],MRT[i,j,k],wspd[i+1,j,k],rh[i+1,j,k])
            UTCIwind[i,j,k] = utci(T[i+1,j,k],T[i+1,j,k],wspd[i+1,j,k],50)
            UTCIrh[i,j,k] = utci(T[i+1,j,k],T[i+1,j,k],0.5,rh[i+1,j,k])
            UTCImrt[i,j,k] = utci(T[i+1,j,k],MRT[i,j,k],0.5,50)
            
## plot
fig,ax = plt.subplots(1,1,figsize=(10,4))

# LABEL FORMATTING
ax.yaxis.set_major_formatter(FormatStrFormatter('%d'+u'\N{DEGREE SIGN}'+'C'))
ax.xaxis.set_minor_locator(HourLocator(np.arange(6, 25, 6)))    # minor
ax.xaxis.set_minor_formatter(DateFormatter("%Hh"))
ax.get_xaxis().set_tick_params(which='minor', direction='out',pad=2,labelsize=6)
ax.xaxis.set_major_formatter(DateFormatter(" %a\n %d %b"))
plt.setp(ax.get_xticklabels(), ha="left")
ax.get_xaxis().set_tick_params(which='major', direction='out',pad=10,labelsize=10)
ax.grid(alpha=0.2)

# LIMITS
ax.set_xlim(tstart,tend)
ax.set_ylim(2,45)

# TITLES
ax.set_title("Berlin 2009", loc="left",fontweight="bold")
ax.set_ylabel("Temperature")

# DATA PLOTTING
l1, = ax.plot(time,T[:,2,2],"k")
l2, = ax.plot(time[1:],UTCI[:,2,2],"C3",lw=2)
f1 = ax.fill_between(time[1:],UTCIwind[:,2,2],T[1:,2,2],color="C0",alpha=0.3)
f2 = ax.fill_between(time[1:],UTCIrh[:,2,2],T[1:,2,2],color="C2",alpha=0.3)
f3 = ax.fill_between(time[1:],UTCImrt[:,2,2],T[1:,2,2],color="C1",alpha=0.2)
f4 = ax.fill_between(time[1:],np.percentile(UTCI,10,axis=(1,2)),np.percentile(UTCI,90,axis=(1,2)),color="C3",alpha=0.5)

labels = ["Air temperature","UTCI"+10*" "+"UTCI-contribution from:","Wind","Humidity","Radiation"]
ax.legend([l1,(l2,f4),f1,f2,f3],labels,loc=3,ncol=5)

plt.tight_layout()
plt.show()