import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib.dates as mdates
from matplotlib import gridspec
from matplotlib.dates import WeekdayLocator, DayLocator, HourLocator, DateFormatter, drange, date2num, num2date
from matplotlib.ticker import FormatStrFormatter

path = "git/AthleticsChampionshipsHeat/data/berlin2009.nc"

ncfile = nc.Dataset(path)
time_hours = ncfile.variables["time"][:]
T = ncfile.variables["t2m"][:]-273.15
Td = ncfile.variables["d2m"][:]-273.15
u = ncfile.variables["u10"][:]
v = ncfile.variables["v10"][:]
radc = ncfile.variables["ssrd"][:]
ncfile.close()

## convert time
t0 = datetime.datetime(1900,1,1)    # netCDF is in hours since
utc2 = datetime.timedelta(hours=2)  # Berlin summer time is UTC+2
time = [t0+datetime.timedelta(hours=np.float64(i))+utc2 for i in time_hours[:]]

tstart = datetime.datetime(2009,8,15,0)    # start at 8am berlin time on Aug 6
tend = datetime.datetime(2009,8,23,18)    # end at 6pm berlin time on Aug 12

## relative humidity from dew point

# Magnus formula
A1 = 17.625
B1 = 243.04
rh = 100*(np.exp(A1*Td / (B1+Td)) / np.exp(A1*T / (B1+T)))

## wind speed
wspd = np.sqrt(u**2 + v**2)

## radiation
rad = np.diff(radc,axis=0)
rad[rad < 0] = 0
radmax = rad.max()

radr = rad/radmax*100
## plot

fig = plt.figure(figsize=(10,4))
all_ax = gridspec.GridSpec(4, 1, height_ratios=[1,1,6,1],hspace=0)
rad_ax = plt.subplot(all_ax[0])
humid_ax = plt.subplot(all_ax[1])
temp_ax = plt.subplot(all_ax[2])
wind_ax = plt.subplot(all_ax[3])

# TIME FORMATTING
wind_ax.yaxis.set_major_formatter(FormatStrFormatter('%d'+u'\N{DEGREE SIGN}'+'C'))
wind_ax.xaxis.set_minor_locator(HourLocator(np.arange(6, 25, 6)))    # minor
wind_ax.xaxis.set_minor_formatter(DateFormatter("%Hh"))
wind_ax.get_xaxis().set_tick_params(which='minor', direction='out',pad=2,labelsize=6)
#ax.xaxis.set_major_locator(WeekdayLocator(byweekday=range(5)))
wind_ax.xaxis.set_major_formatter(DateFormatter(" %a\n %d %b"))
plt.setp(wind_ax.get_xticklabels(), ha="left")
wind_ax.get_xaxis().set_tick_params(which='major', direction='out',pad=10,labelsize=10)

# TEMP AX FORMATTING
temp_ax.grid(alpha=0.2)

# HUMID AX FORMATTING
humid_ax.yaxis.tick_right()
humid_ax.yaxis.set_major_formatter(FormatStrFormatter(r'%d%%'))
humid_ax.yaxis.set_label_coords(-0.05,0)

# WIND SPEED AXES FORMATTING
wind_ax.yaxis.tick_right()
wind_ax.yaxis.set_major_formatter(FormatStrFormatter(r'%d ms$^{-1}$'))
wind_ax.yaxis.set_label_coords(-0.05,0)
wind_ax.set_yticks([0,6])

# RAD AXES FORMATTING
rad_ax.yaxis.tick_left()
rad_ax.yaxis.set_major_formatter(FormatStrFormatter(r'%d%%'))
rad_ax.yaxis.set_label_coords(1.05,0.2)
rad_ax.set_yticks([100])
#rad_ax.set_yticklabels(r"0%%")

# TICKS AND LABELS ONLY AT BOTTOM PANEL
rad_ax.set_xlabel("")
rad_ax.set_xticks([])
humid_ax.set_xlabel("")
humid_ax.set_xticks([])
temp_ax.set_xlabel("")
temp_ax.set_xticks([])

# AXES LIMITS
for ax in [rad_ax,humid_ax,temp_ax,wind_ax]:
    ax.set_xlim(tstart,tend)

humid_ax.set_ylim(30,100)
wind_ax.set_ylim(0,6)
temp_ax.set_ylim(11,33)
rad_ax.set_ylim(0,100)

# TITLES
rad_ax.set_title("Berlin 2009", loc="left",fontweight="bold")
temp_ax.set_ylabel("Air temperature")
humid_ax.set_ylabel("Relative\nhumidity",rotation=0)
wind_ax.set_ylabel("Wind\nspeed",rotation=0)
rad_ax.set_ylabel("Solar\nradiation",rotation=0)

# DATA PLOTTING
temp_ax.plot(time,T[:,2,2],"C1")
temp_ax.fill_between(time,np.percentile(T,10,axis=(1,2)),np.percentile(T,90,axis=(1,2)),color="C1",alpha=0.5)

humid_ax.fill_between(time,rh[:,2,2],color="C0",alpha=.4)
humid_ax.fill_between(time,np.percentile(rh,10,axis=(1,2)),np.percentile(rh,90,axis=(1,2)),color="C0")
humid_ax.plot(time,rh[:,2,2],"#3333AA")

wind_ax.plot(time,wspd[:,2,2],"C3")
wind_ax.fill_between(time,wspd[:,2,2],color="C3",alpha=.3)

rad_ax.plot(time[1:],radr[:,2,2],"C1")
rad_ax.fill_between(time[1:],np.zeros_like(radr[:,2,2]),radr[:,2,2],color="grey",alpha=.5)
rad_ax.fill_between(time[1:],radr[:,2,2],np.ones_like(radr[:,2,2]),color="#EEEE33",alpha=.8)

plt.tight_layout()
plt.show()