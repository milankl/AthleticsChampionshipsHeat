#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:16:32 2023

@author: navarro
"""

import cdsapi
import pandas as pd
import netCDF4 as nc
import numpy as np
import glob
import datetime
import os
import zipfile

def x_round(x):
    return round(x*4)/4


###########################################
city='daegu'
year='2011'
utc=9
Lat=35.83
Long=128.69
tstart = datetime.datetime(2011,8,27,0)     
tend = datetime.datetime(2011,9,4,0)
###########################################

tab=pd.DataFrame()
Lat=x_round(Lat)
Long=x_round(Long)
time=tstart
while time<=tend:
    c = cdsapi.Client()
    c.retrieve(
        'derived-utci-historical',
        {
            'variable': 'universal_thermal_climate_index',
            'version': '1_1',
            'product_type': 'consolidated_dataset',
            'year': str(time.year),
            'month': str("%02d" % time.month),  
            'day': str("%02d" % time.day),
            'area': [Lat-0.25, Long-0.25, Lat+0.25, Long+0.25,],
            'format': 'zip',
        },
        year+'_'+city+'.zip')
    
    with zipfile.ZipFile(year+'_'+city+'.zip',"r") as zip_ref:
        zip_ref.extractall(year+'_'+city)
    
    time=time+datetime.timedelta(hours=24)


files=glob.glob(year+'_'+city+'/*.nc')
files.sort(key=os.path.getmtime)

hour=0+utc
days=[]
for i in range((tend-tstart).days+1):
    days.append((tstart+datetime.timedelta(days=i)).day)

for i in range((tend-tstart).days+1):
    vn = nc.Dataset(files[i])
    utcis = vn.variables["utci"][:]
    utcis=np.array(utcis)
    for h in range(24):
        ligne=[]
        lignetime=tstart+datetime.timedelta(hours=hour)
        ligne.append(str(lignetime.year)+str("%02d" % lignetime.month)+str("%02d" % lignetime.day))
        ligne.append(str("%02d" % lignetime.hour)+':00')
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

tab=tab.rename(columns={0:'Date',1:'Time',2:'point0',3:'point1',4:'point2',5:'point3',
                        6:'point4',7:'point5',8:'point6',9:'point7',10:'point8',
                        11:'Mean',12:'Std'}) 

tab.to_excel(year+'_'+city+'.xlsx', index=False) 

