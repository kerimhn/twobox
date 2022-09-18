#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 09:13:40 2021

@author: elisabethengum
"""

import pandas as pd

from calculate_temp_anomalies import calculate_temp_anomalies
variabel = int(input("Vil du se på historiske data siden 1750, velg 1, hvis du vil se for siste millenium trykk 2\n>> "))
if variabel == 1:
    df=pd.read_csv('historical.csv',index_col = 0,sep=',',encoding = "latin-1")
elif variabel == 2:
    df=pd.read_csv('pmip3.csv',index_col = 0,sep=',',encoding = "latin-1")


# --------------------
# Forcing switches [n=1 off=0]
# --------------------
switch_ghg=1     # Greenhouse gas forcing on=1 off=0
switch_solar=1   # Solar forcing on=1 off=0
switch_volc=0    # Volcanic forcing on=1 off=0
switch_land=1    # Landuse forcing on=1 off=0
switch_aero=0    # Pollution particle forcing on=1 off=0

lambda_planck=-3.2    
lambda_lapse=-0.8     
lambda_water=1.8      
lambda_cloud=0.70     
lambda_albedo=0.30    
lambda_other=0.0

# --------------------
# Deep ocean
# Heat uptake efficency [Wm-2K-1]
# --------------------
# CMIP3
gamma=-0.69  # best guesses [-1 to -0.5]

lambda_sum=sum([lambda_planck,
    lambda_lapse,
    lambda_water,
    lambda_cloud,
    lambda_albedo,
    lambda_other])

df['total_forcing'] = switch_ghg*df['wmghg_data']+switch_solar*df['solar_data']+switch_volc*df['volc_data']+switch_land*df['landuse_data']+switch_aero*df['manaero_data']

radiative_forcing = df['total_forcing'].to_numpy()
idx = df.index

Ts, To = calculate_temp_anomalies(radiative_forcing, lambda_sum, gamma)

temp = pd.DataFrame(index=idx)
temp['Temperatur overflate']=Ts
temp['Temperatur blandingslag']=To
temp.plot()