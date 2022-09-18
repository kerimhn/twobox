#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 10:55:08 2021

@author: elisabethengum
"""
import numpy as np

def calculate_temp_anomalies(radiative_forcing, lambda_sum, gamma):
    # tykkelse av blandingslaget [m]
    H_MIX= 100
    # rest of the ocean [m]     
    H_DEEP=3700-H_MIX   
    # vannets tetthet (kg m-3)
    RHO = 1000
    # spesifikk varmekapasitet for vann(J kg-1 K-1)
    CPO = 4200       
    # andel av jordens overflate dekket av vann
    f_o=0.7
    # effektiv varmekapasitet for atmosfære-hav-systemet [J m-2 K-1]
    CEFF_M=f_o*H_MIX*CPO*RHO
    CEFF_D=f_o*H_DEEP*CPO*RHO
    
    Dt=365*24*60*60
    Nt = len(radiative_forcing)
    Ts_init=0
    To_init=0
    Ts=np.array(()) 
    To=np.array(())
    for t in range(Nt):
        # --------------
        # Temperature tendencies [K/s]
        #     dTs/dt, dTs_dt
        # --------------
        
        if t==0:
            dTs_dt=(radiative_forcing[t]+(lambda_sum*Ts_init)+(gamma*(Ts_init-To_init)))/CEFF_M
            dTo_dt=-gamma*(Ts_init-To_init)/CEFF_D
        else:
            dTs_dt=(radiative_forcing[t]+(lambda_sum*Ts[t-1])+(gamma*(Ts[t-1]-To[t-1])))/CEFF_M
            dTo_dt=-gamma*(Ts[t-1]-To[t-1])/CEFF_D
        
        # ----------------------------------------------------------------------
        # Step temperature forward in time (by "Euler forward"
        # or "forward-in-time" method)
        #----------------------------------------------------------------------
        if t==0:
            Ts=np.append(Ts, Ts_init+dTs_dt*Dt)
            To=np.append(To, To_init+dTo_dt*Dt)
        else:
            Ts=np.append(Ts, Ts[t-1]+dTs_dt*Dt)
            To=np.append(To, To[t-1]+dTo_dt*Dt)
    return Ts, To