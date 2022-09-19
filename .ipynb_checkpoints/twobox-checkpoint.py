#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#####################################################################
The TwoBoxmodel is a very simple box model calculating the changes in 
global SST temperature with time based on prescribed feedback factors
and heat capacity with interactions to the deep ocean.

Author: Asgeir Sorteberg,

Geophysical Institute, University of Bergen.

email: asgeir.sorteberg@gfi.uib.no

Feb. 2011

#####################################################################
"""

import numpy as np

# Inputs are the total radiative forcing, sum of feedback factors and deep ocean heat uptake
# Output is mixed layer temperature (Ts) and deep ocean temperature (To)
def calculate_temp_anomalies(radiative_forcing, lambda_sum, gamma):
    # Thickness of the mixed layer [m]
    H_MIX= 100
    # Thickness of the deep ocean [m]     
    H_DEEP=3700-H_MIX   
    # Density of water (kg m-3)
    RHO = 1000
    # Heat capacity of water (J kg-1 K-1)
    CPO = 4200       
    # Proportion of Earth's surface covered by water
    f_o=0.7
    # Effectiv heat capacity [J m-2 K-1]
    CEFF_M=f_o*H_MIX*CPO*RHO
    CEFF_D=f_o*H_DEEP*CPO*RHO
    
    Dt=365*24*60*60 # Number of seconds in a year
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
