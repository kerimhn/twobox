#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy import matlib
import math
import scipy.io
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from Orbit_funcs import *
from Orbit91 import Orbit91

def inso(t, day, lat=65, qinso=-2, author=1, res=1):
    
    '''
    t      = kilo-years before present
    day    = day of year from Jan 1st, this fraction of a KY is added to t.
             give either tuple (e.g. (1,365) for every day of the year) to compute
             a range of values, or give one specific day
             negative values or day gives output for equinoxes (not supported yet).
    lat    = lattitude in degrees [-90:90].  Would like to add lon. 
    qinso  = -2 is mean daily insolation, -1 is max daily insolation,
              else hours above input threshold. 
    author = 1 gives  Berger and Loutre (1992) solution.
             2 gives  Quinn, Tremaine, and Duncan (1991) solution,
    res    = by default 365 points are calculated per year, if different
             than one, this multiplies the resolution.
             
    returns: 
    Insolation in Watts/m^2 at top of the atmosphere; 
    Time
    

    Example: I=inso([20, 19.999],(1,365),65,-2,1);
    Gives mean daily insolation for a 2 year period from 20 to 19.999 KY BP
    at 65 degrees N using the solution of Berger and Loutre (1991). 
    
    Refs:

    Berger, A. and Loutre, M. F., "Astronomical solutions for paleoclimate
    studies over the last 3 million years", Earth Planet. Sci. Lett., v111,
    p369-382, 1992.

    Quinn, T.R. et al. "A Three Million Year Integration of the Earth's
    Orbit"  The Astronomical Journal v101, p2287-2305, 1991.

    Jonathan Levine (2001), UC Berkeley, wrote the original version of this code
    and I made some modifications.
    '''
    
    # some checks
    if type(day) not in [list,tuple,int]:
        raise ValueError('day should be of type list, tuple or int')
    elif type(day)==tuple:
        day=np.arange(day[0]-1, day[1])
    elif type(day)==int:
        day=[day]
    if type(t) == int:
        t=[t]
    if type(t) not in [list, int]:
        raise ValueError('t should be of type list or int')
    if type(lat) in [int, float]:
        lat=[lat]
    if type(lat) not in [list, int, float]:
        raise ValueError('lat should be of type list or int')
        
        
    if author==2:
        mat = scipy.io.loadmat('orbitparameters.mat')
        ECC=mat['ecc']
        oblE=mat['oblE']
        periE=mat['periE']
        Otime=mat['time'].flatten()
        w=mat['w']  
    else:
        Otime, ECC, oblE, w = Orbit91()  # gets values from Orbit91.txt
        
    I=[]
    T=[]
    L = 1365   #the solar constant, in W/m2 at 1 AU

    #interpolate to times t
    etf = CubicSpline(Otime,ECC) #eccentricity(t)
    et=etf(t)
    otf = CubicSpline(Otime,oblE) #obliquity(t)
    ot=otf(t)
    wtf = CubicSpline(Otime,w) #longitude of perihelion
    wt=wtf(t)

    for tnx in range(len(t)):   #loop over times in history
        tt=t[tnx]               #time for this iteration
        etx=et[tnx]             #eccentricity at tt
        otx=ot[tnx]             #obliquity at tt
        wtx=wt[tnx]*np.pi/180   #argument of precession at tt
        
        #model of the orbit is an ellipse with 365 points, one for
        #each day of a non-leap year.  The eccentricity, obliquity, and
        #precessional parameter are all
        #taken from data of Quinn, Tremaine, and Duncan (1991) or Berger
        #and Loutre (1991).  The x,y,z
        #coordinate system is: x points to perihelion, and z is normal to
        #the modern ecliptic
        #Note using four times the spatial resolution.
        
        [xorb, yorb]=ellipse(etx, 365*res)
        zorb=np.zeros(len(xorb))
        dist=np.sqrt(xorb**2+yorb**2+zorb**2)
        unit=np.array([xorb/dist, yorb/dist, zorb/dist]).T 
        #from perihelion, rotate counterclockwise by w degrees
        vtry=[np.cos(wtx)*unit[:,0][0]-np.sin(wtx)*unit[:,1][0], np.cos(wtx)*unit[:,1][0]+np.sin(wtx)*unit[:,0][0], unit[:,2][0]]
        #closest match is the vernal equinox day
        #dif = (unit-repmat(vtry,length(unit),1)).^2; %difference between the unit vector to each day and the vernal equinox
        dif = (unit-np.tile(vtry,(np.shape(unit)[0],1)))**2
        veq=np.where((sum(dif.T).T)==np.min((sum(dif.T)).T))[0][0]
        #now veq holds the number of days between perihelion and vernal equinox

        vernal = res*80  #vernal equinox, March 20, is 80th day of year (if not a leap year)
        #the following lines create variables with the index corresponding to day of the year,
        #not days since perihelion, as ellipse function generates
    
        pl=np.arange(len(xorb))
        if veq>vernal:
            pl[0:len(xorb)-(veq-vernal)]=np.arange(veq-vernal, len(xorb))
            pl[len(xorb)-(veq-vernal):len(xorb)]=np.arange(veq-vernal)
        if veq<vernal:
            pl[0:vernal-veq]=np.arange(len(xorb)+(veq-vernal), len(xorb))
            pl[vernal-veq:len(xorb)]=np.arange(len(xorb)+(veq-vernal))  
    
        xorb=xorb[pl]
        yorb=yorb[pl]
        zorb=zorb[pl]
        dist=dist[pl]
    
        #Calculate the day of spring equinox
    
        ## skipped for now!!##

        #get the direction of the north pole from the vernal equinox and the obliquity
        #the vector north3 is perpendicular to the vector to the Sun on March 20
        #and makes an angle of oblE with the normal to the orbit.  Further, it gives
        #the correct sense of summer (ie north3 points to Sun on June 20).
        firstpt = [-xorb[vernal],-yorb[vernal],-zorb[vernal]] #from Earth to Sun on vernal equinox
        firstpt=firstpt/np.linalg.norm(firstpt)
               
        #another vector in the orbit plane
        vec2=[xorb[vernal]-xorb[vernal+90*res], yorb[vernal]-yorb[vernal+90*res], zorb[vernal]-zorb[vernal+90*res]]
        vec2=vec2/np.linalg.norm(vec2)
        #normal to plane is given by cross(firstpt,vec2)
        angmom=np.cross(firstpt,vec2)
        angmom=angmom/np.linalg.norm(angmom)
    
        #need the vector cross(angmom,firstpt)
        vec3=np.cross(angmom,firstpt)
        north3=np.cos(otx)*angmom+np.sin(otx)*vec3
        northproj=math.atan2(north3[1], north3[0])
    
        #lattitudinal loop  
        sun=np.zeros((len(lat), len(day)))
        for lnx in range(len(lat)):
            latitude=lat[lnx]
    
            longit=np.arange(0,360/res, 1/res)
    
            #latitiude ring with 360 degreees longitude
            xring=np.cos(latitude*np.pi/180)*np.cos(longit*np.pi/180)
            yring = np.cos(latitude*np.pi/180)*np.sin(longit*np.pi/180)
            zring = np.sin(latitude*np.pi/180)*np.ones(np.shape(longit))
    
            #compute orientation of ring for date with respect to the ecliptic
            #coordinate system.  Measures the elevation of the Sun in the sky
            #at all longitudes on that ring.
    
            [xorient,yorient,zorient] = rot(xring,yring,zring,math.acos(north3[2]),np.pi/2+northproj)
            up=np.array([xorient.T, yorient.T, zorient.T])
    
            #day loop
            for dnx in range(len(day)):
                doy=day[dnx]
                if doy>365*res | doy<1:
                    fac=floor(doy/(365*res))
                    t+=fac/1000
                    doy-=fac*365*res
                if doy%1==0:
                    txorb=xorb[doy]
                    tyorb=yorb[doy]
                    tzorb=zorb[doy]
                    tdist=dist[doy]
                else:
                    txorbf=CubicSpline(np.arange(0,365*res), xorb) 
                    txorb=txorbf(doy)
                    tyorbf=CubicSpline(np.arange(0,365*res), yorb) 
                    txorb=tyorbf(doy)
                    tzorbf=CubicSpline(np.arange(0,365*res), zorb) 
                    tzorb=tzorbf(doy)
                    tdistf=CubicSpline(np.arange(0,365*res), dist) 
                    tdist=tdistf(doy)
                sundirect=np.array([-txorb,-tyorb,-tzorb])
                solarelevation=90-np.arccos((up.T).dot(sundirect.T)/tdist)*(180/np.pi)
                sunlight=(solarelevation+abs(solarelevation))/2
                #each time point represents 2 minutes.  Find the daily radiation
                #radky(tnx).d(dnx).ml(:,lnx) = (L/(tdist.^2))*(sin(sunlight*pi/180))
                if qinso==-2:
                    sun[lnx,dnx]=np.mean((L/(tdist**2))*(np.sin(sunlight*np.pi/180)))
                if qinso==-1:
                    sun[lnx,dnx]=np.max((L/(tdist**2))*(np.sin(sunlight*np.pi/180)))
                if qinso>0:
                    pl=np.where((L/(tdist**2))*(np.sin(sunlight*np.pi/180))>qinso)[0]
                    sun[lnx,dnx]=len(pl)/(res*30)  #hours above threshold. 

        I.append(sun[0])
        T.append(t[tnx]-np.array(day)/(res*365*1000))

    return I,T
