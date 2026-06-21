#havel_lake_exchange.py>

import warnings
warnings.simplefilter("ignore", UserWarning)

import numpy as np
import pandas as pd
import math
import datetime
import os
import matplotlib
#matplotlib.use('Agg') #-- active when on cluster
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import colors
import os
from pandas import Timestamp
import time
import glob
import calendar
import scipy
from scipy import stats
from scipy.interpolate import griddata
import multiprocessing as mp
from multiprocessing.pool import ThreadPool as Pool
import re
import matplotlib.dates as dates
import datetime as dt
from datetime import datetime as dtm
import xlsxwriter
import openpyxl
from time import sleep
import random
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.dates import DateFormatter
from pathlib import Path
import sys
from tqdm import tqdm
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
pd.options.mode.chained_assignment = None  # default='warn'

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

lake_file='data/data_lake/pegel_taeglich_ifb_2003_2023.txt'
havel_file='data/data_havel/havel_wasserstand_tw_01_11_1961.csv'
h_v_file='data/data_lake/h_v.txt' #-height-volume file

#- import model parameters from havel_lake_params.py
import havel_lake_params as params

pi=math.pi

#####################################################################################################
#####################################################################################################
def SCHI_OFF(date_calc,lake_mASL):
    return 0



#####################################################################################################
#####################################################################################################
def SCHI_DAM_REFURB2015(date_calc,lake_mASL):

    Qin_day=0
    Qout_day=0

    ##-import params from havel_lake_params.py
    havel_gauge_zero=params.havel_gauge_zero
    lake_gauge_zero=params.lake_gauge_zero
    lake_bot_elev=params.lake_bot_elev
    mud_thresh=params.mud_thresh #-mud threshold to subtract from bed levels
    schifgr_havel_bed=params.schifgr_havel_bed #- mASL western side of dam 
    schifgr_highest_elev=params.schifgr_highest_elev #- mASL highest point (midpoint) of schiffgraben (where flow into lake/havel begins) 
    schifgr_lake_bed=params.schifgr_lake_bed #-mASL mouth to lake
    d_pipe_dam=params.d_pipe_dam
    l_pipe_dam=params.l_pipe_dam
    f_pipe_dam=params.f_pipe_dam
    n_pipe_dam=params.n_pipe_dam
    schifgr_dam_bed=params.schifgr_dam_bed
    flow_occur_thresh=params.flow_occur_thresh #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 
    #-channel geometry at control section (Sec EE)
    l=params.l #- total channel length (m)
    lx=params.lx #-length of channel at control section (m) Sec EE in dwg file; distance from lake side  
    b=params.b #-channel bed width at control section (m) pipe under dam considered (rough value)
    alpha=params.alpha #-side angles of walls (degrees) Section EE = 75deg (in dwg file)
    n=params.n #-channel roughness (From Chow's book for earth channels)
    #-channel geometry at mid point (sec CC)
    l_mid=params.l_mid #- channel length of schifgr_highest_elev (m) at Sec CC (dwg file) dist. from lake side 
    alpha_mid=params.alpha_mid #side angles of walls (degrees) at lmid Section CC = 47deg (in dwg file)
    b_mid=params.b_mid #- channel bed width at midpoint (sec CC)
    l_pipe=params.l_pipe #-outflow pipe length m (hypolimnetic withdrawal)
    d_pipe=params.d_pipe #-outflow pipe diameter m
    f_pipe=params.f_pipe #-pipe friction factor
    #-settings for SCHI-PU scenario
    inflow_pump=params.inflow_pump #-inflow pump rate l/s

    pipe_flow_threshold=schifgr_dam_bed+d_pipe_dam #-pipe free/pressurized flow threshold 
    schifgr_highest_elev+=flow_occur_thresh
    schifgr_havel_bed-=mud_thresh
    schifgr_highest_elev-=mud_thresh
    schifgr_lake_bed-=mud_thresh

    if int(date_calc.date().year)>2015:
        schifgr_highest_elev-=0.15
        schifgr_lake_bed-=0.15

    df=pd.read_csv(havel_file, encoding='latin', delimiter=';')
    df['date']=pd.to_datetime(df['Datum'], format='%d.%m.%Y')
    df=df[['date','Tagesmittelwert']]
    df=df.replace(-777,np.nan)
    df=df.set_index('date')[['Tagesmittelwert']].interpolate().dropna().reset_index()
    df['havel_mASL']=df['Tagesmittelwert']/100.000+havel_gauge_zero

    ####-function to calculate flow of channel under different flow regimes:
    #-when flow regime 4 or -4 (flooded and pressurized flow regime)
    #-- 4: inflow, Y1(Havel) and Y2(Lake) constant for day; Y is water surface elev mASL
    #-- -4: outflow, Y1(Lake) and Y2(Havel) constant for day; Y is water surface elev mASL
    def flow_reg_4(y1,y2,flow_regime):
        ###-discharge calculation between two connected reservoirs via pressurized pipe
        ###- method from 
        ###- bottom slope of pipe assumed to be zero; energy line slope assumed identical to water surface diff
        ###-using Energy equation and Darcy_Wisbach
        if y1>y2:
            V=math.sqrt(\
                ((y1-y2)*d_pipe_dam*2*9.81)\
                /(f_pipe_dam*l_pipe_dam)
                )
            
            Q=V*pi*(d_pipe/2)**2
            if math.isnan(Q):
                return 0
            else:
                return 2*Q*86400 #- factor of 2 corresponds to two pipes
        else:
            return 0
    
    #-when flow regime 3 or -3 (flooded and open flow regime)
    #-- 3: inflow, Y1(Havel) and Y2(Lake) constant for day; Y is water surface elev
    #-- -3: outflow, Y1(Lake) and Y2(Havel) constant for day; Y is water surface elev        
    def flow_reg_3(y1,y2,flow_regime):
        ###- discharge calculation between two connected reservoirs through free flow in pipe under dam
        ###- Manning flow in a circular pipe; 
        ###- bottom slope assumed to be zero; energy line slope assumed identical to water surface slope
        ###- sources:
        ###---- https://pdhlibrary.com/sites/default/files/2018021-Partially%20Full%20Pipe%20Flow%20Calculations.pdf
        ###---- pdf in docs directory: Partially Full Pipe Flow Calculations.pdf
        
        #-compute hydraulic radius
        r=d_pipe_dam/2

        ##- hydraulic Radius Less than Half Full Flow (page 9 of pdf)
        if y1<=schifgr_dam_bed+d_pipe_dam/2:
            #-theta: angle of sector opening to open surface
            theta=2*math.acos((r-(y1-schifgr_dam_bed))/r)
            A=(r**2*(theta-math.sin(theta)))/2
            P=r*theta #-perimeter
            R=A/P #-hydraulic radius
        
        else: #-Hydraulic Radius More than Half Full Flow (page 13 of pdf)
            theta=2*math.acos((r-(y1-schifgr_dam_bed))/r)
            A=pi*r**2-((r**2*(theta-math.sin(theta)))/2)
            P=2*pi*r-r*theta #-perimeter
            R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        S=(y1-y2)/l_pipe_dam
        n=1+(y1/d_pipe_dam)**0.54-(y1/d_pipe_dam)**1.2 #-Variable Manning Roughness Coeff (pahe 18 of pdf)        
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400#- factor of 2 corresponds to two pipes

    #-when flow regime 1 or -1 (dry flow regime)
    #-- 1: inflow, Y1(Havel) and Y2(Lake); y is water level
    #-- -1: outflow, Y1(Lake) and Y2(Havel); y is water level
    def flow_reg_1(y1,flow_regime):
        if flow_regime==1:
            y2_bed=schifgr_lake_bed+flow_occur_thresh
        elif flow_regime==-1:
            y2_bed=schifgr_havel_bed+flow_occur_thresh

        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- l_mid distance of midpoint from lake side of channel (to calculate long. slope)
        ###- daily flow assumed uniform, energy line slope identical to bottom slope
        b=b_mid
        alpha=alpha_mid
        #- find flow depth at control section (overall uniform)
        y=y1-schifgr_highest_elev

        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        #-b_elev:bed elev at control section
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        S=(schifgr_highest_elev-y2_bed)/l
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400#- factor of 2 corresponds to two pipes

    df=df.reset_index()

    
    havel_mASL=df[df.date == date_calc]['havel_mASL'].iloc[0]
    
    ####- channel flooded
    if lake_mASL>schifgr_highest_elev and havel_mASL>schifgr_highest_elev:
        #-pressurized flow under dam
        if lake_mASL>pipe_flow_threshold or havel_mASL>pipe_flow_threshold:
            #- -4 outflow, both sides>pipe_flow_threshold, lake>havel => channel flooded, pressurized dam outflow
            if lake_mASL>havel_mASL:
                flow_regime=-4
                Qout_day=flow_reg_4(lake_mASL, havel_mASL, flow_regime)

            #- 4 outflow, both sides>pipe_flow_threshold, lake<havel => channel flooded, pressurized dam inflow
            elif lake_mASL<havel_mASL:
                flow_regime=4
                Qin_day=flow_reg_4(havel_mASL, lake_mASL, flow_regime)

        #- free flow under dam
        elif lake_mASL>schifgr_highest_elev and havel_mASL<schifgr_highest_elev:
            #- -3 outflow, both sides<=pipe_flow_threshold, lake>havel => channel flooded, free pipe outflow
            if lake_mASL>havel_mASL:
                flow_regime=-3
                Qout_day=flow_reg_3(lake_mASL, havel_mASL, flow_regime)
            elif havel_mASL>lake_mASL:
                flow_regime=3
                Qin_day=flow_reg_3(havel_mASL, lake_mASL, flow_regime)

    #- -1 outflow, lake>havel & lake>schifgr_lake_bed => channel dry, outflow
    elif lake_mASL>schifgr_highest_elev and havel_mASL<schifgr_highest_elev:
        flow_regime=-1
        Qout_day=flow_reg_1(lake_mASL, flow_regime)

    #- 1 inflow, lake<havel & Havel>schifgr_lake_bed => channel dry, inflow
    elif lake_mASL<schifgr_highest_elev and havel_mASL>schifgr_highest_elev:
        flow_regime=1
        Qin_day=flow_reg_1(havel_mASL, flow_regime)

    return Qin_day,Qout_day



#####################################################################################################
#####################################################################################################
#-existing Schiffgraben, with the dam and pipe underneath
def SCHI_DAM_SIPH_EXST(date_calc,lake_mASL): 

    Qin_day=0
    Qout_day=0

    ##-import params from havel_lake_params.py
    havel_gauge_zero=params.havel_gauge_zero
    lake_gauge_zero=params.lake_gauge_zero
    lake_bot_elev=params.lake_bot_elev
    mud_thresh=params.mud_thresh #-mud threshold to subtract from bed levels
    schifgr_havel_bed=params.schifgr_havel_bed #- mASL western side of dam 
    schifgr_highest_elev=params.schifgr_highest_elev #- mASL highest point (midpoint) of schiffgraben (where flow into lake/havel begins) 
    schifgr_lake_bed=params.schifgr_lake_bed #-mASL mouth to lake
    flow_occur_thresh=params.flow_occur_thresh #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 
    #-channel geometry at control section (Sec EE)
    l=params.l #- total channel length (m)
    lx=params.lx #-length of channel at control section (m) Sec EE in dwg file; distance from lake side  
    b=params.b #-channel bed width at control section (m) pipe under dam considered (rough value)
    alpha=params.alpha #-side angles of walls (degrees) Section EE = 75deg (in dwg file)
    n=params.n #-channel roughness (From Chow's book for earth channels)
    #-channel geometry at mid point (sec CC)
    l_mid=params.l_mid #- channel length of schifgr_highest_elev (m) at Sec CC (dwg file) dist. from lake side 
    alpha_mid=params.alpha_mid #side angles of walls (degrees) at lmid Section CC = 47deg (in dwg file)
    b_mid=params.b_mid #- channel bed width at midpoint (sec CC)
    l_pipe=params.l_pipe #-outflow pipe length m (hypolimnetic withdrawal)
    d_pipe=params.d_pipe #-outflow pipe diameter m
    f_pipe=params.f_pipe #-pipe friction factor
    #-settings for SCHI-PU scenario
    inflow_pump=params.inflow_pump #-inflow pump rate l/s

    schifgr_highest_elev+=flow_occur_thresh
    schifgr_havel_bed-=mud_thresh
    schifgr_highest_elev-=mud_thresh
    schifgr_lake_bed-=mud_thresh

    df=pd.read_csv(havel_file, encoding='latin', delimiter=';')
    df['date']=pd.to_datetime(df['Datum'], format='%d.%m.%Y')
    df=df[['date','Tagesmittelwert']]
    df=df.replace(-777,np.nan)
    df=df.set_index('date')[['Tagesmittelwert']].interpolate().dropna().reset_index()
    df['havel_mASL']=df['Tagesmittelwert']/100.000+havel_gauge_zero

    ##### Outflow pipe (two connected reservoirs) ####
    #- used Energy equation and Darcy_Weisbach equation
    #- for flow calculations see example of URL below:
    #- https://www.youtube.com/watch?v=lYpKPl2zP24
    #- for specific calculations see docs/V_Darcy_Weisbach.jpg
    def PIPE_flow(higher_mASL,lower_mASL):
        #-using Energy equation and Darcy_Wisbach
        V=math.sqrt(\
            ((higher_mASL-lower_mASL)*d_pipe*2*9.81)\
            /(f_pipe*l_pipe)
            )
        Q=-V*pi*(d_pipe/2)**2
        if math.isnan(Q):
            return 0
        else:
            return Q*86400

    #### pipe flow #####
    df = df.set_index('date').resample('D').asfreq().reset_index().interpolate()

    if lake_mASL>pipe_flow_threshold and havel_mASL>pipe_flow_threshold:
        #- 2 inflow, both sides>schifgr_lake_bed, lake<havel => channel flooded,inflow
        if lake_mASL<havel_mASL:
            df.loc[i, 'outflow_regime']=20
            Q_in_day=PIPE_flow(havel_mASL, lake_mASL)
        else:
            df.loc[i, 'outflow_regime']=-20
            Qout_day=PIPE_flow(lake_mASL,havel_mASL)

    return Qin_day,Qout_day
    



#####################################################################################################
#####################################################################################################

def SCHI_OP_EXST(date_calc,lake_mASL):

    Qin_day=0
    Qout_day=0

    ##-import params from havel_lake_params.py
    havel_gauge_zero=params.havel_gauge_zero
    lake_gauge_zero=params.lake_gauge_zero
    lake_bot_elev=params.lake_bot_elev
    mud_thresh=params.mud_thresh #-mud threshold to subtract from bed levels
    schifgr_havel_bed=params.schifgr_havel_bed #- mASL western side of dam 
    schifgr_highest_elev=params.schifgr_highest_elev #- mASL highest point (midpoint) of schiffgraben (where flow into lake/havel begins) 
    schifgr_lake_bed=params.schifgr_lake_bed #-mASL mouth to lake
    flow_occur_thresh=params.flow_occur_thresh #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 
    #-channel geometry at control section (Sec EE)
    l=params.l #- total channel length (m)
    lx=params.lx #-length of channel at control section (m) Sec EE in dwg file; distance from lake side  
    b=params.b #-channel bed width at control section (m) pipe under dam considered (rough value)
    alpha=params.alpha #-side angles of walls (degrees) Section EE = 75deg (in dwg file)
    n=params.n #-channel roughness (From Chow's book for earth channels)
    #-channel geometry at mid point (sec CC)
    l_mid=params.l_mid #- channel length of schifgr_highest_elev (m) at Sec CC (dwg file) dist. from lake side 
    alpha_mid=params.alpha_mid #side angles of walls (degrees) at lmid Section CC = 47deg (in dwg file)
    b_mid=params.b_mid #- channel bed width at midpoint (sec CC)
    l_pipe=params.l_pipe #-outflow pipe length m (hypolimnetic withdrawal)
    d_pipe=params.d_pipe #-outflow pipe diameter m
    f_pipe=params.f_pipe #-pipe friction factor
    #-settings for SCHI-PU scenario
    inflow_pump=params.inflow_pump #-inflow pump rate l/s

    schifgr_highest_elev+=flow_occur_thresh
    schifgr_havel_bed-=mud_thresh
    schifgr_highest_elev-=mud_thresh
    schifgr_lake_bed-=mud_thresh

    df=pd.read_csv(havel_file, encoding='latin', delimiter=';')
    df['date']=pd.to_datetime(df['Datum'], format='%d.%m.%Y')
    df=df[['date','Tagesmittelwert']]
    df=df.replace(-777,np.nan)
    df=df.set_index('date')[['Tagesmittelwert']].interpolate().dropna().reset_index()
    df['havel_mASL']=df['Tagesmittelwert']/100.000+havel_gauge_zero

    ####-function to calculate Manning flow of channel under different flow regimes:
    #-when flow regime 2 or -2 (flooded flow regime)
    #-- 2: inflow, Y1(Havel) and Y2(Lake) constant for day; Y is water surface elev
    #-- -2: outflow, Y1(Lake) and Y2(Havel) constant for day; Y is water surface elev
    def flow_reg_2(y1,y2,flow_regime):
        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- lx distance of control section from lake side of channel (to calculate flow height h)
        ###- bottom slope assumed to be zero; energy line slope assumed identical to water surface slope
        
        #- find flow depth at control section
        y=lx*(y1-y2)/l

        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        #-b_elev:bed elev at control section
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        S=(y1-y2)/l
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400


    #-when flow regime 1 or -1 (dry flow regime)
    #-- 1: inflow, Y1(Havel) and Y2(Lake); y is water level
    #-- -1: outflow, Y1(Lake) and Y2(Havel); y is water level
    def flow_reg_1(y1,flow_regime):
        if flow_regime==1:
            y2_bed=schifgr_lake_bed+flow_occur_thresh
        elif flow_regime==-1:
            y2_bed=schifgr_havel_bed+flow_occur_thresh

        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- l_mid distance of midpoint from lake side of channel (to calculate long. slope)
        ###- daily flow assumed uniform, energy line slope identical to bottom slope
        b=b_mid
        alpha=alpha_mid
        #- find flow depth at control section (overall uniform)
        y=y1-schifgr_highest_elev

        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        #-b_elev:bed elev at control section
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        S=(schifgr_highest_elev-y2_bed)/l
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400

    df=df.reset_index()

    ####- channel flooded
    havel_mASL=df[df.date == date_calc]['havel_mASL'].iloc[0]
    if lake_mASL>schifgr_highest_elev and havel_mASL>schifgr_highest_elev:
        #- -2 outflow, both sides>schifgr_lake_bed, lake>havel => channel flooded,outflow
        if lake_mASL>havel_mASL:
            flow_regime=-2
            Qout_day=flow_reg_2(lake_mASL, havel_mASL, flow_regime)

        #- 2 inflow, both sides>schifgr_lake_bed, lake<havel => channel flooded,inflow
        elif lake_mASL<havel_mASL:
            flow_regime=2
            Qin_day=flow_reg_2(havel_mASL, lake_mASL, flow_regime)

    #- -1 outflow, lake>havel & lake>schifgr_lake_bed => channel dry, outflow
    elif lake_mASL>schifgr_highest_elev and havel_mASL<schifgr_highest_elev:
        flow_regime=-1
        Qout_day=flow_reg_1(lake_mASL, flow_regime)

    #- 1 inflow, lake<havel & Havel>schifgr_lake_bed => channel dry, inflow
    elif lake_mASL<schifgr_highest_elev and havel_mASL>schifgr_highest_elev:
        flow_regime=1
        Qin_day=flow_reg_1(havel_mASL, flow_regime)
    

    return Qin_day,Qout_day




#####################################################################################################
#####################################################################################################

def SCHI_FLI_EXST(date_calc,lake_mASL):

    Qin_day=0
    Qout_day=0

    ##-import params from havel_lake_params.py
    havel_gauge_zero=params.havel_gauge_zero
    lake_gauge_zero=params.lake_gauge_zero
    lake_bot_elev=params.lake_bot_elev
    mud_thresh=params.mud_thresh #-mud threshold to subtract from bed levels
    schifgr_havel_bed=params.schifgr_havel_bed #- mASL western side of dam 
    schifgr_highest_elev=params.schifgr_highest_elev #- mASL highest point (midpoint) of schiffgraben (where flow into lake/havel begins) 
    schifgr_lake_bed=params.schifgr_lake_bed #-mASL mouth to lake
    flow_occur_thresh=params.flow_occur_thresh #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 
    #-channel geometry at control section (Sec EE)
    l=params.l #- total channel length (m)
    lx=params.lx #-length of channel at control section (m) Sec EE in dwg file; distance from lake side  
    b=params.b #-channel bed width at control section (m) pipe under dam considered (rough value)
    alpha=params.alpha #-side angles of walls (degrees) Section EE = 75deg (in dwg file)
    n=params.n #-channel roughness (From Chow's book for earth channels)
    #-channel geometry at mid point (sec CC)
    l_mid=params.l_mid #- channel length of schifgr_highest_elev (m) at Sec CC (dwg file) dist. from lake side 
    alpha_mid=params.alpha_mid #side angles of walls (degrees) at lmid Section CC = 47deg (in dwg file)
    b_mid=params.b_mid #- channel bed width at midpoint (sec CC)
    l_pipe=params.l_pipe #-outflow pipe length m (hypolimnetic withdrawal)
    d_pipe=params.d_pipe #-outflow pipe diameter m
    f_pipe=params.f_pipe #-pipe friction factor
    #-settings for SCHI-PU scenario
    inflow_pump=params.inflow_pump #-inflow pump rate l/s

    schifgr_highest_elev+=flow_occur_thresh
    schifgr_havel_bed-=mud_thresh
    schifgr_highest_elev-=mud_thresh
    schifgr_lake_bed-=mud_thresh

    df=pd.read_csv(havel_file, encoding='latin', delimiter=';')
    df['date']=pd.to_datetime(df['Datum'], format='%d.%m.%Y')
    df=df[['date','Tagesmittelwert']]
    df=df.replace(-777,np.nan)
    df=df.set_index('date')[['Tagesmittelwert']].interpolate().dropna().reset_index()
    df['havel_mASL']=df['Tagesmittelwert']/100.000+havel_gauge_zero

    ####-function to calculate Manning flow of channel under different flow regimes:
    #-when flow regime 2 or -2 (flooded flow regime)
    #-- 2: inflow, Y1(Havel) and Y2(Lake) constant for day; Y is water surface elev
    #-- -2: outflow, Y1(Lake) and Y2(Havel) constant for day; Y is water surface elev
    def flow_reg_2(y1,y2,flow_regime):
        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- lx distance of control section from lake side of channel (to calculate flow height h)
        ###- bottom slope assumed to be zero; energy line slope assumed identical to water surface slope
        
        #- find flow depth at control section
        y=lx*(y1-y2)/l

        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        #-b_elev:bed elev at control section
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        S=(y1-y2)/l
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400


    #-when flow regime 1 or -1 (dry flow regime)
    #-- 1: inflow, Y1(Havel) and Y2(Lake); y is water level
    #-- -1: outflow, Y1(Lake) and Y2(Havel); y is water level
    def flow_reg_1(y1,flow_regime):
        if flow_regime==1:
            y2_bed=schifgr_lake_bed+flow_occur_thresh
        elif flow_regime==-1:
            y2_bed=schifgr_havel_bed+flow_occur_thresh

        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- l_mid distance of midpoint from lake side of channel (to calculate long. slope)
        ###- daily flow assumed uniform, energy line slope identical to bottom slope
        b=b_mid
        alpha=alpha_mid
        #- find flow depth at control section (overall uniform)
        y=y1-schifgr_highest_elev

        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        #-b_elev:bed elev at control section
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        S=(schifgr_highest_elev-y2_bed)/l
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400

    df=df.reset_index()
    df = df.set_index('date').resample('D').asfreq().reset_index().interpolate()

    ####- channel flooded
    havel_mASL=df[df.date == date_calc]['havel_mASL'].iloc[0]
    if lake_mASL>schifgr_highest_elev and havel_mASL>schifgr_highest_elev:
        if lake_mASL<havel_mASL:
            flow_regime=2
            Qin_day=flow_reg_2(havel_mASL, lake_mASL, flow_regime)

    #- 1 inflow, lake<havel & Havel>schifgr_lake_bed => channel dry, inflow
    elif lake_mASL<schifgr_highest_elev and havel_mASL>schifgr_highest_elev:
        flow_regime=1
        Qin_day=flow_reg_1(havel_mASL, flow_regime)

    return Qin_day, Qout_day
        


#####################################################################################################
#####################################################################################################

def SCHI_OP_CLN(date_calc,lake_mASL,channel_bed_mASL):

    Qin_day=0
    Qout_day=0

    ##-import params from havel_lake_params.py
    havel_gauge_zero=params.havel_gauge_zero
    lake_gauge_zero=params.lake_gauge_zero
    lake_bot_elev=params.lake_bot_elev
    mud_thresh=params.mud_thresh #-mud threshold to subtract from bed levels
    schifgr_havel_bed=params.schifgr_havel_bed #- mASL western side of dam 
    schifgr_highest_elev=params.schifgr_highest_elev #- mASL highest point (midpoint) of schiffgraben (where flow into lake/havel begins) 
    schifgr_lake_bed=params.schifgr_lake_bed #-mASL mouth to lake
    flow_occur_thresh=params.flow_occur_thresh #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 
    #-channel geometry at control section (Sec EE)
    l=params.l #- total channel length (m)
    lx=params.lx #-length of channel at control section (m) Sec EE in dwg file; distance from lake side  
    b=params.b #-channel bed width at control section (m) pipe under dam considered (rough value)
    alpha=params.alpha #-side angles of walls (degrees) Section EE = 75deg (in dwg file)
    n=params.n #-channel roughness (From Chow's book for earth channels)
    #-channel geometry at mid point (sec CC)
    l_mid=params.l_mid #- channel length of schifgr_highest_elev (m) at Sec CC (dwg file) dist. from lake side 
    alpha_mid=params.alpha_mid #side angles of walls (degrees) at lmid Section CC = 47deg (in dwg file)
    b_mid=params.b_mid #- channel bed width at midpoint (sec CC)
    l_pipe=params.l_pipe #-outflow pipe length m (hypolimnetic withdrawal)
    d_pipe=params.d_pipe #-outflow pipe diameter m
    f_pipe=params.f_pipe #-pipe friction factor
    #-settings for SCHI-PU scenario
    inflow_pump=params.inflow_pump #-inflow pump rate l/s

    schifgr_highest_elev+=flow_occur_thresh
    schifgr_havel_bed-=mud_thresh
    schifgr_highest_elev-=mud_thresh
    schifgr_lake_bed-=mud_thresh

    #- schiffgraben bed is assumed clean and horizontal
    schifgr_bed=channel_bed_mASL

    #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 
    flow_occur_thresh=schifgr_bed+0.01

    df=pd.read_csv(havel_file, encoding='latin', delimiter=';')
    df['date']=pd.to_datetime(df['Datum'], format='%d.%m.%Y')
    df=df[['date','Tagesmittelwert']]
    df=df.replace(-777,np.nan)
    df=df.set_index('date')[['Tagesmittelwert']].interpolate().dropna().reset_index()
    df['havel_mASL']=df['Tagesmittelwert']/100.000+havel_gauge_zero

    #-when flow regime 1 or -1 
    #-- 1: inflow, Y1(Havel) and Y2(Lake); y is water level
    #-- -1: outflow, Y1(Lake) and Y2(Havel); y is water level
    def flow_reg_1(y1):
        y1-=schifgr_bed #-covert elev to hight
        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- lx distance of control section from lake side of channel (to calculate flow height h)
        ###- daily flow assumed uniform, energy line slope identical to bottom slope
        #- find flow depth at control section (overall uniform)
        y=(l-lx)*(y1)/l

        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        S=(y1)/l #y2 is assumed zero
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400


    #-when flow regime 2 or -2 
    #-- 2: inflow, Y1(Havel) and Y2(Lake) constant for day; Y is water surface elev
    #-- -2: outflow, Y1(Lake) and Y2(Havel) constant for day; Y is water surface elev
    def flow_reg_2(y1,y2):
        y1-=schifgr_bed; y2-=schifgr_bed #-covert elev to hight
        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- lx distance of control section from lake side of channel (to calculate flow height h)
        ###- bottom slope assumed to be zero; energy line slope assumed identical to water surface slope
        
        #- find flow depth at control section
        y=lx*y1/l
        
        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        #-b_elev:bed elev at control section
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius
        
        ###- energy line slope S identical to water surface slope (may change later)
        S=(y1-y2)/l
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400

    df=df.reset_index()
    df = df.set_index('date').resample('D').asfreq().reset_index().interpolate()

    ####- channel flooded
    havel_mASL=df[df.date == date_calc]['havel_mASL'].iloc[0]
    if lake_mASL>flow_occur_thresh and havel_mASL>flow_occur_thresh:
        #- -2 outflow, both sides>schifgr_lake_bed, lake>havel => channel flooded,outflow
        if lake_mASL>havel_mASL:
            Qout_day=flow_reg_2(lake_mASL, havel_mASL)

        #- 2 inflow, both sides>schifgr_lake_bed, lake<havel => channel flooded,inflow
        elif lake_mASL<havel_mASL:
            Qin_day=flow_reg_2(havel_mASL, lake_mASL)
            

    #- -1 outflow, lake>havel & lake>schifgr_lake_bed => channel dry, outflow
    elif lake_mASL>flow_occur_thresh and havel_mASL<flow_occur_thresh:
        Qout_day=flow_reg_1(lake_mASL)

    #- 1 inflow, lake<havel & Havel>schifgr_lake_bed => channel dry, inflow
    elif lake_mASL<flow_occur_thresh and havel_mASL>flow_occur_thresh:
        Qin_day=flow_reg_1(havel_mASL)

    return Qin_day,Qout_day



#####################################################################################################
#####################################################################################################

def SCHI_FLIO_PIPE(date_calc,lake_mASL,d_pipe):

    Qin_day=0
    Qout_day=0

    ##-import params from havel_lake_params.py
    havel_gauge_zero=params.havel_gauge_zero
    lake_gauge_zero=params.lake_gauge_zero
    lake_bot_elev=params.lake_bot_elev
    mud_thresh=params.mud_thresh #-mud threshold to subtract from bed levels
    schifgr_havel_bed=params.schifgr_havel_bed #- mASL western side of dam 
    schifgr_highest_elev=params.schifgr_highest_elev #- mASL highest point (midpoint) of schiffgraben (where flow into lake/havel begins) 
    schifgr_lake_bed=params.schifgr_lake_bed #-mASL mouth to lake
    flow_occur_thresh=params.flow_occur_thresh #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 
    #-channel geometry at control section (Sec EE)
    l=params.l #- total channel length (m)
    lx=params.lx #-length of channel at control section (m) Sec EE in dwg file; distance from lake side  
    b=params.b #-channel bed width at control section (m) pipe under dam considered (rough value)
    alpha=params.alpha #-side angles of walls (degrees) Section EE = 75deg (in dwg file)
    n=params.n #-channel roughness (From Chow's book for earth channels)
    #-channel geometry at mid point (sec CC)
    l_mid=params.l_mid #- channel length of schifgr_highest_elev (m) at Sec CC (dwg file) dist. from lake side 
    alpha_mid=params.alpha_mid #side angles of walls (degrees) at lmid Section CC = 47deg (in dwg file)
    b_mid=params.b_mid #- channel bed width at midpoint (sec CC)
    l_pipe=params.l_pipe #-outflow pipe length m (hypolimnetic withdrawal)
    d_pipe=params.d_pipe #-outflow pipe diameter m
    f_pipe=params.f_pipe #-pipe friction factor
    #-settings for SCHI-PU scenario
    inflow_pump=params.inflow_pump #-inflow pump rate l/s

    schifgr_highest_elev+=flow_occur_thresh
    schifgr_havel_bed-=mud_thresh
    schifgr_highest_elev-=mud_thresh
    schifgr_lake_bed-=mud_thresh

    df=pd.read_csv(havel_file, encoding='latin', delimiter=';')
    df['date']=pd.to_datetime(df['Datum'], format='%d.%m.%Y')
    df=df[['date','Tagesmittelwert']]
    df=df.replace(-777,np.nan)
    df=df.set_index('date')[['Tagesmittelwert']].interpolate().dropna().reset_index()
    df['havel_mASL']=df['Tagesmittelwert']/100.000+havel_gauge_zero

    havel_mASL=df[df.date == date_calc]['havel_mASL'].iloc[0]

    ##### Outflow pipe (two connected reservoirs) ####
    #- used Energy equation and Darcy_Wisbach equation
    #- for flow calculations see example of URL below:
    #- https://www.youtube.com/watch?v=lYpKPl2zP24
    #- for specific calculations see docs/V_Darcy_Wisbach.jpg
    def PIPE_outflow(lake_mASL,havel_mASL):
        if lake_mASL>havel_mASL:
            #-using Energy equation and Darcy_Wisbach
            V=math.sqrt(\
                ((lake_mASL-havel_mASL)*d_pipe*2*9.81)\
                /(f_pipe*l_pipe)
                )
            
            Q=V*pi*(d_pipe/2)**2
            if math.isnan(Q):
                return 0
            else:
                return Q*86400
        else:
            return 0

    ####-function to calculate Manning flow of channel under different flow regimes:
    #-when flow regime 2 or -2 (flooded flow regime)
    #-- 2: inflow, Y1(Havel) and Y2(Lake) constant for day; Y is water surface elev
    #-- -2: outflow, Y1(Lake) and Y2(Havel) constant for day; Y is water surface elev
    def flow_reg_2(y1,y2,flow_regime):
        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- lx distance of control section from lake side of channel (to calculate flow height h)
        ###- bottom slope assumed to be zero; energy line slope assumed identical to water surface slope
        
        #- find flow depth at control section
        y=lx*(y1-y2)/l

        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        #-b_elev:bed elev at control section
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        S=(y1-y2)/l
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400

    #-when flow regime 1 or -1 (dry flow regime)
    #-- 1: inflow, Y1(Havel) and Y2(Lake); y is water level
    #-- -1: outflow, Y1(Lake) and Y2(Havel); y is water level
    def flow_reg_1(y1,flow_regime):
        if flow_regime==1:
            y2_bed=schifgr_lake_bed+flow_occur_thresh
        elif flow_regime==-1:
            y2_bed=schifgr_havel_bed+flow_occur_thresh

        ###-discharge calculation between two connected reservoirs
        ###- method from Open-Channel Hydraulics by Chow
        ###- l_mid distance of midpoint from lake side of channel (to calculate long. slope)
        ###- daily flow assumed uniform, energy line slope identical to bottom slope
        b=b_mid
        alpha=alpha_mid
        #- find flow depth at control section (overall uniform)
        y=y1-schifgr_highest_elev

        #-compute hydraulic radius
        #-alpha: angle of side walls, b: bed width of section, h: water depth (flow height) 
        #-b_elev:bed elev at control section
        m=1/math.tan(alpha*pi/180)#-side slope (1:m; 1 vertical, m horizontal)
        B=b+m*y #-channel width at flow surface
        A=y*(b+B)/2 #-area of flow
        P=b+2*y*math.sqrt(1+m**2) #-wetted perimeter
        R=A/P #-hydraulic radius

        ###- energy line slope S identical to water surface slope (may change later)
        if flow_regime==1:
            S=(schifgr_highest_elev-y2_bed)/l
        elif flow_regime==-1:
            S=(schifgr_highest_elev-y2_bed)/(l-l_mid)
        Q=A*(1/n)*R**(2/3)*math.sqrt(S)
        if math.isnan(Q):
            return 0
        else:
            return Q*86400

    df=df.reset_index()

    #### channel inflow #####
    #- channel flooded
    df = df.set_index('date').resample('D').asfreq().reset_index().interpolate()

    havel_mASL=df[df.date == date_calc]['havel_mASL'].iloc[0]
    if lake_mASL>schifgr_highest_elev and havel_mASL>schifgr_highest_elev:
        #- 2 inflow, both sides>schifgr_lake_bed, lake<havel => channel flooded,inflow
        if lake_mASL<havel_mASL:
            inflow_regime=2
            Qin_day=flow_reg_2(havel_mASL, lake_mASL, inflow_regime)
            
    #-channel dry
    #- 1 inflow, lake<havel & Havel>schifgr_lake_bed => channel dry, inflow
    elif lake_mASL<schifgr_highest_elev and havel_mASL>schifgr_highest_elev:
        inflow_regime=1
        Qin_day=flow_reg_1(havel_mASL, inflow_regime)

    if lake_mASL>havel_mASL:
        Qout_day=PIPE_outflow(lake_mASL,havel_mASL)

    return Qin_day,Qout_day



#####################################################################################################
#####################################################################################################

def SCHI_PU(date_calc,lake_mASL,inflow_pump_ls,d_pipe):

    Qout_day=0

    ##-import params from havel_lake_params.py
    havel_gauge_zero=params.havel_gauge_zero
    lake_gauge_zero=params.lake_gauge_zero
    lake_bot_elev=params.lake_bot_elev
    mud_thresh=params.mud_thresh #-mud threshold to subtract from bed levels
    schifgr_havel_bed=params.schifgr_havel_bed #- mASL western side of dam 
    schifgr_highest_elev=params.schifgr_highest_elev #- mASL highest point (midpoint) of schiffgraben (where flow into lake/havel begins) 
    schifgr_lake_bed=params.schifgr_lake_bed #-mASL mouth to lake
    flow_occur_thresh=params.flow_occur_thresh #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 
    #-channel geometry at control section (Sec EE)
    l=params.l #- total channel length (m)
    lx=params.lx #-length of channel at control section (m) Sec EE in dwg file; distance from lake side  
    b=params.b #-channel bed width at control section (m) pipe under dam considered (rough value)
    alpha=params.alpha #-side angles of walls (degrees) Section EE = 75deg (in dwg file)
    n=params.n #-channel roughness (From Chow's book for earth channels)
    #-channel geometry at mid point (sec CC)
    l_mid=params.l_mid #- channel length of schifgr_highest_elev (m) at Sec CC (dwg file) dist. from lake side 
    alpha_mid=params.alpha_mid #side angles of walls (degrees) at lmid Section CC = 47deg (in dwg file)
    b_mid=params.b_mid #- channel bed width at midpoint (sec CC)
    l_pipe=params.l_pipe #-outflow pipe length m (hypolimnetic withdrawal)
    d_pipe=params.d_pipe #-outflow pipe diameter m
    f_pipe=params.f_pipe #-pipe friction factor
    #-settings for SCHI-PU scenario
    inflow_pump=params.inflow_pump #-inflow pump rate l/s

    schifgr_highest_elev+=flow_occur_thresh
    schifgr_havel_bed-=mud_thresh
    schifgr_highest_elev-=mud_thresh
    schifgr_lake_bed-=mud_thresh

    inflow_pump=inflow_pump_ls*86400/1000 #-convert l/s to m3/d

    df=pd.read_csv(havel_file, encoding='latin', delimiter=';')
    df['date']=pd.to_datetime(df['Datum'], format='%d.%m.%Y')
    df=df[['date','Tagesmittelwert']]
    df=df.replace(-777,np.nan)
    df=df.set_index('date')[['Tagesmittelwert']].interpolate().dropna().reset_index()
    df['havel_mASL']=df['Tagesmittelwert']/100.000+havel_gauge_zero

    havel_mASL=df[df.date == date_calc]['havel_mASL'].iloc[0]

    ##### Outflow pipe (two connected reservoirs) ####
    #- used Energy equation and Darcy_Wisbach equation
    #- for flow calculations see example of URL below:
    #- https://www.youtube.com/watch?v=lYpKPl2zP24
    #- for specific calculations see docs/V_Darcy_Wisbach.jpg
    def PIPE_outflow(lake_mASL,havel_mASL):
        #-using Energy equation and Darcy_Wisbach
        V=math.sqrt(\
            ((lake_mASL-havel_mASL)*d_pipe*2*9.81)\
            /(f_pipe*l_pipe)
            )
        
        Q=V*pi*(d_pipe/2)**2
        if math.isnan(Q):
            return 0
        else:
            return Q*86400

    df=df.reset_index()

    #### channel inflow #####
    #- channel flooded
    df = df.set_index('date').resample('D').asfreq().reset_index().interpolate()

    if lake_mASL>havel_mASL:
        Qout_day=PIPE_outflow(lake_mASL,havel_mASL)

    return inflow_pump,Qout_day

    
    