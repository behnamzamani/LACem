"""
    Python script for LACem model
	Behnam Zamani, TU Berlin
    
    Version 1.15
    predicts prediction of lake surface
    predicts groundwater levels (using catchment surface hydrology)
    New >> predicts time overlap of wells (in addition to lake level and met data)

    Version 1.30
    predicts surface water exchange with nearby water bodies through channels (using Manning's formula)
    For Sacrower See case: inflow/outflow via Schiffgraben + Flap gate cases + hypolimnetic withdrawal

    Version 1.31
    Evaporation mode Priestley-Taylor method (Chow's book eq. 3.5.27) added 
    Evap. model can be chosen between Penman and Priestley-Taylor

    """
import sys
import numpy as np
import pandas as pd
import math
import datetime
import os
from pandas import Timestamp
import time
import glob
import calendar
from statsmodels.tsa.seasonal import seasonal_decompose
from scipy import stats
from scipy import signal
import pymannkendall as mk
import multiprocessing as mp
import re
from time import sleep
import itertools
from tqdm import tqdm
import random
import fileinput
import urllib.request
from progress.bar import Bar
import subprocess
import shutil
import netCDF4 as nc
import matplotlib
#matplotlib.use('Agg') #-- active when on cluster
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import colors
import os
from pandas import Timestamp
import matplotlib.dates as dates
import pyny3d.geoms as pyny
import warnings
import random
import decimal
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore",category=matplotlib.MatplotlibDeprecationWarning)
import scipy
import sys

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

version=os.path.splitext(os.path.basename(__file__))[0]

######---- model pre-settings:
ignore_snow=0 #- 1=ignore snow contribution to surface runoff in SCS method, 0=do nothing
SCS_surf_flow=1 #-1: include surface flow into lake (SCS method); 0: exclude surface flow into lake (SCS method)


###- surface Manning flow scenario, to calculate inflow/outflow between lake and nearby water bodies through canals
###- see function manning_surf_flow for more information
###- Available scenarios:
###-  0 SCHI_OFF
###-  1 SCHI_DAM_REFURB2015
###-- 2 SCHI_OP_EXST
###-- 3 SCHI_FLI_EXST
###-- 4 SCHI_FLIO_PIPE
###-- 5 SCHI_PU


#-call external modules/functions
# functions defined in modules/havel_lake_exchange.py
from havel_lake_exchange import *
from havel_lake_params import inflow_pump,d_pipe

#-create dictionary of Manning surface flow scenarios (to choose automatically by model depending on scenario)
surf_flow_dict={
'SCHI_OFF':SCHI_OFF,
'SCHI_DAM_REFURB2015':SCHI_DAM_REFURB2015,
'SCHI_DAM_SIPH_EXST':SCHI_DAM_SIPH_EXST,
'SCHI_OP_EXST':SCHI_OP_EXST,
'SCHI_FLI_EXST':SCHI_FLI_EXST,
'SCHI_FLIO_PIPE':SCHI_FLIO_PIPE,
'SCHI_PU':SCHI_PU
}

surf_scenarios=[
'SCHI_DAM_REFURB2015',
'SCHI_OP_EXST',
'SCHI_FLI_EXST',
'SCHI_FLIO_PIPE',
'SCHI_PU',
   ]

def run_0():    return model(surf_scenarios[0])
def run_1():    return model(surf_scenarios[1])
def run_2():    return model(surf_scenarios[2])
def run_3():    return model(surf_scenarios[3])
def run_4():    return model(surf_scenarios[4])

all_processes=[

    mp.Process(target=run_0),
    mp.Process(target=run_1),
    mp.Process(target=run_2),
    mp.Process(target=run_3),
    mp.Process(target=run_4),
]



#######-----------------------------------------------------------------------###############
##- area correction factor (subsurface_flow_coeff): to correct effective area of surface flow
###- effective area of catchment is the portion of catchment where most of water directly flows on
###-- surface into the lake before penetrating into groundwater
####--- using SCS method (the rest area is ignored as it penetrates to growundwater)
####--  area correction factor needs calibration 
subsurface_flow_coeff=0.735 #- for from 2012: calibrated value for Sacrower See: 0.65; then 0.668 

#-adjust inflow outflow (not realistic, only to adjust for AEM3D)
xadjust_inflow=1.0
xadjust_outflow=1.0


#-input data for AEM3D boundary condition files
aem3d_boundaries=[] #- remains empty, will be filled automatically

# future changes:
#-later define a function: lake surface area as function of h, replace wherever lake surface area
#- model settings and input files to read from separate config files

################################ meteorological settings and files ######################################
##-evaporation model settings
##---- 1: Penman equation
##---- 2: Priestley-Taylor method
evaporation_method=2

#-Priestley-Taylor method (Chow's book eq. 3.5.27)
alpha_evap=0.8

##-SCS land use data
file_built_20_25='data/data_landuse/blt_20_25.csv'
file_built_25_30='data/data_landuse/blt_25_30.csv'
file_built_30_38='data/data_landuse/blt_30_38.csv'
file_built_38_65='data/data_landuse/blt_38_65.csv'
file_built_65_100='data/data_landuse/blt_65_100.csv'
file_nonwoody_50_75='data/data_landuse/nwd_50_75.csv'
file_nonwoody_75_100='data/data_landuse/nwd_75_100.csv'
file_woody_50_100='data/data_landuse/wd_50_100.csv'


#####################################################################################################
##### Model Function ################################################################################ 
#####################################################################################################

def model(surf_scenario):
    print('modeling',surf_scenario)

    if surf_scenario=='SCHI_PU':
        surf_scenario_plot='SCHI_PU_%sLS_d%scm'%(inflow_pump,d_pipe*100)
    else:
        surf_scenario_plot=surf_scenario

    if not os.path.isdir('results_%s_%s' %(version,surf_scenario_plot)):
        os.mkdir('results_%s_%s' %(version,surf_scenario_plot), 0o666)

    ##- met data file (source data produkt_klima_tag*.txt can be downloaded directly from DWD)
    df_dwd=pd.read_csv('data/data_weather/produkt_klima_tag_18930101_20231231_03987.txt',delimiter=';')
    df_dwd['date']=pd.to_datetime(df_dwd['MESS_DATUM'], format='%Y%m%d')
    df_dwd.columns = df_dwd.columns.str.strip() #-remove initial white space from headers
    df_dwd=df_dwd[['date','FM','RSK','NM','VPM','PM','UPM','TXK','TMK','TNK']]
 
    ##- lake height-volume relationship file
    lake_h_v='data/data_lake/h_v.txt' #-orig: 'data/data_lake/h_v.txt'

    ##- lake level data file:
    lake_level_file='data/data_lake/pegel_taeglich_ifb_2003_2023.txt'
    met_file_dwd_solar='data/data_weather_solar/produkt_st_tag_19451231_20240630_03987.txt'

    ##-Havel water level data:
    havel_file='data/data_havel/havel_wasserstand_tw_01_11_1961.csv'

    ######################### Plot module ############################
    plot_module=1 #- 0: no plot, 1: plots output


    ######################### input data/Lake info/model settings ####################################
    elev_type=2 #- 1 = height, normal height over lake bottom
                #- 2 = height in aem3d format, extracted h-v from Hydrohub lite
                # (height over lake bottom, needs correction)

    bot_depth=36 #-depth (m) of lake's deepest point (from absolute zero surface of bathymetry)


    lon = 13.062     #-longitude in degrees
    lat = 52.382     #-lattitude in degrees
    elev_met= 81        #- elevation of met station in mASL

    lake_bot_elev=-6.11 #- lake bottom elevation in mASL (deepest)
    elev_indic=28.15 #- elevation of the bottom of lake level staff gauge in mASL (pegel PNP NHN) 
    havel_gauge_zero=28 #-mASL

    gamma=0.067 #-KPa/degC psycometric constant for elevation 0-100m (FAO Annex2 table 2.2) 

    r=0.08  #- Albedo for water

    #catchment_area=9304928.52 #- m2 total area of catchment calculated from GIS_SCS including builtup

    #-total area of catchment calculated from GIS_SCS, builtup areas excluded (CNs weighted from CLC method below)
    #- lake area will be subtracted automatically by the model
    catchment_area=6183076.421 #- m2

    #-reduced area of catchment calculated from GIS; only takes the area of manipulated area as 
    #-- area of runoff directly effecting the lake (just to test)
    catchment_area=2358860.847 #- m2

    lake_surf_area=1017154.45 #- m2 surface area of lake (to be defined as a function of elevation later)

    area_reed=75782.81  #- measured area of reed zone in m2 (extracted from GIS layer)

    z0=0.03 #- cm roughness of water surface for  (0.01-0.06) from Chow's hydrology Table 2.8.2
    z2=2 #- height of windspeed measurement


    ###################### Groundwater information ####################################
    #- Based on the segmented approach compute the groundwate inflow/outflow.
    #- Reference: Rosenberry, D.O., LaBaugh, J.W., 2008. https://doi.org/10.3133/tm4D2.
    #- Well IDs are used to call groundwater data files stored in data_gw directory.
    #- gw_data files must be .dat files in ascii format (see gw_data directory for more info).
    #- Well IDs will be recognized automatically from filenames with the format well_ID.dat.  
    #- Thus the filenames have suffix of well IDs.
    #- The bathymetric segments (bathymetry-area data for each segment) in data_bathy folder
    #-  will be stored in bathy_area_segments.csv. Each segment corresponds to one well in 
    #-  data_gw folder, thus the bathymetry-area data for each segment bathy_area_segments.csv
    #-  should be in one collumn whose header is the same well ID given to well_ID.dat file

    segmented_approach=1


    fmt = '%Y-%m-%d' # date format
    fmt_de = '%d.%m.%Y' # date format German
    #-- crop coefficients

    rho_w=997 #-density of water

    Kc_ini_grass=0.35   #-- initial phase (Bermuda grass)
    Kc_mid_grass=0.9    #-- middle phase
    Kc_end_grass=0.65   #-- final phase
    Dini_grass=61
    Lini_grass=10
    Ldev_grass=25       #-development phase
    Lmid_grass=35
    Lend_grass=105

    Kc_ini_reed=1.0     #-- Reed swamp standing water
    Kc_mid_reed=1.2
    Kc_end_reed=1.0
    Dini_reed=122
    Lini_reed=10
    Ldev_reed=30
    Lmid_reed=80
    Lend_reed=20

    Kc_ini_tree=0.5     #-- wallnut tree
    Kc_mid_tree=1.0
    Kc_end_tree=0.6518
    Dini_tree=92
    Lini_tree=20
    Ldev_tree=10
    Lmid_tree=130
    Lend_tree=30

    ################################################################################


    pi=math.pi


    df_dwd_solar=pd.read_csv(met_file_dwd_solar, delimiter=';',  skipinitialspace=True)

    #- fill missing data with NaN and then interpolate the missing values
    df_dwd=df_dwd.fillna(-999); df_dwd=df_dwd.interpolate()
    df_dwd['date']=pd.to_datetime(df_dwd['date'])

    #- read landuse data
    df_built_20_25=pd.read_csv(file_built_20_25)
    df_built_25_30=pd.read_csv(file_built_25_30)
    df_built_30_38=pd.read_csv(file_built_30_38)
    df_built_38_65=pd.read_csv(file_built_38_65)
    df_built_65_100=pd.read_csv(file_built_65_100)
    df_nonwoody_50_75=pd.read_csv(file_nonwoody_50_75)
    df_nonwoody_75_100=pd.read_csv(file_nonwoody_75_100)
    df_woody_50_100=pd.read_csv(file_woody_50_100)

    area_grass=df_nonwoody_50_75['AREA'].sum()+df_nonwoody_75_100['AREA'].sum()-area_reed
    area_tree=df_woody_50_100['AREA'].sum()


    ###################################################################################
    #- compute SCS runoff from daily rainfall
    ####################################################################################
    #- CNs for built up areas (not used)
    #- CNs of the total watershed are extracted from GIS as one final weighted average CN value (see GIS data)
    #- soil group was changed from B (my initial guess) to C based on GIS data (SCS.mxd)

    CN_B_20=68
    CN_B_25=70
    CN_B_30=72
    CN_B_38=75
    CN_B_65=85

    #-CNs for woody and nowoody
    CN_B_nwd_50=69
    CN_B_nwd_75=61
    CN_B_wd=55

    ###-- Read and average rainfall data of two stations
    #df_rain1=pd.read_excel(xls, 'Einzeldaten Chemie')

    #def scs():
    soil_group='B' #-B: shallow loess, sandy loam
    CN_wood_good_cover=55

    #- calculate average imperviousness from landuse tables
    #- total area of each file
    #- (ignored until surface drainage of built-up areas is clarified)
    area_built_20_25=df_built_20_25['AREA'].sum()
    area_built_25_30=df_built_25_30['AREA'].sum()
    area_built_30_38=df_built_30_38['AREA'].sum()
    area_built_38_65=df_built_38_65['AREA'].sum()
    area_built_65_100=df_built_65_100['AREA'].sum()

    area_nonwoody_50_75=df_nonwoody_50_75['AREA'].sum()
    area_nonwoody_75_100=df_nonwoody_75_100['AREA'].sum()
    area_woody_50_100=df_woody_50_100['AREA'].sum()

    #- total area of each landuse
    area_built_total=area_built_20_25+area_built_25_30+area_built_30_38+area_built_38_65+\
    area_built_65_100
    area_nonwoody_total=area_nonwoody_50_75+area_nonwoody_75_100
    area_woody_total=area_woody_50_100

    area_total=area_built_total+area_nonwoody_total+area_woody_total # including built up
    #area_total=area_nonwoody_total+area_woody_total # excluding built up (to include built up areas deactivate this)

    #- weighted CN
    CN_wtd_built=((area_built_20_25/area_built_total)*CN_B_20)+\
        ((area_built_25_30/area_built_total)*CN_B_25)+\
        ((area_built_30_38/area_built_total)*CN_B_30)+\
        ((area_built_38_65/area_built_total)*CN_B_38)+\
        ((area_built_65_100/area_built_total)*CN_B_65)

    CN_wtd_built=0 #- to ignore builtup areas (disable to compute weighted CN above)


    CN_wtd_nowoody=((area_nonwoody_50_75/area_nonwoody_total)*CN_B_nwd_50)+\
            ((area_nonwoody_75_100/area_nonwoody_total)*CN_B_nwd_75)

    CN_wtd_woody=((area_woody_50_100/area_woody_total)*CN_B_wd)

    #- CN including built up
    CN_wtd=((CN_wtd_built*area_built_total)+(CN_wtd_nowoody*area_nonwoody_total)+\
        (CN_wtd_woody*area_woody_total))/area_total

    #- CN excluding built up (to include built up areas deactivate this)
    #CN_wtd=((CN_wtd_nowoody*area_nonwoody_total)+\
    #    (CN_wtd_woody*area_woody_total))/area_total

    #- CN weighted from GIS data_landuse/CN_table_GIS.csv (to include default built-in CN approach deactivate this)
    CN_wtd=80.6 #- from GIS

    S=(25400-254*CN_wtd)/CN_wtd #- S in mm

    #- filter and ignore snowy days to calculate only rain in SCS

    if ignore_snow==1:
        df_dwd.loc[df_dwd['RSKF'] == 7,'RSK']=0 #- excluding snow from rainfall

    df_dwd['SCS_R']=((df_dwd['RSK']-0.2*S)**2)/(df_dwd['RSK']+0.8*S) #- Runoff in mm/d

    #df_dwd['SCS_R']=df_dwd['RSK']

    ##########################################################################
    #- compute evapotranspiration and open water evaporation
    ##########################################################################

    #-- input met data in DWD (Deutsche Wetter Dienst) format
    ##- crop type: reed, tree, grass
    ##- areas in m2
    ##- source: Fao website chapters 3 and 4
    ##- https://www.fao.org/3/X0490E/x0490e07.htm#solar%20radiation
    ##- https://www.fao.org/3/X0490E/x0490e08.htm#chapter%204%20%20%20determination%20of%20eto

    G=0 #-soil heat flux, zero (ignored) for daily calculations

    Lz=360-(int(lon/15.0)*15)      #- longitude of the center of local time zone [degrees west of Greenwich]
    Lm = 360-lon     #- longitude of the measurement site [degrees west of Greenwich]

    def Kc(yday,crop_type): #- function to compute crop coefficient using day of year
        if crop_type=='grass': 
            Kc_ini=Kc_ini_grass; Kc_mid=Kc_mid_grass; Kc_end=Kc_end_grass; Dini=Dini_grass; \
            Lini=Lini_grass; Ldev=Ldev_grass; Lmid=Lmid_grass; Lend=Lend_grass
        if crop_type=='reed':
            Kc_ini=Kc_ini_reed; Kc_mid=Kc_mid_reed; Kc_end=Kc_end_reed; Dini=Dini_reed; \
            Lini=Lini_reed; Ldev=Ldev_reed; Lmid=Lmid_reed; Lend=Lend_reed
        if crop_type=='tree': 
            Kc_ini=Kc_ini_tree; Kc_mid=Kc_mid_tree; Kc_end=Kc_end_tree; Dini=Dini_tree; \
            Lini=Lini_tree; Ldev=Ldev_tree; Lmid=Lmid_tree; Lend=Lend_tree 

        if yday>=Dini and yday<=Dini+Lini:
            return Kc_ini
        elif yday>Dini+Lini and yday<Dini+Lini+Ldev:
            return Kc_ini+(yday-(Dini+Lini))*((Kc_mid-Kc_ini)/(Ldev)) #- Linear interpolation between Kc_ini and Kc_mid for development period
        elif yday>=Dini+Lini+Ldev and yday<=Dini+Lini+Ldev+Lmid:
            return Kc_mid
        elif yday>Dini+Lini+Ldev+Lmid and yday<=Dini+Lini+Ldev+Lmid+Lend:
            return Kc_mid+(yday-(Dini+Lini+Ldev+Lmid))*((Kc_end-Kc_mid)/(Lend)) #- Linear interpolation between Kc_mid and Kc_end for development period
        else:
            return 0

    #-convert date column and select same dates as df_solar
    df_dwd_solar = pd.read_csv(met_file_dwd_solar, delimiter=';',  skipinitialspace=True)
    df_dwd_solar['date']=pd.to_datetime(df_dwd_solar['MESS_DATUM'].astype(str), format='%Y%m%d')
    #df_dwd_solar=df_dwd_solar[df_dwd_solar['FG_STRAHL'] !=-999.0]
    df_dwd_solar=df_dwd_solar.replace(-999.0,np.nan)

    df_dwd_solar.loc[:, ['FG_STRAHL']] = df_dwd_solar.loc[:, ['FG_STRAHL']]\
        .interpolate(method='linear', limit_direction='both', limit=100)

    df_dwd_solar=df_dwd_solar.reset_index().drop(['eor'], axis=1)

    df_dwd=df_dwd[df_dwd['date'].isin(df_dwd_solar['date'])]
    df_dwd=df_dwd.reset_index()


    #-calculate vapour pressure deficit (es - ea) (chapter 3, box 7)
    df_dwd['eTmax']=0.6108*np.exp((17.27*df_dwd['TXK'])/(df_dwd['TXK']+237.3)) #-kPa
    df_dwd['eTmin']=0.6108*np.exp((17.27*df_dwd['TNK'])/(df_dwd['TNK']+237.3)) #-kPa
    #df_dwd['es']=(df_dwd['eTmax']+df_dwd['eTmin'])/2 #-saturation vapour pressure kPa (FAO method)
    df_dwd['es']=0.1*df_dwd['VPM']/(df_dwd['UPM']/100) #-saturation vapour pressure kPa 
    df_dwd['vap_deficit']=df_dwd['es']-df_dwd['VPM']*0.1 #-calculate vapour pressure deficit (es - ea) kPa
    df_dwd['delta']=(4098*(0.6108*np.exp(17.27*df_dwd['TMK'])))\
        /((df_dwd['TMK']+237.3)**2) #-Slope of saturation vapour pressure curve (Delta capital)

    #- calculate solar radiation
    df_dwd.date=pd.to_datetime(df_dwd.date)
    df_dwd['doy']=df_dwd['date'].dt.dayofyear # day of year
    df_dwd['dr']=1+0.033*np.cos(2*pi*df_dwd['doy']/365) # inverse relative distance Earth-Sun
    df_dwd['dlt']=0.0409*np.sin(2*pi*df_dwd['doy']/365-1.39) # solar declination, delta (small)
    phi = lat*pi/180 #- lat in radians
    df_dwd['omega_s']=np.arccos(-np.tan(phi)*np.tan(df_dwd['dlt'])) # sunset hour angle
    df_dwd['Ra']=(24*60)*0.0820*df_dwd['dr']*(df_dwd['omega_s']*np.sin(phi)*np.sin(df_dwd['dlt'])+
        np.cos(phi)*np.cos(df_dwd['dlt'])*np.sin(df_dwd['omega_s']))/math.pi # extraterrestrial radiation
    df_dwd['FG_STRAHL']=df_dwd_solar['FG_STRAHL']
    df_dwd['Rs']=df_dwd_solar['FG_STRAHL']/100 # shortwave radiation in MJ/m2/day
    df_dwd['Rns']=(1-0.23)*df_dwd['Rs'] # Net solar or net shortwave radiation 
    # 0.23=alpha albedo or canopy reflection coefficient,grass reference crop
    df_dwd['Rso']=(0.75+2e-5*elev_met)*df_dwd['Ra'] #clear-sky radiation [MJ m-2 day-1]
    df_dwd['Rnl']=4.903e-9*(((df_dwd['TXK']+273.16)+(df_dwd['TNK']+273.16))/2)*\
    (0.34-0.14*np.sqrt(df_dwd['VPM']*0.1))*((1.35*df_dwd['Rs']/df_dwd['Rso'])-0.35) #- net longwave radiation
    df_dwd['Rn']=df_dwd['Rns']-df_dwd['Rnl'] #- Net radiation

    #-compute ETc (Penman-Montith eq.)
    df_dwd['ET0']=\
        (0.408*df_dwd['delta']*(df_dwd['Rn']-G)+gamma*(900/(df_dwd['TMK']+273))*df_dwd['FM']*df_dwd['vap_deficit'])\
        /(df_dwd['delta']+gamma*(1+0.34*df_dwd['FM'])) #- Penman-Monteith equation, ET for reference crop

    #- calculate crop coefficients Kc for each crop type
    df_dwd['Kc_grass']=df_dwd['doy'].apply(Kc, crop_type='grass')
    df_dwd['Kc_reed']=df_dwd['doy'].apply(Kc, crop_type='reed')
    df_dwd['Kc_tree']=df_dwd['doy'].apply(Kc, crop_type='tree')

    area_total=area_reed+area_tree+area_grass
    fract_reed=area_reed/area_total; fract_tree=area_tree/area_total; fract_grass=area_grass/area_total

    df_dwd['ET_reed']= df_dwd['Kc_reed']*df_dwd['ET0']*fract_reed
    df_dwd['ET_grass']= df_dwd['Kc_grass']*df_dwd['ET0']*fract_grass
    df_dwd['ET_tree']= df_dwd['Kc_tree']*df_dwd['ET0']*fract_tree

    #- compute actual evapotranspiration ETc
    df_dwd['ET']=(df_dwd['ET_reed']+df_dwd['ET_grass']+df_dwd['ET_tree']) #-mm/d
    df_dwd['ET']=df_dwd['ET'] #- mm/d

    #- compute Eo 
    #-Penman equation from Chow's book eqs. 3.5.10-3.5.26
    #- latent heat (lambda)=0.40
    k=0.4 #- von Karman constant (Chow's book, eq. 3.5.14) 
    df_dwd['rho_a']=(df_dwd['PM']*100*0.0289652+df_dwd['VPM']*100*0.018016)/\
        (8.31446*(df_dwd['TMK']+273)) #- kg/m3 density of humid air as function of pressure (Pa) and temp.
    df_dwd['B']=(0.622*(k**2)*df_dwd['rho_a']*df_dwd['FM'])/(df_dwd['PM']*100*rho_w*(np.log(z2*100/z0))**2) #- Chow's book, eq. 3.5.18
    df_dwd['Ea']=df_dwd['B']*(df_dwd['es']*1000-df_dwd['VPM']*100) #- Chow's book, eq. 3.5.17
    df_dwd['Er']=0.0353*df_dwd['Rn']

    df_dwd['Eo_penman']= (df_dwd['delta']/(df_dwd['delta']+gamma))*df_dwd['Er']+\
        (gamma/(df_dwd['delta']+gamma))*df_dwd['Ea']#- mm/d Penman Eq. from Chow's book eq. 3.5.26
    df_dwd.Eo_penman=df_dwd.Eo_penman.mask(df_dwd.Eo_penman.lt(0),0) #- set very small negative evaporations to 0

    #-Priestley-Taylor method (Chow's book eq. 3.5.27)
    df_dwd['Eo_prst']=alpha_evap*(df_dwd['delta']/(df_dwd['delta']+gamma))*df_dwd['Er']
    df_dwd.Eo_prst=df_dwd.Eo_prst.mask(df_dwd.Eo_prst.lt(0),0) #- set very small negative evaporations to 0


    if evaporation_method==1:
      df_dwd['Eo']=df_dwd['Eo_penman']
    elif evaporation_method==2:
        df_dwd['Eo']=df_dwd['Eo_prst']


    ####################################################################################
    #-function for lake volume interpolation
    ####################################################################################
    #- convert direct values of lake level indicator to heights above bottom
    df_h_v=pd.read_csv(lake_h_v, delimiter='\t')

    def lake_volume_mASL(h):
        #- interpolates and calculates lake volume from lake level for each time step
        #- flexible function: lake elevation date doesn't need to have same data resolution with other data (will be interpolated)
        #- requires v-h data of lake as a dataframe: columns: 'height' (height above lake bottom) and 'volume'
        
        #- find closest index of height to h
        closest_i=df_h_v.iloc[(df_h_v['elev']-h).abs().argsort()[:1]].index.values.astype(int)[0]
        x = df_h_v['elev'].to_list()
        y = df_h_v['volume'].to_list()
        interp = scipy.interpolate.interp1d(x, y)

        #return closest_i
        return interp(h)


    ####################################################################################
    #-function for calculating surface flow regimes
    ####################################################################################

    #- function to compute and interpolate lake height (over bottom) from lake volume 
    #- invert function of function lake_volume
    def lake_height_mASL(vol):
        #- find closest index of volume to vol
        closest_i=df_h_v.iloc[(df_h_v['volume']-vol).abs().argsort()[:1]].index.values.astype(int)[0]
        #if closest_i<1 or closest_i>df_h_v.index[-1]-1: #- ignore calculating first and last row (python index error)
        #    return np.nan

        x = df_h_v['volume'].to_list()
        y = df_h_v['elev'].to_list()
        interp = scipy.interpolate.interp1d(x, y)
        return interp(vol)

        #x=[float(df_h_v['volume'][closest_i]), float(df_h_v['volume'][closest_i+1])]
        #y=[float(df_h_v['height'][closest_i]), float(df_h_v['height'][closest_i+1])]   
        #return round(np.interp(vol, x, y),4)



    #################################################################################
    #- compute lake water balance
    #################################################################################
    #- initial estimates:
    #-- lake ref height (mean 0 depth)=29.52 mASL is read from LFU data aver. of max. min. heights
    #-- bottom elevation: -6.11 mASL 
    #-  these values must be updated after percise depth/pegel measurements from lake

    #-find time overlap of groundwater data
    import sys
    files=glob.glob('%s/data/data_gw/well_*[!_profile].dat'%dname)
    earliest=[];latest=[]
    for file in files:
        df=pd.read_csv(file, skiprows=7, delimiter='\t')
        df['date']=pd.to_datetime(df['date'], format='%d.%m.%Y')
        earliest.append(df['date'].iloc[0])
        latest.append(df['date'].iloc[-1])

    wells_latest_start=max(earliest)
    wells_earliest_end=min(latest)

    df_level=pd.read_csv(lake_level_file, delimiter='\t')
    df_level.date=pd.to_datetime(df_level.date)
    df_level = df_level.set_index('date')
    df_level['pegel']=df_level['pegel'].interpolate(limit_direction='both')
    df_level = df_level.reset_index()


    start_date=df_level['date'].iloc[0]
    end_date=df_level['date'].iloc[-1]

    #-- find time overlap of all dataframes
    latest_start=max(df_dwd['date'].iloc[0],df_level['date'].iloc[0],wells_latest_start)
    earliest_end=min(df_dwd['date'].iloc[-1],df_level['date'].iloc[-1],wells_earliest_end)

    df_level=df_level[df_level['date'].between(latest_start, earliest_end)]
    df_dwd=df_dwd[df_dwd['date'].between(latest_start, earliest_end)]

    df_dwd=df_dwd.reset_index().drop(['index'], axis=1)
    df_level=df_level.reset_index().drop(['index'], axis=1)
    df_out=df_dwd.merge(df_level, on='date', how='left').fillna(np.nan)

    ##- interpolate missing lake surface
    df_out['actual_mASL'] = df_out['pegel']\
        .astype(float).interpolate(method='linear', limit_direction='both', limit=500)

    ##- lake level: interpolate missing data, and merge to df_out
    df_out['lake_vol_measured']=df_out['actual_mASL'].astype(float)\
    .apply(lambda x: lake_volume_mASL(x) if pd.notnull(x) else x)

    df_out=df_out.interpolate()


    ############################################################################################
    #- compute groundwater inflow/outflow into/from lake 
    #- the segmented approach is used to compute the groundwate inflow/outflow
    #- reference: Rosenberry, D.O., LaBaugh, J.W., 2008. https://doi.org/10.3133/tm4D2.
    #-
    #--- (unsed) For this module, the bathymetry data of inflow/outflow sides are required in one csv file
    ##-- (unsed) This file must be stored in data_bathy directory with name bathy_area_segments.csv
    ##-- (unsed) including areas of each bathymetry's height range with height of range's bottom 
    ##-- (unsed) corresponding to each segment. 
    #-
    #- There must be at least two segments at two sides of lake.
    ##-- header of each segment must be the well ID of the corresponding segment
    #     (to calculate inflow/outflow surface areas as a function of surface elevation).
    #- The model recognizes inflow/outflow automatically with respect to the segment's flow direction
    ##-- which is calculated using the head difference between lake and well for each segment
    ##-- using Darcy's law.
    #- Shoreline length and well-lake shoreline distance are calculated from GIS 
    #- Bottom elevation of aquifer for each segment is imported from borehole logs of corresponding wells 

    df_bathy_inflow=pd.read_csv('data/data_bathy/bathy_inflow.csv')
    df_bathy_outflow=pd.read_csv('data/data_bathy/bathy_outflow.csv')

    depth_max=max(max(df_bathy_inflow['Z'].abs()),max(df_bathy_outflow['Z'].abs()))
    hmax_mASL=lake_bot_elev+depth_max

    #- function to compute surface area of lake's inflow/outflow faces from bathymetry for groundwater in/outflow
    #--- at each segment (not used)
    def bathy_surf_area(h, well_id): #- h (cm) daily measured surface elevation direct from level indicator (Pegel/pegel)
        sgmnt_id=str(well_id) #- segment ID is same with ID of corresponding well
        h_day=elev_indic-lake_bot_elev+(h/100) #- daily height over bottom, computed from available data
        df=pd.read_csv('data/data_bathy/bathy_area_segments.csv')
        df_smaller=df[df['height']<h_day] #- choose depths smaller than daily height
        h_x=df[df.height > h_day]['height'].iloc[0]   #- dx: height of top of area segment corresponding to h_day
        h_xm1=df_smaller['height'].iloc[-1]     #- dx-1: height of bottom of area segment corresponding to h_day
        area_x=df[df.height > h_day][sgmnt_id].iloc[0]  #-area of area segment corresponding to h_day
        area_day=df_smaller[sgmnt_id].sum()+(area_x*(((h_x-h_xm1)-(h_x-h_day))/(h_x-h_xm1)))
        return area_day



    #-read groundwater files and compute groundwater exchange with lake (segmented approach)
    #- see comments on top of script for more details

    #- create a library of segments to insert wells and segment data
    well_IDs = []
    segs={} #- segments library
    files = [f for f in os.listdir('%s/data/data_gw' %(dname)) if f.startswith('well_') and f.endswith('.dat')]
    for file in files:    
        with open('%s/data/data_gw/%s' %(dname,file), 'r') as fin:
            for line in fin.readlines():
                if 'well_id' in line:
                    well_id=int(re.findall(r'-?\d+\.?\d*', line)[0])   #- read well ID
                    well_IDs.append(well_id)
    segs = {a : [] for a in well_IDs}

    #- read well data and recognize well IDs and other specifications from data_gw directory
    # each segment number is well_ID: well_ID[well_ID, K_md, well_distance_lake, segment_shoreline_length]
    for seg in segs:
        segs[seg].append(seg) #- 1st item well_ID
        with open('%s/data/data_gw/well_%s.dat' %(dname,seg), 'r') as fin:
            for line in fin.readlines():
                if 'remark' in line: #- remark of well (to perform dummy calculations, if necessary)
                    segs[seg].append(re.findall(r'\w[-\w]*', line)[1]) #- 2nd item remark
                if 'K_md_mean' in line:  #- hydraulic conductivity of segment (Kf)
                    segs[seg].append(float(re.findall(r'-?\d+\.?\d*', line)[0]) ) #- 3rd item K_md
                if 'well_distance_lake' in line: #- distance between well and lake shoreline 
                    segs[seg].append(float(re.findall(r'-?\d+\.?\d*', line)[0])) #- 4th item well_dist_lake
                if 'segment_shoreline' in line: #- shoreline length corresponding to segment
                    segs[seg].append(float(re.findall(r'-?\d+\.?\d*', line)[0])) #- 5th item sgmnt_shore
                if 'aquif_bot_elev' in line: #- bottom of aquifer (to compute aquifer effective thickness)
                    segs[seg].append(float(re.findall(r'-?\d+\.?\d*', line)[0])) #- 6th item aquif_bot_elev


    #- compute groundwater inflow/outflow for each segment
    dummy_wells=[]
    for well_id in segs:
        if segs[well_id][1] == 'dummy_mean': #- check well remark to consider in dummy wells
            dummy_wells.append(well_id)

        df_gw=pd.read_csv('data/data_gw/well_%s.dat' %well_id, skiprows=7, delimiter='\t')
        df_gw.date=pd.to_datetime(df_gw.date, format='%d.%m.%Y')
        df_gw=df_gw[df_gw['date'].between(latest_start, earliest_end)]
        df_gw=df_gw.reset_index().drop(['index'], axis=1)
        df_gw=df_gw.rename(columns={'water_table': 'wtr_table_%s'%well_id})
        df_out=df_out.merge(df_gw, on='date', how='left').fillna(np.nan)
        ##- interpolate missing groundwater table (to daily values)
        df_out['wtr_table_%s'%well_id] = df_out['wtr_table_%s'%well_id]\
            .astype(float).interpolate(method='linear', limit_direction='both', limit=2000)
        #if df_out['wtr_table_%s'%well_id]

        if segs[well_id][1] == 'dummy_mean': #- wells mit "dummy_mean" remark are not taken into account here
            #- with respect to lake level before (water table as mean of lake level day before and Havel level)
            #- "dummy_mean" wells will be computed in the daily lake water balance computation
            #- with respect to lake level before (water table as mean of lake level day before and Havel level) 
            df_out['wtr_table_dummy_%s'%well_id]=0

    orig_wells = [i for i in well_IDs if i not in dummy_wells]

    df_out['deltaV_lake']=df_out['lake_vol_measured'].diff() #- volume difference with previous day
    df_out['SCS_R_m3']=df_out['SCS_R']*0.001*(catchment_area*subsurface_flow_coeff-lake_surf_area)
    df_out['evap_m3']=df_out['Eo']*0.001*lake_surf_area  
    df_out['evapotransp_m3']=df_out['ET']*0.001*(catchment_area-lake_surf_area)
    df_out['rain_direct_in_m3']=df_out['RSK']*0.001*lake_surf_area

    #-delta of groundwater (contrib. of groundwater water balance without segmented approach)
    df_out['deltaGW_wb_m3day']=df_out['deltaV_lake']+df_out['evap_m3'].abs()+df_out['evapotransp_m3'].abs()\
        -df_out['SCS_R_m3'].abs()-df_out['rain_direct_in_m3'].abs()


    ############################################################################################
    ############ compute lake surface elevation using water balance model:
    df_out['lake_vol_calc']=0 #- lake volume using water balance
    #df_out['vol_calc_WB_model']=0 #- lake volume using water balance

    #df_out['h_calc_WB_model']=0 #- lake surface height using volume

    df_out['Q_gw_bal_m3day']=0
    df_out['Q_gw_in_m3day']=0
    df_out['Q_gw_out_m3day']=0
    df_out['lake_mASL_calc']=0
    df_out['manning_inflow']=0
    df_out['manning_outflow']=0
    df_out['pump_inflow']=0
    df_out['pipe_outflow']=0

    df['Q_surf1']=0 #-surface flow boundary 1
    df['flow_reg_surf1']=0 #-flow regime of surface flow boundary 1


    #- prepare dummy well water table columns

    #df_out['b_sgmnt_%s' %well_id]=df_out['actual_mASL']-aquif_bot_elev #- effective thickness of aquifer
    #df_out['head_diff_%s' %well_id]=df_out['wtr_table_%s'%well_id]-df_out['actual_mASL']
    #df_out['Q_sgmnt_%s' %well_id]=k_md*df_out['b_sgmnt_%s' %well_id]*sgmnt_shore*df_out['head_diff_%s' %well_id]/well_dist_lake



    for well_id in segs:
        df_out['Q_sgmnt_%s' %well_id]=0


    #- function to compute dynamic mean K (Kf) of aquifer as a function of daily water table
    #- dynamic mean K changes with respect to water table and K of corresponding layers
    def K_dyn(well_id): #-K_dyn will be calculated as m/d
        df=pd.read_csv('data/data_gw/well_%s_profile.dat' %well_id, delimiter='\t')
        df['thickness_cum']=df['thickness'].cumsum() #- cumulative depth of each layer
        depth_tot=df['thickness'].sum() #- total depth of well profile (dry+wet)
        wtr_table_thickness=row['wtr_table_%s'%well_id]-aquif_bot_elev #- thickness of water table (wet)
        wtr_table_depth=depth_tot-(row['wtr_table_%s'%well_id]-aquif_bot_elev) #- depth of water table from well top edge
        #print(well_id,depth_tot,wtr_table_depth,wtr_table_thickness,row['wtr_table_%s'%well_id],aquif_bot_elev)
        wtr_table_i=df[df['depth']>= wtr_table_depth].index.values.astype(int)[0]#- water table layer index
        h_wet_wtr_layer=df['thickness_cum'][wtr_table_i]-wtr_table_depth #- thickness of wet part of water table layer
        K_prod_table_wet=df['k_cmd'][wtr_table_i]*h_wet_wtr_layer #- K x thickness of wet part of water table layer
        df_wet=df[(df['thickness_cum'] >= wtr_table_depth)].copy() #- other layers than water table layer
        df_wet['k_cmd*thickness']=df_wet['thickness']*df_wet['k_cmd'] #-K x depth of other wet layers         
        K_cm_day=(df_wet['k_cmd*thickness'].sum()+K_prod_table_wet)/wtr_table_thickness #- daily K (Kf) value based on water table
        return K_cm_day/100




    ##################################################################################################
    ###########- compute daily volume from water balance and corresponding lake height
    ##################################################################################################
    nrows = df_out.shape[0]
    nrows-=1

    with tqdm(total=df_out.shape[0], ascii=' █') as pbar:
        for i, row in df_out.iloc[1:nrows].iterrows():  # Iter begins from second row of dataframe (first row as initial conditions)
            pbar.update(1)
            pbar.set_description('%s' %row['date'])
            #print('%s' %row['date'], end='\r') 
            if SCS_surf_flow==0:
                row['SCS_R_m3']=0
            #- GW in/out: segmented model
            for well_id in well_IDs:
                remark=segs[well_id][1]
                K_md_mean=segs[well_id][2] #- mean K (Kf) value of aquifer (will be used in case no profile available)
                well_dist_lake=segs[well_id][3]
                sgmnt_shore=segs[well_id][4]
                aquif_bot_elev=segs[well_id][5]

                #- compute K_md of aquifer for each day, using water table
                #- check if well profile exists
                if os.path.isfile('%s/data/data_gw/well_%s_profile.dat' %(dname,well_id)):
                    K_md_day=K_dyn(well_id)
                else:
                    K_md_day=K_md_mean #- if no well profile, then mean K value is considered
            
                #- effective thickness of aquifer
                if i==1:  #-First day of calculation
                    df_out.loc[i,'b_sgmnt_%s' %well_id]=df_out['actual_mASL'].iloc[i-1]-aquif_bot_elev 
                    df_out.loc[i,'head_diff_%s' %well_id]=df_out['wtr_table_%s'%well_id].iloc[i-1]-df_out['actual_mASL'].iloc[i-1] #-delta_h
                else: #-Other days    
                    df_out.loc[i,'b_sgmnt_%s' %well_id]=df_out['lake_mASL_calc'].iloc[i-1]-aquif_bot_elev #- effect. thick. aquif.
                    df_out.loc[i,'head_diff_%s' %well_id]=df_out['wtr_table_%s'%well_id].iloc[i-1]-df_out['lake_mASL_calc'].iloc[i-1]
                
                #-calculate Q_segment of day
                df_out.loc[i,'Q_sgmnt_%s' %well_id]=K_md_day*df_out['b_sgmnt_%s' %well_id].iloc[i]*sgmnt_shore*\
                    df_out['head_diff_%s' %well_id].iloc[i]/well_dist_lake


            #- compute daily total GW in/out/balance in separate columns
            for well_id in segs:
                if df_out['Q_sgmnt_%s' %well_id][i]>0:
                    df_out.loc[i, 'Q_gw_in_m3day'] += df_out['Q_sgmnt_%s' %well_id][i]*xadjust_inflow
                elif df_out['Q_sgmnt_%s' %well_id][i]<0:
                    df_out.loc[i, 'Q_gw_out_m3day'] += df_out['Q_sgmnt_%s' %well_id][i]*xadjust_outflow
                
                df_out.loc[i, 'Q_gw_bal_m3day'] = df_out['Q_gw_in_m3day'][i]+df_out['Q_gw_out_m3day'][i]

            ###- compute lake volume ----#######################       
            #df_out.loc[0, 'lake_vol_calc']=lake_volume(df_out['actual_mASL'][0])
            if i==1:  #-First day of calculation
                df_out.loc[i,'lake_vol_calc']=lake_volume_mASL(df_out['actual_mASL'].iloc[i-1])\
                    +df_out['Q_gw_in_m3day'].iloc[i-1]+df_out['Q_gw_out_m3day'].iloc[i-1]\
                    +df_out['SCS_R_m3'].iloc[i-1]+df_out['rain_direct_in_m3'].iloc[i-1]\
                    -df_out['evap_m3'].iloc[i-1]-df_out['evapotransp_m3'].iloc[i-1]

                #-- add surface flow based on scenario:
                if surf_scenario=='SCHI_PU':
                    df_out.loc[i,'lake_vol_calc']=df_out['lake_vol_calc'].iloc[i]\
                        +surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1], inflow_pump, d_pipe)[0]\
                        -surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1], inflow_pump, d_pipe)[1]
                    df_out.loc[i,'pump_inflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1], inflow_pump, d_pipe)[0]
                    df_out.loc[i,'pipe_outflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1], inflow_pump, d_pipe)[1]

                elif surf_scenario=='SCHI_FLIO_PIPE':
                    df_out.loc[i,'lake_vol_calc']=df_out['lake_vol_calc'].iloc[i]\
                        +surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1], d_pipe)[0]\
                        -surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1], d_pipe)[1]
                    df_out.loc[i,'manning_inflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1], d_pipe)[0]
                    df_out.loc[i,'pipe_outflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1], d_pipe)[1]

                else:
                    df_out.loc[i,'lake_vol_calc']=df_out['lake_vol_calc'].iloc[i]\
                        +surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1])[0]\
                        -surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1])[1]
                    df_out.loc[i,'manning_inflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1])[0]
                    df_out.loc[i,'manning_outflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['actual_mASL'].iloc[i-1])[1]        
        
                #- calculate lake height corresponding to computed lake volume
                df_out.loc[i, 'lake_mASL_calc']=lake_height_mASL(df_out['lake_vol_calc'].iloc[i])

            else: #-Other days
                df_out.loc[i,'lake_vol_calc']=df_out['lake_vol_calc'].iloc[i-1]\
                    +df_out['Q_gw_in_m3day'].iloc[i-1]+df_out['Q_gw_out_m3day'].iloc[i-1]\
                    +df_out['SCS_R_m3'].iloc[i-1]+df_out['rain_direct_in_m3'].iloc[i-1]\
                    -df_out['evap_m3'].iloc[i-1]-df_out['evapotransp_m3'].iloc[i-1]

                #-- add surface flow based on scenario:
                if surf_scenario=='SCHI_PU':
                    df_out.loc[i,'lake_vol_calc']=df_out['lake_vol_calc'].iloc[i]\
                        +surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1], inflow_pump, d_pipe)[0]\
                        -surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1], inflow_pump, d_pipe)[1]
                    df_out.loc[i,'pump_inflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1], inflow_pump, d_pipe)[0]
                    df_out.loc[i,'pipe_outflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1], inflow_pump, d_pipe)[1]

                elif surf_scenario=='SCHI_FLIO_PIPE':
                    df_out.loc[i,'lake_vol_calc']=df_out['lake_vol_calc'].iloc[i]\
                        +surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1], d_pipe)[0]\
                        -surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1], d_pipe)[1]
                    df_out.loc[i,'manning_inflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1], d_pipe)[0]
                    df_out.loc[i,'pipe_outflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1], d_pipe)[1]

                else:
                    df_out.loc[i,'lake_vol_calc']=df_out['lake_vol_calc'].iloc[i]\
                        +surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1])[0]\
                        -surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1])[1]
                    df_out.loc[i,'manning_inflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1])[0]
                    df_out.loc[i,'manning_outflow']=\
                        surf_flow_dict[surf_scenario](df_out['date'].iloc[i-1],df_out['lake_mASL_calc'].iloc[i-1])[1]  
        
                #- calculate lake height corresponding to computed lake volume
                df_out.loc[i, 'lake_mASL_calc']=lake_height_mASL(df_out['lake_vol_calc'].iloc[i])
                
    ##- write output
    df_out=df_out.sort_values('date').drop_duplicates('date',keep='last')

    if not os.path.isdir('results_%s_%s' %(version,surf_scenario_plot)):
        os.mkdir('results_%s_%s' %(version,surf_scenario_plot), 0o666)

    with open('results_%s_%s/water_balance.csv' %(version,surf_scenario_plot),'w') as fout:
        df_out.to_csv(fout, sep=',', lineterminator='\n')

    #-remove last row of df_out (returns zero as the last day of calculation)
    df_out=df_out[:-1]

    ###########################################################################################
    ##### plot module
    ###########################################################################################
    if plot_module==1:
        print('\nPlotting results...')
        decades = mdates.YearLocator(10)   # every decade
        half_decades = mdates.YearLocator(5)   # every decade
        years = mdates.YearLocator(1)   # every year
        two_years = mdates.YearLocator(2)   # every year
        months = mdates.MonthLocator(6)  # every month
        yearsFmt = mdates.DateFormatter('%Y')

        #- exclude incomplete years of beginning and end (to compute annual means)
        if df_out['date'].dt.month[0] != 1:
            df_out=df_out[df_out['date'].dt.year > df_out['date'].dt.year[0]]
        if int(df_out.iloc[[-1]]['date'].dt.month) != 12:
            df_out=df_out[df_out['date'].dt.year < int(df_out.iloc[[-1]]['date'].dt.year)]

        df_out['date']=pd.to_datetime(df_out['date'])
        df_out_annu=df_out.resample('YE', on='date', label='left').mean(numeric_only=True).reset_index()
        df_out_sum_annu=df_out.resample('YE', on='date', label='left').sum(numeric_only=True).reset_index()

        df_out_mon=df_out.resample('ME', on='date', label='left').mean(numeric_only=True).reset_index()
        df_out_sum_mon=df_out.resample('ME', on='date', label='left').sum(numeric_only=True).reset_index()

        #####################################################################################
        #-- plot annual means
        #######################################
        fig, ax = plt.subplots(figsize=(10,10))
        plt.subplots_adjust(hspace=0.001)
        plt.suptitle('Annual means', fontsize=12)

        #-met data ##########################
        ax1=plt.subplot(511)
        ax=ax1
        ax.set_title('Met. data', y=1.0, pad=-14, fontsize=10)
        #line_L1=ax.plot(df_out_annu['date'],df_out_annu['TXK'], color='red',linewidth=0.7, label='max')
        line_L2=ax.plot(df_out_annu['date'],df_out_annu['TMK'], color='green',linewidth=0.7, label='Air temp. mean')
        #line_L3=ax.plot(df_out_annu['date'],df_out_annu['TNK'], color='deepskyblue',linewidth=0.7, label='min')
        #line4=ax.plot(df_out_annu['date'],df_out_annu['TGK'] color='green',linewidth=0.7, label='surf')
        ax.set_ylabel(r'Temp.($\degree$C)', fontsize=10)
        
        #- plot right axis
        ax_right = ax.twinx()
        line_R1=ax_right.plot(df_out_annu['date'],df_out_annu['UPM'], color='magenta',linewidth=0.7, label='rel. hum.')

        ax_right.set_ylabel('humid. (%)', fontsize=10)

        #- prepare legend
        
        lines=line_L2+line_R1
        #lines=line1+line2

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        #plt.title('met data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.axes.xaxis.set_ticklabels([])
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')



        #-hydro data #######################
        ax2=plt.subplot(512)
        ax=ax2
        ax.set_title('Rain and evap.', y=1.0, pad=-14, fontsize=10)
        line_L1=ax.plot(df_out_sum_annu['date'],df_out_sum_annu['SCS_R_m3']/1e6, color='blue',linewidth=0.7, label='Runoff')
        line_L2=ax.plot(df_out_sum_annu['date'],df_out_sum_annu['evapotransp_m3']/1e6, color='salmon', linewidth=0.5,\
         linestyle='-.', label='Evapotransp.')
        line_L3=ax.plot(df_out_sum_annu['date'],df_out_sum_annu['evap_m3']/1e6, color='chocolate',linewidth=0.7, \
            label='Evap. lake.')
        line_L4=ax.plot(df_out_sum_annu['date'],(0.001*df_out_sum_annu['RSK']*catchment_area)/1e6, color='turquoise',\
            linewidth=0.5, label='Rain')

        ax.set_ylabel('volume (x$10^{6}$ $m^{3}$ $d^{-1}$)', fontsize=10)
        #ax.tick_params(axis='y', colors='blue')

        #- plot right axis
        #ax_right = ax.twinx()
        
        #ax_right.set_ylabel('lake height (m)', fontsize=10)

        #- prepare legend
        
        lines=line_L1+line_L2+line_L3+line_L4
        #lines=line1+line2

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        #plt.title('hydro data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.axes.xaxis.set_ticklabels([])
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')


        #-level and volume data ##############
        ax3=plt.subplot(513)
        ax=ax3
        ax.set_title('Groundwater balance', y=1.0, pad=-14, fontsize=10)
        df_out_noNaN=df_out[df_out['pegel'].notnull()].copy()
        #df_out_noNaN=df_out.dropna()
        #print(df_out_noNaN); sleep(100)
        df_out_noNaN['date']=pd.to_datetime(df_out_noNaN['date'])
        df_out_noNaN_mean=df_out_noNaN.resample('YE', on='date', label='left').mean(numeric_only=True).reset_index()
        
        #- plot left axis
        #line_L1=ax.plot(df_out_annu['date'], df_out_annu['deltaGW_m3day']/1e6, color='teal', linewidth=1, label='GW bal. WB_bal')
        line_L2=ax.plot(df_out_annu['date'], df_out_annu['Q_gw_bal_m3day'], color='indigo', linewidth=1, label='GW balanance')
        ax.set_ylabel('$\Delta$GW ($m^{3}$ $d^{-1}$)', fontsize=9)
        ax.tick_params(axis='y')
        #ax.invert_yaxis()
        
        #- plot right axis
        ax_right = ax.twinx()
        line_R1=ax_right.plot(df_out_noNaN_mean['date'], df_out_noNaN_mean['actual_mASL'], color='plum', \
            linestyle='--', dashes=(7, 5), linewidth=0.5, label='lake surf measured.')
        #line_R2=ax_right.plot(df_out_noNaN_mean['date'], df_out_noNaN_mean['gw_elev_mASL'], color='brown', \
        #    linestyle='--',dashes=(10, 5), linewidth=0.5, label='GW surf.')

        ax_right.set_ylabel('Elevation (mASL)', fontsize=9)
        ax_right.tick_params(axis='y', colors='plum')

        #- prepare legend
        lines=line_L1+line_L2+line_R1

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        
        #plt.title('level data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')


        #- Groundwater data ###############
        ax3=plt.subplot(514)
        ax=ax3
        ax.set_title('Groundwater in/out', y=1.0, pad=-14, fontsize=10)
        df_out_noNaN=df_out[df_out['pegel'].notnull()].copy()
        #df_out_noNaN=df_out.dropna()
        #print(df_out_noNaN); sleep(100)
        df_out_noNaN['date']=pd.to_datetime(df_out_noNaN['date'])
        df_out_noNaN_mean=df_out_noNaN.resample('YE', on='date', label='left').mean(numeric_only=True).reset_index()
        
        #- plot left axis
        line_L1=ax.plot(df_out_annu['date'], df_out_annu['Q_gw_in_m3day']/1e6, color='green', linewidth=1, label='GW inflow')
        line_L2=ax.plot(df_out_annu['date'], df_out_annu['Q_gw_out_m3day']/-1e6, color='salmon', linewidth=1, label='GW out (GW model)')
        #line_L3=ax.plot(df_out_annu['date'], (df_out_annu['Q_gw_in_m3day']-df_out_annu['deltaGW_wb_m3day'])/-1e6\
        #    , color='salmon', linewidth=1, label='GW out (GW model)')
        ax.set_ylabel('GW flow (x$10^{6}$ $m^{3}$ $d^{-1}$)', fontsize=9)
        ax.tick_params(axis='y')
        #ax.invert_yaxis()

        #- plot right axis
        ax_right = ax.twinx()
        line_R1=ax_right.plot(df_out_noNaN_mean['date'], df_out_noNaN_mean['actual_mASL'], color='plum', \
           linestyle='--', dashes=(7, 5), linewidth=0.5, label='lake surf. measured')
        #line_R2=ax_right.plot(df_out_noNaN_mean['date'], df_out_noNaN_mean['gw_elev_mASL'], color='brown', \
        #    linestyle='--',dashes=(10, 5), linewidth=0.5, label='GW surf.')

        ax_right.set_ylabel('Elevation (mASL)', fontsize=9)
        ax_right.tick_params(axis='y', colors='plum')

        #- prepare legend
        lines=line_L1+line_L2+line_R1

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        
        #plt.title('level data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')



        #- plot well elevations ###############
        ax3=plt.subplot(515)
        ax=ax3

        #- create a random colormap
        color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                 for i in range(len(well_IDs))]
        
        #- create a random list of linewidths
        width = np.random.uniform(low=0.3, high=1, size=(len(well_IDs),))

        #- plot left axis
        for well_id, i in zip(well_IDs, range(len(well_IDs))):
            line=ax.plot(df_out_annu['date'], df_out_annu['wtr_table_%s' %well_id], color=color[i], linewidth=width[i], label='%s' %well_id)

        ax.set_title('Water table', y=1.0, pad=-14, fontsize=10)
        ax.set_ylabel('Water table (mASL)', fontsize=9)
        
        ax.legend(loc=1, frameon=True, fontsize=7)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')
        
        plt.savefig('results_%s_%s/wb_plot_ann.png' %(version,surf_scenario_plot), dpi=600)



        #####################################################################################
        #-- plot monthly means
        #######################################
        fig, ax = plt.subplots(figsize=(10,10))
        plt.subplots_adjust(hspace=0.001)
        plt.suptitle('monthly means', fontsize=12)

        #-met data ##########################
        ax1=plt.subplot(511)
        ax=ax1
        ax.set_title('Met. data', y=1.0, pad=-14, fontsize=10)
        line_L2=ax.plot(df_out_mon['date'],df_out_mon['TMK'], color='green',linewidth=0.7, label='Air temp. mean')

        ax.set_ylabel(r'Temp.($\degree$C)', fontsize=10)
        
        #- plot right axis
        ax_right = ax.twinx()
        line_R1=ax_right.plot(df_out_mon['date'],df_out_mon['UPM'], color='magenta',linewidth=0.7, label='rel. hum.')

        ax_right.set_ylabel('humid. (%)', fontsize=10)

        #- prepare legend
        
        lines=line_L2+line_R1
        #lines=line1+line2

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        #plt.title('met data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.axes.xaxis.set_ticklabels([])
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')



        #-hydro data #######################
        ax2=plt.subplot(512)
        ax=ax2
        ax.set_title('Rain and evap.', y=1.0, pad=-14, fontsize=10)
        line_L1=ax.plot(df_out_sum_mon['date'],df_out_sum_mon['SCS_R_m3']/1e6, color='blue',linewidth=0.7, label='Runoff')
        line_L2=ax.plot(df_out_sum_mon['date'],df_out_sum_mon['evapotransp_m3']/1e6, color='salmon', linewidth=0.5,\
         linestyle='-.', label='Evapotransp.')
        line_L3=ax.plot(df_out_sum_mon['date'],df_out_sum_mon['evap_m3']/1e6, color='chocolate',linewidth=0.7, label='Evap. lake.')
        line_L4=ax.plot(df_out_sum_mon['date'],df_out_sum_mon['RSK']*catchment_area/(1e6*1000), color='turquoise',\
            linewidth=0.5, label='Rain')

        ax.set_ylabel('volume (x$10^{6}$ $m^{3}$ $d^{-1}$)', fontsize=10)
        #ax.tick_params(axis='y', colors='blue')

        #- plot right axis
        #ax_right = ax.twinx()
        
        #ax_right.set_ylabel('lake height (m)', fontsize=10)

        #- prepare legend
        
        lines=line_L1+line_L2+line_L3+line_L4
        #lines=line1+line2

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        #plt.title('hydro data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.axes.xaxis.set_ticklabels([])
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')



        #-level and volume data ##############
        ax3=plt.subplot(513)
        ax=ax3
        ax.set_title('Groundwater balance', y=1.0, pad=-14, fontsize=10)
        df_out_noNaN=df_out[df_out['pegel'].notnull()].copy()
        #df_out_noNaN=df_out.dropna()
        #print(df_out_noNaN); sleep(100)
        df_out_noNaN['date']=pd.to_datetime(df_out_noNaN['date'])
        df_out_noNaN_mean=df_out_noNaN.resample('ME', on='date', label='left').mean(numeric_only=True).reset_index()
        
        #- plot left axis
        line_L1=ax.plot(df_out_mon['date'], df_out_mon['Q_gw_bal_m3day'], color='indigo', linewidth=1, label='GW balance')
        ax.set_ylabel('$\Delta$GW ($m^{3}$ $d^{-1}$)', fontsize=9)
        ax.tick_params(axis='y')
        ax.invert_yaxis()
        
        #- plot right axis
        ax_right = ax.twinx()
        line_R1=ax_right.plot(df_out_noNaN_mean['date'], df_out_noNaN_mean['actual_mASL'], color='plum', \
            linestyle='--', dashes=(7, 5), linewidth=0.5, label='lake surf.')
        #line_R2=ax_right.plot(df_out_noNaN_mean['date'], df_out_noNaN_mean['gw_elev_mASL'], color='brown', \
        #    linestyle='--',dashes=(10, 5), linewidth=0.5, label='GW surf.')

        ax_right.set_ylabel('Elevation (mASL)', fontsize=9)
        ax_right.tick_params(axis='y', colors='plum')

        #- prepare legend
        lines=line_L1+line_R1

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        
        #plt.title('level data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')


        #- Groundwater data ###############
        ax3=plt.subplot(514)
        ax=ax3
        ax.set_title('Groundwater in/out', y=1.0, pad=-14, fontsize=10)
        df_out_noNaN=df_out[df_out['pegel'].notnull()].copy()
        #df_out_noNaN=df_out.dropna()
        #print(df_out_noNaN); sleep(100)
        df_out_noNaN['date']=pd.to_datetime(df_out_noNaN['date'])
        df_out_noNaN_mean=df_out_noNaN.resample('ME', on='date', label='left').mean(numeric_only=True).reset_index()
        
        #- plot left axis
        line_L1=ax.plot(df_out_mon['date'], df_out_mon['Q_gw_in_m3day'], color='green', linewidth=1, label='GW inflow')
        line_L2=ax.plot(df_out_mon['date'], df_out_mon['Q_gw_out_m3day']*-1, color='salmon', linewidth=1, label='GW outflow')
        ax.set_ylabel('GW flow ($m^{3}$ $d^{-1}$)', fontsize=9)
        ax.tick_params(axis='y')
        #ax.invert_yaxis()
        
        #- plot right axis
        ax_right = ax.twinx()
        line_R1=ax_right.plot(df_out_noNaN_mean['date'], df_out_noNaN_mean['actual_mASL'], color='plum', \
            linestyle='--', dashes=(7, 5), linewidth=0.5, label='lake surf.')
        #line_R2=ax_right.plot(df_out_noNaN_mean['date'], df_out_noNaN_mean['gw_elev_mASL'], color='brown', \
        #    linestyle='--',dashes=(10, 5), linewidth=0.5, label='GW surf.')

        ax_right.set_ylabel('Elevation (mASL)', fontsize=9)
        ax_right.tick_params(axis='y', colors='plum')

        #- prepare legend
        lines=line_L1+line_L2+line_R1

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        
        #plt.title('level data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')


        #- plot well elevations ###############
        ax3=plt.subplot(515)
        ax=ax3

        #- create a random colormap
        color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                 for i in range(len(well_IDs))]
        
        #- create a random list of linewidths
        width = np.random.uniform(low=0.3, high=1, size=(len(well_IDs),))

        #- plot left axis
        for well_id, i in zip(well_IDs, range(len(well_IDs))):
            line=ax.plot(df_out_mon['date'], df_out_mon['wtr_table_%s' %well_id], color=color[i], linewidth=width[i], label='%s' %well_id)

        #- plot actual lake surface
        line_obs=ax.plot(df_out_mon['date'], df_out_mon['actual_mASL'], color='turquoise', linewidth=1, linestyle='-.', label='Lake act.')

        ax.set_title('Water table', y=1.0, pad=-14, fontsize=10)
        ax.set_ylabel('Water table (mASL)', fontsize=9)
        
        ax.legend(loc=1, frameon=True, fontsize=7)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')
        
        plt.savefig('results_%s_%s/wb_plot_mon.png' %(version,surf_scenario_plot), dpi=600)


        #####################################################################################
        #-- surf elev
        #######################################

        df_havel=pd.read_csv(havel_file, encoding='latin', delimiter=';')
        df_havel['date']=pd.to_datetime(df_havel['Datum'], format='%d.%m.%Y')
        df_havel=df_havel[['date','Tagesmittelwert']]
        df_havel=df_havel.replace(-777,np.nan)
        df_havel=df_havel.set_index('date')[['Tagesmittelwert']].interpolate().dropna().reset_index()
        df_havel['havel_mASL']=havel_gauge_zero+df_havel['Tagesmittelwert']/100
        #df_havel=df_havel.loc[(df_havel['date'])]

        fig, ax = plt.subplots(figsize=(10,4))
        plt.subplots_adjust(hspace=0.001)
        plt.suptitle('scenario: %s'%surf_scenario_plot, fontsize=10)

        #-level data ##############
        ax1=plt.subplot(211)
        ax=ax1
        ax.set_title('Surface height', y=1.0, pad=-14, fontsize=10)

        #- plot left axis
        line_L1=ax.plot(df_out['date'], df_out['lake_mASL_calc'], color='indigo', linewidth=1.3, label='Model (segmented)')
        #line_L2=ax.plot(df_out['date'], df_out['h_calc_WB_model'], color='blue', linewidth=1.5, label='Model (GWout: WB)')
        line_obs=ax.plot(df_out['date'], df_out['actual_mASL'], color='fuchsia', linewidth=0.6, linestyle='-', label='measured')
        
        #line_havel=ax.plot(df_havel['date'],df_havel['havel_mASL'],
        # color='violet', linestyle='-', linewidth=0.8, label='Havel')

        #-add standard deviation band:
        std=df_out['actual_mASL'].std()
        line_std=plt.fill_between(df_out['date'], df_out['actual_mASL']-std, df_out['actual_mASL']+std\
            , color='violet', alpha=0.2)


        ax.set_ylabel('surface height (m)', fontsize=9)
        ax.set_ylim(28.5,31.0)
        #- prepare legend
        lines=line_L1+line_obs

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)
        
        #plt.title('level data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')

        #-volume data ##############
        ax2=plt.subplot(212)
        ax=ax2
        ax.set_title('Lake volume', y=1.0, pad=-14, fontsize=10)

        #- plot left axis
        line_L1=ax.plot(df_out['date'], df_out['lake_vol_calc'], color='green', linewidth=1, label='Model (GWout: segmented)')
        #line_L2=ax.plot(df_out['date'], df_out['vol_calc_WB_model'], color='lime', linewidth=1.5, label='Model (GWout: WB)')
        line_obs=ax.plot(df_out['date'], df_out['lake_vol_measured'], color='tomato', linewidth=0.6, linestyle='-', label='measured')
        ax.set_ylabel('volume (m3)', fontsize=9)

        #- prepare legend
        lines=line_L1+line_obs

        labs=[l.get_label() for l in lines]
        ax.legend(lines, labs, loc=1, frameon=True, fontsize=7)

        ax.set_ylim(1.3e7,1.6e7)
        
        #plt.title('level data', fontsize=10)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(two_years)
        ax.xaxis.set_major_formatter(yearsFmt)
        
        ax.format_xdata = mdates.DateFormatter('%Y')

        plt.savefig('results_%s_%s/wb_plot_height_vol.png' %(version,surf_scenario_plot), dpi=600)


     #####################################################################################
        #-- plot monthly means
        #######################################
        fig, ax = plt.subplots(figsize=(10,4))
        plt.subplots_adjust(hspace=0.001)
        plt.suptitle('monthly means', fontsize=12)

        #-met data ##########################
        ax1=plt.subplot(111)
        ax=ax1

        #- create a library of segments to insert wells and segment data
        well_IDs_Havel = []
        files = [f for f in os.listdir('%s/data/data_gw/Havel_wells_orig_levels' %(dname)) if f.endswith('.dat')]
        for file in files:    
            with open('%s/data/data_gw/Havel_wells_orig_levels/%s' %(dname,file), 'r') as fin:
                for line in fin.readlines():
                    if 'well_id' in line:
                        well_id=int(re.findall(r'-?\d+\.?\d*', line)[0])   #- read well ID
                        well_IDs_Havel.append(well_id)

        #- create a random colormap
        color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                 for i in range(len(well_IDs_Havel))]
        
        #- create a random list of linewidths
        width = np.random.uniform(low=0.3, high=1, size=(len(well_IDs_Havel),))

        #- plot left axis
        for well_id, i in zip(well_IDs_Havel, range(len(well_IDs_Havel))):
            df_well=pd.read_csv('data/data_gw/Havel_wells_orig_levels/well_%s.dat' %well_id, delimiter='\t', skiprows=7) 
            df_well['date']=pd.to_datetime(df_well['date'].astype(str), format=fmt_de)
            line_well=ax.plot(df_well['date'], df_well['water_table'], color=color[i], linewidth=width[i], label='%s' %well_id)
        
        #- plot level of Havel
        df_Havel=pd.read_csv('data/data_gw/Havel_wells_dummy_levels/well_%s.dat' %well_id, delimiter='\t', skiprows=7)
        df_Havel['date']=pd.to_datetime(df_Havel['date'].astype(str), format=fmt_de)
        df_Hvl_2 = df_Havel[df_Havel.set_index(['date']).index.isin(df_well.set_index(['date']).index)]
        line_Havel=ax.plot(df_Hvl_2['date'], df_Hvl_2['water_table'], color='red', linewidth=1, label='Havel')
        
        #- plot actual lake surface
        #line_obs=ax.plot(df_out['date'], df_out['actual_mASL'], color='turquoise', linewidth=1, linestyle='-', label='Lake act.')

        ax.set_title('Water table , Havel', y=1.0, pad=-14, fontsize=10)
        ax.set_ylabel('Elevation (mASL)', fontsize=9)
        
        ax.legend(loc=1, frameon=True, fontsize=7)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(decades)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(half_decades)
        ax.format_xdata = mdates.DateFormatter('%Y')
        

        plt.savefig('results_%s_%s/Havel_wells.png' %(version,surf_scenario_plot), dpi=600)


    #####################################################################################
        #-- plot havel vs lake level
        #######################################
        fig, ax = plt.subplots(figsize=(10,4))
        plt.subplots_adjust(hspace=0.001)
        plt.suptitle('monthly means', fontsize=12)

        #-met data ##########################
        ax1=plt.subplot(111)
        ax=ax1

        #- create a library of segments to insert wells and segment data
        well_IDs_Havel = []
        files = [f for f in os.listdir('%s/data/data_gw/Havel_wells_orig_levels' %(dname)) if f.endswith('.dat')]
        for file in files:    
            with open('%s/data/data_gw/Havel_wells_orig_levels/%s' %(dname,file), 'r') as fin:
                for line in fin.readlines():
                    if 'well_id' in line:
                        well_id=int(re.findall(r'-?\d+\.?\d*', line)[0])   #- read well ID
                        well_IDs_Havel.append(well_id)

        #- plot level of Havel
        df_Havel=pd.read_csv('data/data_gw/Havel_wells_dummy_levels/well_%s.dat' %well_id, delimiter='\t', skiprows=7)
        df_Havel['date']=pd.to_datetime(df_Havel['date'].astype(str), format=fmt_de)
        df_Hvl_2 = df_Havel[df_Havel.set_index(['date']).index.isin(df_well.set_index(['date']).index)]
        line_Havel=ax.plot(df_Hvl_2['date'], df_Hvl_2['water_table'], color='red', linewidth=0.5, label='Havel')
        
        #- plot actual lake surface
        line_obs=ax.plot(df_out['date'], df_out['actual_mASL'], color='blue', linewidth=1, linestyle='-', label='Lake act.')

        ax.set_title('Lake , Havel', y=1.0, pad=-14, fontsize=10)
        ax.set_ylabel('Elevation (mASL)', fontsize=9)
        
        ax.legend(loc=1, frameon=True, fontsize=7)
        ax.grid(color='grey',linestyle=':', linewidth=0.2)
        ax.xaxis.set_major_locator(decades)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(half_decades)
        ax.format_xdata = mdates.DateFormatter('%Y')
        

        plt.savefig('results_%s_%s/Havel_lake.png' %(version,surf_scenario_plot), dpi=600)


    #- prepare aem3d boundary coditions

    #- groundwater

    df_gw=df_out.filter(regex='Q_sgmnt_')

    with open('results_%s_%s/gw_daily.out' %(version,surf_scenario_plot),'w') as fout:
        df_gw.to_csv(fout, sep='\t', index=False, lineterminator='\n')



if __name__ == '__main__':
    for p in all_processes:
      p.start()

    for p in all_processes:
      p.join()