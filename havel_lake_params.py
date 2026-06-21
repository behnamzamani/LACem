#-script to read canal data
#- data are my original values (Schiffgraben data rom 1m DEM)

##############-settings for Schiffgraben / see havel_lake_exchange.py for info--------####
##-elevations all in mASL
havel_gauge_zero=28
lake_gauge_zero=28.149
lake_bot_elev=-7

mud_thresh=0.0 #-mud threshold to subtract from bed levels (orig 0.0)
schifgr_havel_bed=29.57 #- mASL western side of dam 
schifgr_dam_bed=29.68 #- mASL canal bed at dam (pipe bed, Sec FF)
schifgr_highest_elev=29.90 #- mASL highest point (midpoint) of schiffgraben (where flow into lake/havel begins) was 29.96 
schifgr_lake_bed=29.8 #-mASL mouth to lake (was 29.84)

d_pipe_dam=0.3 #-diameter of each pipe under dam (there are two pipes)
l_pipe_dam=15 #-pipe length m (under dam)
f_pipe_dam=0.02 #-pipe friction factor
n_pipe_dam=0.025 #-pipe roughness

flow_occur_thresh=0.01 #-threshold of flow occurrence at schifgr_bed (due to surface tension/waves) 

#-channel geometry at control section (Sec EE)
l=640 #- total channel length (m)
lx=493 #-length of channel at control section (m) Sec EE in dwg file; distance from lake side  
b=0.5 #-channel bed width at control section (m) pipe under dam considered (rough value)

alpha=75 #-side angles of walls (degrees) Section EE = 75deg (in dwg file)
n=0.025 #-channel roughness (From Chow's book for earth channels)

#-channel geometry at mid point (sec CC)
l_mid=211 #- channel length of schifgr_highest_elev (m) at Sec CC (dwg file) dist. from lake side 
alpha_mid=36 #side angles of walls (degrees) at lmid Section CC = 47deg (in dwg file)
b_mid=53 #- channel bed width at midpoint (sec CC)


##############-settings for hypolimnetic withdrawal pipe -----------------------------####
l_pipe=2800 #-outflow pipe length m (hypolimnetic withdrawal)
d_pipe=0.5 #-outflow pipe diameter m
f_pipe=0.02 #-pipe friction factor

#-settings for SCHI-PU scenario
inflow_pump=60 #-inflow pump rate l/s