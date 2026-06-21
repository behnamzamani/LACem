LACem model description:
------------------------
The Lake-Aquifer-Catchment exchange model (LACem) is a hydro(geo)logical model, which calculates the hydrological exchange between the lake and its surrounding hydrological system, which includes surface and subsurface hydrological exchanges, atmosphere-cahtchment exchanges and lake-atmosphere exchanges. 
This model is, in an initial step, written specifically for Lake Sacrower See in Berlin/Potsdam area in Germany. However, it may be used for other lakes. In Lake Sacrower See, the Schiffgraben canal, connecting the lake with the nearby Havel river, is an ungauged and rarely flooded canal. It is mostly flooded only when extremely high water levels occur in the Havel occur, or the water level of the lake is very high. The side walls and bottom of the lake are the interface of the lake water body with the surrounding aquifer. Hence, the solid boundary conditions of the model (side walls/bottom) are defined as the groundwater inflow/outflow boundary conditions. 
LACem represents different interface zones with respect to the available groundwater monitoring wells around the lake. The model takes into account and calculates the different components of the lake water balance equation at each time step (daily), where the lake-groundwater exchange is modeled actively by the USGS-segmented approach by (Rosenberry and LaBaugh, 2008), and each segment represents an area of the lake-aquifer interface affected by its corresponding observation well. The groundwater exchange is calculated using Darcy’s law, for which LACem calculates the horizontal hydraulic conductivity K of each segment as a weighted mean of the wet layers of the well profile for each time step (i.e. K varies dynamically as a function of time and groundwater level). Therefore, the model may be considered a semi-3D model, as it integrates several 1D components (segments) of lake-aquifer exchange. 
The surface exchanges include Evapotranspiration of the catchment (through Penman-Monteith approach), Evaporation from open surface (selectable between the two Penman and Priestley-Taylor methods) and surface water inflow/outflows via open channels (using Manning’s equation) and pipes (Darcy-Weisbach equation). Having added all the components mentioned in each time step, the Groundwater inflow/outflow direction at each segment is then calculated automatically by LACem at the lake-aquifer interface with respect to the water level difference between the lake and the corresponding well. The lake water level, as a control variable of the model, is then calculated from the lake volume of each time step. 
==============================================================================================================

Model structure:

Main Script:
-------------
There is a main script  "LACem_v1.31_surf_flow_2EvapModels_parallel.py", which is the 
only script to be run. The other 3 python scripts will be automatically run by the main script 
based on the model requirements. Settings of the model can be changed in the field of "model pre-settings" and "model function" fields. In a future version these model settings and data filenames will be optimized to be given in a separate file named "lacem_params.py".

Surface flow: The surface runoff can be calculated either using  There is a general calibration parameter "subsurface_flow_coeff" which adjusts the fraction of runoff which will not penetrate into groundwater and flow directly into the lake. 

Havel-Lake exchange settings :
------------------------------
the file "havel_lake_exchange.py" includes the hydraulic calculations of the scenarios designed to control and maintain the water lake in Lake Sacrower See using the water exchange between the lake and the Havel river via Schiffgraben canal. Do NOT change anything within this file!
The settings of the hydraulic exchange e.g. settings of the canal, hypolimnetic withdrawal pipe, pump etc. can be changed within the file havel_lake_params.py.

Data_structure:
---------------------------
The data must be stored within the subfolderst  of the "data" folder with respect to their data type:

data_lake: includes one file for the lake water level time series and one file for height-volume relationship (Hypsographic Curve).
 
data_weather: One file for daily weather data with the standard format of German Weather Service (DWD) for daily historical weather data. The solar radiation data will be stored in a separate file.

data_gw: groundwater data, includes the groundwater data of each segment. each segment corresponds to each observation well adjacent or in the neighborhood of the lake. For each segment there will be two files defined:
1. well_X.dat in ASCII format: X in the filename designates the well ID, the first 6 lines designate well ID (well_id, same with X), any remark, mean hydraulic conductivity (K_md_mean), distance of the lake to the shoreline (well_distance_lake), the length of the shoreline corresponding to this segment (segment_shoreline), elevation of the well bottom (aquif_bot_elev) in mASL. The following lines (beginning from line 8) are the time series of the water table elevation (in mASL) of the well. 

2. well_X_profile.dat: the well profile, which is a table including 4 headers of depth of the layer (depth), geologic type (type),	 thickness of the layer (thickness), hydraulic conductivity of the layer	(k_cmd).

The hydraulic conductivity of the entire aquifer at each segment is calculated in a dynamic way with respect to the water table elevation, so that only the k_cmd value of the wet layers below water table will be considered in  computing the mean hydraulic conductivity of the segment. In case no K_cmd value is available, the model calculates the mean conductivity given in as k_md in well_X.dat file.

data_hydraulic: includes the water levels of the nearby water bodies which are have water exchange with the modeled lake (e.g. other lakes, ponds or ungauged rivers).

data_landuse: land use data to be used in SCS rainfall-runoff model.
