#!/usr/bin/perl

#### so that there is not confusion over the years averaged

use Cwd;
$myDir = getcwd;

##### ----------Model-------- #####

$ACCDataDir = "/discover/nobackup/aromanou/TEST2/";      #TO BE CHANGED BY USER
$DataDir = "/discover/nobackup/aromanou/TESTNCO/";       #TO BE CHANGED BY USER
$RUN = E186f9aF40oQ32cc_2;                               #TO BE CHANGED BY USER
$yrini = 1995;                                           #TO BE CHANGED BY USER
$yrend = 2004;                                           #TO BE CHANGED BY USER

#### pot_temp
$variable = pot_temp;                                  #TO BE CHANGED BY USER
$nctag = "oijl";                                       #TO BE CHANGED BY USER
$ilev = 1;                                             #define here level_index to do 1=surface    
$lat = lato;                                           #TO BE CHANGED BY USER
$lon = lono;                                           #TO BE CHANGED BY USER
$area = oxyp;                                          #TO BE CHANGED BY USER
$depth = zoc;                                          #TO BE CHANGED BY USER

#### co2 flux
#$variable = Gas_Exchange_CO2n;                        #TO BE CHANGED BY USER
#$nctag = "taij";                                      #TO BE CHANGED BY USER
#$ilev = 1;                                            #define here level_index to do 1=surface    
#$lat = lat;                                           #TO BE CHANGED BY USER
#$lon = lon;                                           #TO BE CHANGED BY USER
#$area = axyp;                                         #TO BE CHANGED BY USER



##### -------Observations/Another Run----- #####

#### potential temperature
$ObsDir = $DataDir;                                    #TO BE CHANGED BY USER
$RUN2 = "WOA13_AnnMon_onEgrid";                        #TO BE CHANGED BY USER
$ObsFilename = "$RUN2.nc";
$variable_obs = "temp_mon";                            #TO BE CHANGED BY USER
$lon_obs = lon;                                        #TO BE CHANGED BY USER
$lat_obs = lat;                                        #TO BE CHANGED BY USER

#### co2 flux
#$ObsDir = $DataDir;                                   #TO BE CHANGED BY USER
#$RUN2 = "Takahashi_onEgrid";                          #TO BE CHANGED BY USER
#$ObsFilename = "$RUN2.nc";
#$variable_obs = "flux_molCm2yr_tak";                  #TO BE CHANGED BY USER
#$lon_obs = lona;                                      #TO BE CHANGED BY USER
#$lat_obs = lata;                                      #TO BE CHANGED BY USER
#$mon = mon;

$underscore = "_";

