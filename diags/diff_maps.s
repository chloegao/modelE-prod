#!/usr/bin/perl

# creates a 2D plot difference of model - gridded observations
# from a specified level

###### ---------Script--------- ########
print "Doing diff_maps. Subtract model 2D climatology map from 2D obs map.\n";

$FirstFileName = "$variable.ANN$yrini-$yrend.map_lev$ilev.$RUN.nc";
print "First file $FirstFileName \n";

$SecondFileName = "$variable_obs.map_lev$ilev.$RUN2.nc";
print "Second file $SecondFileName \n";


chdir $DataDir;

$OutputFileName = "$variable.$yrini-$yrend.diffMap_lev$ilev.$RUN$underscore$RUN2.nc";

# Make all variables have the same name to use ncdiff
system "ncrename -v $variable_obs,diff $SecondFileName dummy2.nc";  
system "ncrename -v $variable,diff $FirstFileName dummy1.nc";                      

#rename lat/lon from lato/lono
if ($nctag eq "taij") {
  system "ncrename -O -v $lat_obs,lat -d $lat_obs,lat dummy2.nc";  
  system "ncrename -O -v $lon_obs,lon -d $lon_obs,lon dummy2.nc";
}else{
  system "ncrename -O -v $lat,lat -d $lat,lat dummy1.nc";  
  system "ncrename -O -v $lon,lon -d $lon,lon dummy1.nc";  
}

system "ncdiff -v diff dummy1.nc dummy2.nc $OutputFileName";
#system "ncatted -O -a long_name,V,o,c,$variable 'DIFFERENCE MODEL MINUS WOA13' $OutputFileName";
system "rm dummy*";

