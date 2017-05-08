#!/usr/bin/perl

# Computes the global average seasonal cycle

##### -------- Inputs ------- #####
do 'user_input.s';

###### ---------Script--------- ########
print "Doing SeaonsalCycleObs. Calculates the lat/lon average at a specifc depth per month. \n";

if ($nctag eq "taij"){
print ("We use flux in mol, (C/m2/yr)\n");
}

chdir $DataDir;
$OutputFileName = "$variable_obs.seasonalCycle_lev$ilev.$ObsFilename";
print "$OutputFileName\n";

system "ncks -O -v $variable_obs $ObsDir$ObsFilename dummy1.nc";

#############
# add oxyp to file
print "Add $area to file, for weighting\n";
system "ncks -A -a -v $area $area$underscore$RUN.nc dummy1.nc";
system "ncwa -O -a $lat_obs,$lon_obs -w $area -d $lat_obs,-90.0,90.0 -d $lon_obs,-180.0,180.0  dummy1.nc dummy.nc";

# Remove unwanted terms
system "ncks -O -x -v $lon_obs dummy.nc dummy.nc";
system "ncks -O -x -v $lat_obs dummy.nc dummy.nc";

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable_obs dummy.nc $OutputFileName";

}else{
  system "ncks -O -v $variable_obs dummy.nc $OutputFileName";
}

system "ncrename -O -h -v $variable_obs,seasonalCycle_ts $OutputFileName";
system "rm -f dummy*.nc";


