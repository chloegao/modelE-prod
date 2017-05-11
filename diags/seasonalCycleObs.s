#!/usr/bin/perl

# computess the global average seasonal cycle

###### ---------Script--------- ########
print "Doing SeaonsalCycleObs. Calculates the lat/lon average at a specifc depth per month. \n";
do 'user_input.s';

if ($nctag eq "taij"){
print ("We use flux in mol, (C/m2/yr)\n");
}

$computes = "seasonalCycle_ts";
$new_variablename = "$variable$underscore$computes";
chdir $DataDir;
$OutputFileName = "$variable_obs.$computes.lev$ilev.$ObsFilename";
print "$OutputFileName\n";

system "ncks -O -v $area $area$underscore$RUN.nc dummy1.nc";
system "ncrename -O -h -d $lat,$lat_obs -v $lat,$lat_obs dummy1.nc";
system "ncrename -O -h -d $lon,$lon_obs -v $lon,$lon_obs dummy1.nc";


if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -A -F -d $depth,$ilev,$ilev,1 -v $variable_obs $ObsDir$ObsFilename dummy1.nc";
}else{
  system "ncks -A -v $variable_obs $ObsDir$ObsFilename dummy1.nc";
}

system "ncwa -h -O -v $variable_obs -w $area -a $lat_obs,$lon_obs dummy1.nc $OutputFileName";

# Remove unwanted terms
system "ncks -O -x -v $lon_obs $OutputFileName $OutputFileName";
system "ncks -O -x -v $lat_obs $OutputFileName $OutputFileName";

system "ncrename -O -h -v zoc,dep -d zoc,dep $OutputFileName";
system "ncrename -O -h -v $variable_obs,$new_variablename $OutputFileName";
system "rm -f dummy*.nc";



