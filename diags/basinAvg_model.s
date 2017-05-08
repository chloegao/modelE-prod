#!/usr/bin/perl

# creates a lat/dep plot from 3D model data
# of lat plot from 2D model data
# in different ocean basins, Pacific and Atlantic

###### ---------Script--------- ########

chdir $DataDir; 

$InputFileName = "ANN$yrini-$yrend.$variable$RUN.nc";
print "$InputFileName\n";
$OutputFileName = "$variable.ANN$yrini-$yrend.Basin.$RUN.nc";
print "$OutputFileName\n";

# If files don't exist, need to scale and pull out variable
if (! -e "ANN$yrini-$yrend.$nctag$RUN.nc"){
  system "scaleacc ANN$yrini-$yrend.acc$RUN.nc $nctag";
}
if (! -e "ANN$yrini-$yrend.$variable$RUN.nc"){
system "ncks -v $variable ANN$yrini-$yrend.$nctag$RUN.nc $InputFileName";
}


#match variable names from model and Mask
if ($lat eq "lat"){
system "ncks -v $variable $InputFileName modDummy.nc";
}else{
 system "ncrename -O -v $lat,lat -d $lat,lat $InputFileName dummy5.nc";
 system "ncrename -O -v $lon,lon -d $lon,lon dummy5.nc modDummy.nc";
}

##### ----- ATL BASIN ----- ########

print "Atlantic Basin\n";
$basin = "Atl";

# Rename variable in model so we can multiply the two files
system "ncrename -v $variable,AtlMask modDummy.nc";
if (defined($depth)){
  system "ncbo --op_typ=multiply AtlMaskO.nc modDummy.nc dummy.nc";
}else{
  system "ncbo --op_typ=multiply AtlMaskA.nc modDummy.nc dummy.nc";
}
# Average over longitude
system "ncwa -v AtlMask -a lon dummy.nc dummy1.nc";

# Rename variable from AtlMask
$variable_new = "$variable$basin";
print "$variable_new \n";
system "ncrename -v AtlMask,$variable_new dummy1.nc";

# Remove unwanted terms
system "ncks -O -x -v lon dummy1.nc $OutputFileName";

system "rm -R -f *dummy*";


##### ----- PAC BASIN ----- ########

print "Pacific Basin\n";
$basin = "Pac";

system "ncrename -v AtlMask,PacMask modDummy.nc";

if (defined($depth)) {
  system "ncbo --op_typ=multiply PacMaskO.nc modDummy.nc dummy.nc";
}else{
  system "ncbo --op_typ=multiply PacMaskA.nc modDummy.nc dummy.nc";
}

system "ncwa -v PacMask -a lon dummy.nc dummy1.nc";

# rename new variable
$variable_new = "$variable$basin";
print "$variable_new \n";
system "ncrename -v PacMask,$variable_new dummy1.nc";

# Remove unwanted terms
system "ncks -O -x -v lon dummy1.nc dummyPac.nc";

#combine
system "ncks -A -v $variable_new dummyPac.nc $OutputFileName";

system "rm -R -f *dummy* *modDummy*";

