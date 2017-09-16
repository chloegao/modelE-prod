#!/usr/bin/perl

# creates a lat/dep difference plot from 3D model data - obs
# in different ocean basins, Pacific and Atlantic

##### -------- Inputs ------- #####
do 'user_input.s';

$Atl = "Atl";
$Pac = "Pac";
$ann = "ann";

if (index($variable_obs, "_mon") != -1) {
   print "'$variable_obs' contains mon.\n";
  $variable_obs =~ s/mon/ann/g; 
}

$FirstInput = "$variable.ANN$yrini-$yrend.Basin.$RUN.nc";
$varAtl1 = "$variable$Atl";
$varPac1 = "$variable$Pac";

$SecondInput = "$variable_obs$underscore$ann.Basin.$RUN2.nc";
$varAtl2 = "$variable_obs$underscore$ann$Atl";
$varPac2 = "$variable_obs$underscore$ann$Pac";
$OutputFileName = "$variable.ANN$yrini-$yrend.BasinDiff.$RUN$underscore$RUN2.nc";


###### ---------Script--------- ########
chdir $DataDir;

system "ncrename -v $varAtl2,$varAtl1 -v $varPac2,$varPac1 $SecondInput dummy.nc";
system "ncdiff -v $varAtl1 $FirstInput dummy.nc dummy1.nc";
system "ncdiff -v $varPac1 $FirstInput dummy.nc $OutputFileName";

system "ncks -A -a -v $varAtl1 dummy1.nc $OutputFileName";

system "rm -R -f dummy*";
