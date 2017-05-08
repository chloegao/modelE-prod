#!/usr/bin/perl

# lat/dep plot from 3D model data form ojl
# in different ocean basins, Pacific and Atlantic

$variable_amoc= "sf";
$nctag = ojl;
$Atl = "Atl";
$Pac = "Pac";
$underscore = "_";
$varAtl = "$variable_amoc$underscore$Atl";
$varPac = "$variable_amoc$underscore$Pac";

$lon = lono2;

###### ---------Script--------- ########

chdir $DataDir; 

$FileNameAtl = "ANN$yrini-$yrend.$varAtl$RUN.nc";
print "$FileNameAtl\n";
$FileNamePac = "ANN$yrini-$yrend.$varPac$RUN.nc";
print "$FileNamePac\n";
$OutputFileName = "$variable_amoc.ANN$yrini-$yrend.AMOCBasin.$RUN.nc";
print "$OutputFileName\n";

# If files don't exist, need to scale and pull out variable
if (! -e "ANN$yrini-$yrend.$nctag$RUN.nc"){system "scaleacc ANN$yrini-$yrend.acc$RUN.nc $nctag";}
if (! -e "ANN$yrini-$yrend.$varAtl$RUN.nc"){system "ncks -v $varAtl ANN$yrini-$yrend.$nctag$RUN.nc $FileNameAtl";}
if (! -e "ANN$yrini-$yrend.$varPac$RUN.nc"){system "ncks -v $varPac ANN$yrini-$yrend.$nctag$RUN.nc $FileNamePac";}

 system "rm -f dummy*.nc";
