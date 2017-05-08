#!/usr/bin/perl

# creates a 2D plot difference of model - gridded observations
# from a specified level


$var = "seasonalCycle_ts";

$FirstFileName = "$variable.$yrini-$yrend.seasonalCycle_lev$ilev.$RUN.nc";
print "First file $FirstFileName \n";

$SecondFileName = "$variable_obs.seasonalCycle_lev$ilev.$RUN2.nc";
print "Second file $SecondFileName \n";

###### ---------Script--------- ########
chdir $DataDir;

$OutputFileName = "$variable.$yrini-$yrend.seasonalCycleDiff_lev$ilev$RUN$underscore$RUN2.nc";

# Make all variables have the same name to use ncdiff
system "ncrename -v $var,diff $SecondFileName dummy2.nc";  
system "ncrename -v $var,diff $FirstFileName dummy.nc";                      

# Need to add month dimension to model
system "ncrename -d record,mon dummy.nc dummy1.nc";
system "ncks -v mon dummy2.nc dummy_month.nc";
system "ncks -A -a -v mon dummy_month.nc dummy1.nc";


system "ncdiff -v diff dummy1.nc dummy2.nc $OutputFileName";
#system "ncatted -O -a long_name,V,o,c,$variable 'DIFFERENCE MODEL MINUS WOA13' $OutputFileName";
system "rm dummy*";

