#!/usr/bin/perl

# creates one 2D map climatology at specified level
# for annual/DJF/JJA data

###### ---------Script--------- ########
do 'user_input.s';

print "Doing clim_maps to make a climatology 2D map (lat/lon) \n";
chdir $DataDir;
@data_array = ("ANN","DJF","JJA");

#NOTE: DJF starts at yrini+1;

for my $months (@data_array){
$yrin = $yrini;
if ($months eq "DJF"){
  $yrin = $yrini+1;
}

$InputFileName = "$months$yrin-$yrend.acc$RUN.nc";
print "$InputFileName \n";

system "scaleacc $InputFileName $nctag";

$OutputFileName = "$variable.$months$yrin-$yrend.map_lev$ilev.$RUN.nc";
print "$OutputFileName.\n";

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable $months$yrin-$yrend.$nctag$RUN.nc $OutputFileName";
}else{
  system "ncks -O -v $variable $months$yrin-$yrend.$nctag$RUN.nc $OutputFileName";
}

}
