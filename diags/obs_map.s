#!/usr/bin/perl

# creates one 2D map  at specified level

###### ---------Script--------- ########

chdir $DataDir;

if (index($variable_obs, "_mon") != -1) {
   print "'$variable_obs' contains mon.\n";
  $variable_obs =~ s/mon/ann/g; 
}

$InputFileName = "$RUN2.nc";
print "$InputFileName \n";

$OutputFileName = "$variable_obs.map_lev$ilev.$RUN2.nc";
print "$OutputFileName.\n";

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable_obs $ObsDir$InputFileName $OutputFileName";
}elsif (defined($mon)){
  system "ncks -O -v $variable_obs $ObsDir$InputFileName dummy.nc";
  system "ncwa -O -v $variable_obs -a mon dummy.nc dummy1.nc";
  system "ncks -O -x -v mon dummy1.nc $OutputFileName";
}else{
  system "ncks -O -v $variable_obs $ObsDir$InputFileName $OutputFileName";
}


