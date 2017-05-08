#!/usr/bin/perl

# creates one 2D map climatology at specified level

##### -------- Inputs ------- #####
=pod
#$ACCDataDir = "/discover/nobackup/kflynn2/Ekfobio_test2/perlTest_Ekfobio_test2_2400-2415/";      #TO BE CHANGED BY USER
#$DataDir = "/discover/nobackup/kflynn2/FluxNCO/";       #TO BE CHANGED BY USER
#$RUN = "Ekfobio_test2";                               #TO BE CHANGED BY USER
#$yrini = 2400;                                           #TO BE CHANGED BY USER
#$yrend = 2402;                                           #TO BE CHANGED BY USER

#$variable = Gas_Exchange_CO2n;                                  #TO BE CHANGED BY USER
#$nctag = taij;

#$ilev = 1;   #define here level_index to do

#if ($nctag eq "taij")
#{
#$lat = lat;
#$lon = lon;
#$area = axyp;
#} elsif ($nctag eq "oijl")
#{
#$lat = lato;
#$lon = lono;
#$area = oxyp;
#$depth = zoc;
#}
#print "Your lat is $lat\n";

=cut
###### ---------Script--------- ########
print "Doing clim_maps to make a climatology 2D map (lat/lon) \n";
chdir $DataDir;

$InputFileName = "ANN$yrini-$yrend.acc$RUN.nc";
print "$InputFileName \n";

system "scaleacc $InputFileName $nctag";

$OutputFileName = "$variable.ANN$yrini-$yrend.map_lev$ilev.$RUN.nc";
print "$OutputFileName.\n";

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable ANN$yrini-$yrend.$nctag$RUN.nc $OutputFileName";
}else{
  system "ncks -O -v $variable ANN$yrini-$yrend.$nctag$RUN.nc $OutputFileName";
}


