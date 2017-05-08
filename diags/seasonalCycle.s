#!/usr/bin/perl

print "seasonalCycle.s: computes the models global average seasonal cycle at a specified depth\n\n";

###### ---------Script--------- ########

chdir $DataDir;
$OutputFileName = "$variable.$yrini-$yrend.seasonalCycle_lev$ilev.$RUN.nc";
print "Output Name $OutputFileName\n";

@months = qw(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC);
print "Averaging annual cycle for years $yrini to $yrend \n";

foreach $month (@months) 
{
   system "scaleacc $month$yrini-$yrend.acc$RUN.nc $nctag";
   system "ncks -O -v $variable $month$yrini-$yrend.$nctag$RUN.nc dummy_$month.nc";
}
system "ncecat dummy_{JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC}.nc concat_dummy.nc";
system "rm -f dummy_*.nc";

#Rename missing_value to FillValue before computation
system "ncatted -O -a missing_value,$variable,d,f, concat_dummy.nc";
system "ncatted -O -a _FillValue,$variable,c,f,-1e30 concat_dummy.nc";   


#############
if (! -e "$area$underscore$RUN.nc") {
system "ncks -O -v $area JAN$yrini-$yrend.$nctag$RUN.nc $area$underscore$RUN.nc";
}
system "ncks -A -a -v $area $area$underscore$RUN.nc concat_dummy.nc";
system "ncwa -O -a $lat,$lon -w $area -d $lat,-90.0,90.0 -d $lon,-180.0,180.0 concat_dummy.nc dummy.nc";

# Remove unwanted terms
system "ncks -O -x -v $lon dummy.nc dummy.nc";
system "ncks -O -x -v $lat dummy.nc dummy.nc";

if (defined($depth)){  
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable dummy.nc $OutputFileName";
} else  {
  system "ncks -O -v $variable dummy.nc $OutputFileName";
}

system "ncrename -O -h -v $variable,seasonalCycle_ts $OutputFileName";
system "rm -f *dummy*.nc";



