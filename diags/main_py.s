#!/usr/bin/perl

#### so that there is not confusion over the years averaged

use Cwd;
$myDir = getcwd;
print "current directory $myDir \n";

##### -------- INPUTS ------- #####
chdir $myDir;
do 'user_input.s';
print "\n";

##### -------- Script ------- #####
##### -------- averaging files, annual means, climatologies, seasonal means  ------- #####
# copies over accfiles, sumfiles for annual means, climaological means
# and seasonal means, seasonal cycles and then deletes accfiles

=pod
chdir $myDir;
print "Doing avgACC.s \n";
do 'avgACC.s';
print "\n";

##### -------- mean annual cycle at a certain level  ------- #####
chdir $myDir;
print "nctag is $nctag\n";
print "Doing seasonalCycle.s \n";
do 'seasonalCycle.s';
print "\n";
$OutputFileName1 = "$OutputFileName";

##### -------- mean annual cycle diff from obs    ------- #####
chdir $myDir;
print "Doing seasonalCycleObs.s \n";
print "$ObsDir$ObsFilename \n";
do 'seasonalCycleObs.s';
print "\n";
$OutputFileName2 = "$OutputFileName";

chdir $myDir;
print "Doing diff_seasons.s \n";
do 'diff_seasons.s';
print "\n";

print "datadir =  $DataDir \n";
print "outputfilename1 =  $OutputFileName1 \n";
print "outputfilename2 =  $OutputFileName2 \n";

###### ---------PYTHON Script--------- ########
##invoke the python script to plot model seasonal cycle
chdir $myDir;
system "python3 plot_line.py $DataDir$OutputFileName1 $DataDir$OutputFileName2";
###### ------------------------------- ########


##### -------- climatology maps at certain level  ------- #####
chdir $myDir;
print "Doing clim_maps.s \n";
do 'clim_map.s';
print "\n";
$OutputFileName1 = "$OutputFileName";

##### -------- difference maps from observations   ------- #####
chdir $myDir;
print "Doing obs_map.s \n";
do 'obs_map.s';
print "\n";
$OutputFileName2 = "$OutputFileName";

chdir $myDir;
print "Doing diff_maps.s \n";
do 'diff_maps.s';
print "\n";
$OutputFileName3 = "$OutputFileName";

###### ---------PYTHON Script--------- ########
##invoke the python script
chdir $myDir;
system "python3 plot_map.py $DataDir$OutputFileName1 $DataDir$OutputFileName2 $DataDir$OutputFileName3";
###### ------------------------------- ########


##### -------- global averaged timeseries at a certain level   ------- #####
chdir $myDir;
print "Doing map_glbavg_ts.s \n";
do 'map_glbavg_ts.s';
print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script to plot model seasonal cycle
chdir $myDir;
system "python3 plot_linets.py glbAvg_ts $DataDir$OutputFileName";
###### ------------------------------- ########

##### -------- vertical sections in different basins   ------- #####
##### -------- diff basin vertical sections from obs   ------- #####
chdir $myDir;
print "Doing basinAvg_obs.s \n";
do 'basinAvg_obs.s';
print "\n";

chdir $myDir;
print "Doing basinAvg_model.s \n";
do 'basinAvg_model.s';
print "\n";

chdir $myDir;
print "Doing diffBasinAvg.s \n";
do 'diffBasinAvg.s';
print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script
chdir $myDir;
system "python3 plot_section.py $DataDir$OutputFileName";
###### ------------------------------- ########

##### -------- AMOC vertical sections                  ------- #####
  chdir $myDir;
  print "Doing basinAvg_AMOC.s \n";
  do 'basinAvg_AMOC.s';
  print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script
chdir $myDir;
print "$DataDir$FileNameAtl \n";
print "$DataDir$FileNamePac \n";
system "python3 plot_moc.py $DataDir$FileNameAtl $DataDir$FileNamePac";
###### ------------------------------- ########

##### -------- AMOC max at 26N timeseries              ------- #####
  chdir $myDir;
  print "Doing AMOCvertMax.s \n";
  do 'AMOCvertMax.s';
  print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script to plot model seasonal cycle
chdir $myDir;
system "python3 plot_linets.py sf_Atl $DataDir$OutputFileName";

###### ------------------------------- ########

=cut
##### -------- Current Transport timeseries ------- #####
## transports for Kuroshio, Gulf Stream and ACC
  chdir $myDir;
  print "Doing currentTrasp.s \n";
  do 'currentTrasp.s';
  print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script to plot model seasonal cycle

chdir $myDir;
$GS = GulfStream;
$KS = Kuroshio;
$DP = DrakesPassage;
system "python3 plot_text.py $DataDir$GS$yrini-$yrend.txt $DataDir$KS$yrini-$yrend.txt $DataDir$DP$yrini-$yrend.txt";
###### ------------------------------- ########

