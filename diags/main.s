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

##### -------- mean annual cycle diff from obs    ------- #####
chdir $myDir;
print "Doing seasonalCycleObs.s \n";
print "$ObsFilename \n";
do 'seasonalCycleObs.s';
print "\n";


chdir $myDir;
print "Doing diff_seasons.s \n";
do 'diff_seasons.s';
print "\n";

##### -------- climatology maps at certain level  ------- #####
chdir $myDir;
print "Doing clim_maps.s \n";
do 'clim_map.s';
print "\n";

##### -------- global averaged timeseries at a certain level   ------- #####
chdir $myDir;
print "Doing map_glbavg_ts.s \n";
do 'map_glbavg_ts.s';
print "\n";

##### -------- difference maps from observations   ------- #####
chdir $myDir;
print "Doing obs_map.s \n";
do 'obs_map.s';
print "\n";

chdir $myDir;
print "Doing diff_maps.s \n";
do 'diff_maps.s';
print "\n";

##### -------- vertical sections in different basins   ------- #####
chdir $myDir;
print "Doing basinAvg_obs.s \n";
do 'basinAvg_obs.s';
print "\n";

##### -------- diff basin vertical sections from obs   ------- #####

chdir $myDir;
print "Doing basinAvg_model.s \n";
do 'basinAvg_model.s';
print "\n";

chdir $myDir;
print "Doing diffBasinAvg.s \n";
do 'diffBasinAvg.s';
print "\n";

##### -------- AMOC vertical sections                  ------- #####
  chdir $myDir;
  print "Doing basinAvg_AMOC.s \n";
  do 'basinAvg_AMOC.s';
  print "\n";

##### -------- AMOC max at 26N timeseries              ------- #####
  chdir $myDir;
  print "Doing AMOCvertMax.s \n";
  do 'AMOCvertMax.s';
  print "\n";
