use Env;

# -----------------------------------------------------------------------------
sub getEnvironment 
{
  my $env = shift;
  my $compiler = shift;
  my $branch = shift;

  if ($compiler eq "intel") {
    return getIntelEnvironment($env->{$compiler}, $branch);
  } elsif ($compiler eq "gfortran") {
    return getGfortranEnvironment($env->{$compiler}, $branch);
  } elsif ($compiler eq "nag") {
    return getNagEnvironment($env->{$compiler}, $branch);
  } else {
      print "$compiler is not a valid compiler.";
      exit;    
  }

}

# -----------------------------------------------------------------------------
sub getCommon
{
  my $env     = shift;
   if (defined $ENV{REGWORK}) {
      $env->{NOBACKUP} = $ENV{REGWORK};
   }
   else {
      $env->{NOBACKUP}=$ENV{NOBACKUP};
   }
  if (defined $ENV{MODELBASELINE}) {
      $env->{BASELINE_DIRECTORY} = $ENV{MODELBASELINE};
  }
  else {
      $env->{BASELINE_DIRECTORY} = "/discover/nobackup/modele/modelE_baseline";
  }
   if (defined $ENV{GCMSEARCHPATH}) {
      $env->{GCMSEARCHPATH} = $ENV{GCMSEARCHPATH};
   }
   else {
      $env->{GCMSEARCHPATH} = "/discover/nobackup/projects/giss/prod_input_files";
   }
  $env->{RESULTS_DIRECTORY} = $env->{NOBACKUP} . "/regression_results";
  $env->{SCRATCH_DIRECTORY} = $env->{NOBACKUP} . "/regression_scratch";
  $env->{DECKS_REPOSITORY}=$env->{SCRATCH_DIRECTORY} . "/decks_repository";
  $env->{CMRUNDIR}=$env->{SCRATCH_DIRECTORY} . "/cmrun";
  $env->{EXECDIR}=$env->{SCRATCH_DIRECTORY} . "/exec";
  $env->{SAVEDISK}=$env->{SCRATCH_DIRECTORY} . "/savedisk";
  $env->{MP}="no";
  $env->{OVERWRITE}="YES";
  $env->{OUTPUT_TO_FILES}="YES";
  $env->{VERBOSE_OUTPUT}="YES";
  return $env;
}
# -----------------------------------------------------------------------------
sub getIntelEnvironment
{
  my $env     = shift;
  my $branch     = shift;

  $env=getCommon($env);
  $env->{MPIDISTR}="intel";
  $env->{COMPILER}="intel";
  if ($branch =~ m/AR5/) 
  {
    $env->{BASELIBDIR}="/usr/local/other_old/esmf/2.2.2rp3_intel-10.1.017_impi-3.2.2.006/Linux";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_intel-10.1.013";
    $env->{PNETCDFHOME}="/discover/nobackup/mkelley5/pnetcdf-1.2.0";
  }
  else 
  {
    $env->{BASELIBDIR5}="/usr/local/other/esmf400rp1/intel12_impi32";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_intel-12.0.1.107";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/intel12.0.1.107_impi3.2.2.006";
  }
  $env->{MODELERC}=$env->{SCRATCH_DIRECTORY} . "/intel/modelErc.intel";
  return $env;
}

# -----------------------------------------------------------------------------
sub getGfortranEnvironment 
{
  my $env     = shift;
  my $branch     = shift;
    
  $env=getCommon($env);
  $env->{COMPILER}="gfortran";
  if ($branch =~ m/AR5/) 
  {
    $env->{MPIDISTR}="openmpi";
    $env->{BASELIBDIR}="/usr/local/other_old/esmf/2.2.2rp3_gcc-4.5_openmpi-1.4.2/Linux";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/gcc4.5_openmpi-1.4.2";
    $env->{PNETCDFHOME}="/discover/nobackup/mkelley5/pnetcdf-1.2.0";
  }
  else 
  {
    $env->{MPIDISTR}="openmpi";
    $env->{MPIDIR}="/usr/local/other/SLES11.1/openMpi/1.7.2/gcc-4.8.1-shared";
    $env->{BASELIBDIR5}="/usr/local/other/esmf400rp1/gcc4.7_openmpi-1.4.5";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/gcc-4.8.1_openmpi-1.7.2";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_gcc4.8.1";
  }
  $env->{MODELERC} = $env->{SCRATCH_DIRECTORY} . "/gfortran/modelErc.gfortran";
  return $env;
}

# -----------------------------------------------------------------------------
sub getNagEnvironment 
{
  my $env     = shift;
  my $branch     = shift;
    
  $env=getCommon($env);
  $env->{COMPILER}="nag";
  if ($branch =~ m/AR5/) 
  {
      print "NAG does not work with AR5 branch";
      exit;    
  }
  else 
  {
    $env->{MPIDISTR}="mvapich2";
    $env->{MPIDIR}="/usr/local/other/SLES11.1/mvapich2/1.8.1/nag-5.3-907";
    $env->{BASELIBDIR5}="";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/nag-5.3-886_mvapich2-1.8";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_nag-5.3";
  }
  $env->{MODELERC} = $env->{SCRATCH_DIRECTORY} . "/nag/modelErc.nag";
  return $env;
}

# -----------------------------------------------------------------------------
sub setupENVvariables
{
   my $env = shift;

   if (defined $ENV{MODELROOT}) {
      $env->{GIT_CLONE} = $ENV{MODELROOT};
   }
   else {
      print "ENV variable MODELROOT must be set.\n";
      exit;    
   }
   if (defined $ENV{GCMSEARCHPATH}) {
      $env->{GCMSEARCHPATH} = $ENV{GCMSEARCHPATH};
   }
   else {
      $env->{GCMSEARCHPATH} = "/discover/nobackup/projects/giss/prod_input_files";
   }
   if (defined $ENV{REGWORK}) {
      $env->{NOBACKUP} = $ENV{REGWORK};
   }
   else {
      $env->{NOBACKUP}=$ENV{NOBACKUP};
   }
   if (defined $ENV{GIT_REPOSITORY}) {
     $env->{GIT_REPOSITORY} = $ENV{GIT_REPOSITORY};
   }
   else {
     if (defined $ENV{MOCKMODELE}) {
       $env->{GIT_REPOSITORY} = "/discover/nobackup/ccruz/devel/modelE.mock.git";
     }
     else {
       $env->{GIT_REPOSITORY} = "simplex.giss.nasa.gov:/giss/gitrepo/modelE.git";
     }
   }
   $env->{RESULTS_DIRECTORY} = $env->{NOBACKUP} . "/regression_results";
   $env->{SCRATCH_DIRECTORY} = $env->{NOBACKUP} . "/regression_scratch";
   return $env;

}

# -----------------------------------------------------------------------------
sub saveForDiffreport()
{
   my $env = shift;
   my $cfgFile = shift;
# Save configuration settings in a format that is easily parsed by a bash script
   eval { require "$cfgFile"};
   if ($@) {
      print "Failed to load, because : $@"
   }
# Loop over configurations to get rundecks
  foreach my $deck ( keys %configurations )  {
    if ($deck =~ m/^(nonProduction)/i) {
       $deck = substr $deck, 14;
    }
    push (@rundecks, $deck);
  }

   my $rsize = scalar @rundecks;
   my $csize = scalar @compilers;

   my $file =  $env->{GIT_CLONE} . "/exec/testing/testsOutput/" . "." . "$cfgFile";
   open (FH, "> $file") or die "Can't open $file for write: $!";
   my $i = 0;
   while($i < $rsize)
   {
      print FH "DECK=$rundecks[$i]\n";
      $i++;
   }
   $i = 0;
   while($i < $csize)
   {
      print FH "COMPILER=$compilers[$i]\n";
      $i++;
   }   
   print FH "LEVEL=$level\n";
   print FH "BRANCH=$env->{BRANCH}\n";

   close FH or die "Cannot close $file: $!"; 
}

1;
