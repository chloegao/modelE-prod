use Env;
use CommandPool;
use CommandEntry;
use RegressionTools;
use RegressionEnvs;

# This is the main modelE regression tests script. Its role is to setup the testing
# environment for rundecks specified in the configuration file (required argument).
# It also creates a "pool" of commands that is executed once the setup is done.

$num_args = $#ARGV + 1;
if ($num_args != 1) {
  print "\nUsage: regressionTests.pl config_filename\n";
  exit;
}

$cfgFile = $ARGV[0];
print "Configuration file name: $cfgFile\n";

# Get configurations and other options from configuration file
eval { require "$cfgFile"};
if ($@) {
  print "Failed to load, because : $@"
}

# Create references to local variables
my $branch = $gitbranch;
my $cleanScratch = $doCleanScratch;

# Initialize local data structures
my $env = {};
my $resolutions = {};
my $useCases = {};
my $numProcesses = {};
my @rundecks;

$env->{BRANCH} = $branch;
$env = setupENVvariables($env);

$machine = $ENV{HOST};
if ($machine =~ discover || $machine =~ borg)
{
  # On DISCOVER the default git is quite old so make sure we use the latest
  $ENV{PATH}="/usr/local/other/git/1.7.3.4_GNU/libexec/git-core:".$ENV{PATH};
}

foreach $compiler (@compilers) 
{
  $env->{$compiler} = getEnvironment($env, $compiler, $branch);
  $env->{$compiler}->{BRANCH} = $branch;
}

# Save settings for diffreport
&saveForDiffreport($env, $cfgFile);

# Resolution identifiers for rundeck suite
$resolutions->{nonProduction_E4TcadC12} = "8x10";
$resolutions->{nonProduction_E_AR5_C12} = "8x10";
$resolutions->{EM20}                    = "4x5";
$resolutions->{E1oM20}                  = "4x5";
$resolutions->{E4F40}                   = "2x2.5";
$resolutions->{E4TampF40}               = "2x2.5";
$resolutions->{E4TcadiF40}              = "2x2.5";
$resolutions->{E4arobio_h4c}            = "2x2.5";
$resolutions->{E4arobio_g6c}            = "2x2.5";
$resolutions->{E4TctomasF40}            = "2x2.5";
$resolutions->{E_AR5_CADI}              = "2x2.5";
$resolutions->{SGP4TESTS}               = "0";   # single column model 
$resolutions->{E4C90L40}                = "CS";  # cubed sphere

# numProcesses is resolution based (Note that domain decomposition is along y)  
$numProcesses->{"8x10"}->{GENTLE}      = [1,4];
$numProcesses->{"8x10"}->{AGGRESSIVE}  = [1,4,12];
$numProcesses->{"8x10"}->{INSANE}      = [22];
$numProcesses->{"8x10"}->{POLAR}      = [24];

$numProcesses->{"4x5"}->{GENTLE}       = [1,4];
$numProcesses->{"4x5"}->{AGGRESSIVE}   = [1,4,23];
$numProcesses->{"4x5"}->{INSANE}       = [23,44];
$numProcesses->{"4x5"}->{POLAR}       = [23,46];

$numProcesses->{"2x2.5"}->{GENTLE}     = [1,8];
$numProcesses->{"2x2.5"}->{AGGRESSIVE} = [1,45];
$numProcesses->{"2x2.5"}->{INSANE}     = [45,88];
$numProcesses->{"2x2.5"}->{POLAR}     = [46,90];

$numProcesses->{"0"}->{GENTLE}         = [];
$numProcesses->{"0"}->{AGGRESSIVE}     = [];
$numProcesses->{"0"}->{INSANE}         = [];
$numProcesses->{"CS"}->{GENTLE}        = [6];
$numProcesses->{"CS"}->{AGGRESSIVE}    = [6,48];
$numProcesses->{"CS"}->{INSANE}        = [6,84];
$numProcesses->{"2x2.5"}->{XLDECK}     = [44];  # for extra large (memory) rundecks
$numProcesses->{"2x2.5"}->{XLRUN}      = [88];  # for extra long (>=1mo) runs

# Loop over configurations and copy into local arrays
foreach my $deck ( keys %configurations )  {
  $cnt=0;  
  foreach my $cfg ( @{$configurations{$deck}} )  
  {
    $cnt=$cnt+1;
    if ($cnt == 1) {
      if ($cfg eq 'S') {
        $useCases->{$deck}->{CONFIGURATIONS} = ["SERIAL"];
      }
      elsif ($cfg eq 'M') {
        $useCases->{$deck}->{CONFIGURATIONS} = ["MPI"];
      }
      elsif ($cfg eq 'B') {
        $useCases->{$deck}->{CONFIGURATIONS} = ["SERIAL", "MPI"];
      }
      else {
        $useCases->{$deck}->{CONFIGURATIONS} = ["MPI"];
      }
    }
    elsif ($cnt == 2) {
      $useCases->{$deck}->{DURATION} = $cfg;
    }
    elsif ($cnt == 3) {
      $useCases->{$deck}->{DEBUGFLAGS} = $cfg;
    }
  }   
  $useCases->{$deck}->{NUM_MPI_PROCESSES} = $numProcesses->{$resolutions->{$deck}}->{$level};
  push (@rundecks,$deck);
}

# Create a new command pool:
my $pool = CommandPool->new();

# Create a "clean scratch space" command
# If we are running ALL the tests at the same time we probably want to clean up 
# the scratch space
if ($cleanScratch eq 'YES') {
  my $clean = CommandEntry->new({COMMAND => "rm -rf $env->{SCRATCH_DIRECTORY}/* $env->{RESULTS_DIRECTORY}/*/*;"});
  $pool->add($clean);
}

# git checkout command
my $git = CommandEntry->new(gitCheckout($env)); # Compiler is not important here
$pool->add($git);

# Write .modelErc command
foreach $compiler (@compilers) 
{
  $pool->add(writeModelErcFile($env->{$compiler}));
  unless (-d $env->{$compiler}->{RESULTS_DIRECTORY}."/$compiler") 
  {
    mkdir "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler" or die $!; 
  }
}

my $reference = "$env->{SCRATCH_DIRECTORY}" . "/$branch";
my $shortName;

print "Loop over all configurations...\n";
foreach my $rundeck (@rundecks) 
{ 
  foreach $compiler (@compilers) 
  {
    if ($compiler eq 'nag') 
    {
      $ENV{PATH}="/usr/local/other/SLES11.1/mvapich2/1.8.1/nag-5.3-907/bin:".$ENV{PATH};
      $ENV{LD_LIBRARY_PATH}="/usr/local/other/SLES11.1/mvapich2/1.8.1/nag-5.3-907/lib:".$ENV{LD_LIBRARY_PATH};
    }
    $env->{$compiler}->{RUNDECK} = $rundeck;

    foreach $configuration (@{$useCases->{$rundeck}->{CONFIGURATIONS}}) 
    {

      $env->{$compiler}->{CONFIGURATION} = $configuration;
      $env->{$compiler}->{DEBUGFLAGS} = $useCases->{$rundeck}->{DEBUGFLAGS};
      $env->{$compiler}->{DURATION} = $useCases->{$rundeck}->{DURATION};

      $shortName = $rundeck;
      if ($rundeck =~ m/^(nonProduction)/i) {
         $shortName = substr $rundeck, 14;
      }
      my $tempDir="$env->{SCRATCH_DIRECTORY}/$compiler/$shortName.$configuration";

      my $copy  = createTemporaryCopy($reference, $tempDir, $branch);
      $copy->{STDOUT_LOG_FILE} = "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler/$shortName.$configuration.$compiler.buildlog";
      $pool->add($copy, $git);

      my $build = compileRundeck($env->{$compiler}, $tempDir);

      $pool->add($build, $copy);
	    
      my $previous = $build; # for dependencies
	    
      my @peList = (1);
      if ($configuration eq "MPI") 
      {
	@peList = @{$useCases->{$rundeck}->{NUM_MPI_PROCESSES}};
      }
      foreach my $npes (@peList) 
      {
        my $run = runConfiguration($env->{$compiler}, $tempDir, $npes);
        $pool->add($run, $previous, $compiler);
        $previous = $run;
      }
    }
  }
}

print "\n***************************\n";
print "Starting regression tests.\n";
print "***************************\n\n";
$pool->run(*LOG, 0);
print "Completed modelE regression tests.\n";
close(LOG);

