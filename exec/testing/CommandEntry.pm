#!/usr/bin/perl
# CommandEntry class
package CommandEntry;
use Env;
use Switch;

#  Constructor
sub new 
{
  my $proto = shift;
  my ($args) = @_;

  my $class = ref($proto) || $proto;
  my $self  = {};
  
  $self -> {COMMAND} = " ";
  $self -> {NUM_PROCS} = 1;
  $self -> {QUEUE} = LOCAL; # or BATCH
  $self -> {DEPENDENCY} = ""; # place holder for list of references to other commands
  $self -> {STATUS} = PENDING;
  $self -> {STDOUT_LOG_FILE} = "LOG";
  $self -> {STDERR_LOG_FILE} = "ERR";
  $self -> {COMPILER} = "";
  $self -> {RUNDECK} = "";
  $self -> {BRANCH} = "";

  # override defaults with arguments
  while ( my ($key, $value) = each(%$args) ) 
  {
    $self -> {$key} = $value;
  }
  # Make self an object in CommandEntry class and return reference to it
  bless ($self, $class);
  return $self;
}

sub dependency 
{
  my $self       = shift;
  return ($self -> {DEPENDENCY});
}

sub status 
{
  my $self = shift;
  return ($self -> {STATUS});
}

sub isComplete 
{
  my $self       = shift;
  return ( $self -> {STATUS} eq COMPLETE );
}

sub isRunning 
{
  my $self       = shift;
  return ( $self -> {STATUS} eq RUNNING );
}

sub isPending 
{
  my $self       = shift;
  return ( $self -> {STATUS} eq PENDING );
}

sub isReady 
{
  my $self       = shift;  
  my $pending    = ($self -> isPending());
  my $dependency = ($self -> dependency());
  if ($dependency eq "") { return $pending }
  else {return ($pending and ($self -> dependency() -> isComplete()) )};
}

sub setStatus 
{
  my $self      = shift;
  my $newStatus = shift;
  $self -> {STATUS} = $newStatus;
}

sub setDependency 
{
  my $self = shift;
  my $newDependency = shift;
  $self -> {DEPENDENCY} = $newDependency;
}

sub runInBatch 
{
  my $self          = shift;
  my $commandString = shift;
  my $queue         = shift;
  my $queueString   = " ";

  if ($queue) 
  {
    $queueString = "-q $queue\n";
  }

  my $NODE_SIZE = 12;
  my $nCPUS = $self -> {NUM_PROCS};
  my $nodes = 1 + int(($nCPUS-1)/$NODE_SIZE);
  my $ncpus = $nCPUS;
  $ncpus = $NODE_SIZE if ($ncpus > $NODE_SIZE);


  my $compiler = $self->{COMPILER};
  my $branch   = $self->{BRANCH};
  my $jobname  = $self->{RUNDECK};

  my $onEC2 = $ENV{RUNONEC2};

  my $doMock = $ENV{MOCKMODELE};

  if ($doMock == 1)
  {
    if    ($jobname =~ nonProduction) 
    { $walltime = "00:10:00"; }
    elsif ($jobname =~ M20) 
    { $walltime = "00:15:00"; }
    elsif ($jobname =~ obio || $jobname =~ cad || $jobname =~ amp) 
    { $walltime = "00:20:00"; }
    elsif ($jobname =~ tomas || $jobname =~ AR5_CAD) 
    { $walltime = "00:25:00"; $nodes=2; }
    else 
    { $walltime = "00:05:00"; $nodes=1; }
  }
  else
  {
    if ($jobname =~ M20 || $jobname =~ nonProduction || $jobname =~ C12) 
    { $walltime = "01:00:00"; }
    elsif ($jobname =~ obio) 
    { $walltime = "02:00:00"; }
    elsif ($jobname =~ cadi) 
    { $walltime = "03:00:00"; }
    elsif ($jobname =~ tomas || $jobname =~ AR5_CAD || $jobname =~ amp) 
    { $walltime = "03:00:00"; $nodes=4; }
    else 
    { $walltime = "03:00:00"; $nodes=1; }
  }
  print " runInBatch: doMock=$doMock, onEC2=$onEC2, COMPILER=$compiler\n";
  print " runInBatch: jobname=$jobname, walltime=$walltime, nodes=$nodes\n";

  # since nonProduction rundeck names can be quite long, extract the 
  # nonProduction_ part...
  if ($jobname =~ m/^(nonProduction)/i) {
     $jobname = substr $jobname, 14;
  }
  # ...then make sure PBS job name is only 15 characters long
  my $validPBSname = substr $jobname, 0, 14;

  my $script = <<EOF;
#!/bin/bash
#PBS -l select=$nodes:mpiprocs=12
#PBS -l walltime=$walltime
#PBS -N $validPBSname
#PBS $queueString
#SBATCH -A s1001

cd \$PBS_O_WORKDIR

. /usr/share/modules/init/bash
module purge

EOF

  if ($branch =~ m/AR5/) 
  {
    if ($compiler eq "intel") 
    {

      $script .= <<EOF;
module load comp/intel-10.1.017 mpi/impi-3.2.2.006
EOF

    } 
    elsif ($compiler eq "gfortran")  
    {
      $script .= <<EOF;
module load other/comp/gcc-4.5 other/mpi/openmpi/1.4.2-gcc-4.5
EOF
    }
    else 
    {
       # Nothing to do
    }

  } 
  else 
  {
    if ($compiler eq "intel") 
    {
      $script .= <<EOF;
module load comp/intel-13.1.3.192 mpi/impi-3.2.2.006
EOF
    } 
    elsif ($compiler eq "gfortran")  
    {
      $script .= <<EOF;
module load other/comp/gcc-4.8.1 other/mpi/openmpi/1.7.2-gcc-4.8.1-shared
EOF
    }  
    elsif ($compiler eq "nag")  
    {
      $script .= <<EOF;
module load comp/nag-5.3-907 other/mpi/mvapich2-1.8.1/nag-5.3-907
EOF
    }
    else 
    {
       return;
    }
  }

  $script .= <<EOF;
$commandString
date
env
EOF

  if ($onEC2 == 1)
  {
    `echo '$script' > pbstmp; chmod +x pbstmp; ./pbstmp; rm pbstmp`;
  }
  else
  {
    `echo '$script' > pbstmp; qsub -V pbstmp; rm pbstmp`;
  }
}

sub launch 
{
  my $self = shift;
  my $semaphore = shift;
  my $logFile = $self -> {STDOUT_LOG_FILE};
  my $logErr = $self -> {STDERR_LOG_FILE};

  my $commandString   = "(" . $self -> {COMMAND} . ") 2>&1 >> $logFile; chmod 664 $logFile; /usr/bin/touch $semaphore;";

  my $mode = $self -> {QUEUE};
  print " launch: queue=$mode\n";

  setModuleEnvironment($self->{COMPILER}, $self->{BRANCH});
  $ENV{MODELERC}=$self->{MODELRC};
  if ($self -> {QUEUE} eq LOCAL) 
  {
    `($commandString) &`; # run in background
  }
  else 
  {
    $self -> runInBatch("$commandString", $self -> {QUEUE});
  }
}

sub setModuleEnvironment 
{
    my $compiler = shift;
    my $branch = shift;

    require $ENV{MODELROOT}."/exec/testing/perlreq";
    print " setModuleEnvironment: COMPILER=$compiler, branch=$branch\n";

    module (purge);
    if ($branch =~ m/AR5/) 
    {
      if ($compiler eq "intel") 
      {
        module (load, "comp/intel-10.1.017",  "mpi/impi-3.2.2.006");
      } 
      elsif ($compiler eq "gfortran") 
      {
        module (load, "other/comp/gcc-4.5", "other/mpi/openmpi/1.4.2-gcc-4.5");
      } 
      else 
      {
        # Nothing to do
      }
    } 
    else 
    {

      if ($compiler eq "intel")
      {
        module (load, "comp/intel-13.1.3.192",  "mpi/impi-3.2.2.006");
      }
      elsif ($compiler eq "gfortran")
      {
        module (load, "other/comp/gcc-4.8.1", "other/mpi/openmpi/1.7.2-gcc-4.8.1-shared");
      }
      elsif ($compiler eq "nag") 
      {
        module (load, "comp/nag-5.3-907", "other/mpi/mvapich2-1.8.1/nag-5.3-907");
      } 
      else 
      {
        # Nothing to do
      }
    }
}

1;
