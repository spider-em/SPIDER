#!/usr/bin/perl
#
# <html><head><title>Publish a job</title></head><body><pre>
#
# SOURCE: /usr8/spider/pubsub/publish.perl 
#         New ArDean Leith                                 Feb. 2002  
#         TRACE log added                                  Jan. 2008
#         ($jobidt !=  0)&&($jobidt != -$jobid)            Jan. 2008
#         NFS LOG locking                                  Apr. 2009
#         Remote publish support                           Jul. 2010
#
# PURPOSE: Submit job to PubSub publisher que. 
#
# USAGE EXAMPLE:  publish.perl  "spider pam/acn @apmq 17 grp=17" 
#
# INPUT:
#    -s                        (Optional sync file created when command starts )
#    Command listing           (argument )
#
# OUTPUT: 
#    Success or failure        (returned)  

use Getopt::Std;
use Fcntl;

local($jobid, $runcommand, $n, $machine);

$dir_pubsub = $ENV{'PUBSUB_DIR'};

# Set que location & name
if ($dir_pubsub) 
   {$que = "$dir_pubsub/pubsub.que";}    
else
   {die "Environmental variable: PUBSUB_DIR is not defined \n";}

# Set PubSub usage trace log location & name
$trace = "$dir_pubsub/pubsub.trace";    

# Get que host
if ($ENV{PUBSUB_MASTER}) 
   {$quehost = $ENV{PUBSUB_MASTER};}
else
   {die "Environmental variable: PUBSUB_MASTER is not defined \n";}
#print  "Current que host: $quehost \n";

# Get current directory location
$locnow = $ENV{'PWD'};

# Set user's log location & name (in current directory)
$log = $locnow . "/pubsub.log";

# Get current PID
$jobid  = $$;                    # Current PID is in $$

# Get currentusername
$user   = $ENV{'USER'};          # Current user name from sys. env. 

# Get optional name of sync file to be created when job  runs on subscriber 
getopt("s",);                    # Removes this argument from the line

#Get command to be run from argument line
$runcommand = join(" ",@ARGV);
#print " Runcommand:  $runcommand  \n";

#Get current host computer name
$hostname = $ENV{'HOST'};       # Current host computer from sys. env

if ( $hostname ne $quehost )
   {                            # Host is NOT PubSub master, run it on master
   $do = "ssh $quehost \" '\'cd $locnow ; ./publish.perl -s wait_$jobid $runcommand \" &";
   print " From: $hostname, Publishing: $do \n";
   system($do);
   
   $checkstart = 4;             # Start check time    (seconds)
   $maxstart   = 30 * 60;       # Max start wait time (seconds)
   $waited     = 0;             # Elapsed start waiting time 
   $waitfor    = "wait_$jobid";
   
   #print " From: $hostname, Waiting for : $waitfor \n";
   while ( !( -r $waitfor) && ($waited < $maxstart) )
      {  # Wait till published OK signal file exists or timeout 
      select(undef,undef,undef,$checkstart);    # Waits $checkstart seconds
      $waited = $waited + $checkstart;
      }
      
   if (-r $waitfor) 
      { unlink($waitfor);}      # Remove published OK signal file
   else
      { print " Timeout from remote publish after: $waited sec. \n";}
      
   #print " From: $hostname, Found : $waitfor \n";
   exit;
   }

$checkwait = 3;            # Loop wait time (seconds)
$ntry_mach = 0;            # Number of machines tried so far
$started   = 0;

# Find first OS command in $runcommand
$runcommand1 = $runcommand;
$runcommand1 =~ m/(\w+)/o ;
$run1        = $1;
#print "run1: $run1: \n"; 

# Substitute ';' for any ; in $runcommand to preserve ';' in rsh
$runcommands = $runcommand;
$runcommands =~ s/\;/';'/g; 

# Substitute $jobid for any PUBSUB_JOBID in $runcommand
$runcommands =~ s/PUBSUB_JOBID/$jobid/g; 

while (($started == 0) && ( $ntry_mach < 10 ) )
   { # Loop until started or >= 10 machines tried (arbitrary limit)

   #print "In publish; User: $user Location: $locnow\n";
   # Publish job now by listing in publisher que
   &publish("$que","$jobid","$user","$locnow");

   # Record qued time for use in log line 
   $quedtime = time();
   #system(" echo Waiting for sub to  job: $jobid for: $user \n");

   $ntry_mach++;    # Increment counter for number of machines tried

   # Wait for a subscription, find machine that subscribed
   $machine = '';
   while (!$machine)
      { 
      if ( ! open(PUBQUE, "+< $que")) 
          { system ("echo publish can not open: $que \n") ; return 0; }

      # Shared lock= 1, Exclusive lock=2, Non-blocking request=4,  Free lock=8
      unless (flock(PUBQUE, 1))
         { die " Can not read_lock: $que \n";}

      seek(PUBQUE,0,0) || die "Can't seek to start of: $que \n"; #al oct07

      # Read the  queue and see if this job is assigned to a machine yet
      while (<PUBQUE>)       #read input
         {
         $linet = $_;
         #print  "linet: $linet \n";

         chop($linet);
         ($jobidt,$user,$machval,$rest) = split(/\b\s+/,$linet,4);

         if ($jobidt == -$jobid)
            { # This job is locked, need to start this job on $machval
            $machine = $machval;
            #print "machine: $machine \n"; 
            $line = $linet;
            last;               # No need to read further
            }
         } # End of:  while (<PUBQUE>)   

      # Unlock the que file 
      flock(PUBQUE,8);
      close(PUBQUE);

      if (!$machine) 
         { # Wait $checkwait seconds and loop (will loop forever until $machine)
         select(undef,undef,undef,$checkwait);    # Waits $check seconds
         }

      }  # End of: while ($machine)

   # Machine found -------------------------------------------------------------

   #print "Publish; jobidt: " . $jobidt . "machine:" . $rest . "\n\n";

   # See if $locnow is visable from the subscriber machine
   $tolocrun = "rsh -n $machine \"ls -d $locnow\" 2>&1 ";
   #print "tolocrun:  $tolocrun  \n"; 

   $gotlocrun = `$tolocrun`;     
   $gotlocrun =~ s/^\s+//;    # Delete initial whitespace
   $gotlocrun =~ s/\s+$//;    # Delete initial whitespace

   #print "gotlocrun: $gotlocrun: \n";
   #print "locnow   : $locnow: \n";
   $locnowpsl = $locnow . "/";    # Add final slash (some PWD returns)

   if (($locnow ne $gotlocrun) && ($locnowpsl ne $gotlocrun))
      {
      system("echo ERROR Can not access: $locnow on: $machine");
      #system("echo ERROR :$locnow: NOT SAME AS :$gotlocrun:");
      del_locked_job($jobid);      # Delete this locked job from que 

      # try next machine
      next;     
      }
 
   $torun = "rsh $machine cd $locnow ';' $runcommands ';' rsh -n $quehost \"$dir_pubsub/delete.perl $jobid $locnow/pubsub.log\" &";
   #print "Starting: $torun \n"; 
   system("$torun");
   $torun =~ s/';'/;/g;           # Get rid of ' around ';' 
   #print "Running: $torun \n"; 

   # List starting time in log line 
   $starttime = time();
     ($sec,$min,$hour,$mday,$month,$year,$more) = localtime($starttime); 
   $time = sprintf("%04d-%02d-%02d-%02d:%02d:%02d", 
                       $year+1900,$month+1,$mday, $hour,$min,$sec);

   #system ("echo : locking log file now");
   # Open log file for appending
   if ( ! open(LOG, ">> $log"))
      { system ("echo In publish; can not open: $log \n") ; return 0; }

   # Lock the log file
   ### KLUDGE for ALBANY old kernal on usr11 host not NFS lockable
   $nolock =  ($log =~ m/\/usr11/ );
   if (! $nolock )
      { #print " LOCKING\n"; 
      $status = lock_file(LOG,"w");
      }
   ### End KLUDGE for ALBANY old kernal on usr11 host

   # Copy QUE info to log file and append job-group, starting time & command
   print LOG "$jobid $machine $quedtime $starttime (Start: $time) (Oper: $runcommand) \n";

   # Unlock the log file 
   if (! $nolock )
      { #print " UNLOCKING\n"; 
      $status = unlock_file(LOG);
      }
   close(LOG);

   $ntry        = 0;

   while ( ($started == 0) && ( $ntry < 20 ) )
      #while ( ($started == 0) && ( $ntry < 4 ) )
      {  # Loop while job not running yet, or until ntry >= 20

      #$search = "rsh -n $machine \"ps -u $user -o cmd | grep $run | grep -v tcsh | wc -l\" 2>&1";
      #$procs = `$search` ;

      #see if command is running yet on $machine
      @procst = `rsh -n $machine \"ps -u $user -o cmd\"` ;
      #print "Procst returns: @procst ";

      if ($procst =~ m/[A-Za-z]+/x )
          {  # Something wrong with 'rsh' call
          system("echo ERROR: $procst");
          last;
          }

      @proc = grep { /.$run1.*/ } @procst;
      # print "proc: @proc \n";
      $nprocs = $#proc + 1;
      #print "Found: $nprocs process containing: $run1  on: $machine \n";

      if ($nprocs > 0 )
          {  # runcommand has started, exit loop
          #system("echo runcommand has started with procs: $procs");
          $started = 1;
          }

      else
          {  # See if jobid has already finished and been deleted from que

          #$str = "grep $jobid $log | grep Runtime";
          #print "See if jobid has finished: $str \n";
          $str = "grep $jobid $log | grep Runtime | wc -l 2>&1";
          $started = `$str`;
          #print " Job: $jobid Finished already: $started \n" ;
          }

      if ( $started < 1)
         {
         # Wait $checkwait seconds and loop to see if started after wait
         select(undef,undef,undef,$checkwait); # Waits $check seconds
         $ntry++;                              # Increment starting try counter
         #print "Inner looping: $ntry";
         }

      } # End of: while ( !$started && ( $ntry < 18 ) )
 
   # Delete this job from que now (as started, over, or dead)
   del_locked_job($jobid);      # Delete this job from que 

   #$started=0; # USE THIS AS TEST FOR BADNODES

   if ( $started < 1)
      {  # Try another machine
      #system("echo Outer looping: $ntry_mach");
      print " Job: $jobid Re-published as it failed to run on: $machine\n" ;
      }

   }  # END of: while (($started == 0) && ( $ntry_mach < 20 ) ) xxxxxxxxxx


if  ( ($started > 0) && ( $opt_s) )
   { #print "Creating sync file: $opt_s \n";
   system("touch $opt_s"); }       # Create the requested sync file (empty)


# Open trace file for appending
if ( ! open(TRACE, ">> $trace"))
   { system ("echo In publish; can not open: $trace \n"); }
else
   {
   # Lock the trace file
   #flock(TRACE, 2);
   print TRACE " Starting:   $time $jobid $user $machine $locnow ($runcommand) \n";
   #flock(TRACE,8);
   close(TRACE);
   }


# ---------------------- publish --------------------------------------
sub publish 
    {
    if (@_ > 5)
       {die  "Usage: $0 que jobid user directory  
       Where:
       que is:        publisher que file \n
       jobid is:      job PID \n
       user is:       username \n 
       directory is:  current working directory \n"
       }

    $quefile = $_[0];
    #print "quefile: $quefile \n";

    $jobid = $_[1];
    #print "jobid: $jobid \n";

    $user = $_[2];
    #print "user: $user \n";

    $locnowt = $_[3];
    #print "locnowt is: $locnowt \n";
    
    # Open Publisher que file for appending
    open(PUBQUE, ">> $quefile") ||
        open(PUBQUE, "> $quefile") ||
           die "Publish can not open: $quefile \n";

   # Open trace file for appending
   $pubtime = time();
   ($sec,$min,$hour,$mday,$month,$year,$more) = localtime($pubtime); 
   $ptime = sprintf("%04d-%02d-%02d-%02d:%02d:%02d", 
                    $year+1900,$month+1,$mday, $hour,$min,$sec);
   {$trace = "$dir_pubsub/pubsub.trace";}
   if ( ! open(TRACE1, ">> $trace"))
      { system ("echo In publish(publish); Can not open: $trace \n"); }
   else
      {
      # Append publish  info to trace file
      print TRACE1 " Publishing: $ptime $jobid $user ------ $locnow \n";
      close(TRACE1);
      }

    # No buffering of  output
    my    $old_fh = select(PUBQUE);
    local $|      = 1;
    select($old_fh);
    if ($jobid , 0) 
       {$jobid=-$jobid;}

    # Shared lock= 1, Exclusive lock=2, Non-blocking request=4,  Free lock=8
    unless (flock(PUBQUE, 2 | 4))  #al dec07
       {
       print "Waiting for lock on: $quefile ";
       flock(PUBQUE, 2)  or die "Can't lock: $quefile";
       }
    seek(PUBQUE,0,2) || die "Can't seek to end of: $quefile \n"; #al oct07

    # Put jobid & username in Publisher que file 
    print PUBQUE "$jobid $user $locnowt\n";

    # Close Publisher que file
    flock(PUBQUE, 8);
    close(PUBQUE);

    system(" echo ' 'Publishing jobid: $jobid for: $user\n");

    return 1;
   }


# ---------------------- del_locked_job --------------------------------------

sub del_locked_job 
    {
    if (@_ > 3)
       {die  "Usage: $0 jobid
       Where:
       jobid is:    job PID \n"
       }
    local $linet,$out,$jobidt;

    $jobid = $_[0];
    #print "Deleting locked jobid: $jobid \n";

    # Open QUE file
    if ( ! open(PUBQUE, "+< $que")) 
          { system ("echo del_locked_job can not open: $que \n") ; return 0; }

    #No buffering of  output
    my    $old_fh = select(PUBQUE);
    local $|      = 1;
    select($old_fh);

    # Shared lock= 1, Exclusive lock=2, Non-blocking request=4,  Free lock=8
    unless (flock(PUBQUE, 2))
       { die " Can not write_lock: $que \n";}
    seek(PUBQUE,0,0) || die "Can't seek to start of: $que \n"; #al oct07

    # Read the whole queue and rewrite the queue without the deleted job
    $out = '';
    while (<PUBQUE>)       # Read que input
       {
       $linet = $_;
       ($jobidt) = split(/\b\s+/,$linet,1);  # First field in line
       #print  "jobidt NOW: $jobidt";   

        if (($jobidt !=  0) && ($jobidt != -$jobid))
          { # Not for current job, keep line in que file
          $out .= $linet;
          #print  "jobidt: $jobidt jobid: $jobid : $out :";   
          #print  "out NOW: $out";   
          }
       }

   #print  "QUE NOW: $out";   

   seek(PUBQUE,0,0);      # Goto que beginning again
   print PUBQUE $out;     # Put all output lines in que file
   truncate(PUBQUE, tell(PUBQUE));

   # Unlock the que file 
   flock(PUBQUE,8);
   close(PUBQUE);

   }        


sub lock_file  # ------------------- lock_file --------------------
   {
   my($file, $type) = @_;
   my($flag, $struct);

   $flag = ($type eq "r") ? F_RDLCK : F_WRLCK;
   $struct = pack("ssx32", $flag, 0);
   return(fcntl($file, F_SETLKW, $struct));
   }

sub unlock_file # -------------------- unlock_file ----------------
   {
   my($file) = @_;
   my($struct, $status);

   $struct = pack("ssx32", F_UNLCK, 0);
   $status = fcntl($file,F_SETLKW, $struct);
   return($status);
   }


# </pre></body></html>

