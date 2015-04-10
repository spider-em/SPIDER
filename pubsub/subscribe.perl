#!/usr/bin/perl
#
# <html><head><title>Subscribe to a job from publisher que</title></head><body><pre>
# SOURCE: /usr8/spider/pubsub/subscribe.perl 
#
# HISTORY:  New                                     Sep. 2002 ArDean Leith
#           Added local permit files                Mar. 2005 ArDean Leith
#           Local permit fallthru                   Apr. 2005 ArDean Leith
#           Local & master permit checked           Sep. 2005 ArDean Leith
#           Local permit subsidiary                 Jun. 2007 ArDean Leith
#           Partition                               Jun. 2007 ArDean Leith
#           Seek, removed demon, removed openque    Jan. 2008 ArDean Leith
#           In que locked msg., removed /n          Jan. 2008 ArDean Leith
#           Timeout added to rsh                    Jul. 2010 ArDean Leith
#           localpermit job skip                    Oct. 2010 ArDean Leith
#
# PURPOSE:       Subscribe to job from pubsub que. Left running like a deamon
#
# USAGE EXAMPLE: subscribe.perl pubsub.que 
#
# INPUT:         queue      (argument #1) OPTIONAL (Has hardwired default)
#
# OUTPUT:        command    (returned)  
#
#
# ------------------ subscribe command --------------------------------

#use POSIX qw(uname);
use File::Basename;

# String containing 'ps' identifier for SPIDER job

if (! $ARGV[0]) 
   {die "Usage: $0 [publish_que]  
         Purpose: Subscribe to job from publisher que
         Argument:
           publish_que is: file containing publisher que to be watched.\n"; 
   }

local($locnow);

$master    = $ENV{'PUBSUB_MASTER'};    # Site specific name of master node
# Site specific name cluster partition (can usually be left undefined)
$partition = $ENV{'PUBSUB_PARTITION'}; # Site specific name cluster partition 

system("echo  "); 

$dir_pubsub = $ENV{'PUBSUB_DIR'};

$que = "$dir_pubsub/$ARGV[0]$partition";
system("echo Subscribe que: $que  \n"); 

$permitfile = "$dir_pubsub/pubsub.permit$partition";
system ("echo Permitfile: $permitfile \n");

# Converting this process to a 'daemon' not implemented as it fails 

$checkwait    = 2;     # Wait time before checking QUE again
$afterjob     = 0;     # Do not skip any job

$time_to_die  = 0;     # RUNS FOREVER!!

until ($time_to_die)
   {   
   # Find a job
   
   #system(" echo Calling findajob, afterjob:  $afterjob \n");
   #($jobid,$locnow) = &findajob("$que",$afterjob);
   ($jobid,$locnow) = &findajob("$que");
   #system(" echo After findajob using afterjob: $afterjob, jobid: $jobid \n");
 
   $locnow =~ s/\n$//;
   $locnow =~ s/ //;
   #system(" echo subscribe found jobid:  $jobid locnow: $locnow \n");
   $started = 0;
   
   if ($jobid > 0) 
      {  # Got a job from the QUE, find a subscriber for this job

      #system(" echo After findajob using afterjob: $afterjob, jobid: $jobid \n");

      # Find a node for the job to run on
      $found_a_node = 0;

      until ($found_a_node)
         {
         #system(" echo Hunting node for jobid: $jobid \n");

         $base  = basename($permitfile);
         $localpermitfile = "$locnow/$base";

         #system(" echo localpermitfile is: $localpermitfile \n");
         if ( -r $localpermitfile ) 
            {
            #system("echo Using localpermitfile: $localpermitfile");
            open(PERMIT, "< $localpermitfile") ||
                die "Can not open local permissions file: $localpermitfile \n";
            $localpermit = 1;
            }
         else
            {
            #system(" echo Using: $permitfile");
            open(PERMIT, "< $permitfile") ||
                die "Can not open permissions file: $permitfile \n";
            $localpermit = 0;
            }

         $out     = ' ';
         $n       = 0;
         $foundid = 0;

         while (<PERMIT>)   # Read permit file, skip comment lines (#)
            {
            # Read input line from permit file, skip comment lines (anything containing #) --------
            #system("echo  raw permit input: $_ ");
            if (/\s*#/)
               { #system("echo  Skipping: $_ ");   
               next;  
               }

            #system("echo  permit string:       $_ ");
            ($machine,$limit,$days,$time1,$time2,$check,$comment) = /(\w+)\s+(\d+)\s+(\d+)\s+(\d+:\d\d)\s+(\d+:\d\d)\s+(\d+)\s+(.*)/;    # split at spaces
            #system(" echo limit,days,time1,time2,check: $limit,$days,$time1,$time2,$checkok \n");
            #system(" echo limit: $limit  \n");

            if ( $localpermit )
               { # Check to see if machine OK by master permit file --------------------------

               open(MASTERPERMIT, "< $permitfile") ||
                   die "Can not open master permissions file: $permitfile \n";
               @master_line = <MASTERPERMIT>;    # Read whole file
               close(MASTERPERMIT); 
 
               $okinmaster = 0;
               foreach(@master_line)
                  {
                  if (m/$machine/ && ! m/#+.*$machine/)    
                     {  # Matching machine in master permit     
                     # print "match and not comment : " . $_ .  "\n";
                     $okinmaster = 1;
                     last;
                     } 
                  }
               if (! $okinmaster)
                  {next;}    # Machine OK in local permit but not permitted in master .permit
               } # End of: if ( $localpermit )

            # Find how many SPIDER pubsub jobs are already running on this machine --------------
            $locked = `grep '$machine ' $que`;
            if ($locked)   
               {
               system("echo On $machine: a different job is already starting");
               next;
               }   
            #system("echo On $machine: No jobs starting");

            $SIG{ALRM} = sub {die "timeout" };
            $SIG{CHLD} = 'IGNORE';
	    
	    eval 
	       {
	       alarm(90);
	       
               # SEARCHES FOR A STRING SPECIFIC TO COMMAND USAGE %%%% (spider)
               $procs = `rsh -n $machine "ps -u root -N -o cmd | grep 'spider' | grep -v grep  | wc -l" 2>&1`;
               #$procs = `rsh -n $machine "ps -ef | grep spider | grep -v grep | wc -l" 2>&1`;  
               #$procs = `rsh -n $machine "ps -u root -N -o cmd | grep 'tcsh.*spider' | grep -v grep  | wc -l" 2>&1`;
               #$procs = `rsh -n $machine ps -ef | grep spider | wc -l`; 
	       alarm(0);
	       }; # End of: eval 
	       
             if ($@) 
	        {
	        if ($@ =~ /timeout/)
	           { $procs = "Connnection to: $machine timeout"; } # Timed out	
                else
	           { alarm(0); die; }      # Propagate unexpected other exception
	        }
	            	       
            #system("echo Number of Pubsub jobs on $machine: $procs");
            if ($procs =~ m/[A-Za-z]+/x)
                {  # Something wrong with rsh call
                system("echo FAILURE: $procs");
                next;
                }
            chop($procs);

            if ($procs >= $limit)
               {next;}    # This machine has too many jobs running

            # Check to see if time OK for running on this machine ---------------------------
            # Get current time
            ($sec,$min,$hour,$mday,$month,$year,$dayofweek,$line) = localtime(time());
            ($hour1,$min1) = split(/:/,$time1);
            ($hour2,$min2) = split(/:/,$time2);

            $mins0 = $hour  * 60 + $min ;
            $mins1 = $hour1 * 60 + $min1;
            $mins2 = $hour2 * 60 + $min2;

            #system("echo days: $days      day of week: $dayofweek  \n");
            #system("echo mins0: $mins0    mins1: $mins1     mins2: $mins2   \n");
            #system("echo hour: $hour      hour1: $hour1 hour2: $hour2 \n");

            if  (($mins0 < $mins1) || ($mins0 > $mins2))
                 { 
                 #system("echo ($mins0 < $mins1) || ($mins0 > $mins2) \n");
                 next;
                 }
 
            if (($days == 5) && (($dayofweek == 0) || ($dayofweek == 6)))
                 { 
                 #system("echo ($days == 5) && (($dayofweek != 0) || ($dayofweek != 6)) \n");
                 next;
                 }
 
            # OK time to start a job on this machine
            # system("echo On: $machine  time OK to start a job \n");

            # OK to start a job on this machine -------------------------------------------------
            system("echo On: $machine  Current jobs: $procs OK to start: $jobid \n");

            # Lock the job in the que by giving it negative number
            $gotid = &lockajob($que,$jobid,$machine);
            if ($gotid == $jobid)
               {  # Locked OK
               #system("echo Locked: $gotid \n");
               $found_a_node = 1;

               # Start the job
               #system("echo Can start job: $jobid  on: $machine \n");
               $started = 1;
               last;
               }        # End of  if     ($gotid == $jobid)
	       
            elsif ($gotid < 0)
               {
               system("echo Error -- lockajob could not find job: $jobid  in: $que!!");
               $found_a_node = 1; 
               last;
               }        # End of:  if     ($gotid < 0)
            }           # End of:  while  (PERMIT)


         close(PERMIT);
	    
         $found_a_node = 1;  # Force it to check another job
	 
	 if (!$started && 
	      $localpermit && 
	      $jobid != $afterjob)  # Want to skip jobs till after this one
	    {
	    $afterjob = $jobid ;
            #system("echo Set afterjob: $afterjob \n"); 
            last;            # Find another job
	    }

	 else   # debug info
	    {
	    $afterjob     = 0;
            #system("echo Finished permitfile, started: $gotid \n"); 
	    }

        }              # End of:  until ($found_a_node) --------------
	              
     }                 # End of:  if ($jobid > 0)                  

     if (!$started && ($afterjob == 0 )) 
        {select(undef,undef,undef,$checkwait);}  # Waits $check seconds and loops

  }                    # End of:  until ($time_to_die)

exit;



# ----------------------- findajob --------------------------------- 
# USAGE:  findajob ("pubsub.que") 
#
# INPUT:
#    Que name & location (argument #1)
#    Jobid of last job in loacal permit file (argument #2)
#
# OUTPUT: 
#    Found job id      (Returned)  


sub findajob 
   {  # Finds next job from que, returns job info
   if (@_ < 1)
      {die  "Usage: $0 que  \n
      que is:      file containing publisher que. \n "
      }
   #   afterjob is: number of last jobid used from que if not started. \n"
   local($rest,$jobid,$user,$n,$locnow);

   $que      = $_[0];

   #print  "Opening: " . $que . " afterjob: ".$afterjob . "\n"; 
   if ( ! open(PUBQUE, "+< $que"))
      { system ("echo Subscriber can not open: $que \n") ; return 0; }

   $foundid = 0;
   $rest    = ' ';
   $n       = 0;
   $usenext = (afterjob > 0); # Only use next job if > afterjob
   
   #   {$usenext = 0;}   # Only use next job if > afterjob

   # Shared lock= 1, Exclusive lock=2, Non-blocking request=4,  Free lock=8
   unless (flock(PUBQUE, 1))
      { die " Can not read_lock: $que \n";}
   seek(PUBQUE,0,0) || die "Can't seek to start of: $que \n"; #al jan08

   while (<PUBQUE>)       # Read input from que file
      {
      ($jobidt, $user, $locnow, $rest) = split(/\b\s+/,$_,3);

      $n++;
      #print "input[".$n."]:    " . $_."\n";
      #print "Found jobidt[".$n."]:   " . $jobidt . "  user: " . $user . "\n";
      #print "locnow[".$n."]:   " . $locnow . "  rest: " . $rest . "\n";

      if ($jobidt > 0)
         {   # Not locked job, need to start this job
         #print " Found qued jobidt: $jobidt  locnow: $locnow \n";

         if ($jobidt == $afterjob) 
            {
	    $usenext = 1;   # Can use next job found now
            #system("echo Can use job after: $afterjob at : $jobidt \n"); 
            next;
            }
	    	    
         if ($afterjob > 0 && ! $usenext)
	    {
            #system("echo Can not use this job: $jobidt \n"); 
            next; 
	    }   # Can not use this job
	 
	 $beforeid = 0;  # Can use next job found
	 
         $foundid = $jobidt;
         last;
         }    # End of: if ($jobidt > 0)
      }       # End of: while (<PUBQUE>)
      
   flock(PUBQUE,8);    # Unlock que file

   #system("echo Endof findajob: $foundid $afterjob $usenext \n");
       
   if ($usenext && !$foundid )
      {  # No job to run after: afterjob so reset afterjob
      $afterjob = 0; 
      #system("echo Set afterjob: $afterjob \n"); 
      }
      
   #system ("echo foundid,usenext:  $foundid,$usenext,$jobidt\n") ;
   
   return $foundid, $locnow;
   }


# ----------------------- lockajob --------------------------------- 
# USAGE example:  lockajob ("pubsub.que",1,node1) 
#
# INPUT:
#    queue            (argument #1)
#    job number       (argument #2)  
#    machine          (argument #3)  
#
# OUTPUT: 
#    locked number    (returned)  
#
 
sub lockajob 
   {  # Finds specific job in que, locks job, returns job id
   if (@_ < 1)
      {die  "Usage: $0 que jobid machine\n
      que is:      Publisher que.\n
      jobid is:    Job PID to start.\n
      machine is:  Machine to start this job on. \n"
      }
   local($que,$jobid,$machine,$user,$rest);

   $que      = $_[0];
   $jobid    = $_[1];
   $machine  = $_[2];
   #print  "lockajob  que: $que  jobid: $jobid  machine: $machine\n"; 

   $lockedid = 0;
   if ( ! open(PUBQUE, "+< $que"))
       { system ("echo lockajob can not open: $que \n") ; return $lockedid; }

   # No buffering of  output
   my    $old_fh  = select(PUBQUE);
   local $|       = 1;
   select($old_fh);

   $out      = '';
   $n        = 0;
   $lockedid = -1;

   # Shared lock= 1, Exclusive lock=2, Non-blocking request=4,  Free lock=8
   unless (flock(PUBQUE, 2))
      { die " Can not write_lock: $que \n";}
   seek(PUBQUE,0,0) || die "Can't seek to start of: $que \n"; #al jan08

   while (<PUBQUE>)       # Read input lines from que file
      {
      ($idnow,$user,$rest) = split(/\b\s+/,$_,3);
      $n++;
      #print "idnow[ $n ]: $idnow $user \n";

      if ($idnow == $jobid)
         {   # Want to lock this job now
         #print "idnow: " . $idnow . "rest:" . $rest . "\n\n";

         $out .= "-$jobid $user $machine $rest";  # alter lock & copy line

         $lockedid = $jobid;
         #system ("echo lockedid,rest:  $lockedid, $rest \n") ;

         while (<PUBQUE>)         # Read remaining input, copy to OUTPUT
            {$out .=  $_; }
         last;
         }
      else
         {   # Another job, just copy to OUTPUT
         $out .=  $_;
         # print "copy: " . $_ ;
         }
      }   # END of: while (<PUBQUE>)      

   if ($lockedid >= 0)
      {
      # Copy new publisher que to PUBQUE file 
      seek(PUBQUE,0,0) || die "Can't seek to start of $que";
      print PUBQUE $out;
      truncate(PUBQUE, tell(PUBQUE));
      }

   flock(PUBQUE, 8);
   close(PUBQUE);

   return $lockedid;
   }

# </body></pre></html>
