#!/usr/bin/perl
#
# SOURCE: /usr8/spider/pubsub/delete.perl    
#         new                     ArDean Leith  Jun 2001
#         Msg. formatting         ArDean Leith  Aug 2010
#
# PURPOSE: Put final statistics for job in PubSub log file
#
# USAGE:  delete.perl 13681 pubsub.log
#
# INPUT:
#    job id    (argument #1)
#    log file  (argument #2)
#
# OUTPUT: 
#    job id          (returned)  

use Fcntl;

# Get the jobid from arguments 
$deljobnum      = @ARGV[0];

# Get the log file name from arguments 
$log       = @ARGV[1];
#print  "log: $log \n"; 

# Terminate this jobid in the publisher log

# Open LOG file
open(LOG, "+< $log") ||
   die "Delete.perl can not open: $log \n";;
    
# Lock the opened log file

#%%%% KLUDGE for ALBANY old kernal on usr11 host not NFS lockable
$nolock =  ($log =~ m/\/usr11/ );
if (! $nolock )
   { #print " LOCKING\n"; 
   $status = lock_file(LOG,"w");
   }
#%%%% End KLUDGE for ALBANY old kernal on usr11 host

    
# Read the log and rewrite the log for terminating job
$out = '';
while (<LOG>)
   {
   $line = $_;
   #print  "line: $line \n";

   chop($line);
   ($jobnumber,$remaining) = split(/\b\s+/,$line,1);

   $gotruntime =  ($line =~ m/.*Runtime.*/ ); # Do not alter this line if already done

   if ( ($jobnumber == $deljobnum) && (!$gotruntime) )
      {  # Found this job, list ending time in log 
      ($jobnumber,$machine,$quedtime,$starttime,$remaining) = split(/\b\s+/,$line,5);

      # Get current time
      $endtime = time();
      #($sec,$min,$hour,$mday,$month,$year,$more) = localtime($endtime);
      #$time = sprintf("%04d-%02d-%02d  %02d:%02d:%02d", 
      #                $year+1900,$month+1,$mday, $hour,$min,$sec);

      # Find elapsed run time
      $runtime = $endtime - $starttime;
 
      #$out .= "$deljobnum $machine $quedtime $starttime $remaining (Ended: $time) (Runtime: $runtime) \n";
      $out .= "$deljobnum $machine $quedtime $starttime $remaining (Runtime: $runtime) \n";

      while (<LOG>)       # Read remining lines from QUE file 
         { $out .= $_; }  # Keep line in LOG file
      last;               # Finished reading now
      }
   else
      { $out .= $_; }   # Different job or not a locked job
   }
    
seek(LOG,0,0);    # Goto log beginning again
print LOG $out;   # Put all output lines in log file
truncate(LOG, tell(LOG));

# Unlock the log file 
if (! $nolock )
   { #print " UNLOCKING\n"; 
   $status = unlock_file(LOG);
   }
close(LOG);

print  " Finished job: $deljobnum \n"; 
$ret = 0;
exit $ret;


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
