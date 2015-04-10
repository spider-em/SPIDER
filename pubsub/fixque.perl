#!/usr/bin/perl
#
# <html><head><title>Publish a job</title></head><body><pre>
#
# Que location : fixque.perl                     Nov. 2003 ArDean Leith    
#                              Can do +jobid now Mar. 2005 ArDean Leith
# LOCATI)ON: /net/bali/usr1/spider/pubsub/ 
#
# PURPOSE: Remove zombie job from publisher que. 
#
# USAGE:  fixque.perl  [pubsub.que] jobid 
#
# INPUT:
#        que              (optional argument)
#        jobid number     (argument)         

#Que location & name
$que = "pubsub.que";
if ($dir_pubsub) 
   {$que = "$dir_pubsub/pubsub.que";}
if ($ARGV[0] && $ARGV[1])
   { #Que location & jobid
   $que   = $ARGV[0];
   $jobid = $ARGV[1];   
   }
else
   {
   # Get jobid
   $jobid = $ARGV[0];
   }
   
#print  "jobid: $jobid \n";

if ( ! open(PUBQUE, "+< $que")) 
    { system ("echo Can not open: $que \n") ; exit 0; }

#Exclusive lock=2, Non-blocking request=4,  Free lock=8
unless (flock(PUBQUE, 2))
   { die " Can not write_lock: $que \n";}

# Read the  queue and rewrite the queue after deleting job
$out   = '';
$line  = '';
while (<PUBQUE>)       #read input
   {
   $linet = $_;
   #print  "linet: $linet \n";

   chop($linet);
   ($jobidt,$user,$machval,$rest) = split(/\b\s+/,$linet,4);
   #print  "jobidt: $jobidt \n";

   if ($jobidt == $jobid)
      { # This job is to be deleted from que
      $line = $linet;
      print "Removing: " . $line  . "\n\n";
      }
   else
      {
      $out .= $_;
      }
   }

if( $line)
   {   # Removed something
   seek(PUBQUE,0,0);    # Goto que beginning again
   print PUBQUE $out;   # Put all output lines in que file
   truncate(PUBQUE, tell(PUBQUE));
   }

# Unlock the que file 
flock(PUBQUE,8);
close(PUBQUE);

exit;

