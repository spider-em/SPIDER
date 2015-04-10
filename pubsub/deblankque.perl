#!/usr/bin/perl
#
# <html><head><title>Remove blank lines from que</title></head><body><pre>
#
# Que location : deblankque.perl     ArDean Leith    Nov. 2003
#
# LOCATI)ON: /net/bali/usr1/spider/pubsub/ 
#
# PURPOSE: Remove blank lines from publisher que. 
#
# USAGE:  deblankque.perl  
#
# INPUT:
#        que              (optional argument)

#Que location & name
$que = "pubsub.que";
if ($dir_pubsub) 
   {$que = "$dir_pubsub/pubsub.que";}
if ($ARGV[0] && $ARGV[1])
   { #Que location
   $que   = $ARGV[0];
   }
   
if ( ! open(PUBQUE, "+< $que")) 
    { system ("echo fixque can not open: $que \n") ; exit 0; }

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

   if ($jobidt =~ m/\d+/)
      { # This is a job line
      $out .= $_;
      #print "Found jobid: $jobidt \n"; 
      }
   else
      {
      $line = $linet;
      }
   }

#print "In publish; jobidt: " . $jobidt . "machine:" . $rest . "\n\n";
print  $out;   # Put all output lines on terminal

if( $line)
   {
   seek(PUBQUE,0,0);    # Goto que beginning again
   print PUBQUE $out;   # Put all output lines in que file
   truncate(PUBQUE, tell(PUBQUE));
   }

# Unlock the que file 
flock(PUBQUE,8);
close(PUBQUE);

exit;

