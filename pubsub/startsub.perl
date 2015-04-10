#!/usr/bin/perl
## <html><head><title>Start a subsciber running</title></head><body><pre>
#
# File: startsub.perl     ArDean Leith Sep 02
#
# Source location: /net/bali/usr1/spider/pubsub
#
# Purpose: Start subscriber process (if not running) on this machine
#
# Usage example: startsub.perl 
# Test  example: startsub.perl 'subscribe.perl pubsub.que' 
#
# INPUT:
# $command      (subscriber startup string)

# OUTPUT: 
# success or failure     (returned)  

#print "arg0: ".$ARGV[0] . "\n";

if (@ARGV[0] eq "-h"  || @ARGV[0] eq "?")
   {    #help requested
   @op = split(m;/;,$0);     # drop directory list
   $op =$op[@op-1];
   printf(" Purpose: Start subscriber process. \n");
   printf("  \n");
   printf(" Usage: $op [que]  \n");
   printf("        optional argument is: que. \n");
   printf("  \n");
   printf(" Usage example: $op pubsub.que \n"); 
   printf("  \n");
   exit;
   }

# que is in first argument
$que = $ARGV[0];
if (!$que)
   {$que = "pubsub.que";}
#print "que: ". "$que" . "\n";

$search = "subscribe.perl";      # name of subscriber daemon 

$master = $ENV{'PUBSUB_MASTER'}; # name of node where subscriber is run
#print "master: ". "$master" . "\n";

($host) = split(/\./,$ENV{'HOST'});
if ($host != $master)
  { die " Subscribe must only be run on $master\n";}

# find out if "subscribe" is already running here using howmany
#print "search: " . $search . "\n";
$already = &howmany($search);

if ($already < 1)
   {  # subscribe not running, so start it up now
   $command = "$search $que";

   # get pubsub directory name
   $dir_pubsub = $ENV{'PUBSUB_DIR'};
   if ($dir_pubsub) 
      {$command = "$dir_pubsub/$command";}

   system("$command &");
   print "Started: $command \n";
   }
else
   {
   print "$search already running on: $master \n";
   }

exit 1;

# ----------------------- howmany --------------------------------- 

# Purpose: find out how specified processes are running
#
# Usage example:  howmany ("spider.*\@pub");
#
# INPUT:
# Process identifier for 'ps'           (argument #1)
#
# OUTPUT: 
# Number of processes running           (returned)  
#

sub howmany 
   {
   if (@_ < 1)
      {die  "Usage: $0 ps_string 
       ps_string is string to be searched for in ps. \n"
      }
   #print  "hunting for: " . $_[0] . "\n"; 

   open(PS, "ps -ef |")   || die "cannot fork: $!";
   $psstring = <PS>;                  # skip ps header line

   $num = 0;
   while ($psstring = <PS>)
      {
      #print  $psstring . "\n"; 
      if ($psstring =~ m#$_[0]#)           #count matching processes
         {
         $num++;
         # print  $num .  ": " . $psstring ; 
         }
      }

      #print  "# of $_[0] processes: " . $num . "\n"; 

   return $num;
   }
1;    #needed for require


# </body></pre></html>
