#!/usr/bin/perl
#
# <html><head><title>List Pubsub que </title></head><body><pre>
#
# New                           Nov. 2007 ArDean Leith    
#                             
# SOURCE: /net/bali/usr1/spider/pubsub/ 
#
# PURPOSE: List contents of publisher que. 
#
# USAGE:  listque.perl  [pubsub.que]  
#
# INPUT:
#        que              (optional argument)

#Que location & name
$que = "pubsub.que";
if ($dir_pubsub) 
   {$que = "$dir_pubsub/pubsub.que";}
if ($ARGV[0])
   { #Que location 
   $que   = $ARGV[0];
   }
   
if ( ! open(PUBQUE, "< $que")) 
    { die "Can not open: $que \n"; }

#Exclusive lock=2, Non-blocking request=4,  Free lock=8
unless (flock(PUBQUE, 2))
   { die " Can not write_lock: $que \n";}

# Read the  queue and rewrite the queue after deleting job
$line  = 0;
 while (<PUBQUE>)       # Read input
   {
   $line++;
   print  " $_";
   }

# Unlock the que file 
flock(PUBQUE,8);
close(PUBQUE);
if ($line < 1) { print " Q empty \n";}

exit;

