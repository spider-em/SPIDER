#!/usr/local/bin/perl
#
# File: distribute.perl     ArDean Leith         Mar. 2001
#
# Source location: /net/bali/usr1/spider/bin/
#
# Purpose: distribute SPIDER jobs over specified computers
#
# Master task. Starts on one node only. 
# Coordinates and synchronizes all tasks. 
# Usage example:  distribute spi/acn @apmq_master
#
# INPUT:
# Project code/data code               (argument #1)
# Name of SPIDER procedure distributed (argument #2) 
# Micrograph numbers                   (queried) 
# List of computers to be used         (queried) 
#
#OUTPUT: 
#leave_while_running.jnk  

$tmpfile = "leave_while_running.jnk";

if ($#ARGV < 1) 
    {print ("\n distribute needs two command line arguments!\n\n");}

if (@ARGV[0] eq "-h"  || @ARGV[0] eq "?" || @ARGV[0] eq "")
   {    #help requested
   @op = split(m;/;,$0);     # drop directory list
   $op =$op[@op-1];
   printf(" Purpose: Run SPIDER procedure on specified machines. \n");
   printf("  \n");
   printf(" Usage: $op PRJ/DAT  procedure  \n");
   printf("        PRJ/DAT   is: SPIDER PROJECT/DATACODE. \n");
   printf("        Procedure is: SPIDER procedure to be distributed is in second argument \n");    printf("        (will only work if computers accept: rsh) \n");
   printf("  \n");
   printf(" Usage example: $op  PRJ/DAT @apmq_master \n"); 
   printf("  \n");
   exit;
   }

#extensions are in first argument
$EXT = $ARGV[0];

#SPIDER procedure to be distributed is in second argument
$RUN = $ARGV[1];


print "\n";
print "Faster computers: sicily,bali,crete,baffin,cuba,barra,ganga,nevis,ceram";
print "\n";
print "Older compters:  saba,lucia,penang,borneo,kodiak";
print "\n";
print "\n";

print "Enter computers (separate with commas): \n";
$machines = <STDIN>;                         # get input
@machinesa = eval qq(\($machines\));         #convert string to list

#sylt,bali,sicily,crete,flores,barra

print "\n";
print "Enter micrograph numbers (e.g. 1..3,5..10): "; 
$nums = <STDIN>;
print "\n";

@fields = split(/\,| /,$nums);

foreach $val (@fields)
   {
   @list= eval ($val);
   # print  "list: " . @list . "\n";
   push(@numa,@list);     # can read in 2...5,6 etc.
   }

#open the temp file that is being created
open (OUT, "> $tmpfile") ||
    die "Cannot create:  $tmpfile in this directory." ;
chmod (0770,$tmpfile);                      #make temp file executable

$locnow = $ENV{'PWD'};
$locnow =~ s#/net/\w*/##,$ENV{'PWD'};
#print " locnow: " . $locnow ." \n";  

print OUT "# remotely runs processes for SPIDER -- AP MQ.\n"; 

#create task for each micrograph
$i = 0; 
foreach $num (@numa)
   { 
    print OUT "rsh " . $machinesa[$i] . " -n " ;  
    print OUT " cd /net/" . $ENV{'HOST'} . '/'. $locnow." ";  
    print OUT "\';\' ";                     # 
    print OUT "spider $EXT ";
    print OUT $RUN . " " . $num . " X77=" . $num . " \& \n";

    $i++;          #increment machine selector
    $i = $i % @machinesa; 
    }

#execute the command file
close(OUT);
system  "./$tmpfile"; 

exit;

