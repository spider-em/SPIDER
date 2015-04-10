#!/usr/local/bin/perl
#
# SOURCE:  ren5.perl
# PURPOSE: Renames file series
#
#USAGE:
# \ls *.html >! junk.lis
# perl  /home/leith/pgm/perl/ren5.perl < junk.lis >! junk

#$\ = "\n";              # set output record separator

line: while (<>) 
    {
    chop;       # strip record separator

    $string = $_;

#   Remove _ 
    $string =~ s/_//;
    $string =  lc $string;
 
#    print $string; 

     print "mv $_  $string \n";
    }
