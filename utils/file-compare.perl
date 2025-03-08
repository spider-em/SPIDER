#!/usr/bin/env perl
#
# SOURCE: /home/dean/Software/spider-2020-2025/utils/file-compare.perl
#
# PURPOSE:  Compare all files in a directory 
#
# CHANGES:     Author:    ArDean Leith Feb 2025
# USAGE:       /home/dean/Software/spider-github-2025/utils/file-compare.perl 

$edit_dir1 = qq(/home/dean/Software/spider-github-2025/src/);# Files to compared 
$edit_dir2 = qq(/home/dean/Software/spider-2020-2025/src/);  # Files to be edited 
$new_dir   = qq($edit_dir/Attic-jnk);                        # Dir for edited files
       	
# Actual program begins here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use strict; use warnings;

use Text::Diff;

#foreach $file (@files)
#  {
#  $oldfile = "$edit_dir/$file";
#  $newfile = "$new_dir/$file";
 
#  print " File: $file \n";
  # Check for .spi files, do not edit them
#  if ($file !~ /$\.html/ )    
#     {  
#3     #print  " Not html: $file\n";
#3     next;
#3     }
#3  }

#my $diffs = diff 'file1' => 'file2';

print $diffs;

$file = "redlin.f";

$file1 = $edit_dir1/$file; 
$file2 = $edit_dir2/$file;

my $diffs = diff $file1 => $file2;

print  " File: $file \n";
print $diffs;

