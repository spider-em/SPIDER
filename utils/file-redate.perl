#!/usr/bin/env perl
#
# SOURCE: /home/dean/Software/spider-2020-2025/utils/file-redate.perl
#
# PURPOSE:  Redate all files in a directory 
#
# CHANGES:     Author:    ArDean Leith Feb 2025
# USAGE:       /home/dean/Software/spider-github-2025/utils/file-compare.perl 

#$bad-edit-dir = qq(/home/dean/Software/spider-github-2025/src/);  # Files to compared

 
#$gud-edit-dir = qq(/home/dean/Software/spider-2020-2025/src/);    # Files to be edited 

$bad_date_dir = qq(/spider-2020-2025/src/);    # ok  date Files
$gud_date_dir = qq(/spider-github-2025/src/);  # bad date Files

#------- active
$bad_date_dir = qq(diff-ed-git);      # Bad  date dir
$gud_date_dir = qq(diff-ed-spi);      # Good date dir

$gud_out_dir  = qq(diff-ed-git);      # Good output dir
       	
# Actual program begins here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   my $fnum = 0;
     
   opendir(GUD_OUT_DIR,$gud_out_dir) || 
      die "Can't open directory: $gud_out_dir \n";
   print "Opened dir: $gud_out_dir \n";
   
   @all_files = readdir(GUD_OUT_DIR);    # Read in file names

allf:
   foreach $file (@all_files)
      {
      # Check for non src files, do not edit them
      #print "\nFile: $file \n";

      if ($file !~ /$\.f/ )    
         { next ; }    #print "Not src: $file\n"; 
  
      $output = `ls -l $gud_out_dir/$file`;
      #print "Updating file: $gud_out_dir/$file \n";

      ($atime_b, $mtime_b) = (stat("$bad_date_dir/$file"))[8,9];
      #print "Bad time:  $atime_b  $mtime_b \n";

      # Record good input file's atime & mtime
      ($atime_g, $mtime_g) = (stat("$gud_date_dir/$file"))[8,9];
      #print "Good time: $atime_g  $mtime_g \n";
    
      # Set times on output file
      utime($atime_g, $mtime_g, "$gud_out_dir/$file");
       
      $output = `ls -l $gud_out_dir/$file`;
      print "Updated: $gud_out_dir/$file : $output";

      $fnum++;
      ###if ($fnum > 4) {exit;}
    
      }   # End of: foreach ...
 
#-rw-r--r-- 1 dean users 3330 Apr 27  2009 diff-ed-git/acrs_2.f
#-rw-r--r-- 1 dean users 3330 Apr 27  2009 diff-ed-spi/acrs_2.f
