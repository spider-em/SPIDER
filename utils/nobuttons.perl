#!/usr/bin/perl
#
# SOURCE: /usr8/spider/utils/raw2docs.perl
#       exa buttons and links                  ArDean Leith  Jan. 2013
#       Nobuttons flag supported               ArDean Leith  Mar. 2015
#
# PURPOSE: Puts buttons, header, and trailer on html documents. 
#          If -wadsworth copies files to external Web site
#          Most HTML files (except if links) under the directory docs have 
#          their replacements done.
# 
# NOTE:    The $rel_path to buttons file is suitably modified enroute.

# Previously in: htdefs.ph;
$spider_root   = qq($ENV{'SPIDER_ROOT'});
$spider_root   = "/usr8/spider";

$rawdocs_dir   = qq($spider_root/rawdocs);
$manhtml_dir   = qq($spider_root/docs/man);
$docs_dir      = qq($spider_root/docs);
$butt_path     = "buttons";
$man_path      = "man";                # This is HTML man - not ASCII man!
$exa_path      = "exa";                # Examples dir

#$wads_docs_dir = qq(nnewton:/usr3/WWW/wwwinternal/spider_doc/spider/docs);
$wads_docs_dir = 'spider-stage:/export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider/docs/';    

# For Wadsworth copy command compress, update, preserve Exec., preserve time,recursive,verbose
$sendit = qq(rsync -zuEtrv);

# --------------------------- Main routine ---------------------------

#print "sendit: $sendit \n";

$fnum      = 0;
$wadsworth = 0; 

$rel_path  = ".";
while (<$rawdocs_dir/*.html>)
   {
   s/$rawdocs_dir\///e;
   $fn = $_;

   #print " Opening fn: $rawdocs_dir/$fn \n";
   
   open(fro, "< $rawdocs_dir/");

   # Search for 'nobuttons' flag in source html comment
   system("grep '<!-- NOBUTTONS -->' $rawdocs_dir/$fn > /dev/null");
   #system("grep 'garbage' $rawdocs_dir/$fn");
   
   $wantbuttons = 0;
   if ($? == -1) 
      { print " Failed to execute: $! \n"; } 
   elsif ($? & 127) 
      { printf " Grep died with signal: %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without'; } 
   else 
      { 
      #printf " Child exited with value: %d \n", $? >> 8; 
      $wantbuttons = $? >> 8;
      }
      
   #print " Grep of: $rawdocs_dir/$fn Returned: $?  Wantbuttons: $wantbuttons \n";
   
   if ($wantbuttons)
      { print " Want to add buttons to: $fn \n";}
   else 
      { print " Do not want to add buttons to: $fn\n";}

   
    
   # Record input file's atime & mtime
   ($atime, $mtime) = (stat("$rawdocs_dir/$fn"))[8,9];

   print " Opening fn: $docs_dir/$fn \n";
exit;
#   open(too, "> $docs_dir/$fn")|| 
#      die "Can't open: $docs_dir/$fn \n";
#
#   if ($fn=~/^[A-Z]/) 
#      {  # Filename starts with uppercase letter
#      while ($linn=<fro>)  {print too $linn;}
#      }
#   else 
#      {&installit;}
   close(fro);
#   close(too);
   #print " Created: $docs_dir/$fn \n";

   # Alter the output files atime & mtime so same as input
   utime($atime, $mtime, "$docs_dir/$fn");
   $fnum++;
   }   

print " Created: $fnum  HTML SPIDER doc files in: $docs_dir \n";
print " ------------------------  \n";

