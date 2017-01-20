#!/usr/bin/perl
#!/usr/local/bin/perl
#
# perl script to create a document file from a particular ascii file.

$key = 0 ;
 while (<>) 
   {
    $key = $key + 1;
    $reg = "7";
    $append_str = $key + $reg;     
    #$new_line = $append_str + $_ ;
    #print $new_line, "\n";
   
    printf "%5d",$key;
    print ' ';
    print $reg;
    print $_ ;
    
   }
