#!/usr/bin/perl
#
# SOURCE: /usr16/software/spider/utils/man2exa.perl
#
# PURPOSE:  Converts ascii manuals of SPIDER in $oldman_dir to exa shell
#
# CHANGES:              Author:                ArDean Leith Apr 2013


$spider_root = qq($ENV{'SPIDER_ROOT'});
$spider_root = "/usr16/software/spider";

$hr          = qq(<p>\n);

$rawdocs_dir = qq($spider_root/rawdocs);
$oldman_dir  = qq($spider_root/man);           # Directory for ascii man pages 
$exa_dir     = qq($spider_root/rawdocs/exa);   # Directory for example man pages 
$htmlext     = qq(html);                       # File extension for html files 


# End of content from htdefs.ph


sub htmlescape  # ----------------- htmlescape ----------------------
   { # Escape html <>& characters where not part of tags
   #  Put this first !!
   $line =~s/([^'])&([^'])/$1&amp;$2/g;    # & (not '&') --> &amp;
   $line =~s/([^'])</$1&lt;/g;		   # < (not '<)  --> &lt;
   $line =~s/>([^'])/&gt;$1/g;             # > (not >')  --> &gt;
   #$line =~s/([^'])"([^'])/$1&quot;$2/g;   # " (not '"') --> &quot;
   #   $line =~s/"/&quot;/g;               # " --> &quot;

   $line =~s/'</</g;                       # Retain '< tag start
   $line =~s/>'/>/g;                       # Retain >' tag end
   $line =~s/'&'/&/g;                      # '&' --> &
   $line =~s/'"'/"/g;                      # '"' --> "
   }

sub lineprepare  # --------------- lineprepare ---------------------
   { # Removes hyphen at end, also calls htmlescape
   if ($line =~s/-\s*\n$//)
       {	# Collapse if hyphen is at the end of the line
       $line .= join(' ',split(' ',shift(@lines)));
       print  "$oldman_dir/$fn.$ext: possible hyphenation: $_ \n";
       }
   $line =~s/\n$//;
   &htmlescape;
   }

sub headerscan  # ---------------- headerscan -----------------------

   {
   # Scan the header. 
   # $oper $description(possibly multiword) $oper
   # (optional line) [optional whitespaces] $update
   # 

   @words = split(' ',$line);
   $oper  = shift(@words); 
   $last  = pop(@words);
   
   if ($oper ne $last) 
      {    # If not 1 word 
      $oper .=" "; $oper .=shift(@words); 
      $llast =pop(@words); $llast .=" "; $llast .=$last;
      if ($oper ne $llast)
         {   # If not 2 words
         $oper .=" "; $oper .=shift(@words); 
         $lllast =pop(@words); $lllast .=" "; $lllast .=$llast;
	 
         if ($oper ne $lllast) 
            {   # If not `3 words
            $oper .=" "; $oper .=shift(@words); 
            $lllast =pop(@words); $lllast .=" "; $lllast .=$llast;
            if ($oper ne $lllast) 
            {print " $fn.$ext:Header Name! \n";}  
            }
         }
      }
   unless ($oper=~/[^a-z]/)
      {print "$oldman_dir.$ext/$fn: strange header\n";};
      
   $defn=join(' ',@words);
   unless ($defn =~/[a-z]/) 
      {print  "$oldman_dir/$fn.$ext: Definition: $defn -- is in caps only\n";}
      
   print TOO qq(<html>\n);
   print TOO qq(<head>\n);
   print TOO qq(<title> $defn </title>\n);
   print TOO qq(   <link rel='stylesheet' href='ex.css' type='text/css' />\n);
   print TOO qq(</head>\n\n);
   
   print TOO qq(<body>\n);
   print TOO qq(<br /> <br />\n\n);

   print TOO qq(<h2>Usage Example - Operation: $oper </h2>\n);
   print TOO qq(<h3>$defn</h3>\n);
   print TOO qq(<h4> Procedure used: <a href="./images/$fn.spi">./images/$fn.spi</a></h4> \n);
   
   print TOO qq(<hr />\n\n);
   
   print TOO qq(<div class="ex">\n);

   while(@lines)
      {
      $line=shift(@lines);
      &lineprepare;
      if ($line=~/\w/
         {
         #if ($line=~/^PURPOSE|^USAGE/) 
         if ($line=~/^USAGE/) 
            {last;}    # Altered for awkward headers
	    
         $outline = join (' ',split(' ',$line));
	 
         #print TOO "<em>($outline)</em><br />\n";
         }
      }
   }


sub usagescan  # ---------------- usagescan -----------------------
   {
   
   
   return unless ($line=~/^USAGE/);
   
   @words = split(' ',$line); shift(@words);
   $line  = join (' ',@words);
   
   print TOO qq(<table class="opl">\n);
   print TOO qq(  <tr> <div class="oplin"> <td class="op">.OPERATION:           </td> <td class="opres"> $line     </td> <td class="com">; &nbsp; $defn </td> </div> </tr> \n); 
   
   while(@lines)
      {
      $line = shift(@lines);
      &lineprepare;
      
      if ($line=~/\w/)
         {	    
         $outline = join (' ',split(' ',$line)); 
	 
         if ($outline=~/^\./) 
            {
            # Line starts with dot e.g. .INPUT FILE: 

	    @Fld    = split(':',$outline);
	    $prompt = $Fld[0];
	    $resp   = $Fld[1];
	    
	    #print qq(Prompt: $prompt:    Response: $resp \n);

	    print TOO qq(  <tr> <div class="in"> <td class="pr"> $prompt: </td> <td class="res"> $resp </td> <td class="com">; &nbsp; xxx </td> </div> </tr>\n);
            }
          }
      }	    
      print TOO qq(</table> \n\n);
   }



# The actual program begins here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



$manfile = $ARGV[0];

$fn = $manfile;

$fn =~s/\.man//;

#print "fn: $fn \n";
 
open(fro, " $oldman_dir/$manfile ") || die "Can't open $oldman_dir/$manfile \n";

#print "opened: $oldman_dir/$manfile \n";

open(TOO, "> $exa_dir/$fn.$htmlext") || die "Can't write to $exa_dir/$fn.$htmlext \n";

print "opened: $exa_dir/$fn.$htmlext \n";


   @lines = <fro>;

   LINE_MAIN: 
   while(@lines)
      {
      $line=shift(@lines);
      &lineprepare;              # Handles - at end and '< for html
      if ($line=~/^Filename:/)
         {next LINE_MAIN;}
	 
      elsif ($line=~/\w/)
         {
         if ($line=~/^\S/) {last;}
         $outline=join (' ',split(' ',$line));
         print TOO qq($outline\n);
         print "$oldman_dir/$fn.$ext: strange beginning\n";
         } 
      }

   &headerscan; 
   
   &usagescan;
   
   # &purposescan;
   # &warningscan;
   # &seealsoscan;
   # &notesscan;
   # &subscan;
   # &callerscan; 

   print TOO qq(</div>\n\n);

   print TOO qq(<br /> <hr /> <br />\n\n);

   print TOO qq(<table class="verybigimgs">\n);  
   print TOO qq(  <tr> <th class="lab_in">FIRST INPUT IMAGE <br />  ()         </th> </tr>\n);
   print TOO qq(  <tr> <td><img  class="img_in"  src="./images/.jpg">          </td> </tr>\n);
   print TOO qq(  <tr> <td class="nam_in">                                     </td> </tr>\n);
   print TOO qq(</table>\n\n);

   print TOO qq(<br /> <br />\n\n);

   print TOO qq(<table class="verybigimgs">\n);  
   print TOO qq(  <tr> <th class="lab_in">IMAGE <br /> ( )                     </th> </tr>\n);
   print TOO qq(  <tr> <td><img  class="img_in" src="./images/ .jpg">          </td> </tr>\n);
   print TOO qq(  <tr> <td class="nam_in">                                     </td> </tr>\n);
   print TOO qq(</table>\n\n);

   print TOO qq(<br /> <hr /> <br />\n\n);

   print TOO qq(<table class="verybigimgs">\n);  
   print TOO qq(  <tr> <th class="lab_out">OUTPUT IMAGE                        </th> </tr>\n);
   print TOO qq(  <tr> <td><img  class="img_out" src="./images/$fn.jpg">       </td> </tr>\n);
   print TOO qq(  <tr> <td class="nam_out">$fn                                 </td> </tr>\n);
   print TOO qq(</table>\n\n);

   print TOO qq(</body>\n);
   print TOO qq(</html>\n);
   
   
   close(fro);
   close(TOO);






