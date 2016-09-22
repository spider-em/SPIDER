#!/usr/bin/perl
#
# SOURCE: /usr8/spider/utils/old2raw.perl
# CHANGES:              Author:                Sandeep Sibal Nov. 1994
#       Changed index to table format          ArDean Leith  Jan. 1997 
#       Adds "Eulerian angles" link            ArDean Leith  Feb. 2002
#       Cosmetic changes                       ArDean Leith  Feb. 2006
#       Do not menu obsolete 'AP' operations   ArDean Leith  Sep. 2006
#       MPI                                    ArDean Leith  Oct. 2006
#       Altered for Linux                      ArDean Leith  Mar. 2009
#       File modification times                ArDean Leith  Mar. 2009
#       usr8                                   ArDean Leith  Jul. 2010
#       addexamples                            ArDean Leith  Aug. 2011
#       Do not menu obsolete 'RT' operations   ArDean Leith  Oct. 2011
#       Alpahabet links                        ArDean Leith  Oct. 2011
#       Do not menu obsolete 'AS' operations   ArDean Leith  Mar. 2012
#       Do not menu obsolete 'PK C' operations ArDean Leith  Dec. 2012
#       Do not menu obsolete 'PK 3' operation  ArDean Leith  Jun. 2013
#       Do not menu obsolete 'HI *' operations ArDean Leith  Aug. 2013
#       Do not menu obsolete 'RM *' operations ArDean Leith  Aug. 2013
#       Do not menu obsolete 'PJ *' operations ArDean Leith  Feb. 2014
#       Removed obsolete 'RM *' operations     ArDean Leith  Sep. 2014
#       Moved raw man dir                      ArDean Leith  Mar. 2015
#
# PURPOSE:  Converts ascii text manuals of SPIDER in $text_mans_dir to HTML 
#           (only files with .man  extensions are processed), and 
#           places the HTML files in $raw_docs_dir.
#           Also creates an index of the operations.


# Locations and definitions previously in: htdefs.ph;
$LOCK_EX       = 1;
$LOCK_UN       = 8;
$spider_root   = qq($ENV{'SPIDER_ROOT'});
$spider_root   = "/usr8/spider";

$hr            = qq(<p>\n);
$mans_dir      = qq(man);                     # Directory for headerized html man pages wrt docs
$text_mans_dir = qq($spider_root/man);        # Directory for ascii text man pages 
$raw_mans_dir  = qq(jnk_raw_mans);            # Directory for pre-header html man pages wrt utils
$raw_docs_dir  = qq($spider_root/rawdocs);    # Directory for the pre-header html docs 
$src_dir       = qq($spider_root/src);        # Directory for tSPIDER source code
$exa_dir       = qq($spider_root/docs/exa);   # Directory for example html pages 
$htmlext       = qq(html);                    # File extension for html files 
$index_dir     = $raw_docs_dir;               # Directory for the pre-header html operation index 
$ferr          = "htdocs_warnings.LOG";       # Error File in utils directory


sub euler  # -------------------------- euler -------------------------
   {    # Add link for "Euler.. angle."  al Feb 2002
   $line =~s#(Euler\w* +angles*)#<a href="../euler.html">$1</a>#g;	
   }

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
       print ferr "$text_mans_dir/$fn.$ext: possible hyphenation: $_ \n";
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
      $llast =pop(@words);$llast .=" "; $llast .=$last;
      if ($oper ne $llast)
         {   # If not 2 words
         $oper .=" "; $oper .=shift(@words); 
         $lllast =pop(@words);$lllast .=" "; $lllast .=$llast;
         if ($oper ne $lllast) 
            {   # If not 3 words
            $oper .=" "; $oper .=shift(@words); 
            $lllast =pop(@words);$lllast .=" "; $lllast .=$llast;
            if ($oper ne $lllast) 
            {print " $fn.$ext:Header Name! \n";}  
            }
         }
      }
   unless ($oper=~/[^a-z]/)
      {print ferr "$text_mans_dir.$ext/$fn: strange header\n";};
      
   $defn = join(' ',@words);
   unless ($defn =~/[a-z]/) 
      {print ferr "$text_mans_dir/$fn.$ext: Definition: $defn -- is in caps only\n";}
    
   $scripting = $defn =~ /Scripting/;  # For scripting example page 
   #print "Scripting: $scripting \n";
   
     
   $definition{$oper} = $defn;				# This is an associative array!
   $htmlurl{$oper}    = "$mans_dir/$fn.$htmlext";	# Also an associative array!
   $localtype{$oper}  = "$ext";		                # Also an associative array! (al)
   
   print TOO qq(<head>\n);
   print TOO qq(<title>SPIDER: $oper \($defn\)</title>\n);
   print TOO qq(   <link rel='stylesheet' href='niceman.css' type='text/css' />\n);
   print TOO qq(</head>\n\n);
   print TOO qq(<body>\n);
   print TOO qq(<h2>$oper - $defn</h2>\n);
   while(@lines)
      {
      $line=shift(@lines);
      &lineprepare;
      if ($line=~/\w/)
         {
         if ($line=~/^PURPOSE|^USAGE/) 
            {last;}    # Altered for awkward headers
         $outline = join (' ',split(' ',$line));
         print TOO "<em>($outline)</em><br />\n";
         }
      }
   }


sub purposescan # ---------------- purposescan -----------------------
   {
   return unless ($line=~/^PURPOSE/);
   $first="yes";
   @words=split(' ',$line);shift(@words);
   $line=join (' ',@words);$line=" ".$line; unshift(@lines,$line);

   print TOO $hr;
   print TOO qq(<dl><dt><strong>PURPOSE</strong>\n);
   print TOO qq(<p>\n<dd>);
   while(@lines)
      {
      $line=shift(@lines);
      if ($first eq "yes") 
         {$first="no";}
      else {&lineprepare;}
         if ($line=~/\w/)
            {
            if ($line=~/^\S/)
               {
               &addexample;    # Add example link
               last;
	       }
            &euler;  # Add Eulerian angles link al Feb 2002
            $outline=join (' ',split(' ',$line));
            print TOO qq($outline\n);
            }
         }
      print TOO qq(</dl>\n);
      }

sub warningscan  # ---------------- warningscan ----------------- (UNUSED)
   {
   return unless ($line=~/^WARNING/);
   $first="yes";
   @words=split(' ',$line);
      shift(@words);$line=join (' ',@words);
      $line=" ".$line;unshift(@lines,$line);
   print TOO $hr;
   print TOO qq(<dl><dt><strong>WARNING</strong>\n);
   print TOO qq(<p>\n<dd>);
   while(@lines) 
      {
      $line=shift(@lines);
      if ($first eq "yes") {$first="no";}
      else{&lineprepare;}
      if ($line=~/\w/) 
         {
         if ($line=~/^\S/) {last;}
         $outline = join (' ',split(' ',$line));
         print TOO qq($outline\n);
         }
      else {print TOO qq(<p>\n);}
      }
   print TOO qq(</dl>\n);
   }
   
sub seealsoscan # ---------------- seealsoscan -----------------------
   {
   return unless ($line=~/^SEE ALSO|^EE ALSO/);
   $first="yes";
   @words=split(' ',$line);
   shift(@words);
   shift(@words);
   $line=join (' ',@words);
   $line=" ".$line;
   unshift(@lines,$line);
   print TOO $hr;
   print TOO qq(<dl><dt><strong>SEE ALSO</strong>\n);
   print TOO qq(<p>\n<dd>);
   print TOO qq(<table>);

   while(@lines)
      {
      $line=shift(@lines);
      if ($first eq "yes") 
         {$first="no";}
      else 
         {&lineprepare;}
      if ($line=~/\w/)
         {
         if ($line=~/^\S/) 
            {last;}
         $outline=join (' ',split(' ',$line));
         ($operr,$defnn)=split(/ \(/,$outline);
         $defnn=~s/\)//;
         $loperr=$operr; 	
         $loperr=~tr/[A-Z]/[a-z]/;
         $loperr=~s/\s//g;
         if (-T "$text_mans_dir/$fn.man")
            {    # man chapter 
            &getdescscan;
            }
         else
            {   # Chapter not found
            print TOO qq(<strong>$operr</strong> [$defnn]<br />\n);
            print ferr "$text_mans_dir/$fn.$ext: Reference $text_mans_dir/$fn.* not found!\n";
            }
         }
      }	
      print TOO qq(</table>);
      print TOO qq(</dl>\n);			
   }

sub getdescscan # ---------------- getdescscan -----------------------
   {
   # SCAN THE HEADER OF THE SEE ALSO FILE FOR DESCRIPTION. 
   # $oper $description(possibly multiword) $oper

   unless (open(FROE, "< $text_mans_dir/$loperr.man"))               
      {
      print ferr "From:  $fn.$ext   Can not find See also: $loperr.man\n";
      return;
      }
   @linesd=<FROE>;
   $n=0;
   while (@linesd[$n] =~ /^\s$/)
      {$n++;}
   @wordsd=split(' ',@linesd[$n]);

   $operd=shift(@wordsd); 
   $lastd=pop(@wordsd);
   if ($operd ne $lastd) 
      {    # If not 1 word in OP eat other words
      $operd  .=" "; $operd  .=shift(@wordsd); 
      $llastd =pop(@wordsd);
      $llastd .=" "; 
      $llastd .=$lastd;
      if ($operd ne $llastd)
         {   #if not 2 words
         $operd .=" ";  $operd .=shift(@wordsd); 
         $lllastd =pop(@wordsd);
         $lllastd .=" "; 
         $lllastd .=$llastd;
         if ($operd ne $lllastd) 
            {   #if not 3 words
            $operd .=" "; $operd .=shift(@wordsd); 
            $lllastd =pop(@wordsd);
            $lllastd .=" "; 
            $lllastd .=$llastd;
            if ($operd ne $lllastd) 
               {print "$fn.$ext   Desc Header Name: $operd! \n";}  
            }
         }
      }
   $defnd=join(' ',@wordsd);
   unless ($defnd =~/[a-z]/) 
      {print ferr "From /$fn.$ext: Descr: $defnd  is in caps only\n";}
   # print qq(defnd: $defnd \n);
   print TOO qq(<tr><td><a href="$loperr.$htmlext"><strong>$operr</strong></a></td><td> [$defnd]</td></tr>\n);
   close(FROE);
   }

sub includescan # ---------------- includescan -----------------------
   {
   # Include text from another file #INCLUDE filename.  Note that if filename
   # starts with 'cp' it has already been included using 'cpp'

   return unless ($line=~/^#INCLUDE/);

   @words = split(' ',$line) ;
   $incfile = pop(@words);
   #print " Including file: $incfile  \n";

   unless (open(FROI, "< $text_mans_dir/$incfile"))               
      {
      print      "From: $fn.$ext  Can not include: $text_mans_dir/$incfile\n";
      print ferr "From: $fn.$ext  Can not include: $text_mans_dir/$incfile\n";
      return;
      }

   while (<FROI>)       # Insert whole included file contents
      { print  TOO $_;   }	
   close(FROI);

   $line = "  ";        # Blank the line for next subroutine
   }





sub usagescan  # ---------------- usagescan -----------------------
   {
   return unless ($line=~/^USAGE/);
   @words=split(' ',$line);shift(@words);
   print TOO $hr;
   print TOO qq(<dl><dt><strong>USAGE</strong></dt>\n);
   $line=join (' ',@words);
   # if ($oper ne $line) {print ferr "$text_mans_dir/$fn.$ext: Strange .OPERATION\n";}
   if ($oper ne $line) {print ferr " HEADER CMD: $oper           USAGE: $line  \n";}
   print TOO qq(<p>\n<dd>.OPERATION: $line<br />\n);

   $previous="short";  # Short margin length
   $PRE = '';
   $GOT = '';
   while(@lines)
      {
      $line=shift(@lines);
      &lineprepare;
      
      &includescan;         # Put in any included files in here
      if ($line=~/\w/)
         {
         if ($line=~/^\S/) 
            {last;}
         &euler;    # Add Eulerian angles link al Feb 2002

         #  Do not alter preformatted text
         if ($line=~/<PRE>/)   { $PRE=1};
         if ($line=~/<\/PRE>/) { $PRE=''};
         if ($PRE)
            {  # This is preformatted text, do not squish
            print TOO qq($line \n);
            }
         else
            {
            $outline=join (' ',split(' ',$line)); 
            if ($outline=~/^\./) 
               {
               # Line starts with . e.g. .INPUT FILE: 
               if ($previous eq "short") 
	           {print TOO qq(<dd>);}
		   
               if (! $GOT) 
                  {print TOO qq(<br />);} 
               print TOO qq($outline<br />\n);
               $previous="long";
               $GOT='';
               }
            elsif ($line=~/^\s\s\s\s\s\s\s\s/)
               { # A long left margin (normal)
               if ($previous eq "short") 
                  {print TOO qq(</dd> <dd>);}
               print TOO qq($outline\n);
               $previous = "long";
               $GOT='';
               }
            else 
               { # A short left margin 
               if ($previous eq "long") 
                  {print TOO qq(</dd> <dt>);}
               print TOO qq($outline\n);
               $previous = "short";
               $GOT='';
               }
            }
         }
      else
        {print TOO qq(<p>\n); $GOT=1;}
      }
   print TOO qq(</dl>\n);		
   }



sub addexample # ---------------- addexample -----------------------
   {
   #print  qq(fn: $fn.\n);
   
   if ($fn=~/^iq/)
      {
      if (-T "$exa_dir/$fn.html")
         {   # Specific IQ example page needed 
         print TOO qq(&nbsp; <a href="../exa/$fn.$htmlext">Example</a>.\n);
         return;
         }
      else
         {   # General IQ Example page needed 
         #print  qq(fn: $fn.\n);
         print TOO qq(&nbsp; <a href="../exa/inquiry.$htmlext">Example</a>.\n);
         return;
         }
      }
      
   if ($fn=~/^sdic|^udic/) 
      {   # 'SD IC' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/sdic.$htmlext">Example</a>.\n);
      return;
      }
   if ($fn=~/^udfind/) 
      {   # 'UD FIND' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/udfind.$htmlext">Example</a>.\n);
      }
   if ($fn=~/^udnext/) 
      {   # 'UD NEXT' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/udnext.$htmlext">Example</a>.\n);
      return;
      }
   if ($fn=~/^hit|^hir/) 
      {   # 'HI' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/hi.$htmlext">Example</a>.\n);
      return;
      }
   if ($fn=~/^hidr/) 
      {   # 'HID' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/hid.$htmlext">Example</a>.\n);
      }
   if ($fn=~/^ec/) 
      {   # 'EC' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/ec.$htmlext">Example</a>.\n);
      return;
      }
   if ($fn=~/^ro$|^roi|^rosd/) 
      {   # 'RO' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/ro.$htmlext">Example</a>.\n);
      return;
      }
   if ($fn=~/^msi|^msf|^ms3/) 
      {   # 'MS' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/ms.$htmlext">Example</a>.\n);
      return;
      }
   if ($fn=~/^sa3|^sae|^sap/) 
      {   # 'SA' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/sa.$htmlext">Example</a>.\n);
      return;
      }
   if ($fn=~/^ins/) 
      {   # 'IN' Example page needed 
      print TOO qq(&nbsp; <a href="../exa/in.$htmlext">Example</a>.\n);
      return;
      }

   if ($scripting > 0) 
      {   # Scripting Example page needed 
      print TOO qq(&nbsp; <a href="../exa/scripting.$htmlext">Example</a>.\n);
      }
   elsif (-T "$exa_dir/$fn.html")
      {    # Example page present 
      print TOO qq(&nbsp; <a href="../exa/$fn.$htmlext">Example</a>.\n);
      }
   }




sub notesscan  # ---------------- notesscan -----------------------
        {   # Processing NOTES:
	
	return unless ($line=~/^NOTES|^NOTE|^Note/);
	$first="yes";
	@words=split(' ',$line);shift(@words);
           $line=join (' ',@words);$line=" ".$line;unshift(@lines,$line);
	print TOO $hr;
	print TOO qq(<strong>NOTES</strong>\n);
	print TOO qq(<ol>\n);
	if ($line=~/^(\s+)(\d+)/) 
           {$numbered_notes="yes";}
	else 
           {$numbered_notes="no"; print TOO qq(<p>\n<li>);}
	$prev="full";
	notes_line:
        $PRE = '';
	while(@lines)
            {
            $line=shift(@lines);
            if ($first eq "yes")
               {$first="no";}
            else 
               {&lineprepare;}

            &includescan;         # Put in any included files in notes

            if ($line=~/\w/)
                {
                if ($line=~/^\S/)  # Stop if line does not start with  whitespace  
                   {last;}        

                 #-------------------- Do not alter preformatted text
                if ($line=~/<PRE>/)
                   { $PRE=1};
                if ($line=~/<pre>/)
                   { $PRE=1};
                if ($line=~/<\/PRE>/)
                   { $PRE=''};
                if ($line=~/<\/pre>/)
                   { $PRE=''};
                if ((! $PRE) && $line=~s/^(\s+)(\d+)\.//) 
                   {print TOO qq(<p>\n<li>); $numbered_notes="yes";}
                elsif ((! $PRE) &&$line=~s/^(\s+)(\d+)\)//) 
                   {print TOO qq(<p>\n<li>); $numbered_notes="yes";}
                elsif ((! $PRE) &&($numbered_notes eq "no")&&($prev eq "blank"))
                   {print TOO qq(<p>\n<li>);}
                   &euler;  #add Eulerian angles link  al Feb 2002
                #-------------------- Do not alter preformatted text
                if ($PRE)
                   {  # This is preformatted text,do not squish
		   print TOO qq($line \n);
                   }
                else
                   {  # Not preformatted, OK to squish
		   $outline=join (' ',split(' ',$line));
		   print TOO qq($outline\n);
                   }
                $prev="full";
		}
	   else{$prev="blank";}
	}
	print TOO qq(</ol>\n);
}

sub subscan  # -------------------------------- subscan --------------------------------
   {

   return unless ($line=~/^SUBROUTINE|^\(SUBROUTINE|^Subroutine|\(Subroutine|^subroutine|^\(subroutine/);
   $first="yes";
   @words=split(' ',$line);shift(@words);
      $line=join (' ',@words);$line=" ".$line;unshift(@lines,$line);
   @words=();
   print TOO $hr;
   print TOO qq(<dl><dt><strong>SUBROUTINES</strong>\n);
   print TOO qq(<p>\n<dd>);
   while(@lines)
      {
      $line=shift(@lines);
      if ($first eq "yes") 
         {$first="no";}
      else
         {&lineprepare;}

      if ($line=~/\w/) 
         {
         if ($line=~/^\S/) 
            {last;}
         @phrases=split(' ',$line);
         while(@phrases) 
            {@wordss=split(/\W/,shift(@phrases)); push(@words,@wordss);}
         }
      }
   $first_sub="yes";
   foreach $word (@words)
      {
      $wordd = $word;
      $word=~tr/[A-Z]/[a-z]/;
      if ($first_sub eq "no") 
         {print TOO qq(, );}
      else 
         {$first_sub = "no";}

      if ($nosource)
         { print TOO qq(<a href="../nosource.html">$wordd</a> \n);}
      elsif ( -e "$src_dir/$word.f") 
         { print TOO qq(<a href="../../src/$word.f">$wordd</a>); }
      else 
         { 
         print ferr "$text_mans_dir/$fn.$ext: $src_dir/$word.f not found\n";
         print TOO qq($wordd);
         }
      }
   print TOO qq(</dl>\n);		
   }


sub callerscan # ---------------- callerscan -----------------------
   {
   #print " nosource: $nosource \n";
   #if ($nosource) {print " got nosource: $nosource \n"};

   return unless 
      ($line=~/^Caller|^CALLER/);
   $first="yes";
   @words=split(' ',$line);shift(@words);
      $line=join (' ',@words);$line=" ".$line;unshift(@lines,$line);
   @words=();
   print TOO $hr;
   print TOO qq(<dl><dt><strong>CALLER</strong>\n);
   print TOO qq(<p>\n<dd>);
   while(@lines)
      {
      $line=shift(@lines);
      if ($first eq "yes")
         {$first="no";}
      else
         {&lineprepare;}
      if ($line=~/\w/)
         {
         if ($line=~/^\S/) {last;}
         @phrases=split(' ',$line);
         while(@phrases) 
            {@wordss=split(/\W/,shift(@phrases)); push(@words,@wordss);}
         }
      }
   $first_caller="yes";
   foreach $word (@words)
      {
      $wordd = $word;
      $word=~tr/[A-Z]/[a-z]/;
      if ($first_caller eq "no") 
         {print TOO qq( , );}
      else 
         {$first_caller = "no";}

      if ($nosource)
         { 
         print TOO qq(<a href="../nosource.html">$wordd</a> \n);
         #print " sorry: $nosource \n";
         }
      elsif ( -e "$src_dir/$word.f") 
         { print TOO qq(<a href="../../src/$word.f">$wordd</a> \n);}
      else 
         { 
         print ferr "$text_mans_dir/$fn.$ext: $src_dir/$word.f not found\n";
         print TOO qq($wordd \n);
         }
      }
      print TOO qq(</dl>\n);		
   }	


# xxxxxxxxxxxxxxxxxxxxxxx  The main script begins here xxxxxxxxxxxxxxxxxxxxxx

$nosource = 0;
if ($ARGV[0])
   { $nosource = 1; }    # WANT SOURCE LINK TO NULL jan 05 al   
   
open (ferr, "> $ferr");  # Error output file

system("mkdir -p $raw_mans_dir") ;
   
opendir(oldman_dir,$text_mans_dir) 
   || die "Can't open: $text_mans_dir\n";


@files = readdir(oldman_dir);	# Read in ascii file names
closedir(oldman_dir);
$fnum = 0;			# Number of html files created
f:	
foreach $fn (@files)
   {
   if ($fn =~s/\.man//) 
      {$ext="man";} 
   else 
      {next f;}
   if ($fn eq "pocketcalc") 
      {next f;}
   unless (-T "$text_mans_dir/$fn.$ext") 
      {next f;}	# Not a valid file

   # Run cpp to add include files, option -P strips first line == cpp id #
   if ($fn =~/^cp/)
      {open(fro, "/usr/bin/cpp -P -traditional $text_mans_dir/$fn.$ext |") || die "Can't open: $text_mans_dir/$fn.$ext\n";}
   else
      {open(fro, "< $text_mans_dir/$fn.$ext") || die "Can't open: $text_mans_dir/$fn.$ext\n";}

   open(TOO, "> $raw_mans_dir/$fn.$htmlext") || die "Can't write to: $raw_mans_dir\n";
  
   # Find input file's atime & mtime
   ($atime, $mtime) = (stat("$text_mans_dir/$fn.$ext"))[8,9];

   @lines=<fro>;

   $remain="no";

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
         print ferr "$text_mans_dir/$fn.$ext: strange beginning\n";
         } 
      }

   &headerscan; 
   &purposescan;
   &warningscan;
   &seealsoscan;
   &usagescan;
   &notesscan;
   &subscan;
   &callerscan; 

   if ($#lines >0) 
      {$remain="yes";print TOO qq(<pre>$line\n);}

   while(@lines)
      {     #Any other stuff that might have remained
      $line=shift(@lines);
      &lineprepare;
      print TOO qq($line\n);
      }
   if ($remain eq "yes") {print TOO qq(</pre>\n);}

   print TOO qq(</body>\n);
   $fnum++;
   close(fro);
   close(TOO);

   # Alter the output files atime & mtime
   utime($atime, $mtime, "$raw_mans_dir/$fn.$htmlext");

   #print " Created: $raw_mans_dir/$fn.$htmlext \n";
   }
print " \n";
print " Created: $fnum  HTML SPIDER manual files in: $raw_mans_dir \n";
print " ------------------------\n";

# Create index of operations for local and sendaway use -------------------------------
$local_index="yes";
for ($i = 1; $i <=2; $i++)
    {
    open(indexfile,">$index_dir/operations_doc.html") 
        || die "Can't write operations_doc.html into $index_dir\n";

    print indexfile qq(<head>\n);
    print indexfile qq(<title>SPIDER: An Index of SPIDER Operations</title>\n);
    print indexfile qq(</head>\n);
    print indexfile qq(<body>\n);

    # put a NAME link in for the alphabet selector
    print indexfile qq(<a name="alpha-sel"></a>\n);
    print indexfile qq(<h2 align="center">MANUAL OF SPIDER OPERATIONS</h2>\n);

    print indexfile qq(<p>Use the <a href="crossref.html">Cross-Reference Index</a>
                       to find operation for a desired task. \n);
    print indexfile qq(A <strong>||</strong> symbol indicates the operation contains OMP parallel directives. \n);
    print indexfile qq(A <strong>||</strong>MPI symbol indicates the operation contains both OMP and MPI parallel directives. \n);

    print indexfile qq(<em>See also: SPIDER's <a href="user_doc.html#pocketcalc">Pocket Calculator</a> </em><br />\n);
    print indexfile qq(<dl>\n);

    # Make an alphbetical indexing hyperlist to go to operations quickly (al)
    $oper1old = 'A';
    @operindex = ('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z');

    print indexfile qq(<h2>\n);
    foreach $letter (@operindex)
      { print indexfile qq(<a href="#$letter">$letter</a> ); }

    print indexfile qq(\n</h2>\n);
    $oper1old = 'A';
    print indexfile qq(<table border=0>\n);
    print indexfile qq(<a name="$oper1old"></a>\n);

    # Put in index of operations 
    foreach $oper (sort keys(%definition))
        {
        $oper1 = substr($oper,0,1);
	
        if ($oper1 gt $oper1old ) 
            {  # Put in hypertag for this letter of alphabet (al)
            #print qq($oper1old $oper1 \n);
            $oper1old = $oper1;
	    
            if ($oper1 ne "B")
	       {	    
               # Put alphabet link line
               print indexfile qq(<h3>\n);
               foreach $letter (@operindex)
                  { print indexfile qq(<a href="#$letter">$letter</a> ); }
               print indexfile qq(\n</h3>\n);
               }
	    
            # Put link to the alphabet selector
            print indexfile qq(<tr><td><b> <a href="#alpha-sel"> </a> </b></td><tr></table>\n);

            print indexfile qq(<a name="$oper1"></a><table border=0>);
            }

        # Leave out some obsolete 'AP' operations from index
        if ($oper =~ m/AP MQ|AP NQ|AP RQ|AP MD|AP RD|AP RN|AP MH/)
            {next;}
        # Leave out some obsolete 'OR' operations from index
        if ($oper =~ m/OR MQ|OR NQ|OR 2|OR 2/)
            {next;}
        # Leave out some obsolete 'AD'.. operations from index
        if ($oper =~ m/^AD$|^AD F$|^AD R$|^SU$|^MU$|^MU M$|^MU O$|^MU D$/)
            {next;}
        # Leave out some obsolete 'BP' operations from index
        if ($oper =~ m/BP 32F OLD|BP 3F OLD/)
            {next;}
        # Leave out some obsolete 'RT' operations from index
        if ($oper =~ m/RT M|RT B|RT 3|RT 3D|RT 3DS|RT C|RT 3A|RT 3AS|RT 3L|RT 3LS/)
            {next;}
        # Leave out some obsolete 'AS' operations from index
        if ($oper =~ m/AS AD|^AS$/)
            {next;}
        # Leave out some obsolete 'RTD' operations from index
        if ($oper =~ m/^RTD/)
            {next;}
        # Leave out obsolete 'RT' operation from index
        if ($oper =~ m/^RT$/)
            {next;}
        # Leave out obsolete 'PK C' operations from index
        if ($oper =~ m/^PK C$|^PK DC$|^PK 3$/)
            {next;}
        # Leave out obsolete 'RO' operations from index (use 'RO SD')
        if ($oper =~ m/^RO$/)
            {next;}
        # Leave out obsolete 'FI X' operation from index (use 'FI H')
        if ($oper =~ m/^FI X$/)
            {next;}
        # Leave out obsolete 'HI *' operations from index (use 'HIS')
        if ($oper =~ m/^HI$|^HI D|^HI M|^HI R|^HI T/)
            {next;}
        # Leave out some obsolete 'PJ *' operations from index
        if ($oper =~ m/^PJ AT|^PJ COL|^PJ SHAD|^PJ ST/)
            {next;}
        # Leave out obsolete 'HD D' operation from index
        if ($oper =~ m/^HD D/)
            {next;}


        if (($local_index eq "yes") || ($localtype{$oper} eq "man"))
            {    
            # Print the index line for the operation
            print indexfile qq(<tr><td><a href="$htmlurl{$oper}"> <strong>$oper</strong></a> <td> $definition{$oper}</td></tr>\n);
            }
        }
    print indexfile qq(</table>\n);
    print indexfile qq(</dl>\n);
    print indexfile qq(</body>\n);
    close(indexfile);
    $local_index = "no";
    }
print " Created operation index: $index_dir/operations_doc.html \n";
print " ------------------------\n";

close(ferr);




