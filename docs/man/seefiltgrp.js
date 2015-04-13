/*
 seefiltgrp : print the links for the filter operations.
       Parameter OP is the operation from the calling page, to
       prevent it from appearing in it's own SEE ALSO list.
*/

function seefiltgrp(op) {
   document.write("<DL><DT><STRONG>SEE ALSO</STRONG><P>")
   document.write("<DD>")

   if (op != "FD") { 
      document.write('<A HREF="fd.html"><STRONG>FD</STRONG></A> [Filter according to a Document file]<BR>') }
   if (op != "FF") { 
      document.write('<A HREF="ff.html"><STRONG>FF</STRONG></A> [Fourier Filter]<BR>') }
   if (op != "FF PLOT") { 
      document.write('<A HREF="ffplot.html"><STRONG>FF PLOT</STRONG></A> [Fourier Filter - design filter]<BR>') }
   if (op != "FP") { 
      document.write('<A HREF="fp.html"><STRONG>FP</STRONG></A> [Fourier interpolation]<BR>') }
   if (op != "FQ") { 
      document.write('<A HREF="fq.html"><STRONG>FQ</STRONG></A> [Filter Quick] <BR>') }
   if (op != "FQ NP") { 
      document.write('<A HREF="fqnp.html"><STRONG>FQ NP</STRONG></A> [Filter Quick - No Padding] <BR>') }
   if (op != "FT") { 
      document.write('<A HREF="ft.html"><STRONG>FT</STRONG></A> [Fourier Transformation]<BR>') }

   document.write("</DL><P>")
}
