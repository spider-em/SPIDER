 set encoding iso_8859_1
 set xlabel "Eigenvalue number"
 set ylabel "%"
 set xrange [0.2:08+.5]
 set boxwidth 0.5
 set style fill solid
 plot \
  "data/cas_ca_eigpct_doc.dat" using 1:4 title 'eigenvalues' with boxes
 set terminal postscript eps enhanced monochrome 'Times-Roman' 24
 set output "data/cas_ca_eigpct.eps"
 replot
