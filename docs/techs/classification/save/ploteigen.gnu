# SET NUMBER OF FACTORS HERE
LAST = 8

# SET POSTSCRIPT OUTPUT FILENAME HERE
#set output 'EIGENVALUES.eps'

#set terminal postscript eps enhanced monochrome 'Times-Roman' 24
set encoding iso_8859_1
set xlabel 'Eigenvalue number'
set ylabel '%'
set boxwidth 0.5
LAST = LAST + 0.5
#set xrange [0.2:LAST]
plot 'jnk-eigpct_doc.dat' using 1:4 title 'eigenvalues' with boxes


#plot 'eigenvalues.dat' using 1:4 title 'eigenvalues' with boxes

