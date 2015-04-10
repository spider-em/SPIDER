set title " Correlation"
set xlabel "Angle"
set ylabel "CC"
set xrange [0:   359]
set yrange [-1.0:1.0]
plot "jnkgp_DATA.tst" using 1:  2 title "Radius:   5" , \
 "" using 1:  3 title "Radius:   7" , \
 "" using 1:  4 title "Radius:   9" , \
 "" using 1:  5 title "Radius:  11" , \
 "" using 1:  6 title "Radius:  13" , \
 "" using 1:  7 title "Radius:  15" , \
 "" using 1:  8 title "Radius:  17" , \
 "" using 1:  9 title "Radius:  19" , \
 "" using 1: 10 title "Radius:  21" , \
 "" using 1: 11 title "Radius:  23" , \
 "" using 1: 12 title "Radius:  25" , \
 "" using 1: 13 title "Radius:  27"
