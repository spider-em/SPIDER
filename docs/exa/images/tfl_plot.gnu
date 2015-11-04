set terminal jpeg color enhanced "Helvetica" 20
set output "tfl_plot.jpg"
ylabel "Amplitude"
set xlabel "Radius (pixels)"
set title "Amplitude vs Radius"
plot 'tfl_doc.dat' using 1:3 notitle with lines
