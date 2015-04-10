# To store plot as a postscipt file (to print it, for example)
# remove "#" from the following lines
#set term post  
#landscape color solid "Times-Roman" 24
#set output 'plotpower1.ps'

set xrange [60:200]
set yrange [.0002:.0007]
plot 'roo001.ext' using 1:3 title "mic001" with lines,\
'roo003.ext' using 1:3 title "mic003" with lines,\
'roo013.ext' using 1:3 title "mic013" with lines,\
'roo014.ext' using 1:3 title "mic014" with lines,\
'roo020.ext' using 1:3 title "mic020" with lines
