max_height_bad = 7000
average_bad    = 1000
std_dev_bad    = 300

max_height_good = 500
average_good    = 1200
std_dev_good    = 300

histbad(x)=max_height_bad*exp(-(x-average_bad)**2/std_dev_bad**2)
histgood(x)=max_height_good*exp(-(x-average_good)**2/std_dev_good**2)

plot  [0:3000][] \
'histcccbad.dat'  using 3:4 title 'bad'  with boxes, \
'histcccgood.dat' using 3:4 title 'good' with boxes, \
histbad(x) title 'fitted bad', \
histgood(x) title 'fitted good'

