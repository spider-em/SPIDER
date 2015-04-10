#!/usr/bin/env python

import os, sys
import spyder

# this definition wraps Spider's FS operation in a Python function.
def file_statistics(file):
    mySpider.toSpider("FS [max],[min],[avg],[sig]", file)
    a1 = mySpider.getreg('max')
    a2 = mySpider.getreg('min')
    a3 = mySpider.getreg('avg')
    a4 = mySpider.getreg('sig')
    return a1,a2,a3,a4



# you need to give it the full path to the SPIDER executable
spider_executable = "/net/bali/usr1/spider/bin/spider"

# start a Spider session
mySpider = spyder.SpiderSession(spiderexec=spider_executable, dataext="dat")

# Send commands to Spider via the toSpider method.
# All arguments are strings (in quotes), or variables that evaluate to strings.
mySpider.toSpider("MO")
mySpider.toSpider("tmp001")
mySpider.toSpider("64,64")   # each line can be sent to Spider individually
mySpider.toSpider("T")

filename = "tmp002"

# or use one line, separated by commas. 
mySpider.toSpider("MO",filename, "64,64", "R", "N")

# Get register values from Spider: first send variable names to Spider,
mySpider.toSpider("FS [max],[min],[avg],[sig]", filename)

# then use getreg() to get the values.
a1 = mySpider.getreg('max')
a2 = mySpider.getreg('min')
a3 = mySpider.getreg('avg')
a4 = mySpider.getreg('sig')
print "file statistics for %s: max %f min %f avg %f std %f" % (filename,a1,a2,a3,a4)

# You can 'wrap' a Spider operation to behave as a Python function (defined above)
max, min, avg, sig = file_statistics(filename)

print "file statistics for %s:" % (filename)
print "max %f, min %f, avg %f, std %f" % (max, min, avg, sig)

# How about a shortened version of toSpider, to save typing.
spi = mySpider.toSpider
spi("DC S", filename, "_1", "2,2")


# end the Spider session
mySpider.close()




