#!/usr/bin/env python

import os, sys
from struct import unpack
import re

"""
There are 2 streams:
 The Python program sends commands to Spider as if they were typed at the
 .OPERATION: prompt.

 The only information Python gets from Spider are register values, via
 an external fifo pipe.

 The spider session is started by creating an instance of the SpiderSession
 class:

      sp = SpiderSession(dataext='dat')

 Then you use the instance methods (functions of sp)
 - send commands to Spider with sp.toSpider("op", "infile","outfile","args")
 - get register values from Spider w/ sp.getreg("[var]")
"""

external_pipename = "TMP_SPIDER_PIPE.pipe"
      
isSpiderReg = re.compile("[xX]\d\d")

def isSpiderRegister(s):
    " returns 1 if string matches pattern: 'X' or 'x' followed by 2 digits"
    if isSpiderReg.match(s): return 1
    else: return 0 

class SpiderSession:
    def __init__(self, spiderexec=None, dataext='.dat'):
        # spider executable
        if spiderexec == None:
            if os.environ.has_key('SPIDER_LOC'):
                self.spiderexec = os.path.join(os.environ['SPIDER_LOC'],'spider')
            else:
                self.spiderexec = '/usr/local/spider/bin/spider'
        else:
            self.spiderexec = spiderexec

        self.dataext = dataext
        if dataext[0] == '.': self.dataext = dataext[1:]
        
        # create the external fifo pipe to get [registers] from Spider
        self.pipename = external_pipename
        if os.path.exists(self.pipename):
            try: os.remove(self.pipename)
            except: pass
        os.mkfifo(self.pipename)

        # Start spider process, initialize with some MD commands.
        self.spider = os.popen(self.spiderexec, 'w')
        self.toSpider(self.dataext)
        self.toSpider("MD", "TERM OFF")
        self.toSpider("MD", "RESULTS OFF")
        self.toSpider("MD", "PIPE", self.pipename)
        # open the pipe as a readable file (to get regs from Spider)
        self.fromSpider = open(self.pipename, 'r')
        
        # In Linux, the first register from Spider is piped in a different
        # format than subsequent calls to getreg().
        # The next 4 lines initialize the register pipeline with the first call.
        if sys.platform[:5] == 'linux':
            self.toSpider("[i] = 3.241","PI REG", "[i]")
            rp = ""
            while len(rp) < 9:
                rp += self.fromSpider.readline()
            a,regval,c = unpack('ffc',rp)

    def toSpider(self, *args):
        " each item is a line sent to Spider"
        for item in args:
            self.spider.write(str(item) + '\n')
        self.spider.flush()

    def getreg(self, varname):
        varname = varname.strip()
        if varname[0] != '[' and not isSpiderRegister(varname):
            varname = '[' + varname + ']'
        self.toSpider("PI REG", varname)

        rp = ""   # sometimes it takes multiple readlines to get all 13 bytes
        while len(rp) < 13:
            rp += self.fromSpider.readline()
        
        if sys.platform[:5] == 'linux':
            a,b,regval,c = unpack('fffc',rp)            
        elif sys.platform[:4] == 'irix':
            regval,c = unpack('fc',rp)
            
        return regval
    
    def close(self, delturds=1):
        self.toSpider("en d")          # end the spider process,
        self.fromSpider.close()        # close the fifo pipe
        try: os.remove(self.pipename)
        except: pass
        if delturds:
            for file in ['fort.1', 'jnkASSIGN1', 'LOG.'+self.dataext]:
                if os.path.exists(file):
                    try: os.remove(file)
                    except: pass

    
# --------------------------------------------------------------
if __name__ == '__main__':

    sp = SpiderSession(dataext='dat')

    sp.toSpider("[size]=117")
    s = sp.getreg('size')
    print "---------------------- size =  %f" % s
    sp.toSpider("x11=7.7")
    s = sp.getreg('x11')
    print "---------------------- x11 =  %f" % s
    sp.close()
