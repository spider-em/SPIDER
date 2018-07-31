#!/usr/bin/env python

# SOURCE:  spider/tools/bin/pytest.py
# PURPOSE: Testing imports

"""
test.py : 
pytest.py 
"""

import os, sys
import time

from   tkMessageBox        import askyesno

print ' '                                           
print '\n'.join(sys.path)                           
print ' '                                           

print  "pytest.py           import numpy"               
import numpy                                          

print  "pytest.py,          import SpiderImageSeries" 
from   Spider               import SpiderImageSeries  

print  "pytest.py           import SpiderImagePlugin" 
import SpiderImagePlugin                              
   
print  "pytest.py,          import Tkinter *"         
from   Tkinter import *

print  "pytest.py,          import ImageTk"           
from   PIL                  import ImageTk

print  "pytest.py           import Image"             
from   PIL import Image                              

print  "pytest.py           import Pmw"               
import Pmw                                           

#print '\n'.join(sys.path)                            
print ' '                                           

EXIT = '_exit_'

