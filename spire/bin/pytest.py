#!/usr/bin/env python
#
# SOURCE:  spider/spire/bin/pytest.py
#
# PURPOSE: Testing imports
#
# Spider Python Library
# Copyright (C) 2006-2018  Health Research Inc., Menands, NY
# Email:    spider@health.ny.gov

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
# import SpiderImagePlugin  # 2018 al                             
from   Spider               import SpiderImagePlugin # 2018 al                             
   
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

