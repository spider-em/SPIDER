# GG.py : GLOBAL VARIABLES - GUI
"""
SPIRE - The SPIDER Reconstruction Engine
Copyright (C) 2006-2012  Health Research Inc.

HEALTH RESEARCH INCORPORATED (HRI),
ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455

Email:  spider@wadsworth.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
"""

import os, sys
from Tkinter import *
import Queue

#import Spire  #al Sept 2010
import LocalVars

applicationname = "SPIRE: SPIDER Reconstruction Engine"

#applicationversion = Spire.__version__  #al
__version__ = '1.5.5'#al oct12
applicationversion = __version__#al

bigLabeltxt = "SPIRE: Spider Reconstruction Engine"
mainImage = "pics/ribo1.gif" #"pics/goat.gif"   #"atom.gif")

defaultconfigdir = LocalVars.defaultconfigdir
#Spire.__path__ = [os.path.join(os.path.dirname(defaultconfigdir), 'Spire')] #al
__path__ = [os.path.join(os.path.dirname(defaultconfigdir), 'Spire')] #al
#SPIREHOME = Spire.__path__[0] #al
SPIREHOME = __path__[0]#al

defaultconfig = os.path.join(defaultconfigdir, "SingleParticle.xml")

def setMainImage(img):  # Finds name of main image
    global mainImage
    mimg = img
    if os.path.exists(img):
        mimg = img
    else:
        im = os.path.join(SPIREHOME, img)
        if os.path.exists(im):
            mimg = im
        else:
            im = os.path.join(SPIREHOME, "pics", img)
            if os.path.exists(im):
                mimg = im
                
    if os.path.exists(mimg):
        mainImage = mimg
        return 1
    else:
        return 0

setMainImage(mainImage)

hostlist = []
ssh_hosts = []
webBasedir = ""

# get spider, jweb commands from LocalVars.py
spiderProgram = LocalVars.spidercmd
jwebcmd = LocalVars.jwebcmd
if hasattr(LocalVars,'hostlist'):
    hostlist = LocalVars.hostlist
if hasattr(LocalVars,'ssh_hosts'):
    ssh_hosts = LocalVars.ssh_hosts
if hasattr(LocalVars,'webBasedir'):
    webBasedir = LocalVars.webBasedir
if hasattr(LocalVars,'defaultconfig'):
    defaultconfig = LocalVars.defaultconfig

MAXDIRSIZE = 17000 

# platform-dependent variables
if sys.platform[:4] == 'irix':
    platform = "irix"
elif sys.platform[:5] == 'linux':
    platform = "linux"
elif os.name == 'nt':
    platform = 'nt'
else:
    print "UNSUPPORTED PLATFORM: %s, %s" % (sys.platform, os.name)

if platform == 'irix':
    if webBasedir == "":
        #webBasedir = 'file://localhost/net/bali/usr1/spider/docs/'
        webBasedir = 'http://www.wadsworth.org/spider_doc/spider/docs/'
    defaultEditor = "nedit"
    
elif platform == 'linux':
    if webBasedir == "":
        webBasedir = 'http://www.wadsworth.org/spider_doc/spider/docs/'
        #webBasedir = 'file:///usr/local/spider/docs/'
    defaultEditor = "xedit"  # ..for KDE. otherwise use "emacs"?
    
elif os.name == 'nt':
    if webBasedir == "":
        webBasedir = "C:/Program Files/SPIDER/docs/"
    defaultEditor = "Wordpad.exe"

HelpURL = os.path.join(webBasedir, "index.html")
SpireHelpURL = os.path.join(webBasedir, "spire/index.html")

# dictionary that links extensions to external applications
if hasattr(LocalVars, 'AppDict'):
    AppDict = LocalVars.AppDict
else:
    AppDict = {".gnu" : 'gnuplot',".gp"  : 'gnuplot',".pdf" : 'acroread',
               ".ps"  : 'gs',".tif" : 'jweb', ".gif" : 'jweb',".jpg" : 'jweb' }

# careful that this is accessed only by the main thread
disp_fnums = StringVar()

topwindow = ""
msg = ""
havebusy = 1
balloon = 1

class GUIOutput:
    " implements a thread-safe text box. ie, threads can write to it "
    def __init__(self, msgbox, logfile=None):
        self.msgbox = msgbox
        self.logfile = logfile
        self.uselogfile = 0
        self.queue = Queue.Queue() # write commands put data on the queue
        self.update_que()
        
    def __call__(self, string, color=None):
        self.write(str(string) + '\n', color)
        
    def write(self, txtstring, color=None):
        # item in the queue is a list: [text, color]
        self.queue.put([txtstring, color])

    def update_que(self):
        try:
            while 1:
                item = self.queue.get_nowait()
                if item == None:
                    pass
                else:
                    if len(item) < 2:
                        print "Error: text sent to Main window in incorrect format:"
                        print str(item)
                    else:
                        txtstring = item[0]
                        color = item[1]
                        self.msgbox.component('text').configure(state='normal')
                        if color == None:
                            self.msgbox.component('text').insert('end', txtstring)
                        else:
                            self.msgbox.component('text').insert('end', txtstring, color)
                        self.msgbox.component('text').see('end')
                        self.msgbox.component('text').update_idletasks()
                        self.msgbox.component('text').configure(state='disabled')
                        if self.uselogfile:
                            fp = open(self.logfile, 'a')
                            fp.write(txtstring)
                            fp.close()
        except Queue.Empty:
            pass
        
        self.msgbox.after(100, self.update_que)

# End GUIOutput -------------------------------

class SysFonts:
    """ sets up default fonts """
    def __init__(self):
        self.font = {'small' : ('Helvetica',11,'normal'),
                     #'medium': ('Arial', 12, 'normal'),
                     #'large' : ('Times', 16, 'bold'),
                     #'output': ('Courier', 12, 'normal')}
                     'medium': ('Arial', 11, 'bold'),
                     'large' : ('Times', 14, 'bold'),
                     'output': ('Courier', 11, 'normal')}
        self.font['outbold'] = (self.font['output'][0],
                                self.font['output'][1], 'bold')
    def setsmallfont(self):
        self.font = {'small' : ('Helvetica',10,'normal'),
                     'medium': ('Arial', 10, 'bold'),
                     'large' : ('Times', 12, 'bold'),
                     'output': ('Courier', 10, 'normal')}
        self.font['outbold'] = (self.font['output'][0],
                                self.font['output'][1], 'bold')

    def stringify(self):
        return str(self.font)
        
    def setFont(self,whichfont,newfont):
        self.font[whichfont] = newfont
        if whichfont == 'output':
            self.font['outbold'] = (self.font['output'][0],
                                    self.font['output'][1], 'bold')

    def show(self):
        print self.font['small']
        print self.font['medium']
        print self.font['large']
        print self.font['output']

class SysColors:
    """ defaults system colors """
    def __init__(self):
        self.background = 'gray90'
        self.bgd01 = '#aaaaee'  # '#99ccff'
        self.bgd02 = '#ddddff'  # '#3399ff' 
        self.bgd03 = '#222266'  # '#0000cc' 
        self.entrybgd = 'white'
        self.foreground = 'black'
        self.attribs = ['background', 'bgd01', 'bgd02', '.bgd03',
                        'entrybgd', 'foreground']

    def stringify(self):
        "return self as stringified dictionary"
        d = {}
        for a in self.attribs:
            if hasattr(self, a):
                d[a] = eval('self.' + a)
        return str(d)

    def setAttributes(self,dict):
        " given a dictionary, resets all its attributes "
        keys = dict.keys()
        for key in keys:
            if hasattr(self,key):
                value = dict[key]
                s = "self.%s = '%s'" % (key, value)
                exec(s)
            

editorlist = ['jot', 'nedit','gedit', 'xedit', 'emacs', 'kwrite']

"""
Graphical preferences are UserPrefs, while nongraphical ones
are in SystemPrefs. Both SystemPrefs and UserPrefs put in
SysInterface and saved in $HOME/.spidui
"""
class UserPrefs:
    """ user preferences for fonts, colors, etc. """
    def __init__(self):
        self.balloonHelp = 1
        self.useDatabase = 1
        self.delAsk = 1      # 1 = ask before deleting
        self.saveResults = 0   # refers to the spireout file (Sep,05)
        self.MaxResultLines = 100
        self.MaxDisplayFiles = 100  # for tables of output files
        self.editor = defaultEditor
        self.savelogfile = 0
        self.useFilenumsEntry = 1  # entry takes precedence over Filenums file
        self.startupJweb = 0
        self.displayProgram = 'jweb' # or qview or xv
        self.useLocalPrefsFile = 0  # 0 = /home/user/.spire, 1 = projdir/.spire
        
        self.font = SysFonts()
        self.colors = SysColors()

remotecmdlist = ['rsh', 'ssh']
netexeclist   = ['local', 'remote']
shellcmdlist  = ['csh', 'sh']

class SystemPrefs:
    def __init__(self):
        self.MaxDirsize = MAXDIRSIZE # size of stat[directory], not no.of files
        self.filenumFile =   "" #'filenums'
        self.filenumSymbol = '[FILENUMS]'
        self.paramFile = 'params'
        self.configFile = defaultconfig
        self.configdir = defaultconfigdir # directory of default configurations
        self.logfile = 'spirelog'
        
        #self.host = os.environ["HOST"]
        self.remotecmd = 'rsh'     # rsh | ssh
        self.netexec = 'local'     # local | remote
        self.shellcmd = 'csh'      # sh | csh
        
        self.checkPrior = 0
        self.helpURL = HelpURL
        # How the spider command is called on this system
        self.spider = spiderProgram
        self.spireout = 'spireout'  # outputs from results file
        self.AppDict = AppDict

class SysInterface:
    """ a container for both user & system classes, so both can
        be pickled into a single file """
    def __init__(self):
        self.user = UserPrefs()
        self.system = SystemPrefs()
        self.file = "" # os.path.join(os.environ['HOME'],".spire")


# Then create shorthand links

interface = SysInterface()  # instance
prefs = interface.user
sysprefs = interface.system

if prefs.useLocalPrefsFile == 0 or not hasProject():
    interface.file = os.path.join(os.environ['HOME'],".spire")
else:
    interface.file = os.path.join(GB.P.projdir, '.spire')

spid = prefs.font  # instance with defaults
colors = prefs.colors

sysbgd = colors.background
bgd01 = colors.bgd01
#bgd01 = '#F9C396'
bgd02 = colors.bgd02
bgd03 = colors.bgd03

brdr = 2
brelief = 'raised'   # buttons
frelief = 'raised' # frames
padx=5
pady=5
pad2x = 2 * padx
pad2y = 2 * pady
