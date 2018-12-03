#!/usr/bin/env python
#
# SOURCE:   verifybyview.py 
#
# MODIFICATIONS:
#    TO DO: when changing #columns, link max to current
#    TO DO: checkbutton to enable/disable paint-select
#    TO DO: use arbitrary label for class montage
#    TO DO: add arbitrary label to particle montage
#    TO DO: in createMontage, use unsaved select values
#    2016-07-05 -- import montage
#    2012-04-25 -- can search for BYHAND files (rather than reopen classes to select them)
#    2012-04-24 -- if initial filename field is blank, will use default
#    2012-04-24 -- if particle list not found, will pop up error window
#    2012-04-24 -- added exit button for initial-filename window
#    2009-07-28 -- different class colors noted in good-class doc
#    2009-05-29 -- bound carriage-return to update directory
#    2009-05-14 -- set maximum #columns
#    2009-05-08 -- if BYHAND file exists, colors deselected particles also
#    2009-04-01 -- added shortcut for filename toggle
#    2009-02-26 -- only class averages need to be in initial_directory
#               -- integer labels no longer displayed as float
#               -- no longer asks for image range for large montages
#    2008-11-09 -- "Templates" window works more nicely
#    2007-12-20 -- attached GNU GPL
#    2007-11-27 -- distinguished preferences for size in class-average & particle montages
#    2007-11-27 -- bug fix in selectAll and clearAll
#    2007-10-05 -- added changelog
#    2007-10-02 -- bug fix in reading existing good-classes file
#    2007-09-04 -- can read stacks
#               -- carriage-return bound to "OK" buttons
#    2007-09-03 -- open classes now colored yellow
#               -- text label in montage now user-defined
#    2007-05-30 -- converts SPIDER images internally to "luminosity"
#    2006-12-20 -- window position saved
#    2006-10-20 -- "invert" button colored obnoxiously for visibility
#    2006-08-29 -- saves settings to ".verifybyview"
#    2006-03-17 -- uses new "Spiderutils" library instead of older "spiderutils"
#    2005-05-12 -- right click now selects whole class
#    2005-04-19 -- added sliders for contrast and brightness (slow)
#    2005-04-12 -- bug fix: red, deselected particles had been inadvertently saved to output doc
#    2005-04-08 -- backs up pre-existing output files, rather than overwriting them
#    2005-04-07 -- reads in pre-existing good-particle lists
#               -- program remembers colors and display size for particle montages
#    print "verifybyview.py, Modified 2016 July 5"  
#
# Spider Python Library
# Copyright (C) 2006-2018  Health Research Inc., Menands, NY
# Email:    spider@health.ny.gov

import os,sys,re
import Pmw
import montage 

from   Tkinter           import *
from   tkFileDialog      import askopenfilename, askdirectory
from   tkMessageBox      import showerror
from   PIL               import Image, ImageTk
from   PIL               import ImagePalette
from   math              import sqrt
from   Spider            import SpiderImagePlugin #2018 import *
from   Spider            import Spiderutils

def backup(filename):
    if os.path.exists(filename):
        found_vacancy = 0
        tiebreaker = 0
        shortdir = os.path.basename(os.path.dirname(filename))
        while not found_vacancy:
            test_filename = filename + '_' + str(tiebreaker)
            if os.path.exists(test_filename):
                tiebreaker += 1
            else:
                found_vacancy = 1
                short_old = os.path.join(shortdir, os.path.basename(filename))
                short_new = os.path.join(shortdir, os.path.basename(test_filename))
                print 'Renamed', short_old, 'to', short_new
                os.rename(filename, test_filename)

def parseGeometry(geometry):
    # Adapted from http://effbot.org/tkinterbook/wm.htm#wm.Wm.wm_geometry-method

    m = re.match("(\d+)x(\d+)([-+]\d+)([-+]\d+)", geometry)
    return map(int, m.groups())

def clickOK(parent,win):  
    win.destroy()

class selectionClass:
    def __init__(self, value, key, color, label, activecolor=None):
        self.value = value # int
        self.key = key     # string of a number
        self.color = color
        self.label = label # for radiobutton menu
        self.active = ""
        if activecolor != None:
            self.active = activecolor

# ******************************************************************************
# class pmontage
#       input image list can be a list of filenames or a previously loaded
#       list of Image.images
#       Create the toplevel window outside and pass it in
# ******************************************************************************

class pmontage:
    def __init__(self, master, title=None,
                 partdoc=None, filerange=None, parttemplate=None, 
                 ncol=None, maxcols=12, useLabels=0, savefilename=None, reverse=0, 
                 sizeVar=1, selectedColor='1', deselectedColor='Red', 
                 contrast=0, brightness=127, xcoord=10, ycoord=32):

        self.partdoc = partdoc
        self.filerange = filerange
        self.parttemplate = parttemplate

        # load imagelist
        self.loadImagelist()
        if self.error:
#            print "loadImagelist: self.error=", self.error
            return

        if title == None: title = "Images"
        self.top = master   #Toplevel(master)
        self.top.title(title)
        self.partLabels = useLabels
        # set some size variables
        self.xmax, self.ymax = self.top.winfo_screenwidth(), self.top.winfo_screenheight()

        fn = Spiderutils.template2filename(self.parttemplate,self.imagelist[0])
        im = Image.open(fn)
        if im.istack:
            self.isStack = True
#            print "Particles are stacks:", fn, self.parttemplate,self.imagelist[0]
        else:
            self.isStack = False
#            print "Particles are unstacked"
        
        # get size of a single image
        self.xsize,self.ysize = im.size

        self.maxcols = maxcols
	self.maxcolVar = IntVar()
	self.maxcolVar.set(maxcols)
#	print "maxcolVar",self.maxcolVar.get(),type(self.maxcolVar.get())

        if ncol == None:
            self.ncol = self.montagesize()
        else:
            self.ncol = ncol

	if savefilename != None:
            self.savefilename = savefilename
        else:
            self.savefilename = ""

        # some Tk Variables
        self.ncolVar = StringVar()

        self.sizeVar = IntVar()
        self.sizeVar.set(sizeVar)

        self.showVar = IntVar()
        self.showVar.set(self.partLabels)
        self.bd = 2  # border around images
        sysbgd = self.systembackground = "#d9d9d9"
        self.actbgd = "#ececec"

        # a selection class = ( number, 'value (name)', 'key (desc)', 'label (highlight)')
        self.selectClasses = {}
        self.selectClasses['0'] = selectionClass(0,'0', sysbgd, 'Deselect',self.actbgd)
        self.selectClasses['1'] = selectionClass(1,'1','green', 'Green',  'light green')
        self.selectClasses['2'] = selectionClass(2,'2','black', 'Black',  'gray50')
        self.selectClasses['3'] = selectionClass(3,'3','blue',  'Blue',   'light blue')
        self.selectClasses['4'] = selectionClass(4,'4','gray75','Lt gray','gray50')
        self.selectedColor = StringVar()
        self.selectedColor.set(selectedColor)
        self.deselectedColor = StringVar()
        self.deselectedColor.set(deselectedColor)
        self.selectedColor.trace_variable('w', self.selectcallback)

        self.numgood = 0  # number of good particles
        self.last = len(self.imagelist)

        self.reverse = reverse
        if reverse == 1: self.revOrder()

        # bind shortcuts
        master.bind('<Control-s>', self.saveSelections)
        master.bind('<Control-r>', self.readSelections)
        master.bind('<Control-w>', self.closeWindow)
        master.bind('<Control-a>', self.selectAll)
        master.bind('<Control-i>', self.invertSelect)
        master.bind('<Control-v>', self.invertSelect)
        master.bind('<Control-f>', self.displayFilenames)
        master.bind('<Control-t>', self.test)
        master.bind('<Delete>'   , self.clearAll)

        # contrast stuff
        self.contrast = DoubleVar()
        self.brightness = DoubleVar()
        self.maxcontrast = 127
        self.maxbrightness = 255
        self.contrast.set(contrast)
        self.brightness.set(brightness)

        # window positions
        self.xcoord = xcoord
        self.ycoord = ycoord

        # draw stuff
        self.makeMenus()
        self.createMontage()

    def test(self, event=None):
        print 'test winfo:', self.top.winfo_rootx(), self.top.winfo_rooty()
        geom = self.top.geometry()
        print 'winfo_geometry:', geom
        parsed = parseGeometry(geom)
        print 'parsed:', parsed[2], parsed[3]

    def closeWindow(self, event=None):
        geom = self.top.geometry()
        parsed = parseGeometry(geom)

        offset_x = self.top.winfo_rootx() - parsed[2]
        offset_y = self.top.winfo_rooty() - parsed[3]

        self.xcoord = self.top.winfo_rootx() - offset_x
        self.ycoord = self.top.winfo_rooty() - offset_y

        self.top.destroy()

    def selectcallback(self, name, index, mode):
        "called whenever self.selectedColor changes "
        key = self.selectedColor.get()
        print 'New color:', self.selectClasses[key].color

    def makeMenus(self):
        # ------- create the menu bar -------
        self.mBar = Frame(self.top, relief='raised', borderwidth=1)
        self.mBar.pack(side='top', fill = 'x')
        self.balloon = Pmw.Balloon(self.top)

        # Make the Display menu
        Dspbtn = Menubutton(self.mBar, text='Display', relief='flat')
        Dspbtn.pack(side=LEFT, padx=5, pady=5)
        Dspbtn.menu = Menu(Dspbtn, tearoff=0)

        Dspbtn.menu.add_command(label='no. columns', underline=0,
                                    command=self.displayParms)
        Dspbtn.menu.add_checkbutton(label='show filenames', underline=0,
                                    variable=self.showVar,
                                    command=self.displayFilenames_1)
        Dspbtn.menu.add_command(label='reverse order', command=self.revRedisplay)

        # 'size' has a submenu of checkbuttons
        Dspbtn.menu.sizes = Menu(Dspbtn.menu, tearoff=0)
        Dspbtn.menu.sizes.add_radiobutton(label='1/2', underline=0,
                                          variable=self.sizeVar, value=0,
                                          command=self.displayParms_2)
        Dspbtn.menu.sizes.add_radiobutton(label='1x', underline=0,
                                          variable=self.sizeVar, value=1,
                                          command=self.displayParms_2)
        Dspbtn.menu.sizes.add_radiobutton(label='2x', underline=0,
                                          variable=self.sizeVar, value=2,
                                          command=self.displayParms_2)
        Dspbtn.menu.add_cascade(label='image size', menu=Dspbtn.menu.sizes)

        # Make the Color submenu
        Dspbtn.menu.colors = Menu(Dspbtn.menu, tearoff=0)
        keys = self.selectClasses.keys()
        keys.sort()
        Dspbtn.menu.colors.good = Menu(Dspbtn.menu.colors, tearoff=0)
        Dspbtn.menu.colors.bad =  Menu(Dspbtn.menu.colors, tearoff=0)
        for key in range(1,4):  # selected colors
            sc = self.selectClasses[str(key)]
            Dspbtn.menu.colors.good.add_radiobutton(label = sc.label, underline=0,
                                        value = sc.value,
                                        foreground = 'white',
                                        background = sc.color,
                                        activebackground= sc.active,
                                        activeforeground= 'black',
                                        variable=self.selectedColor)
        Dspbtn.menu.colors.bad.add_radiobutton(label = 'Deselect', underline=0,
            value = self.systembackground, foreground = 'white', background = self.systembackground, 
            activebackground = self.actbgd, activeforeground= 'black', variable=self.deselectedColor)
        Dspbtn.menu.colors.bad.add_radiobutton(label = 'Red', underline=0,
            value = 'Red', foreground = 'white', background = 'Red', 
            activebackground = 'pink', activeforeground= 'black', variable=self.deselectedColor)
        Dspbtn.menu.colors.add_cascade(label='selected',   menu=Dspbtn.menu.colors.good)
        Dspbtn.menu.colors.add_cascade(label='deselected', menu=Dspbtn.menu.colors.bad)
        Dspbtn.menu.add_cascade(label='colors', menu=Dspbtn.menu.colors)

        Dspbtn.menu.add_separator()
        Dspbtn.menu.add_command(label='Read selection', underline=0,
                                command=self.readSelections)
        Dspbtn.menu.add_command(label='Save selection', underline=0,
                                command=self.saveSelections)
        Dspbtn.menu.add_separator()
        Dspbtn.menu.add_command(label='Close', underline=0,
                                     command=self.closeWindow)
        Dspbtn['menu'] = Dspbtn.menu

        # Make the Select menu
        Selectbtn = Menubutton(self.mBar, text='Select', relief='flat')
        Selectbtn.pack(side=LEFT, padx=5, pady=5)
        Selectbtn.menu = Menu(Selectbtn, tearoff=0)

        Selectbtn.menu.add_command(label='Select all', underline=0,
                                    command=self.selectAll)
        Selectbtn.menu.add_command(label='Clear all', underline=0,
                                    command=self.clearAll)
        Selectbtn.menu.add_command(label='Invert', underline=0,
                                    command=self.invertSelect)

        Selectbtn['menu'] = Selectbtn.menu

        # Make the Save+Close button
        Selbtn = Button(self.mBar, text='Save+Close', command=self.saveClose)
        Selbtn.pack(side=LEFT, padx=5, pady=5)

        # Make the Close button
        Clsbtn = Button(self.mBar, text='Close', command=self.closeWindow)
        Clsbtn.pack(side=LEFT, padx=5, pady=5)

        # Make the Invert button
        key = self.selectedColor.get()
        sc = self.selectClasses[key]
        badcolor = self.deselectedColor.get()
        Invbtn = Button(self.mBar, text='Invert', command=self.invertSelect, 
            background = sc.color, foreground = badcolor, 
            activeforeground = sc.color, activebackground = badcolor)
        Invbtn.pack(side=LEFT, padx=5, pady=5)

        # Make the Help menu
        Helpbtn = Menubutton(self.mBar, text='Help', underline=0, relief='flat')
        Helpbtn.pack(side=RIGHT, padx=5, pady=5)
        Helpbtn.menu = Menu(Helpbtn, tearoff=0)
        Helpbtn.menu.add_command(label='Keyboard shortcuts', underline=0, 
                                 command=self.shortcuts)
        Helpbtn['menu'] = Helpbtn.menu

        # Make the contrast slidebars
        f2 = Frame(self.top)
        scon = Scale(f2, label='contrast', orient=HORIZONTAL,
                     from_=0, to=127, #sliderlength=width,
                     variable=self.contrast, command=self.update)
        sbri = Scale(f2, label='brightness', orient=HORIZONTAL,
                     from_=0, to=255, # resolution=0.01, sliderlength=width,
                     variable=self.brightness, command=self.update)
        scon.pack(side='left', padx=5)
        sbri.pack(side='right', padx=5)
        f2.pack(side='bottom', fill='x', expand=0)

    def revOrder(self):
        self.imagelist.reverse()
        self.last = len(self.imagelist)-self.last

    def revRedisplay(self):
        self.revOrder()
        self.displayParms_2()

    def shortcuts(self):
        sc = Toplevel(self.top)
        sc.title("Shortcuts")

        scmain = Frame(sc, borderwidth=2, relief=RIDGE)
        Label(scmain,  text='Keyboard shortcuts:').grid(row=0, sticky=W)

        Label(scmain, text='CTRL-a'     ).grid(row=1, column=0, sticky=W)
        Label(scmain, text='Select all' ).grid(row=1, column=1, sticky=W)
        Label(scmain, text='CTRL-v'     ).grid(row=2, column=0, sticky=W)
        Label(scmain, text='Invert'     ).grid(row=2, column=1, sticky=W)
        Label(scmain, text='CTRL-s'     ).grid(row=3, column=0, sticky=W)
        Label(scmain, text='Save'       ).grid(row=3, column=1, sticky=W)
        Label(scmain, text='CTRL-r'     ).grid(row=4, column=0, sticky=W)
        Label(scmain, text='Read'       ).grid(row=4, column=1, sticky=W)
        Label(scmain, text='CTRL-w'     ).grid(row=5, column=0, sticky=W)
        Label(scmain, text='Close'      ).grid(row=5, column=1, sticky=W)
        Label(scmain, text='Delete'     ).grid(row=6, column=0, sticky=W)
        Label(scmain, text='Clear All'  ).grid(row=6, column=1, sticky=W)
        Label(scmain, text='CTRL-f'     ).grid(row=6, column=0, sticky=W)
        Label(scmain, text='Show labels').grid(row=6, column=1, sticky=W)
        scmain.pack(padx=5, pady=5)

        okframe = Frame(sc)
        Button(okframe, text="OK", command=sc.destroy).pack()  # padx=5, pady=5)
        okframe.pack()

    def montagesize(self):
        " returns number of columns to use "
        nimgs = len(self.imagelist)
        labelheight = 0
        if self.partLabels != 0:
            labelheight = 20 # pixels
        aspect = self.xsize / float(self.ysize + labelheight)
        # approximate a 1.5:1 wd:ht
        nrows = sqrt( 2*nimgs*aspect / 3.0)
        f = nimgs / nrows
        i = int(f)
        if f%i > 0.5:
            i += 1
        ncols = i
        if ncols <  7: ncols =  7  # min num.columns
	if ncols > self.maxcols: ncols = self.maxcols  # max num.columns
        return ncols
    
    def createMontage(self):
        if hasattr(self,'sf'): self.sf.destroy()
        
        self.sf = Pmw.ScrolledFrame(self.top)
        self.fr = self.sf.interior()
        self.display(self.fr)
        self.sf.pack(fill=BOTH, expand=1)
        
        ht = self.fr.winfo_reqheight()
        wd = self.fr.winfo_reqwidth()
        if ht == 1 or wd == 1:
            self.fr.after(500, self.frsize)
        else:
            if ht < self.xmax and wd < self.ymax:
                self.sf.component("clipper").configure(height=ht)
                self.sf.component("clipper").configure(width=wd)
        self.readSelections()

    def frsize(self):
        " try to resize the window after the mainloop starts "
        ht = self.fr.winfo_reqheight()
        wd = self.fr.winfo_reqwidth()
        if ht == 1 or wd == 1:
            self.fr.after(500, self.frsize)
        else:
            if ht < self.xmax and wd < self.ymax:
                self.sf.component("clipper").configure(height=ht)
                self.sf.component("clipper").configure(width=wd)
        
    def checksize(self):
        size = self.sizeVar.get()
        xhalf = int(self.xsize/2)
        xtwice = 2*self.xsize
        yhalf = int(self.ysize/2)
        ytwice = 2*self.ysize
        if size != 1:
            if size == 0:
                x,y = xhalf, yhalf
            elif size == 2:
                x,y = xtwice,ytwice
        else:
            x,y = self.xsize, self.ysize

        return x,y

    def loadImagelist(self):
        # load image list
        partdoc = self.partdoc
        filerange = self.filerange
        if len(partdoc) < 1:
            return
        if not os.path.exists(partdoc):
            print "loadImagelist: Unable to find " + partdoc
            showerror("Error!", "Unable to find " + partdoc)
            self.error = True
            return
        D = Spiderutils.readSpiderDocFile(partdoc)
        self.imagelist = []
        if filerange == None:
            filerange = D.keys()
#        print 'filerange', filerange
        for f in filerange:
            num = int(D[f][0])  # each element is a list (of column data)
##            imgfilename = Spiderutils.template2filename(self.parttemplate, num)
##            imagelist.append(imgfilename)
            self.imagelist.append(num)
        if len(self.imagelist) < 1:
            print "no images loaded"
            master.destroy()
            return
        else: self.error = False
            

    def display(self, parent):
        x,y = self.checksize()

        useLabels = self.partLabels
        self.particle2label = {}
        self.photolist = []
        i = 0
        j = 0
        key = 0

        pal = self.makepalette()

        if self.isStack: im = Image.open(self.parttemplate)

        for num in self.imagelist:
            key += 1  # increment counter

            if self.isStack:
                index = num - 1
                im.seek(index)
            else:
                filename = Spiderutils.template2filename(self.parttemplate, num)
                im = Image.open(filename)

            img = im.convert2byte()
            ix,iy = img.size
            img.convert("L")

            photo = ImageTk.PhotoImage(img, palette=256,
                                       master=parent)

            if ix != x or iy != y:
                icpy = img.copy()  # o.w. orig is resized
                img = icpy.resize( (x,y) )

            img.putpalette(pal)
            photo = ImageTk.PhotoImage(img, palette=256, master=parent)
#            photo.paste(img)

            if hasattr(im,'selectvalue'):
                sc = self.selectClasses[im.selectvalue]
                ib = Label(parent, image=photo, bd=self.bd, borderwidth=3, 
                    background=sc.color)
            else:
                ib = Label(parent, image=photo, bd=self.bd, borderwidth=3)
            ib.photo = photo

            if self.isStack:
                ib.filename = 'slice ' + str(num)  # filename
            else:
                ib.filename = filename

            ib.grid(row=j, column=i, padx=2, pady=2)
            ib.key = key
            ib.bind('<Button-1>', lambda event, w=ib, i=im: self.select(w,i))
            ib.bind('<Control-Button-1>', lambda event, 
                w=ib, i=im: self.shiftSelect(w,i))
            ib.bind('<Shift-Button-1>', lambda event, 
                w=ib, i=im: self.shiftSelect(w,i))
            ib.bind('<Shift-Enter>', lambda event, w=ib, i=im: self.select(w,i))
            self.photolist.append(ib)

            # particle-to-label lookup table
            filenumber = int(Spiderutils.getfilenumber(ib.filename))
            self.particle2label[filenumber] = ib

            if useLabels:
                lb = Label(parent, text=os.path.basename(ib.filename))
                lb.grid(row=j+1, column=i)
            i += 1
            if i > self.ncol-1:
                i = 0
                j += 1
                if useLabels: j += 1 # incr j twice for labels

    def update(self, value):
        p = self.makepalette()
        x,y = self.checksize()

        if self.isStack: im = Image.open(self.parttemplate)

        for counter in range(len(self.photolist)):
            if self.isStack:
                num = self.imagelist[counter]
                index = num - 1
                im.seek(index)
                img = im.convert2byte()
            else:
                filename = Spiderutils.template2filename(self.parttemplate,
                                                self.imagelist[counter])
                img = Image.open(filename).convert2byte()

            # check size
            ix,iy = img.size
            if ix != x or iy != y:
                icpy = img.copy()  # o.w. orig is resized
                img = icpy.resize( (x,y) )

            img.putpalette(p)
            self.photolist[counter].photo.paste(img)

    def makepalette(self):
        " make a 256 palette based on contrast, brightness "
        c = self.contrast.get()
        b = self.brightness.get()
        x1 = b - (self.maxcontrast - c)
        if x1 < 0: x1 = 0
        elif x1 > 255: x1 = 255
            
        x2 = b + (self.maxcontrast - c)
        if x2 < 0: x2 = 0
        elif x2 > 255: x2 = 255
        #print "%f %f %f %f" % (c,b,x1,x2)

        p = range(256)

        if x1 > -1:
            x = int(x1)
            for i in range(0, x): p[i] = 0

        if x2 < 256:
            x = int(x2)
            for i in range(x, 256): p[i] = 255

        if x1 != x2:      
            m = 256 / (x2-x1)
            b = -m * x1
            xx1 = int(x1)
            xx2 = int(x2)
            for i in range(xx1,xx2+1):
                v = int(m*i + b)
                if v < 0: v = 0
                elif v > 255: v = 255
                p[i] = v

        pal = []
        for i in p:
            pal.append(p[i])
            pal.append(p[i])
            pal.append(p[i])
        return pal

    def readSelections(self, event=None):
        if os.path.exists(self.savefilename):
            goodpartdoc = Spiderutils.readSpiderDocFile(self.savefilename)
            partkeys = goodpartdoc.keys()
            for key in partkeys:
                filenumber = int(goodpartdoc[key][0])
                partclass = int(goodpartdoc[key][1])
                sc = self.selectClasses[str(partclass)]
                color = sc.color

                # map particle to label and photo.configure and photo.selectvalue
                photo = self.particle2label[filenumber]
                photo.configure(background=sc.color)
                photo.selectvalue = sc.key

	    # also color deselected particles
	    badcolor = self.deselectedColor.get()
	    sc = self.selectClasses['0']
	    for photo in self.photolist:
	        if not hasattr(photo, 'selectvalue'): 
		    photo.configure(background=badcolor)
		    photo.selectvalue = sc.key
#		filenumber = Spiderutils.getfilenumber(photo.filename)
#		print "readSelections, selectvalue:", photo.selectvalue, ", filenumber:", filenumber

            self.numgood = len(partkeys)
            print 'Read', self.numgood, 'keys from', os.path.basename(self.savefilename)

#        else:
#            print self.savefilename, 'does not exist'

    def invertSelect(self, event=None):
        goodcolor = self.selectedColor.get()
        sc = self.selectClasses[goodcolor]
        badcolor = self.deselectedColor.get()
        newcolor = sc.color
        for photo in self.photolist:
            currentcolor = photo.cget('background')
            if currentcolor == newcolor:              # then deselect
                sc = self.selectClasses['0']
                photo.configure(background=badcolor)
            else:                              # o.w. select with new color
                sc = self.selectClasses[goodcolor]
                photo.configure(background=sc.color)
            photo.selectvalue = sc.key

    def shiftSelect(self, widget, image):
        # widgets are numbered from 0..n-1

        # check if selecting forward or backward
        if widget.key > self.last:
            start = self.last
            end = widget.key
        else:
            start = widget.key-1
            end = self.last

        key = self.selectedColor.get()
        sc = self.selectClasses[key]

        # loop through keys in range
        for counter in range(start, end):
            photo = self.photolist[counter]
            photo.selectvalue = sc.key
            photo.configure(background=sc.color)

        self.last = widget.key

    def select(self, widget, image):
        goodcolor = self.selectedColor.get()
        sc = self.selectClasses[goodcolor]
        badcolor = self.deselectedColor.get()
        newcolor = sc.color
        currentcolor = widget.cget('background')
        #print "current color, new color: %s, %s" % (str(currentcolor),str(newcolor))
        if currentcolor == newcolor:              # then deselect
            sc = self.selectClasses['0']
            widget.configure(background=badcolor)
            widget.selectvalue = sc.key
            image.selectvalue = sc.key
        else:                              # o.w. select with new color
            widget.configure(background=sc.color)
        widget.selectvalue = sc.key
        image.selectvalue = sc.key

        # update last
        self.last = widget.key

    def selectAll(self, event=None):
        goodcolor = self.selectedColor.get()
        sc = self.selectClasses[goodcolor]

        for photo in self.photolist:
            sc = self.selectClasses[goodcolor]
            photo.configure(background=sc.color)
            photo.selectvalue = sc.key

    def clearAll(self, event=None):
        badcolor = self.deselectedColor.get()
        sc = self.selectClasses['0']

        for photo in self.photolist:
            photo.configure(background=badcolor)
            photo.selectvalue = sc.key

    def saveSelections(self, event=None):
        """
    savasfilename can't append to a file - if the file exists an error message pops up
    if output specified, don't ask (put it in window in menubar?)
        """
        if self.savefilename == "":
            self.savefilename = asksaveasfilename(title="Save points into a doc file")
        filename = self.savefilename
        if len(filename) == 0:
            return

        # construct a dictionary to pass to writedocfile
        headers = ['file_number', '    class ']
        F = {}
        key = 0
        for photo in self.photolist:
            filenumber = Spiderutils.getfilenumber(photo.filename)
            try:
                int(filenumber)
                if hasattr(photo, 'selectvalue'):
                    val = photo.selectvalue
                    if int(val) != 0:
                        key = key + 1
                        F[key] = [filenumber, val]
            except:
                print "error getting file number from %s" % photo.filename
                continue

        if key == 0:
            print "No images selected"
            return
        else:
            backup(filename)
            if Spiderutils.writeSpiderDocFile(filename,F, headers=headers, append=1):
                self.numgood = key
                print 'Wrote', self.numgood, 'keys to %s' % os.path.basename(filename)
            else:
                showerror("Error!", "Unable to write to %s" % os.path.basename(filename))

    def saveClose(self, event=None):
        self.saveSelections()
        self.closeWindow()

    def displayFilenames(self, event=None):
        q = self.showVar.get()
#	print "displayFilenames:", self.showVar.get(), self.partLabels
	self.showVar.set(1-q)
#	print "displayFilenames:", self.showVar.get(), self.partLabels
#	print
	self.displayFilenames_1()

    def displayFilenames_1(self):
        q = self.showVar.get()
#	print "displayFilenames:", q, self.partLabels
        if q != self.partLabels:
            self.partLabels = q  
            self.createMontage()

    def displayParms_2(self):
        self.createMontage()
        
    def displayParms(self):
        self.ncolVar.set(str(self.ncol))
        w = Toplevel(self.top)
        w.title("Montage")
        self.parmWindow(w)
        self.top.wait_window(w) # wait for window to be destroyed

        changed = False  # default
	
	# numcol
	try:
            i = int(self.ncolVar.get())
	    j = int(self.maxcolVar.get())
            if i != self.ncol: changed = True
            if j != self.maxcols: changed = True
#                self.ncol = i
#                self.createMontage()
	except:
            print "no. of columns must be an integer"
	    self.maxcolVar.set(self.maxcols)  # I don't know why I need this
            return

        if changed :
	    self.maxcols = j
	    if self.ncol > self.maxcols: 
		self.ncol = self.maxcols
	    else: # if ncol less than maxcols
		if i != self.ncol: 
		    self.ncol = i
		else:
		    self.ncol = self.montagesize()
	    print "displayParms", self.ncol, self.maxcols
	    self.createMontage()

    def parmWindow(self,win):
#	print "maxcolVar",self.maxcolVar.get(),type(self.maxcolVar.get())
        win.bind('<Return>', lambda e, w=win: self.parmquit(w))

	fe = Frame(win)
        Label(fe,text="Number of columns:").grid(row=0, column=0, padx=2, pady=2)
#        Label(fe,text="Number of columns:").pack(side='left', padx=2, pady=2)

        # number of columns
	e = Entry(fe, textvariable=self.ncolVar, width=6)
	e.focus_set()
        e.grid(row=0, column=1, pady=2, padx=2)
#        e.pack(side='left', fill='x', expand=1, pady=2, padx=2)
	
        # number of columns
        Label(fe,text="Maximum columns:").grid(row=1, column=0, padx=2, pady=2)
	maxcolEntry = Entry(fe, textvariable=self.maxcolVar, width=6)
        maxcolEntry.grid(row=1, column=1, pady=2, padx=2)
        fe.pack(side='top',fill='both', expand=1)

       # bottom
        fbut = Frame(win, borderwidth=2) #relief='raised',
        Button(fbut, text='Ok', command=win.destroy).pack(padx=2,pady=2)
        fbut.pack(side='bottom', fill='x', expand=1)

    def parmquit(self,win):
        win.destroy()

# -------------------------------------------------------------------
# classAverages:  inherits from support.montage
#                 overrides makeMenus() and display()
# -------------------------------------------------------------------

class classAverages(montage.montage): 
    
    def __init__(self, master, prefs, 
                 title="verifybyview.py", ncol=None, useLabels=0, reverse=0):

        self.extension = prefs.extension
        if self.extension != "":
            ext = self.extension
            self.filetypes = [ ("", "*"+ext), ("All files", "*") ]
        else:
            self.filetypes = [ ("All files", "*") ]

        self.prefs = prefs
        self.master = master
        self.dir = prefs.dir1
        self.dir_var = StringVar()
        self.dir_var.set(self.dir)
        self.classlist = prefs.class_list_base + self.extension
        self.classlist_var = StringVar()
        self.classlist_var.set(self.classlist)
        self.classavgtemplate = prefs.class_avg_template + self.extension
        self.classavg_var = StringVar()
        self.classavg_var.set(self.classavgtemplate)

        # initialize class-average list, read SPIDER document, and get keys
        imagelist = []
        spi_class_list = Spiderutils.readSpiderDocFile(os.path.join(self.dir, \
                                                                    self.classlist))
        self.class_keys = spi_class_list.keys()

        # loop through classes and substitute in class#
	imagelist.append('-a')
        for class_key in self.class_keys:
            class_num = int(spi_class_list[class_key][0])
            class_avg_name = Spiderutils.template2filename(os.path.join(self.dir, \
                                                           self.classavgtemplate), n=class_num)
            imagelist.append(class_avg_name)
#	    print 'class_avg_name', class_avg_name

        self.reverse = self.prefs.reverse
        if self.reverse == 1: self.revOrder()
        ncol = self.prefs.ncol

        if len(imagelist) < 1: return
        montage.montage.__init__(self, master,
                                imagelist=imagelist,
                                title=title,
                                ncol=ncol,
                                useLabels=self.prefs.caLabels)

        self.top.protocol("WM_DELETE_WINDOW", master.quit)
        
	# set up some templates
        self.partdocfile = ""  
        self.partdoctemplate = StringVar()
        self.partdoctemplate.set(self.partdocfile)   # display only the basename
        self.serfile = ""  
        self.sertemplate = StringVar()
        self.sertemplate.set(self.serfile)        # particle images
        self.goodclasslist = ""  
        self.goodclass_var = StringVar()
        self.goodclass_var.set(self.goodclasslist)
        self.goodparttemplate = ""  
        self.goodpart_var = StringVar()
        self.goodpart_var.set(self.goodparttemplate)

        # a selection class = ( 'value', 'key', 'color', 'label (highlight)')
        self.selectClasses = {}
        self.selectClasses['0'] = selectionClass(0,'0','red'  , 'Red',    'pink')
        self.selectClasses['1'] = selectionClass(1,'1','green', 'Green',  'light green')
        self.selectClasses['2'] = selectionClass(2,'2','black', 'Black',  'gray50')
        self.selectClasses['3'] = selectionClass(3,'3','blue',  'Blue',   'light blue')
        self.selectClasses['4'] = selectionClass(4,'4','gray75','Lt gray','gray50')
        self.goodColor = StringVar()
        self.goodColor.set('1')
        self.badColor = StringVar()
        self.badColor.set('0')
#        self.selectedColor.trace_variable('w', self.selectcallback)


        # CHECK
        self.montageColor = prefs.montageColor          # '1'
        self.badMontageColor = prefs.badMontageColor    # 'Red'

        self.contrast = 0
        self.brightness = 127

        self.avgframe = {}

        self.makeMenus()
        
        self.prefs.caReduction = prefs.caReduction  # 1
        self.sizeVar.set(prefs.caReduction)
        self.partLabels = 0
        
        self.createMontage()

        self.top.bind('<Alt-u>' , self.updateDir)
        self.top.bind('<Return>', self.updateDir)
        self.top.bind('<Alt-s>' , self.saveSelections)
        self.top.bind('<Alt-r>' , self.readSelections)
        self.top.bind('<Alt-p>' , self.prefs.printSettings)
        self.top.bind('<Alt-l>' , self.locateSelections)
        self.top.bind('<Alt-t>' , self.test)

    def test(self, event=None):
#	print "len(imagelist[0])",len(self.imagelist)
#	print "self.dir",self.dir
	print "self.classavgtemplate",self.classavgtemplate
#	print "self.imagelist[0]", self.imagelist[0]
	print "self.imagelist[0].info['filename']", self.imagelist[0].info['filename']
	print "self.imagelist[0].info", self.imagelist[0].info
	print "self.avg_lut", self.avg_lut

    def selectcallback(self, name, index, mode):
        "called whenever self.selectedColor changes "
        key = self.selectedColor.get()
        print 'New color:', self.selectClasses[key].color

    # ------- override the menu bar functions -------
    def makeMenus(self):
        self.mBar = Frame(self.top, relief='raised', borderwidth=1)
        self.mBar.pack(side='top', fill = 'x')
        self.balloon = Pmw.Balloon(self.top)
        
        # Make the File menu
        Filebtn = Menubutton(self.mBar, text='File', underline=0,
                                 relief='flat')
        Filebtn.pack(side=LEFT, padx=5, pady=5)
        Filebtn.menu = Menu(Filebtn, tearoff=0)
        Filebtn.menu.add_command(label='Templates',
                                 command=self.callTemplates)
        Filebtn.menu.add_command(label='Read good-class list',
                                 command=self.readSelections)
        Filebtn.menu.add_command(label='Save selection', underline=0, 
                                 command=self.saveSelections)
        Filebtn.menu.add_command(label='Locate good-particle lists', # underline=0, 
                                 command=self.locateSelections)
        Filebtn.menu.add_command(label='Close all windows',
                                 command=self.closeWindows)
        Filebtn.menu.add_command(label='Quit',
                                 command=self.master.quit)
        Filebtn['menu'] = Filebtn.menu
    
        # Make the Display menu (same as montage.py)
        Dspbtn = Menubutton(self.mBar, text='Display', relief='flat', underline=0)
        Dspbtn.pack(side=LEFT, padx=5, pady=5)
        Dspbtn.menu = Menu(Dspbtn, tearoff=0)

        Dspbtn.menu.add_command(label='no. columns', underline=0,
                                    command=self.displayParms_save)
#        Dspbtn.menu.add_checkbutton(label='show filenames', underline=0,
#                                    variable=self.showVar,
#                                    command=self.displayFilenames)
        Dspbtn.menu.add_command(label='labels', underline=0,
                                    command=self.displayLabels)
        

        # 'size' has a submenu of checkbuttons
        Dspbtn.menu.sizes = Menu(Dspbtn.menu)
        Dspbtn.menu.sizes.add_radiobutton(label='1/2', underline=0,
                                          variable=self.sizeVar, value=0,
                                          command=self.displayParms_2save)
        Dspbtn.menu.sizes.add_radiobutton(label='1x', underline=0,
                                          variable=self.sizeVar, value=1,
                                          command=self.displayParms_2save)
        Dspbtn.menu.sizes.add_radiobutton(label='2x', underline=0,
                                          variable=self.sizeVar, value=2,
                                          command=self.displayParms_2save)
        Dspbtn.menu.add_cascade(label='image size', menu=Dspbtn.menu.sizes)

        Dspbtn.menu.add_command(label='reverse order', command=self.revRedisplay)

        # Make the Color submenu
        Dspbtn.menu.colors = Menu(Dspbtn.menu, tearoff=0)
        keys = self.selectClasses.keys()
        keys.sort()
        Dspbtn.menu.colors.good = Menu(Dspbtn.menu.colors, tearoff=0)
        Dspbtn.menu.colors.bad =  Menu(Dspbtn.menu.colors, tearoff=0)
        for key in keys:
            sc = self.selectClasses[key]
            Dspbtn.menu.colors.good.add_radiobutton(label = sc.label, underline=0,
                value = sc.value, foreground = 'white', background = sc.color,
                activebackground= sc.active, activeforeground= 'black',
                variable=self.goodColor)
            Dspbtn.menu.colors.bad.add_radiobutton(label = sc.label, underline=0,
                value = sc.value, foreground = 'white', background = sc.color,
                activebackground= sc.active, activeforeground= 'black',
                variable=self.badColor)
        Dspbtn.menu.colors.add_cascade(label='selected, good', menu=Dspbtn.menu.colors.good)
        Dspbtn.menu.colors.add_cascade(label='selected, bad',  menu=Dspbtn.menu.colors.bad)
        Dspbtn.menu.add_cascade(label='colors', menu=Dspbtn.menu.colors)

        Dspbtn['menu'] = Dspbtn.menu

        Helpbtn = Menubutton(self.mBar, text='Help', underline=0, relief='flat')
        Helpbtn.pack(side=RIGHT, padx=5, pady=5)
        Helpbtn.menu = Menu(Helpbtn, tearoff=0)
        Helpbtn.menu.add_command(label='Keyboard shortcuts', underline=0, 
                                 command=self.shortcuts)
        Helpbtn['menu'] = Helpbtn.menu

        # add second row for directory field
        self.bBar = Frame(self.top, relief='raised', borderwidth=1)
        self.bBar.pack(side='top', fill='x')

        # add field for directory
        dir_entry = Entry(self.bBar, textvariable=self.dir_var, width=48)
        dir_entry.pack(side='left',padx=2, pady=5)

        dir_button = Button(self.bBar, text='Save+Update', 
            command=self.saveUpdate)
        dir_button.pack(side='left',padx=2, pady=5)
        dir_button = Button(self.bBar, text='Update', underline=0, 
            command=self.updateDir)
        dir_button.pack(side='left',padx=2, pady=5)

    def displayLabels(self):
        beforeText = self.prefs.textLabel
        beforeLabel = self.prefs.caLabels
        labelWindow = Toplevel(self.top)
        self.getLabels(self.top,labelWindow)
        self.top.wait_window(labelWindow)
        
        if self.labelVar.get() != beforeText or self.showVar.get() != beforeLabel:
            print 'label', beforeText, self.showVar.get(), self.labelVar.get()
            self.prefs.textLabel = self.labelVar.get()
            self.prefs.caLabels = self.showVar.get()
            self.createMontage()
            self.prefs.saveSettings()

    def getLabels(self,parent,win):
        labelFrame = Frame(win)
        labelButton = Checkbutton(labelFrame, text='show labels', variable=self.showVar)
        labelButton.grid(row=0, column=0)
        Label(labelFrame, text='text label:').grid(row=1, column=0)
        self.labelVar = StringVar()
        self.labelVar.set(self.prefs.textLabel)
        labelEntry = Entry(labelFrame, textvariable=self.labelVar)
        labelEntry.grid(row=1, column=1)
        labelFrame.pack()
        quitFrame = Frame(win)
        Button(quitFrame, text='OK', command=win.destroy).pack()
        quitFrame.pack()
        win.bind('<Return>', lambda p=parent, w=win: clickOK(p,w)) 

    def labelQuit(self,parent,win):
        win.destroy()

    def displayParms_save(self):
        self.displayParms()
#        print 'ncolVar', self.ncolVar.get()
        self.prefs.ncol = self.ncolVar.get()
        self.prefs.saveSettings()

    def displayParms_2save(self):
        self.prefs.caReduction = self.sizeVar.get()
        self.prefs.saveSettings()
        self.displayParms_2()
        
    def revOrder(self):
        self.class_keys.reverse()

    def revRedisplay(self):
        self.revOrder()
        self.prefs.reverse = self.reverse = 1 - self.reverse
        self.displayParms_2()
        self.prefs.saveSettings()

    def saveUpdate(self):
        self.saveSelections()
        self.updateDir()

    def shortcuts(self):
        sc = Toplevel(self.top)
        sc.title("Shortcuts")

        scmain = Frame(sc, borderwidth=2, relief=RIDGE)
        Label(scmain, text='Shortcuts:').grid(row=0, sticky=W)

        Label(scmain, text='ALT-u').grid(row=1, column=0, sticky=W)
        Label(scmain, text='Update directory'          ).grid(row=1, column=1, sticky=W)
        Label(scmain, text='ALT-r'                     ).grid(row=2, column=0, sticky=W)
        Label(scmain, text='Read good classes'         ).grid(row=2, column=1, sticky=W)
        Label(scmain, text='ALT-s'                     ).grid(row=3, column=0, sticky=W)
        Label(scmain, text='Save good classes'         ).grid(row=3, column=1, sticky=W)
        Label(scmain, text='ALT-l'                     ).grid(row=4, column=0, sticky=W)
        Label(scmain, text='Locate good-particle lists').grid(row=4, column=1, sticky=W)
        Label(scmain, text='ALT-F4'                    ).grid(row=5, column=0, sticky=W)
        Label(scmain, text='Close all'                 ).grid(row=5, column=1, sticky=W)
        scmain.pack(padx=5, pady=5)

        okframe = Frame(sc)
        Button(okframe, text="OK", command=sc.destroy).pack()
        okframe.pack()

    def updateDir(self, event=None):
        self.dir = self.prefs.dir1 = self.dir_var.get()
        print 'New working directory: ', os.path.basename(self.dir)

        # initialize class-average list, read SPIDER document, and get keys
        imagelist = []
        spi_class_list = Spiderutils.readSpiderDocFile(os.path.join(self.dir,self.classlist))
        self.class_keys = spi_class_list.keys()

        # loop through classes and substitute in class#
        for class_key in self.class_keys:
            class_num = int(spi_class_list[class_key][0])
            classavgfullpath = os.path.join(self.dir,self.classavgtemplate)
            class_avg_name = Spiderutils.template2filename(classavgfullpath, class_num)
            imagelist.append(class_avg_name)

        if self.reverse == 1: self.revOrder()

        # create montage of new class-averages
        self.imagelist = loadImageSeries(imagelist)
#        self.ncol = self.montagesize()
        self.createMontage()
        self.prefs.saveSettings()

    def closeWindows(self):
        winlist = self.top.winfo_children()
        for win in winlist:
            if win.winfo_exists():
                if win.winfo_class() == 'Toplevel':
                    win.destroy()        

    def hello(self, im, frame):
        " filename is the classavg image "

        frame.configure(background='yellow')

        if self.partdocfile == "" or self.serfile == "":
            self.callTemplates()

        filename = im.info['filename']
#	print 'filename',filename
#	print "self.dir",self.dir
        filenum = int(Spiderutils.getfilenumber(filename))

        # take filenumber from 'filename' and substitute it into 'tmp' template
        tmp = os.path.join(self.dir,self.partdoctemplate.get())
        fn = Spiderutils.template2filename(tmp, os.path.basename(filename))

        classmontage_win = Toplevel(self.top)
        xcoord = self.prefs.xcoord
        ycoord = self.prefs.ycoord
        classmontage_win.geometry('+%d+%d' % (xcoord, ycoord) )
	filenumber = Spiderutils.getfilenumber(filename)

        goodpartdoc = os.path.join(self.dir, 
            Spiderutils.template2filename(self.goodpart_var.get(),filenumber))

        # if particle directory not listed explicitly, assume same as other files
        partdir = os.path.dirname(self.serfile)
	partbase = os.path.basename(self.serfile)
        if partdir == '':
#            print 'Particles assumed to be in', self.dir
            parttemplate = os.path.join(self.dir,self.serfile)  # + self.extension
        else:
#            print 'Particle directory explicitly named', partdir + '/'
            # if filename (minus directory) has no asterisks, try directory
            if partbase.count("*") == 0:
	        parttemplate = Spiderutils.template2filename(self.serfile,filenumber)
	    else:
		parttemplate = self.serfile
#        print 'PARTTEMPLATE', partdir, partbase.count("*"), partbase, parttemplate
#	print 'GOODPART', self.goodpart_var.get(), goodpartdoc

        title = "images for class: " + os.path.basename(filename)
        m = pmontage(classmontage_win, title=title,
            partdoc=fn, filerange=None, parttemplate=parttemplate, 
            useLabels=self.partLabels, maxcols=self.prefs.maxPartCols, savefilename=goodpartdoc, 
            reverse=self.reverse, sizeVar=self.prefs.montageReduction, 
            selectedColor=self.prefs.montageColor, deselectedColor=self.prefs.badMontageColor, 
            contrast=self.contrast, brightness=self.brightness, 
            xcoord=self.prefs.xcoord, ycoord=self.prefs.ycoord)

        # look for error particle montage
        if m.error:
#            print "function hello: Error"
            # first check path of particles
            absdir = os.path.dirname(os.path.abspath(self.sertemplate.get()))
            if os.path.exists(absdir) == False:
                showerror("Error!", "Unable to find directory %s" % absdir)
            classmontage_win.destroy()
            return

        # draw montage window
#        m.makeMenus()
#        m.createMontage()
#        m.readSelections()
        self.top.wait_window(classmontage_win)

        # remember stuff from particle-montage window
        self.prefs.xcoord = m.xcoord
        self.prefs.ycoord = m.ycoord
        self.prefs.montageReduction = m.sizeVar.get()
        self.prefs.montageColor = m.selectedColor.get()
        self.prefs.badMontageColor = m.deselectedColor.get()
        self.contrast = m.contrast.get()
        self.brightness = m.brightness.get()
        self.prefs.reverse = self.reverse = m.reverse
#        print 'particle labels:', m.partLabels
        self.partLabels = m.partLabels
	self.prefs.maxPartCols = m.maxcols
#	print "self.maxPartCols:",self.prefs.maxPartCols

        self.prefs.saveSettings()

        if m.numgood > 0:
#            print "filenumber", filenumber
#	    print "self.avg_lut[filenum]", self.avg_lut[filenumber]
#	    print "self.imagelist[self.avg_lut[filenum]]", self.imagelist[self.avg_lut[filenumber]]
	    self.imagelist[self.avg_lut[filenum]].selectvalue = self.selectClasses[self.goodColor.get()].value  # WAS 1
#            im.selectvalue = 1  # visited w/good particles
            print "Class", filenum, "selected, category =", self.selectClasses[self.goodColor.get()].value
            frame.configure(background=self.selectClasses[self.goodColor.get()].color)
        else:
            self.imagelist[self.avg_lut[filenum]].selectvalue = 0
#            im.selectvalue = 0  # visited, no good particles
            print "No particles from class", filenum, "selected"
            frame.configure(background=self.selectClasses[self.badColor.get()].color)

    def saveSelections(self, event=None):
#        print 'saveSelections:'

        # construct a dictionary to pass to writedocfile
        headers = ['class_num', '    class ']
        F = {}
        key = 0

        # loop though imagelist
        for classavg in self.imagelist:
            filenumber = Spiderutils.getfilenumber(classavg.info['filename'])
            try:
                int(filenumber)
                if hasattr(classavg, 'selectvalue'):
                    val = classavg.selectvalue
                    if int(val) != 0:
                        key = key + 1
                        F[key] = [filenumber, val]
#                    print 'filenumber:', filenumber, key, F[key]
#                else:
#                    print 'filenumber:', filenumber, 'None'
            except:
                print "error getting file number from %s" % classavg.filename
                continue
#        print 'dictionary:', F

        filename = os.path.join(self.dir,self.goodclass_var.get())
        backup(filename)
        if Spiderutils.writeSpiderDocFile(filename,F, headers=headers, append=1):
            shortdir = os.path.basename(self.dir)
            shortfile = os.path.join(shortdir, self.goodclass_var.get())
            print 'Wrote', key, 'classes to %s' % shortfile
        else:
            showerror("Error!", "Unable to write to %s" % os.path.basename(filename))

    def locateSelections(self, event=None):
        if self.goodparttemplate == "":
            self.callTemplates()
            
        spi_class_list = Spiderutils.readSpiderDocFile(os.path.join(self.dir, self.classlist))
        found_counter = 0
        
        # loop through classes
        for class_key in self.class_keys:
            # get class#
            class_num = int(spi_class_list[class_key][0])
            
            # substitute class# into BYHAND filename template
            goodpartdocbase = Spiderutils.template2filename(self.goodpart_var.get(), class_num)
            goodpartdoc = os.path.join(self.dir, goodpartdocbase)
            
            # check if BYHAND file exists
            if os.path.exists(goodpartdoc) : 
                found_counter += 1
                
                # if it exists, select class
                print "Class", class_num, "selected, category =", self.selectClasses[self.goodColor.get()].value
                
                # update class selection
                class_value = self.selectClasses[self.goodColor.get()].value 
                self.imagelist[self.avg_lut[class_num]].selectvalue =  class_value
                
                # color average
                color = self.selectClasses[str(class_value)].color
                self.avgframe[class_num].configure(background=color)
#            else:
#                print "Class", class_num, "not selected"

        print "Located ", found_counter, "classes"

    def readSelections(self, event=None):
        if self.goodclasslist == "":
            self.callTemplates()

        goodclasspath = os.path.join(self.dir, self.goodclasslist)
        
        if os.path.exists(goodclasspath):
            goodclassdoc = Spiderutils.readSpiderDocFile(goodclasspath)
            try: 
                classkeys = goodclassdoc.keys()
#		goodselectclass = self.selectClasses[self.goodColor.get()]
                for key in classkeys:
                    class_num = int(goodclassdoc[key][0])
		    class_value = int(goodclassdoc[key][1])
#                    color = self.selectClasses[self.goodColor.get()].color
                    color = self.selectClasses[str(class_value)].color
#                    color = background=self.selectClasses[self.goodColor.get()].color
#                    color = goodselectclass.color
                    self.avgframe[class_num].configure(background=color)
                    self.imagelist[self.avg_lut[class_num]].selectvalue =  class_value  # WAS 1
                print len(classkeys), 'classes found'
            except AttributeError:
                print 'AttributeError: 0 keys, maybe?'
        else:
            print self.goodclasslist, 'does not exist'

    def callTemplates(self):
	# call getTemplates
        gtwin = Toplevel(self.top)
        self.getTemplates(gtwin)
        self.top.wait_window(gtwin)  
	# wait until template window gone

        lst = self.classlist_var.get()
        if len(lst) > 0:
            self.classlist = lst = os.path.basename(lst)
            self.prefs.class_list_base = os.path.splitext(lst)[0]

        # if they typed in w/o asterisks..
        avg = self.classavg_var.get()
        if len(avg) > 0:
            self.classavgtemplate = os.path.basename(avg)
            if string.find(avg,"*") < 0:
                avg = os.path.basename(Spiderutils.name2template(avg))
                self.classavg_var.set(avg)
        self.prefs.class_avg_template = os.path.splitext(avg)[0]

        doc = self.partdoctemplate.get()
        if len(doc) > 0:
            self.partdocfile = os.path.basename(doc)
            if string.find(doc,"*") < 0:
                doc = os.path.basename(Spiderutils.name2template(doc))
                self.partdoctemplate.set(doc)
            self.prefs.partdocfile = os.path.splitext(doc)[0]

        ser = self.sertemplate.get()
        if len(ser) > 0:
            self.serfile = ser
	    
#	    # if no asterisk is found
#            if string.find(ser,"*") < 0:
#                ser = Spiderutils.name2template(ser)  # unnecessary?
#                self.sertemplate.set(ser)             # unnecessary?
	    self.prefs.serfile = os.path.splitext(ser)[0]

        gpl = self.goodpart_var.get()
        if len(gpl) > 0:
            self.goodparttemplate = os.path.basename(gpl)
            if string.find(gpl,"*") < 0:
                gpl = Spiderutils.name2template(gpl)
                self.goodpart_var.set(os.path.basename(gpl))
            self.prefs.goodparttemplate = os.path.splitext(gpl)[0]

        gcl = self.goodclass_var.get()
        if len(gcl) > 0:
            self.goodclasslist = gcl = os.path.basename(gcl)
            self.prefs.goodclasslist = os.path.splitext(gcl)[0]

        self.prefs.saveSettings()

    def getTemplates(self, win):
        win.title("Templates")
        extension = self.extension

        if self.partdocfile == "": 
            self.partdoctemplate.set(self.prefs.partdocfile + extension)
        if self.serfile == "": 
            self.sertemplate.set(self.prefs.serfile + extension)
        if self.goodclasslist == "": 
            self.goodclass_var.set(self.prefs.goodclasslist + extension)
        if self.goodparttemplate == "": 
            self.goodpart_var.set(self.prefs.goodparttemplate + extension)

        ftop = Frame(win, padx=5, pady=5)
        ftop.focus_set()
        
        bstat = Button(ftop, text="Class list: ",
                      command = lambda w=win, d='lst': self.setTemplates(w,d))
        estat = Entry(ftop, textvariable = self.classlist_var, width=24)
        bavg = Button(ftop, text="Class average template: ",
                      command = lambda w=win, d='cat': self.setTemplates(w,d))
        eavg = Entry(ftop, textvariable = self.classavg_var, width=24)
        bdoc = Button(ftop, text="Particle doc file template: ",
                      command = lambda w=win, d='doc': self.setTemplates(w,d))
        edoc = Entry(ftop, textvariable = self.partdoctemplate, width=24)
        bser = Button(ftop, text="Particle template: ",
                      command = lambda w=win, d='ser': self.setTemplates(w,d))
        eser = Entry(ftop, textvariable = self.sertemplate, width=24)
        bGC = Button(ftop, text="Good class list: ",
                      command = lambda w=win, d='gcl': self.setTemplates(w,d))
        eGC = Entry(ftop, textvariable = self.goodclass_var, width=24)
        bGP = Button(ftop, text="Good particle list template: ",
                      command = lambda w=win, d='gpl': self.setTemplates(w,d))
        eGP = Entry(ftop, textvariable = self.goodpart_var, width=24)

        # fit stuff onto grid
        Label(ftop, text='INPUTS').grid(row=0, sticky=W)
        bstat.grid(row=1, column=0, sticky='w')
        estat.grid(row=1, column=1, sticky='nsew')
        ftop.columnconfigure(1, weight=1)  # make entry expand
        bavg.grid(row=2, column=0, sticky='w')
        eavg.grid(row=2, column=1, sticky='nsew')
        bdoc.grid(row=3, column=0, sticky='w')
        edoc.grid(row=3, column=1, sticky='nsew')
        bser.grid(row=4, column=0, sticky='w')
        eser.grid(row=4, column=1, sticky='nsew')
        Label(ftop, text='OUTPUTS').grid(row=5, sticky=W)
        bGC.grid(row=6, column=0, sticky='w')
        eGC.grid(row=6, column=1, sticky='nsew')
        bGP.grid(row=7, column=0, sticky='w')
        eGP.grid(row=7, column=1, sticky='nsew')
        ftop.pack(side='top', fill='both', expand=1)
        
        fbut = Frame(win)
        Button(fbut, text="Done", command=win.destroy).pack(padx=5, pady=5)
        fbut.pack()
        win.bind('<Return>', lambda p=self.top, w=win: clickOK(p,w)) 

    def setTemplates(self, parent, which='doc'):
        filename = askopenfilename(parent=parent, filetypes=self.filetypes, initialdir=self.prefs.dir1)
        if len(filename) < 1:
            return

        # particles are the only files in unchanging directory
        if which == 'ser':
            tmp = Spiderutils.name2template(filename)
	    
	    # if filename has asterisks, convert to template
	    if string.find(tmp,"*") > 0:
		self.sertemplate.set(tmp)
	    # if filename doesn't have asterisks, use basename
	    else:
		self.sertemplate.set(os.path.basename(tmp))
		
        elif which == 'lst':  # not a template
            self.classlist_var.set(os.path.basename(filename))
        elif which == 'cat':
            tmp = Spiderutils.name2template(filename)
            self.classavg_var.set(os.path.basename(tmp))
        elif which == 'gcl':  # not a template
            self.goodclass_var.set(os.path.basename(filename))
        elif which == 'gpl':
            tmp = Spiderutils.name2template(filename)
            self.goodpart_var.set(os.path.basename(tmp))
        else:  # file type DOC
            tmp = Spiderutils.name2template(filename)
            self.partdoctemplate.set(os.path.basename(tmp))

    # display here overrides definition in montage class
    def display(self, parent):
        size = self.sizeVar.get()
        xhalf = int(self.xsize/2)
        xtwice = 2*self.xsize
        yhalf = int(self.ysize/2)
        ytwice = 2*self.ysize
        if size != 1:
            if size == 0:
                x,y = xhalf, yhalf
            elif size == 2:
                x,y = xtwice,ytwice
        else:
            x,y = self.xsize, self.ysize
            
        useLabels = self.prefs.caLabels
        i = 0
        j = 0

        self.avg_lut = {}  # initialize class-to-key lookup table
        self.avgframe = {}

        spi_class_list = Spiderutils.readSpiderDocFile(os.path.join(self.dir,self.classlist))

        for class_key in self.class_keys:  # range from 1 to n
            listkey = class_key-1
            im = self.imagelist[listkey]

            ix,iy = im.size
            if ix != x or iy != y:
                icpy = im.copy()  # o.w. orig is resized
                img = icpy.resize( (x,y) )
            else:
                img = im

            photo = ImageTk.PhotoImage(img, palette=256, master=parent)

            class_num = int(spi_class_list[class_key][0])

            # color background
            if hasattr(im,'selectvalue'):
                if im.selectvalue == 1:  # contains good particles
                    f = Frame(parent, padx=2, pady=2, borderwidth=2, 
                        background=self.selectClasses[self.goodColor.get()].color)
                else:  # visited but no good particles
                    f = Frame(parent, padx=2, pady=2, borderwidth=2, 
                        background=self.selectClasses[self.badColor.get()].color)
            else:  # no selectvalue
                f = Frame(parent, padx=2, pady=2, borderwidth=2)

            btn = Button(f,image=photo, command=lambda i=img, f=f: self.hello(i,f))
            btn.bind('<Button-3>', lambda event, i=img, f=f: self.selectWholeClass(i,f))
            btn.pack()

            f.photo = photo
            filename = im.info['filename']
            f.grid(row=j, column=i, padx=2, pady=2)
            if useLabels:
                lb = Label(parent, text=os.path.basename(filename))
                lb.grid(row=j+1, column=i)

#                class_cc = Label(parent, text='CCC=' + str(spi_class_list[class_key][1]))
                class_label_float = float(spi_class_list[class_key][1])

		# check if integer
		if round(class_label_float) == class_label_float: 
		    class_label = str(int(class_label_float))
		else: 
		    class_label = str(class_label_float)
                class_cc = Label(parent, text=self.prefs.textLabel + '=' + class_label)
                class_cc.grid(row=j+2, column=i)
            i += 1
            if i > self.ncol-1:
                i = 0
                j += 1
                if useLabels: j += 2 # incr j twice for labels

            self.avg_lut[class_num]  = listkey  # remember for later
#	    print "display", class_num, listkey
            self.avgframe[class_num] = f

    def selectWholeClass(self, im, frame):
        filename = im.info['filename']

        if self.partdocfile == "" or self.goodparttemplate == "":
            self.callTemplates()

        # load filenumbers
        filename = im.info['filename']
        tmp = os.path.join(self.dir,self.partdoctemplate.get())
        filenum = int(Spiderutils.getfilenumber(filename))
#        print 'filenum:', filenum, type(filenum)
        fn = Spiderutils.template2filename(tmp, filenum)

        if not os.path.exists(fn):
            print "Unable to find " + fn
            return
        F = {}  # initialize output
        D = Spiderutils.readSpiderDocFile(fn)
        files = D.keys()
        for f in files:
            num = int(D[f][0])  # each element is a list (of column data)
            F[f] = [num, 1]

        # write filenumbers
        goodpartdocbase = Spiderutils.template2filename(self.goodpart_var.get(), 
            filenum)
        goodpartdoc = os.path.join(self.dir, goodpartdocbase)
        backup(goodpartdoc)
        headers = ['file_number', '    class ']
        if Spiderutils.writeSpiderDocFile(goodpartdoc, F, headers=headers, append=1):
            print 'Wrote', len(D), "keys to %s" % goodpartdocbase
        else:
            showerror("Error!", "Unable to write to %s" % goodpartdocbase)

        # highlight class average
        self.imagelist[self.avg_lut[filenum]].selectvalue = self.selectClasses[self.goodColor.get()].value  # WAS 1
        print "Class", Spiderutils.getfilenumber(filename), "selected, category =", self.selectClasses[self.goodColor.get()].value
        frame.configure(background=self.selectClasses[self.goodColor.get()].color)

# -------------------------------------------------------------------

def initTemplate(parent, which):
    if which == 'avg' : 
        filename = askopenfilename(parent=parent)
        if len(filename) < 1:return
        tmp = os.path.basename(Spiderutils.name2template(filename))
        class_avg_var.set(os.path.splitext(tmp)[0])
    elif which == 'list' : 
        filename = askopenfilename(parent=parent)
        if len(filename) < 1:return
        tmp = os.path.basename(os.path.splitext(filename)[0])
        class_list_var.set(tmp)
    elif which == 'ext' : 
        filename = askopenfilename(parent=parent, 
           title="Choose file with valid extension")
        if len(filename) < 1:return
        tmp = os.path.splitext(filename)[1]
        datext_var.set(tmp)
    elif which == 'dir' : 
        filename = askdirectory(parent=parent, 
            title="Initial directory")
        if len(filename) < 1:return
        initdir_var.set(filename)

def buttonExit(win):
    win.destroy()
    prefs1.continue_boolean = False

class settings:
    # defaults
    dir1 = 'select/prj001'
    extension = '.dat'
    class_list_base = 'classes_by_ccc'
    class_avg_template = 'classavg***'
    serfile = 'stkfilt'
    partdocfile = 'docclass***'
    goodclasslist = 'goodclasses'
    goodparttemplate = 'byhand***'
    montageReduction = 1
    montageColor = 1
    badMontageColor = 'Red'
    xcoord = 10
    ycoord = 32
    ncol = 10
    maxPartCols = 12
    reverse = 1
    caLabels = 1
    textLabel = 'CCC'
    caReduction = 1
    
    def printSettings(self, event=None):
        print 
        print 'init_directory       ' + self.dir1
        print 'dat_extension        ' + self.extension
        print 'class_list           ' + self.class_list_base
        print 'class_average        ' + self.class_avg_template
        print 'particle_template    ' + self.serfile
        print 'particle_doc_file    ' + self.partdocfile
        print 'good_class_list      ' + self.goodclasslist
        print 'good_particle_list   ' + self.goodparttemplate
        print 'particle_reduction   ' + str(self.montageReduction)
        print 'class_reduction      ' + str(self.caReduction)
        print 'good_particle_color  ' + str(self.montageColor)
        print 'bad_particle_color   ' + self.badMontageColor
        print 'particle_xcoord      ' + str(self.xcoord)
        print 'particle_ycoord      ' + str(self.ycoord)
        print 'num_columns          ' + str(self.ncol)
        print 'max_particle_cols    ' + str(self.maxPartCols)
        print 'reverse_order        ' + str(self.reverse)
        print 'show_filenames       ' + str(self.caLabels)
        print 'text_label           ' + self.textLabel

    def readSettings(self):
        import os, string  # redundant?

        file = '.verifybyview'
        map = {}  # initialize dictionary

        if os.path.exists(file) : 
            input = open(file,'r')
            L = input.readlines()  # read line-by-line
            input.close()

            # read contents
            for line in L:
                try:
                    first  = string.split(line)[0]
                    second = string.split(line)[1]
                    map[first] = second
                except IndexError: pass   # if the field is blank, next steps will supply defaults
        else:
            print "Usage: verifybyview.py {initial_directory or class_list}"

        if map.has_key('init_directory'):      self.dir1 = map['init_directory']
        if map.has_key('class_list'):          self.class_list_base = map['class_list']
        if map.has_key('class_average'):       self.class_avg_template = map['class_average']
        if map.has_key('dat_extension'):       self.extension = map['dat_extension']

        if map.has_key('particle_template'):   self.serfile = map['particle_template']
        if map.has_key('particle_doc_file'):   self.partdocfile = map['particle_doc_file']
        if map.has_key('good_class_list'):     self.goodclasslist = map['good_class_list']
        if map.has_key('good_particle_list'):  self.goodparttemplate = map['good_particle_list']

        if map.has_key('particle_reduction'):  self.montageReduction = int(map['particle_reduction'])
        if map.has_key('class_reduction'):     self.caReduction = int(map['class_reduction'])
        if map.has_key('good_particle_color'): self.montageColor = map['good_particle_color']
        if map.has_key('bad_particle_color'):  self.badMontageColor = map['bad_particle_color']

        if map.has_key('particle_xcoord'):     self.xcoord = int(map['particle_xcoord'])
        if map.has_key('particle_ycoord'):     self.ycoord = int(map['particle_ycoord'])

        if map.has_key('num_columns'):         self.ncol = int(map['num_columns'])
        if map.has_key('max_particle_cols'):   self.maxPartCols = int(map['max_particle_cols'])
        if map.has_key('reverse_order'):       self.reverse = int(map['reverse_order'])
        if map.has_key('show_filenames'):      self.caLabels = int(map['show_filenames'])
        if map.has_key('text_label'):          self.textLabel = map['text_label']

    def saveSettings(self):
#        print 'Saving preferences'
#        self.printSettings()
    
        class_list_base = os.path.splitext(self.class_list_base)[0]
        class_avg_template = os.path.splitext(self.class_avg_template)[0]
    
        list = []
        list.append('init_directory       ' + self.dir1 + '\n')
        list.append('dat_extension        ' + self.extension + '\n')
        list.append('class_list           ' + class_list_base + '\n')
        list.append('class_average        ' + class_avg_template + '\n')

        list.append('particle_template    ' + self.serfile  + '\n')
        list.append('particle_doc_file    ' + self.partdocfile + '\n')
        list.append('good_class_list      ' + self.goodclasslist + '\n')
        list.append('good_particle_list   ' + self.goodparttemplate + '\n')

        list.append('particle_reduction   ' + str(self.montageReduction) + '\n')
        list.append('class_reduction      ' + str(self.caReduction) + '\n')
        list.append('good_particle_color  ' + str(self.montageColor) + '\n')
        list.append('bad_particle_color   ' + self.badMontageColor + '\n')

        list.append('particle_xcoord      ' + str(self.xcoord) + '\n')
        list.append('particle_ycoord      ' + str(self.ycoord) + '\n')

        list.append('num_columns          ' + str(self.ncol) + '\n')
        list.append('max_particle_cols    ' + str(self.maxPartCols) + '\n')
        list.append('reverse_order        ' + str(self.reverse) + '\n')
        list.append('show_filenames       ' + str(self.caLabels) + '\n')
        list.append('text_label           ' + self.textLabel + '\n')

        output = open('.verifybyview','w')
        output.writelines(list)

root = Tk()
root.title('verifybyview.py')
#ncol = None
#root.withdraw()

prefs1 = settings()  # class settings defined above
#prefs1.printSettings()
prefs1.readSettings()
#prefs1.printSettings()

if sys.argv[1:]:
    if os.path.isdir(sys.argv[1]) : 
        prefs1.dir1 = sys.argv[1]
    if os.path.isfile(sys.argv[1]) : 
        prefs1.dir1 = os.path.dirname(sys.argv[1])
        prefs1.class_list_base = os.path.basename(sys.argv[1])

# open new window for class-average template
avg_win = Toplevel(root)
avg_win.transient(root)  # WAS avg_win.lift(root)
avg_win.title("Initial settings")
init_frame = Frame(avg_win)
init_frame.focus_set()

# TO DO: create class for buttons + entries
datext_var = StringVar()
datext_button = Button(init_frame, text='Data extension', 
    command = lambda w=avg_win, d='ext': initTemplate(w,d))
datext_button.grid(row=0, column=0, sticky=W)
datext_var.set(prefs1.extension)
datext_entry = Entry(init_frame, textvariable=datext_var, width=34)
datext_entry.grid(row=0, column=1)

initdir_var = StringVar()
initdir_button = Button(init_frame, text='Initial directory', 
    command = lambda w=avg_win, d='dir': initTemplate(w,d))
initdir_button.grid(row=1, column=0, sticky=W)
initdir_var.set(prefs1.dir1)
initdir_entry = Entry(init_frame, textvariable=initdir_var, width=34)
initdir_entry.grid(row=1, column=1)

class_list_var = StringVar()
class_list_button = Button(init_frame, text='Class list doc file', 
    command = lambda w=avg_win, d='list': initTemplate(w,d))
class_list_button.grid(row=2, column=0, sticky=W)
class_list_var.set(os.path.splitext(prefs1.class_list_base)[0])
class_list_entry = Entry(init_frame, textvariable=class_list_var, width=34)
class_list_entry.grid(row=2, column=1)

class_avg_var = StringVar()
classavg_button = Button(init_frame, text='Class-average template', 
    command = lambda w=avg_win, d='avg': initTemplate(w,d))
classavg_button.grid(row=3, column=0, sticky=W)
class_avg_var.set(os.path.splitext(prefs1.class_avg_template)[0])
class_avg_entry = Entry(init_frame, textvariable=class_avg_var, width=34)
class_avg_entry.grid(row=3, column=1)

init_frame.pack(padx=5, pady=5, ipadx=2, ipady=2)
quit_frame = Frame(avg_win, padx=40)
#Button(quit_frame, text="Done", command=avg_win.destroy).pack(pady=4)
done_button = Button(quit_frame, text="Continue", command=avg_win.destroy)
done_button.pack(side=LEFT)
prefs1.continue_boolean = True
exit_button = Button(quit_frame, text="Exit", command=lambda w=avg_win: buttonExit(w))
exit_button.pack(side=RIGHT)
quit_frame.pack(fill=X)
avg_win.bind('<Return>', lambda p=root, w=avg_win: clickOK(p,w))  # parmquit(p,w))

# wait for initial window
root.wait_window(avg_win)
if not prefs1.continue_boolean : sys.exit()

prefs1.extension = datext_var.get()
prefs1.class_list_base = class_list_var.get() # + prefs1.extension
prefs1.class_avg_template = class_avg_var.get() # + prefs1.extension
prefs1.dir1 = initdir_var.get()

prefs1.saveSettings()

mn = classAverages(root, prefs1)  # , useLabels=1)

root.mainloop()
