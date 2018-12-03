#!/usr/bin/env python
#
# SOURCE:   spider/spire/bin/classavg.py 
#
# PURPOSE:  Tool for creating & viewing class averages
#
# Spider Python Library
# Copyright (C) 2006-2018  Health Research Inc., Menands, NY
# Email:  spider@health.ny.gov

import re,os,sys
import montage
import Pmw

from   Tkinter           import *
from   tkFileDialog      import askopenfilename
from   PIL               import Image, ImageTk

from   Spider            import Spiderutils
from   Spider            import SpiderImagePlugin # 2018 al import *

# -------------------------------------------------------------------
# classAverages:  Inherits from support.montage
#                 overrides makeMenus() and display()

class classAverages(montage.montage): 
    
    def __init__(self, master, imagelist,
                 title="Class Averages",
                 ncol=None, useLabels=0):
        
        self.master = master
        if len(imagelist) < 1: return
        montage.montage.__init__(self, master,
                                imagelist=imagelist,
                                title=title,
                                ncol=ncol,
                                useLabels=useLabels)

        self.top.protocol("WM_DELETE_WINDOW", master.quit)
        # set up some templates
        self.docfile = ""  # full path of docfile - 1st column has particle numbers
        self.doctemplate = StringVar()
        self.doctemplate.set(self.docfile)   # display only the basename
        self.serfile = ""
        self.sertemplate = StringVar()
        self.sertemplate.set(self.serfile)        # particle images
        
        self.extension = os.path.splitext(imagelist[0])[1]
        if self.extension != "":
            ext = self.extension
            self.filetypes = [ ("", "*"+ext), ("All files", "*") ]
        else:
            self.filetypes = [ ("All files", "*") ]

        self.makeMenus()
        self.createMontage()

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
        Filebtn.menu.add_command(label='Class avgs',
                                 command=self.getClassAverages)
        Filebtn.menu.add_command(label='Templates',
                                 command=self.callTemplates)
        Filebtn.menu.add_command(label='Close all windows',
                                 command=self.closeWindows)
        Filebtn.menu.add_separator()
        Filebtn.menu.add_command(label='Quit', underline=0,
                                     command=self.master.quit)
        Filebtn['menu'] = Filebtn.menu
    
        # Make the Display menu (same as montage.py)
        Dspbtn = Menubutton(self.mBar, text='Display', relief='flat')
        Dspbtn.pack(side=LEFT, padx=5, pady=5)
        Dspbtn.menu = Menu(Dspbtn, tearoff=0)

        Dspbtn.menu.add_command(label='no. columns', underline=0,
                                    command=self.displayParms)
        Dspbtn.menu.add_checkbutton(label='show filenames', underline=0,
                                    variable=self.showVar,
                                    command=self.displayParms_1)
        # 'size' has a submenu of checkbuttons
        Dspbtn.menu.sizes = Menu(Dspbtn.menu)
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
        Dspbtn['menu'] = Dspbtn.menu

    def closeWindows(self):
        winlist = self.top.winfo_children()
        for win in winlist:
            if win.winfo_exists():
                if win.winfo_class() == 'Toplevel':
                    win.destroy()        

    def hello(self, im):
        " filename is the classavg image "
        if self.docfile == "" or self.serfile == "":
            self.callTemplates()

        filename = im.info['filename']
            
        #print "hello from " + filename
        tmp = self.doctemplate.get()
        fn = Spiderutils.template2filename(tmp, os.path.basename(filename))
        if len(fn) < 1:
            return
        if not os.path.exists(fn):
            print "Unable to find " + fn
            return
        D = Spiderutils.readdoc(fn, keys='all')
        plist = []
        files = D.keys()
        for f in files:
            num = int(D[f][0])  # each element is a list (of column data)
            imgfilename = Spiderutils.template2filename(self.serfile, num)
            plist.append(imgfilename)
        title = "images for class: " + os.path.basename(filename)
        newtop = Toplevel(self.top)
        m = montage.montage(newtop, imagelist=plist, title=title,
                            useLabels=0)
        m.makeMenus()
        m.createMontage()

    def callTemplates(self):
        gtwin = Toplevel(self.top)
        self.getTemplates(gtwin)
        self.top.wait_window(gtwin)  # wait til template window gone
        # if they typed in w/o asterisks..
        doc = self.doctemplate.get()
        if len(doc) > 0:
            self.docfile = doc
            if string.find(doc,"*") < 0:
                tmp = Spiderutils.name2template(doc)
                self.doctemplate.set(tmp)
        ser = self.sertemplate.get()
        if len(ser) > 0:
            self.serfile = ser
            if string.find(ser,"*") < 0:
                tmp = Spiderutils.name2template(ser)
                self.sertemplate.set(tmp)

    def getTemplates(self, win):
        win.title("Templates")
        ftop = Frame(win)
        bdoc = Button(ftop, text="Doc file template: ",
                      command = lambda w=win, d='doc': self.setTemplates(w,d))
        edoc = Entry(ftop, textvariable = self.doctemplate)
        bser = Button(ftop, text="Particle template: ",
                      command = lambda w=win, d='ser': self.setTemplates(w,d))
        eser = Entry(ftop, textvariable = self.sertemplate)
        bdoc.grid(row=0, column=0, sticky='w')
        edoc.grid(row=0, column=1, sticky='nsew')
        ftop.columnconfigure(1, weight=1)  # make entry expand
        bser.grid(row=1, column=0, sticky='w')
        eser.grid(row=1, column=1, sticky='nsew')
        ftop.pack(side='top', fill='both', expand=1)
        
        fbut = Frame(win)
        Button(fbut, text="Done", command=win.destroy).pack(padx=5, pady=5)
        fbut.pack()

    def setTemplates(self, parent, which='doc'):
        filename = askopenfilename(parent=parent, filetypes=self.filetypes)
        if len(filename) < 1:
            return
        tmp = Spiderutils.name2template(filename)
        if which == 'doc':
            self.doctemplate.set(tmp)
        else:
            self.sertemplate.set(tmp)

    def getClassAverages(self):
        filenames = askopenfilename(filetypes=self.filetypes, multiple=1)
        if len(filenames) < 1:
            return
        self.imagelist = loadImageSeries(filenames)
        self.ncol = self.montagesize()
        self.createMontage()

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
            
        useLabels = self.useLabels
        i = 0
        j = 0
        for im in self.imagelist:
            ix,iy = im.size
            if ix != x or iy != y:
                icpy = im.copy()  # o.w. orig is resized
                img = icpy.resize( (x,y) )
            else:
                img = im

            photo = ImageTk.PhotoImage(img, palette=256, master=parent)
            b = Button(parent,image=photo, command=lambda f=img: self.hello(f))
            b.photo = photo
            filename = im.info['filename']
            b.grid(row=j, column=i)
            if useLabels:
                lb = Label(parent, text=os.path.basename(filename))
                lb.grid(row=j+1, column=i)
            i += 1
            if i > self.ncol-1:
                i = 0
                j += 1
                if useLabels: j += 1 # incr j twice for labels

# -------------------------------------------------------------------
if __name__ == "__main__":

    root = Tk()

    ncol = None
        
    if not sys.argv[1:]:
        filelist = askopenfilename(multiple=1)
        if len(filelist) < 1:
            print "Syntax: python classavg.py imagelist"
            #print "Syntax: python classavg.py -ncol n imagelist"
            sys.exit(1)
    elif sys.argv[1] == "-ncol":
        ncol = int(sys.argv[2])
        filelist = sys.argv[3:]
    else:
        filelist = sys.argv[1:]
    
    #root.withdraw()
    # 2018 al FAILED   Image.register_open("SPIDER", SpiderImageFile)
    mn = classAverages(root, filelist, ncol=ncol, useLabels=1)
    root.mainloop()
    
