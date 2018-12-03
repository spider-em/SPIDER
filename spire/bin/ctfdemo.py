#!/usr/bin/env python
#
# SOURCE:  ctfdemo.py 
#
# PURPOSE: Demonstration of CTF parameters
#
# Spider Python Library
# Copyright (C) 2006-2018  Health Research Inc., Menands, NY
# Email:    spider@health.ny.gov

from   Tkinter import * 
from   math    import *

import Pmw

class CTFplot:
    " default values "
    def __init__(self, master,
                 cs        = 2.0,
                 defocus   = 20000,
                 kev       = 200,
                 pixelsize = 2.82,
                 src       = 0.0,    # Source size
                 spread    = 0.0,    # Defocus spread
                 acr       = 0.0,    # Amplitude contrast ratio
                 gep       = 2):     # Gaussian envelope parameter
        self.top = master
        self.top.title("ctfdemo")
        self.cs      = StringVar();      self.cs.set(cs)
        self.defocus = StringVar(); self.defocus.set(defocus)
        self.kev     = StringVar();     self.kev.set(kev)
        self.pixsize = StringVar(); self.pixsize.set(pixelsize)
        self.src     = StringVar();     self.src.set(src)
        self.spread  = StringVar();  self.spread.set(spread)
        self.acr     = StringVar();     self.acr.set(acr)
        self.gep     = StringVar();     self.gep.set(gep)

        self.vardict = {}

        self.envelopeShow = IntVar()
        self.envelopeShow.set(0)
        self.squared      = IntVar()
        self.squared.set(0)

        self.n             = 250
        self.max_spat_freq = 1.0 / (2.0 * float(self.pixsize.get()))
        self.X             = []
        self.Y             = []
        self.E             = []
        for i in range(self.n):
            self.X.append( (self.max_spat_freq / float(self.n)) * i)
            self.Y.append(0.0)
            self.E.append(0.0)

        self.kappa = -pi**2 / (16.0 * log(2.0))
        self.infinity = 1e50
        self.compute()
        acr = self.acr.get()
        if float(acr) == 0: self.acr.set(0.1)

        # ------- create the menu bar -------
        self.mBar = Frame(master, relief='raised', borderwidth=1)
        self.mBar.pack(side='top', fill = 'x')

        # Make the File menu
        Filebtn = Menubutton(self.mBar, text='File', underline=0,
                                 relief='flat')
        Filebtn.pack(side=LEFT, padx=5, pady=5)
        Filebtn.menu = Menu(Filebtn, tearoff=0)
        Filebtn.menu.add_separator()
        Filebtn.menu.add_command(label='Quit', underline=0,
                                     command=master.quit)
        Filebtn['menu'] = Filebtn.menu

        # Make the Option menu
        Optbtn = Menubutton(self.mBar, text='Options', underline=0,
                                 relief='flat')
        Optbtn.pack(side=LEFT, padx=5, pady=5)
        Optbtn.menu = Menu(Optbtn, tearoff=0)

        Optbtn.menu.add_checkbutton(label='Grid', underline=0,
                                    command=self.showGrid)
        Optbtn.menu.add_checkbutton(label='Show Envelope', underline=0,
                                    command=self.showEnvelope)
        Optbtn.menu.add_checkbutton(label='Squared', underline=0,
                                    command=self.showSquared)
        Optbtn['menu'] = Optbtn.menu
       
        # ------- Widgets start here -------
        ff = Frame(master)
        fg = Frame(ff, relief='raised', borderwidth=2) # Upper left frame for plot
        self.g = Pmw.Blt.Graph(fg) 
        self.curveLine = 'model'
        self.g.line_create(self.curveLine,
                     xdata=tuple(self.X),
                     ydata=tuple(self.Y),
                     color = 'blue',
                     symbol='')

        self.envLine = 'envelope'
        self.g.line_create(self.envLine,
                     xdata=tuple(self.X),
                     ydata=tuple(self.E),
                     color = 'green',
                     symbol='')

        self.g.configure(title='Transfer function demo')
        self.g.legend_configure(hide=1)
        # Fix the limits of the axes (o.w. axes move, not plot)
        ymin,ymax = self.g.axis_limits("y")
        self.g.axis_configure("y", min=ymin, max=ymax)
        xmin,xmax = self.g.axis_limits("x")
        #self.g.axis_configure("x", min=xmin, max=xmax, title="Spatial Frequency")

        self.g_width = 446 #int(self.g.extents("plotwidth"))
        self.g_height = 150 #int(self.g.extents("plotheight"))
        self.showEnvelope()
        self.g.grid(row=0, column=0, sticky='nsew')
        fg.columnconfigure(0, weight=1) 

        fp = Frame(ff, relief='raised', borderwidth=2)  #Upper right frame for pixsize

        plabel = Label(fp,text="pixelsize")
        pslider = Scale(fp, orient='vertical', from_=0, to=6,
                       tickinterval = 1,
                       resolution = 0.01, label ="",
                       variable = self.pixsize,
                       length =  self.g_height,
                       showvalue=0,
                       command=self.xupdate)

        pentry = Entry(fp, textvariable=self.pixsize, width=10, background='white')
        pentry.bind('<Return>', self.xupdate)
        plabel.grid(row=0, column=0)
        pentry.grid(row=1, column=0)
        pslider.grid(row=2, column=0)
   
        fg.grid(row=0, column=0, sticky='nsew')
        fp.grid(row=0, column=1, sticky='ns')
        ff.columnconfigure(0, weight=1)
        ff.rowconfigure(0, weight=1)
        ff.pack(side='top', fill='x') #, expand=1)

        # ------ the set of sliders -------

        self.sf = Pmw.ScrolledFrame(master, horizflex='expand', vertflex='fixed',
                                    vscrollmode='dynamic')
        f = self.sf.interior()
        #f = Frame(master, relief='raised', borderwidth=2)
        self.slider(f, start=0, end=50000, row=0,
                         label='defocus',
                         tickinterval=10000,
                         resolution=50,
                         variable = self.defocus)
        self.slider(f, start=0, end=500, row=1,
                         label='electron\nenergy (kev)',
                         tickinterval=100,
                         variable = self.kev)
        self.slider(f, start=0, end=5.0, row=2,
                         label='spherical\naberration',
                         tickinterval=1,
                         variable = self.cs)
        self.slider(f, start=0, end=0.005, row=3,
                         label='source\nsize',
                         tickinterval=0.001,
                         digits = 2,
                         variable = self.src)
        self.slider(f, start=0, end=500, row=4,
                         label='defocus\nspread',
                         tickinterval=100,
                         variable = self.spread)
        self.slider(f, start=0, end=0.5, row=5,
                         label='amplitude\ncontrast ratio',
                         tickinterval=0.1,resolution=0.01,
                         variable = self.acr)
        self.slider(f, start=0, end=2, row=6,
                         label='Gaussian\nenvelope',
                         tickinterval=0.5,
                         variable = self.gep)
        f.columnconfigure(2, weight=1) # makes column expandable
        self.sf.pack(side='bottom', expand=1, fill='both')


    def slider(self, master, start=0, end=10, row=0, label="",
               tickinterval=1, resolution=None, digits=0,variable = None):

        #variable.trace_variable("w", self.varchange)
        #self.vardict[variable._name] = variable
        
        lab = Label(master, text=label)
        if resolution == None:
            resolution = float(tickinterval) / 50.0
        slider = Scale(master, orient='horizontal', from_=start, to=end,
                       tickinterval = tickinterval,
                       resolution = resolution, label ="",
                       variable = variable,
                       showvalue=0,
                       #length = self.g_width,
                       digits = digits,
                       command=self.update)
        ent = Entry(master, textvariable=variable, width=10, background='white')
        ent.bind('<KeyPress>', self.update)

        lab.grid(row=row, column=0, sticky='ne', padx=5, pady=5)
        ent.grid(row=row, column=1, sticky='n', padx=5, pady=5)
        slider.grid(row=row, column=2, sticky='ew', padx=5, pady=5)
        

    def compute(self):
        cs    = 1e7 * float(self.cs.get())
        kv = float(self.kev.get())
        if kv != 0:
            lmbda = 12.398 / sqrt(kv* (1022+kv))
        else:
            lmbda = self.infinity
        if cs != 0:
            f1    = 1.0 / sqrt(cs*lmbda)
            f2    = sqrt(sqrt(cs*lmbda**3))
        else:
            f1 = self.infinity
            f2 = self.infinity

        pixsize = float(self.pixsize.get())
        if pixsize != 0:
            self.max_spat_freq = 1.0 / (2.0 * pixsize)
        km1   = f2 * self.max_spat_freq
        dk    = km1 / float(self.n)
        q1    = (float(self.src.get()) * f2)**2
        gep = float(self.gep.get())
        if gep != 0:
            env   = 1.0/gep**2
        else:
            env = self.infinity
        env1  = env/f2**2
        f     = -pi**2
        ds1   = f1 * float(self.spread.get())
        kappa = ds1 * self.kappa
        dz1   = f1 * float(self.defocus.get())

        envFlag = self.envelopeShow.get()
        acr = float(self.acr.get())
        squared = self.squared.get()

        for i in range(self.n):
            ak = i * dk
            p  = ak**3 - dz1 * ak
            ch = exp(ak*4 * kappa)
            self.E[i] = (exp(f*q1*p**2)*ch)*2*exp(-env1*ak**2)
            qqt = 2.0*pi*(0.25*ak**4 - 0.5*dz1*ak**2)
            qqt1 = (1.0-acr)*sin(qqt)-acr*cos(qqt)
            self.Y[i] = self.E[i] * qqt1
            if squared:
                self.E[i] = (self.E[i])**2
                self.Y[i] = (self.Y[i])**2

    def update(self, scalevalue=None):
        #if scalevalue != None:
            #print scalevalue
        self.compute()
        self.g.element_configure(self.curveLine, ydata=tuple(self.Y))
        self.g.element_configure(self.envLine, ydata=tuple(self.E))
        self.top.update_idletasks()

    def xupdate(self, scalevalue=None):
        pixsize = float(self.pixsize.get())
        if pixsize != 0:
            self.max_spat_freq = 1.0 / (2.0 * pixsize)
            for i in range(self.n):
                self.X[i] = i* (self.max_spat_freq / float(self.n))
            self.g.element_configure(self.curveLine, xdata=tuple(self.X))
            self.g.element_configure(self.envLine, xdata=tuple(self.X))
            self.g.axis_configure("x", max=self.max_spat_freq)
            self.update()

    def resetYaxis(self):
        ymin = ymax = self.Y[0]
        for i in range(self.n):
            if ymin > self.Y[i]: ymin = self.Y[i]
            if ymax < self.Y[i]: ymax = self.Y[i]
        self.g.axis_configure("y", min=ymin, max=ymax)
       

    def openFile(self):
        return("filename")
    
    def showGrid(self):
        self.g.grid_toggle()

    def showEnvelope(self):
        show = self.envelopeShow.get()
        self.envelopeShow.set(not show)
        if show:
            self.g.element_show([self.curveLine, self.envLine])
        else:
            self.g.element_show([self.curveLine])
        #self.update()

    def showSquared(self):
        self.squared.set(not self.squared.get())
        self.compute()
        self.resetYaxis()
        self.update()

    """ I think the only effect this has is to call update a 2nd time """
    def varchange(self, varName, index, mode):
        variable = self.vardict[varName]
        self.update()
        #print "varchange %s" % variable.get()


# ------- end CTFplot class definition

if __name__ == '__main__':

    master = Tk()
                 
    c = CTFplot(master)
    master.mainloop() 
