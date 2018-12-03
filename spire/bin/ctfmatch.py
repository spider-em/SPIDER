#!/usr/bin/env python
#
# SOURCE:  ctfmatch.py
# 
# PURPOSE: Display CTF doc files 
#
# Spider Python Library
# Copyright (C) 2006-2018  Health Research Inc., Menands, NY
# Email:    spider@health.ny.gov

from   Tkinter        import * 
from   tkFileDialog   import askopenfilename, asksaveasfilename
import tkMessageBox
import Pmw
import string, sys
import os, time
from   math           import *
from   commands       import getoutput

from Spider           import Spiderutils

import webbrowser
webpage = "http://www.wadsworth.org/spider_doc/spider/spire/tools-docs/ctfmatch.html"

def ctfhelp():
    try:
        webbrowser.open(webpage,1)
    except:
        pass

def ctfabout():
    s = "CTFmatch 1.0\n\n" + \
        "A tool for analyzing the output " +\
        "from SPIDER's CTF FIND command."
    tkMessageBox.showinfo("About CTFmatch 1.0", s)

def integer(astring):
    if type(astring) == type(""):
        dot = string.find(astring,".")
        if dot > 0:
            astring = astring[:dot]
    return int(astring)

def writedoc(filename, column1=1, column2=0):
    column1 = integer(column1)
    column2 = integer(column2)
    # If file already exists
    if os.path.exists(filename):
        # if it's a doc file, try to get the last key
        if Spiderutils.isSpiderDocfile(filename):
            lastline = getoutput("tail -1 %s" % filename)
            if len(lastline) > 0:
                key = 1 + int(string.split(lastline)[0])
            else:
                key = 1
        else:
            # if it's not a doc file...
            os.remove(filename)
            key = 1
            
        data = "%5d 2 %11d% 11d\n" % (key, column1, column2)

        try:
            fp = open(filename, 'a')  # append
            fp.write(data)
            fp.close()
        except:
            print "Unable to write to %s" % filename
            return 0
    # if it's a new file
    else:
        try:
            fp    = open(filename, 'w')
            fname = os.path.basename(filename)
            ext   = os.path.splitext(filename)[1]
            ext   = ext[1:] # remove dot
            date,time,id = Spiderutils.nowisthetime()
            h = " ;ctf/%s   %s AT %s   %s\n" % (ext,date,time,fname)
            fp.write(h)
            fp.write(" ; /     MICROGRAPH   DEFOCUS\n")
            key  = 1
            data = "%5d 2 %11d %11d\n" % (key, column1, column2)
            fp.write(data)
            fp.close()
        except:
           print "Unable to create %s" % filename
           return 0
    return 1
            
def readDefocus(filename):
    " returns list of (mic#, defocus) pairs (mic=int defocus=string)"
#    F = spiderutils.readSpiderDocFile(filename, col_list=(1,2))
    F = Spiderutils.readdoc(filename, keys='all')
    if F == None: return []
    keys = F.keys()
    keys.sort()
    
    M = []
    for key in keys:
        mic     = int(F[key][0])
        defocus = str(int(F[key][1]))
        M.append( (mic, defocus) )
    return M
            
    
def readdoc(filename, factor=1.0, squared=1):
    F = Spiderutils.readdoc(filename, keys='all')
    if F == None: return []
    roofile = 0
    A = []; B = []; C = []; D = []; E = []
    keys = F.keys()
    keys.sort()
    
    # get the 1st line of data, test if roo (cols 3 & 4 = 1)
    k = keys[0]
    vals = F[k]
    vlen = len(vals)
    if vlen < 4:
        roofile = 1
    else:  #
        if vals[2] == 1 and vals[3] == 1:
            roofile = 1
            
    # get the data
    for key in keys:
        if roofile:
            a = factor * F[key][0]
            D.append(a)
        else:
            freq = F[key][0] 
            bgd  = factor * F[key][1]
            sub  = factor * F[key][2]
            env  = factor * F[key][3]
            roo  = bgd + sub     # get original spectrum
            if squared == 1:
                k = env + bgd  
                roo = roo * roo  # square the spectrum
                bgd = bgd * bgd  # square the background
                sub = roo - bgd  # new subtracted curve
                env = k*k - bgd  # new envelope
            A.append(freq)
            B.append(bgd)
            C.append(sub)
            D.append(env)
            E.append(roo)

    if roofile:
        return [D]
    else:
        return [A,B,C,D,E]

def getFiles(filetypes=None):
    ft = []
    if filetypes != None:
        ft.append( (filetypes,filetypes )) # "*.dat" --> ("*.dat", "*.dat")
    ft.append(("All files", "*"))
    f = askopenfilename(multiple=1, filetypes=ft)
    "f is a long string under irix, but a tuple of strings in linux. Neat, huh?"
    if type(f) == type("string"):
        return string.split(f)
    else:
        return f  # hopefully a tuple or a list


###################################################################
#
# CTF plot

class CTFplot:
    " default values "
    def __init__(self, master, filename=None, args=None):
        # first set defaults
        self.top = master
        self.cs = StringVar();      self.cs.set(2.0)
        self.defocus = StringVar(); self.defocus.set(20000)
        self.kev = StringVar();     self.kev.set(200)
        self.pixsize = StringVar(); self.pixsize.set(2.82)
        self.src = StringVar();     self.src.set(0.0)
        self.spread = StringVar();  self.spread.set(0.0)
        self.acr = StringVar();     self.acr.set(0.1)
        self.gep = StringVar();     self.gep.set(2.0)
        self.defocusfile = ""
        self.tfed = []
        got_pixsize = 0
        got_kev = 0
        max_defocus = 60000
        # then process args, if any
        if args != None:
            keys = args.keys()
            for key in keys:
                k = key.replace('-','') # delete minus sign
                if k == "defocus":
                    self.defocusfile = args[key]
                elif k == "tfed":
                    self.tfed = args[key]
                else:
                    s = "self." + k + ".set(" + args[key] + ")"
                    try:
                        #print "%s : %s" % (k,s)
                        eval(s)
                        if k == 'pixsize' and float(self.pixsize.get()) != 0:
                            got_pixsize = 1
                        elif k == 'kev' and float(self.kev.get()) != 0:
                            got_kev = 1
                    except:
                        print "key %s not recognized" % key

        self.multfactor = 1.0  #1000.0  # for input data
        self.datamax = 1.0  # used to control height of model
        self.cutoff = 40     # index (integer) of modelmax on x axis
        self.modelmax = 1.0
        self.factor = 1.0
        self.askParms = 1
        # don't start with parameter window if have pixsize and kev
        if got_pixsize and got_kev:
            self.askParms = 0        
        
        self.max_spat_freq = 1.0 / (2.0 * float(self.pixsize.get()))
        self.kappa = -pi**2 / (16.0 * log(2.0))
        self.infinity = 1e50
        self.ymax = StringVar()
        self.xmin = StringVar() ; self.xmin.set(0)
        self.modymax = StringVar()
        self.n = 250
        # arrays for holding data columns
        self.arr = {'frq':[], 'bgd':[], 'sub':[], 'env':[], 'roo':[]}
        self.showlist = {'bgd':1, 'sub':1, 'env':1, 'roo':1, 'model':1}
        # variables for checkbuttons in menus
        self.showRoo = IntVar() ; self.showRoo.set(1)
        self.showBgd = IntVar() ; self.showBgd.set(1)
        self.showSub = IntVar() ; self.showSub.set(1)
        self.showEnv = IntVar() ; self.showEnv.set(1)
        self.showMod = IntVar() ; self.showMod.set(1)
        self.colors = {'bgd':'#00cc00', 'sub':'red', 'env':'#9999ff',
                       'roo':'#ff9900',  'model':'white'}
        # by default, data is squared as it is read in 
        self.squared = IntVar() ; self.squared.set(1)
        # whether to use the empirical envelope
        self.envelope = IntVar() ; self.envelope.set(0)
        self.roofile = ""
        self.tfedfile = ""
        self.savefile = ""
        self.defDict = {}   # [micno] = defocus
        self.fileDict = {}   # [micno] = (basename, fullpath)
            
        # set up model data and read file
        self.X = []
        self.Y = []
        self.E = []

        if filename != None:
            if not self.openFile(filename):
                filename = None

        if len(self.arr['frq']) > 0:
            self.X = self.arr['frq']
        else:
            for i in range(self.n):
                self.X.append((self.max_spat_freq / float(self.n)) * i)
        for i in range(self.n):
            self.Y.append(0.0)
            self.E.append(0.0)

        if filename != None:
            self.newymax()
        else:
            self.ymax.set('1.0')
        self.modymax.set(1.0)
        self.compute() # generate the model

        # ------- create the menu bar -------
        self.mBar = Frame(master, relief='raised', borderwidth=1)
        self.mBar.pack(side='top', fill = 'x')
        self.balloon = Pmw.Balloon(self.top)
        
        # Make the File menu
        Filebtn = Menubutton(self.mBar, text='File', underline=0,
                                 relief='flat')
        Filebtn.pack(side=LEFT, padx=5, pady=5)
        Filebtn.menu = Menu(Filebtn, tearoff=0)
        Filebtn.menu.add_command(label='Open TF ED file',
                                 command=self.callOpenFile)
        Filebtn.menu.add_command(label='Open File series',
                                 command=self.fileSeries)
        Filebtn.menu.add_command(label='Open Defocus file',
                                 command=self.openDefocus)
        Filebtn.menu.add_command(label='Save Defocus as...',
                                 command=self.saveDefocusAs)
        Filebtn.menu.add_separator()
        Filebtn.menu.add_command(label='Quit', underline=0,
                                     command=self.quit)
        Filebtn['menu'] = Filebtn.menu

        
        # Make the Option menu
        Optbtn = Menubutton(self.mBar, text='Options', relief='flat')
        Optbtn.pack(side=LEFT, padx=5, pady=5)
        Optbtn.menu = Menu(Optbtn, tearoff=0)

        Optbtn.menu.add_command(label='Parameters', underline=0,
                                    command=self.callSetParms)
        Optbtn.menu.add_checkbutton(label='Squared data', underline=0,
                                    variable=self.squared,
                                    command=self.showSquared)
        Optbtn.menu.add_checkbutton(label='Grid', underline=0,
                                    command=self.showGrid)
        Optbtn.menu.add_checkbutton(label='Use empirical envelope',
                                    variable=self.envelope,
                                    command=self.replot)
        Optbtn['menu'] = Optbtn.menu
       
        # Parameter menu
        #menuBar.addmenuitem('Model', 'checkbutton', 'reset yscale',
                            #indicatoron=0,
                            #label='Reset Ymax', command=self.showGrid)
        
        # Make the Show menu
        Showbtn = Menubutton(self.mBar, text='Show', relief='flat')
        Showbtn.pack(side=LEFT, padx=5, pady=5)
        Showbtn.menu = Menu(Showbtn, tearoff=0)

        Showbtn.menu.add_checkbutton(label='1D spectrum', 
                                    background = 'black',
                                    foreground = self.colors['roo'],
                                    selectcolor='red',
                                    activeforeground = self.colors['roo'],
                                    variable = self.showRoo,
                                    command=self.replot)
        Showbtn.menu.add_checkbutton(label='Background', 
                                    background = 'black',
                                    foreground = self.colors['bgd'],
                                    selectcolor='red',
                                    activeforeground = self.colors['bgd'],
                                    variable = self.showBgd,
                                    command=self.replot)
        Showbtn.menu.add_checkbutton(label='Subtracted data', 
                                    background = 'black',
                                    foreground = self.colors['sub'],
                                    selectcolor='red',
                                    activeforeground = self.colors['sub'],
                                    variable = self.showSub,
                                    command=self.replot)
        Showbtn.menu.add_checkbutton(label='Envelope', 
                                    background = 'black',
                                    foreground = self.colors['env'],
                                    selectcolor='red',
                                    activeforeground = self.colors['env'],
                                    variable = self.showEnv,
                                    command=self.replot)
        Showbtn.menu.add_checkbutton(label='Model', 
                                    background = 'black',
                                    foreground = self.colors['model'],
                                    selectcolor='red',
                                    activeforeground = self.colors['model'],
                                    variable = self.showMod,
                                    command=self.replot)
        Showbtn['menu'] = Showbtn.menu
        
        # Help menu
        Helpbtn = Menubutton(self.mBar, text='Help', relief='flat')
        Helpbtn.pack(side=RIGHT, padx=5, pady=5)
        Helpbtn.menu = Menu(Helpbtn, tearoff=0)

        Helpbtn.menu.add_command(label='Help', command=ctfhelp)
        Helpbtn.menu.add_command(label='About', command=ctfabout)
        Helpbtn['menu'] = Helpbtn.menu
        
        
        # ------- widgets start here -------
        self.g_width = 446 #int(self.g.extents("plotwidth"))
        self.g_height = 150 #int(self.g.extents("plotheight"))

        ff = Frame(master) # frame that holds everything
        
        " yscale slider "
        fy = Frame(ff, relief='raised', borderwidth=2)
        ylabel = Label(fy,text="y max")
        ymax = string.atof(self.ymax.get())

        self.yslider = Scale(fy, orient='vertical', from_= ymax, to=0.0,
                       tickinterval = ymax/6.0,
                       resolution = 0.01,   #ymax/30.0,
                       label ="",
                       variable = self.ymax,
                       length =  self.g_height,
                       showvalue=0,
                       command=self.yupdate)

        " model height scale "
        mlabel = Label(fy,text="model\nheight")
        mmax = 1.0   #string.atof(self.modymax.get())

        self.mslider = Scale(fy, orient='vertical', from_= mmax, to=0.0,
                       tickinterval = mmax/5.0,
                       resolution = mmax/50.0,
                       label ="",
                       variable = self.modymax,
                       length =  self.g_height,
                       showvalue=0,
                       command=self.update)
        
        ylabel.pack(side='top')
        self.yslider.pack(side='top', padx=5, pady=5)
        self.mslider.pack(side='bottom', padx=5, pady=5)
        mlabel.pack(side='bottom')


        " the main plot "
        fg = Frame(ff, relief='raised', borderwidth=2)
        self.g = Pmw.Blt.Graph(fg, plotbackground="black" ) 
        self.curves = ['bgd','sub','env', 'roo']
        i = 0
        for curve in self.curves:
            self.g.line_create(curve,
                         xdata=tuple(self.X),
                         ydata=tuple(self.arr[curve]),
                         color = self.colors[curve],
                         symbol='')
            i += 1

        # the model curve
        self.g.line_create('model',
                     xdata=tuple(self.X),
                     ydata=tuple(self.Y),
                     color = self.colors['model'],
                     symbol='')
        self.g.legend_configure(hide=1)
        # fix the limits of the axes (o.w. axes move, not plot)
        ymin,ymax = self.g.axis_limits("y")
        self.g.axis_configure("y", min=ymin, max=ymax)
        xmin,xmax = self.g.axis_limits("x")
        self.g.axis_configure("x", min=xmin, max=xmax, title="Spatial Frequency")
        self.xmin.set(xmin)
        xmax = self.max_spat_freq

        xslider = Scale(fg, orient='horizontal', from_= xmin, to=xmax,
                       tickinterval = xmax/4.0,
                       resolution = 0.001,   #xmax/40.0,
                       label ="x min",
                       variable = self.xmin,
                       length =  int(self.g_width) / 2,
                       showvalue=0,
                       command=self.xminupdate)

        xbutton = Button(fg, text="reset ymax", command=self.resetYmax)

        saveBut = Button(fg, text='Save Defocus', command=self.saveDefocus)
        
        fl = Frame(fg, relief='sunken', borderwidth=2)
        #self.roolabel = Label(fl, text=self.roofile, background=filebgdcolor)
        self.tfedlabel = Label(fl, text=os.path.basename(self.tfedfile))
        self.defocuslabel = Label(fl, text=os.path.basename(self.defocusfile))
        self.savelabel = Label(fl, text=os.path.basename(self.savefile))
        self.tfedlabel.pack(side='top', padx=10, pady=5)
        self.defocuslabel.pack(side='top', padx=10, pady=5)
        self.savelabel.pack(side='top', padx=10, pady=5)

        self.g.pack(side='top', fill='both', expand=1)
        xslider.pack(side='left', padx=5, pady=5)
        xbutton.pack(side='left', padx=5, pady=5)
        fl.pack(side='right', padx=10, pady=10)
        saveBut.pack(side='right', padx=10, pady=10)
        
        #fy.pack(side='left', expand=1, fill='y')
        #fg.pack(side='right', expand=1, fill='both')
        #ff.pack(side='top', expand=1, fill='both')

        # ------ the set of sliders -------

        f = Frame(ff, relief='raised', borderwidth=2)
        self.sliderlist = []
        self.slider(f, start=0, end=max_defocus, row=0,
                         label='defocus',
                         tickinterval=10000,
                         resolution = 50,
                         variable = self.defocus)
        self.slider(f, start=0, end=0.005, row=3,
                         label='source\nsize',
                         tickinterval=0.001,
                         digits = 2,
                         variable = self.src)
        self.slider(f, start=0, end=500, row=4,
                         label='defocus\nspread',
                         tickinterval=100,
                         variable = self.spread)
        self.slider(f, start=0, end=2, row=6,
                         label='Gaussian\nenvelope',
                         tickinterval=0.5,resolution=0.01,
                         variable = self.gep)
        f.columnconfigure(2, weight=1) # makes column expandable

        # fy, fg on top; f on bottom
        fy.grid(row=0, column=0, sticky='ns')
        fg.grid(row=0, column=1, sticky='nsew')
        f.grid(row=1, column=0, columnspan=2, sticky='ew')
        ff.columnconfigure(1, weight=1)  # make fg expand
        ff.rowconfigure(1, weight=1)     # make f expand
        
        ff.pack(expand=1, fill='both')
        self.top.update_idletasks()

        if self.askParms == 1:
            self.callSetParms()

        # if called with a list of tfed files
        if len(self.tfed) > 0:
            self.fileSeries(self.tfed)
        if self.defocusfile != "":
            self.openDefocus(self.defocusfile)

        ###### end init ------------------------------------------

    def callSetParms(self):
        w = Toplevel(self.top)
        self.ParmWindow = w
        self.setParms(w)
        self.top.wait_window(w) # wait for window to be destroyed
        self.xupdate()
        self.replot()

    def setParms(self, win=None):
        win.title('Set parameters')
        f = Frame(win)
        labpx = Label(f,text='pixel size(A): ')
        labkv = Label(f,text='electron energy (kev): ')
        labcs = Label(f,text='spherical aberration: ')
        labac = Label(f,text='amplitude contrast ratio: ')
        entpx = Entry(f, textvariable=self.pixsize, width=10, background='white')
        entkv = Entry(f, textvariable=self.kev, width=10, background='white')
        entcs = Entry(f, textvariable=self.cs, width=10, background='white')
        entac = Entry(f, textvariable=self.acr, width=10, background='white')
        labpx.grid(row=0, column=0, sticky='e')
        labkv.grid(row=1, column=0, sticky='e')
        labcs.grid(row=2, column=0, sticky='e')
        labac.grid(row=3, column=0, sticky='e')
        entpx.grid(row=0, column=1)
        entkv.grid(row=1, column=1)
        entcs.grid(row=2, column=1)
        entac.grid(row=3, column=1)
        f.pack(side='top')
        fb = Frame(win)
        b = Button(fb, text='ok', command=win.destroy)
        b.pack(padx=5, pady=5)
        fb.pack()
        self.askParms = 0

    def slider(self, master, start=0, end=10, row=0, label="",
               tickinterval=1, resolution=None, digits=0,variable = None):

        lab = Label(master, text=label)
        if resolution == None:
            resolution = float(tickinterval) / 20.0
        slider = Scale(master, orient='horizontal', from_=start, to=end,
                       tickinterval = tickinterval,
                       resolution = resolution, label ="",
                       variable = variable,
                       showvalue=0,
                       #length = self.g_width,
                       digits = digits,
                       command=self.update)
        self.sliderlist.append(slider)
        ent = Entry(master, textvariable=variable, width=10, background='white')
        ent.bind('<KeyPress>', self.update)

        lab.grid(row=row, column=0, sticky='w', padx=5, pady=5)
        ent.grid(row=row, column=1, sticky='w', padx=5, pady=5)
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

        acr = float(self.acr.get())
        squared = 1
        use_emp_envelope = self.envelope.get()
        if use_emp_envelope:  # use empirical envelope 
            self.modymax.set(1)

        for i in range(self.n):
            ak = i * dk
            p  = ak**3 - dz1 * ak
            ch = exp(ak*4 * kappa)
            self.E[i] = (exp(f*q1*p**2)*ch)*2*exp(-env1*ak**2)
            qqt = 2.0*pi*(0.25*ak**4 - 0.5*dz1*ak**2)
            qqt1 = (1.0-acr)*sin(qqt)-acr*cos(qqt)
            if not use_emp_envelope:
                self.Y[i] = self.E[i] * qqt1
            else:
                self.Y[i] = qqt1
            if squared:
                self.E[i] = (self.E[i])**2
                self.Y[i] = (self.Y[i])**2

        # height of model is always modelheight(0..1) times ymax
        # (unless use_emp_envelope)

        ymax = max(self.Y)
        if use_emp_envelope and len(self.arr['env']) > 0:
            if ymax != 0:
                env = self.arr['env']
                for i in range(self.n):
                    self.Y[i] = (self.Y[i] * env[i]) / ymax
        else:
            if ymax != 0:
                f = float(self.modymax.get()) * float(self.ymax.get())
                self.factor = f / ymax
            for i in range(self.n):
                self.Y[i] = self.Y[i] * self.factor

        if hasattr(self,'g'):
            self.g.element_configure('model', ydata=tuple(self.Y))

    def update(self, scalevalue=None):
        self.compute()
        #for curve in self.curves:
            #self.g.element_configure(curve, ydata=tuple(self.arr[curve]))
        #self.g.element_configure('model', ydata=tuple(self.Y))
        self.top.update_idletasks()

    def replot(self):
        " replots data based on which curves are in display list "
        self.showlist['roo'] = self.showRoo.get()
        self.showlist['bgd'] = self.showBgd.get()
        self.showlist['sub'] = self.showSub.get()
        self.showlist['env'] = self.showEnv.get()
        self.showlist['model'] = self.showMod.get()
        self.newymax()
        self.compute()
        self.showCurves()

    def showCurves(self):
        curves = self.showlist.keys()
        show = []
        for curve in curves:
            if self.showlist[curve]:
                show.append(curve)
        self.g.element_show(show)
        self.top.update_idletasks()

    def newymax(self):
        " compute ymax over all curves except model "
        # get xmin and index into self.X
        xmin = float(self.xmin.get())
        for i in range(self.n):
            if self.X[i] > xmin:
                break
        # i is index
        M = []
        keys = self.showlist.keys()
        for list in keys:
            if list != 'frq' and list != 'model':
                if self.showlist[list] and len(self.arr[list]) > i:
                    newlist = self.arr[list][i:]
                    M.append(max(newlist))
                    #print "%d %s %f" % (i,list,max(self.arr[list]))

        if len(M) == 0: return
        ymax =  max(M)
        if ymax == 0: ymax = 1.0
        #print "ymax %f" % ymax
        self.ymax.set(ymax)

        if hasattr(self,'yslider'):       
            self.yslider.configure(from_= ymax, to=0.0,
                                   tickinterval = ymax/5.0,
                                   resolution = ymax/50.0)
            self.yslider.set(ymax)
            self.g.axis_configure("y", max=ymax, min=0.0)

    def resetYmax(self):
        self.newymax()
        self.compute()
        self.showCurves()

    def xupdate(self, scalevalue=None):
        pixsize = float(self.pixsize.get())
        if pixsize == 0:
            return
        self.max_spat_freq = 1.0 / (2.0 * pixsize)
        for i in range(self.n):
            self.X[i] = i* (self.max_spat_freq / float(self.n))
        for curve in self.curves:
            self.g.element_configure(curve, xdata=tuple(self.X))
        self.g.element_configure('model', xdata=tuple(self.X))
        self.g.axis_configure("x", max=self.max_spat_freq)

    def xminupdate(self, scalevalue=None):
        xmin = float(self.xmin.get())
        if xmin < self.max_spat_freq:
            self.g.axis_configure("x", min=xmin)
        
    def yupdate(self, scalevalue=None):
        ymax = float(self.ymax.get())
        if ymax > 0:
            self.g.axis_configure("y", max=ymax)
            
    #def resetYaxis(self):
        #ymin = ymax = self.Y[0]
        #for i in range(self.n):
         #  if ymin > self.Y[i]: ymin = self.Y[i]
         #   if ymax < self.Y[i]: ymax = self.Y[i]
       #self.g.axis_configure("y", min=ymin, max=ymax)

    def setFileLabels(self):
        if type(self.savefile) == type(""):
            self.savelabel.configure(text=os.path.basename(self.savefile))
        if type(self.tfedfile) == type(""):
            self.tfedlabel.configure(text=os.path.basename(self.tfedfile))
    
    def callOpenFile(self, filename=None):
        if not self.openFile(filename=filename):
            return
        self.xupdate()
        self.replot()
        for curve in ['bgd', 'sub', 'env', 'roo']:
            if len(self.arr[curve]) > 0:
                self.g.element_configure(curve,
                                         ydata=tuple(self.arr[curve]))
        self.setFileLabels()

    def openFile(self, filename=None):
        if filename == None or filename == "":
            filename = askopenfilename()
            if filename == None or filename == "":
                return 0
        A = readdoc(filename, self.multfactor, self.squared.get())
        if len(A) == 0:
            print "openFile: error - fileread returned empty list"
            return 0
        elif len(A) == 1:  # roofile
            self.arr['roo'] = A[0] ; self.showRoo.set(1)
            self.roofile = filename
        elif len(A) == 5:
            self.arr['frq'] = A[0] ; self.X = self.arr['frq']
            self.arr['bgd'] = A[1] #; self.showBgd.set(1)
            self.arr['sub'] = A[2] #; self.showSub.set(1)
            self.arr['env'] = A[3] #; self.showEnv.set(1)
            self.arr['roo'] = A[4] #; self.showRoo.set(1)
            self.tfedfile = filename
        self.n = len(A[0])

        self.E = []
        self.Y = []
        for i in range(self.n):
            self.E.append(0.0)
            self.Y.append(0.0)
        return 1

    def filenumber(self, filename):
        " returns an integer "
        i = string.rfind(filename,'.')
        if i < 0: return -1

        x = i  
        while i > 0:
            i = i-1
            try:
                int(filename[i])
            except:
                break
        try:
            n = int(filename[i+1:x])
        except:
            n = -1
        return n

    def makedisplaylist(self, sort=None):
        if sort == 'defocus' and len(self.defDict) == 0:
            return []
        displaylist = []
        if self.defDict:
            dkeys = self.defDict.keys()
            fkeys = self.fileDict.keys()
            if sort == 'defocus':
                # create a list of (defocus, file) pairs, and sort it
                deflist = []
                for key in dkeys:
                    if key in fkeys:
                        file = self.fileDict[key][0] # base name
                        defocus = self.defDict[key]
                        deflist.append((defocus, file))
                deflist.sort()
                for d in deflist:
                    file = d[1] 
                    defocus = d[0]
                    displaylist.append("%s     %s" % (file, defocus))
            else: # sort == 'file' or None
                for fn in fkeys:
                    if fn in dkeys:
                        file = self.fileDict[fn][0]
                        defocus = self.defDict[fn]
                        displaylist.append("%s     %s" % (file, defocus))
        else:
            # no defocus info, only filenames
            fkeys = self.fileDict.keys()
            fkeys.sort()
            for file in fkeys:
                displaylist.append(self.fileDict[file][0])

        self.displaylist = displaylist
                
        if hasattr(self,'boxwin') and self.boxwin.winfo_exists():
            self.box.setlist(self.displaylist)
            self.boxwin.lift()
            
        return displaylist
            

    def fileSeries(self, flist=None):
        "create a listbox from a user-specified set of files"
        # get the file list
        if flist == None: flist = getFiles()
        if len(flist) < 1: return

        # add basename, fullname to the file dictionary
        self.fileDict = {}
        
        for f in flist:
            basename = os.path.basename(f)
            fn = self.filenumber(basename)
            if fn != -1:
                self.fileDict[fn] = [basename, f]
        self.filelist = flist  # full path

        self.displaylist = self.makedisplaylist()                

        # put the filenames in a scrolled list box
        if hasattr(self,'boxwin') and self.boxwin.winfo_exists():
            self.makedisplaylist()
        else:
            self.boxwin = Toplevel(self.top)
            flabels = Frame(self.boxwin)
            bf = Button(flabels, text='Files',
                        command=lambda self=self,s='files':self.makedisplaylist(sort=s))
            bf.pack(side='left',padx=5,pady=5)
            bd = Button(flabels, text='Defocus',
                        command=self.defocusButtonfunc)
            bd.pack(side='right',padx=5,pady=5)
            flabels.pack(side='top')
            self.box = Pmw.ScrolledListBox(self.boxwin,
                                           items = self.displaylist,
                                           selectioncommand=self.select)
            self.box.pack(side='top', padx=5, pady=5, fill='both', expand=1)
            #d = Button(self.boxwin, text='Read Defocus', command=self.openDefocus)
            #d.pack(side='left', padx=5, pady=5)
            b = Button(self.boxwin, text='Done', command=self.boxwin.destroy)
            b.pack(side='bottom', padx=5, pady=5)
            #self.box.focus_set()


    def select(self):
        sels = self.box.getcurselection()
        if len(sels) < 1: return
        
        f = string.split(sels[0])
        fn = self.filenumber(f[0])  # get number from basename
        filename = self.fileDict[fn][1]  # full path
        if len(f) == 2: # "filename   defocus"
            self.defocus.set(f[1])
        self.callOpenFile(filename=filename)

    def defocusButtonfunc(self):
        " if defocus data loaded, sorts on that column; else loads data"
        if len(self.defDict) == 0:
            self.openDefocus()
        else:
            self.makedisplaylist(sort='defocus')

    def openDefocus(self, filename=None):
        " expects 1st column=mic #, 2nd col=defocus "
        if filename == None:
            filename = askopenfilename(title="Open doc file with defocus values")
        if filename == "" or filename == None or len(filename) == 0:
            return 0
        if not os.path.exists(filename):
            print "Unable to find defocus file: %s" % filename
            return 0

        D = readDefocus(filename) # list of (mic#, defocus) pairs
        if len(D) == 0:
            self.defDict = {}
            return 0
        for item in D:
            micnum = item[0]
            defocus = item[1]
            self.defDict[micnum] = defocus
        self.defocusfile = filename
        self.defocuslabel.configure(text=os.path.basename(filename))

        self.makedisplaylist()
        return 1

    def saveDefocusAs(self):
        filename = asksaveasfilename()
        try:
            os.remove(filename)
        except:
            print "unable to write to %s" % filename
        if filename != "":
            self.saveDefocus(filename=filename)

    def saveDefocus(self, filename=None):
        "save current defocus and file number to a doc file"
        if filename != None:
            self.savefile = filename
        if self.savefile == "":
            filename = asksaveasfilename()
            if filename == "":
                return 
            self.savefile = filename
        self.setFileLabels()

        if self.tfedfile == "":
            print "defocus data will be saved to %s" % self.savefile
            return

        micnum = self.filenumber(os.path.basename(self.tfedfile))
        defocus = int(float(self.defocus.get()))
        outfile = self.savefile
        headers = Spiderutils.getDocfileHeaders(outfile)
        if os.path.exists(outfile):
            # try to replace the line
            d = Spiderutils.readdoc(outfile, keys='all')
            keys = d.keys()
            found = 0
            for k in keys:
                mic = d[k][0]
                if mic == micnum:
                    d[k][1] = defocus
                    found = 1
                    break
            if found:
                Spiderutils.writedoc(outfile,columns=d)
            else:
                Spiderutils.writedoc(outfile,columns=[[micnum],[defocus]],mode='a')
        else:
            Spiderutils.writedoc(outfile,columns=[[micnum],[defocus]],headers=headers)
        #if writedoc(self.savefile, column1=micnum, column2=defocus):
        print "defocus %s saved to %s" % (defocus, outfile)
    
    def showGrid(self):
        self.g.grid_toggle()

    def showSquared(self):
        " variable changed automatically in checkbutton "
        # self.squared.set(self.squared.get())
        if self.tfedfile != "":
            self.callOpenFile(filename=self.tfedfile)

    def quit(self):
        if hasattr(self,'ParmWindow'):
            try:
                self.ParmWindow.destroy()
            except:
                pass
        self.top.quit()

# ------- end CTFplot class definition


def getopts(argv):
    "get command line arguments, return a dictionary"
    opts = {}
    while argv:
        if argv[0][0] == '-':        # find "-name value" pairs
            if argv[0] == "-tfed":
                opts[argv[0]] = argv[1:]
                break  # argv = None
            else:
                opts[argv[0]] = argv[1]  # dict key is "-name" arg
                argv = argv[2:]                    
        else:
            # if there's no flag, assume remainder of args is a tfed file list
            opts["-tfed"] = argv   # argv = argv[1:]
            break
    #print opts
    return opts

def printhelp():
    help = "Usage: ctfmatch.py [-arg value]\n " + \
           "The remaining arguments are (-keyword value) pairs. The following are supported:\n" + \
           "   -cs (spherical aberration)\n" + \
           "   -kev (electron energy)\n" + \
           "   -pixsize (A/pixel)\n" + \
           "   -src (source size)\n" + \
           "   -spread (defocus spread)\n" + \
           "   -acr (amplitude contrast ratio)\n" + \
           "   -gep (Gaussian envelope parameter)\n" + \
           "   -defocus (SPIDER doc file with defocus values)\n" + \
           "   -tfed files* (output from TFED - THIS MUST BE THE LAST ARGUMENT)\n" + \
           "where files* are output doc files from TF ED command with 4 columns:\n" + \
           "   spatial frequency, background, subtracted data, envelope.\n" +\
           "e.g.: ctfmatch.py -cs 2.20 -kev 200 -pixsize 2.82\n" + \
           "e.g.: ctfmatch.py ctf003.dat -pixsize 3.76 -defocus defocus.acn -tfed ctf*\n" + \
           "      ctfmatch.py -help     (prints this message)"
    print help

if __name__ == '__main__':

    nargs = len(sys.argv[1:])

    if nargs > 0:
        h = sys.argv[1]
        if h == '-help' or h == '-h':
            printhelp()
            sys.exit()
        if h[0] != '-':
            argv = sys.argv[1:]
            filename = h
        else:
            argv = sys.argv[1:]
            filename = ""
        args = getopts(argv)
 
    master = Tk()
    master.title("CTF Match")
                 
    if nargs == 0:
        c = CTFplot(master)
    elif nargs > 1 and filename != "":
        c = CTFplot(master, filename=filename, args=args)
    elif nargs > 0 and filename != "":
        c = CTFplot(master, filename=filename)
    elif nargs > 1 and filename == "":
        c = CTFplot(master, args=args)
    else:
        printhelp()
        sys.exit()
    master.mainloop() 
