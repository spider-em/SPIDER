#!/usr/bin/env python

# SOURCE:   ctfgroup.py
# PURPOSE:  Tool for creating defocus groups
#
# Spider Python Library
# Copyright (C) 2006  Health Research Inc. , Menands, NY 
#
# Email:  spider@health.ny.gov
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

import re, os, sys
import Pmw
import pyplot

from   math                import *
from   Tkinter             import *
from   tkMessageBox        import askokcancel
from   tkFileDialog        import askopenfilename, asksaveasfilename
from   commands            import getoutput

from   Spider.Spiderutils  import *

re_nums = re.compile('\d+\D')  # integers followed by one non-int char

"""
def name2template(file):
    " given 'mic021.dat' --> returns mic***.dat "
    if len(file) == 0: return

    m = re_nums.findall(file)  # m = ['021.']
    if len(m) < 1:
        return
    numstr = m[-1]  # if > 1, use last (rightmost) number
        
    stars = "*" * (len(numstr) - 1)
    stars += "."
    a = file.find(numstr)  # find index
    tmp = file.replace(numstr,stars,1)
    return tmp

def template2filename(template, n):
    " given (pic***.dat, 3) --> returns pic003.dat "
    nstars = template.count("*")
    strn = str(n)
    numstr = strn.zfill(nstars)
    sts = "*" * nstars
    filename = template.replace(sts,numstr)
    return filename
"""

# returns "" if user hits 'cancel'
def askfilename(title='Select a file'):
    " if you hit 'cancel', askopenfilename returns an empty tuple "
    f = askopenfilename(title=title)
    if len(f) == 0:
        return ""
    else:
        f = f.replace("/tmp_mnt","")
        return f
    
def asksavefilename():
    f = asksaveasfilename()
    if len(f) == 0:
        return ""
    else:
        f = f.replace("/tmp_mnt","")
        return f

colors = ["#F9C396", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
          "#963939", "#009900", "#990000", "#999999", "#990099", "#999900", "#009999",
          "#CC9360", "#000099", "#F99CF9", "#909000", "#FF6F00", "#C9CFC9", "#FF6F60"]

def getcolor(i):
    if type(i) == type("string"):
        i = int(i)
    if i == 0:
        return "white"
    n = len(colors)
    i = i%n
    return colors[i]

# ########################################################################
#
#
#

class scatter:
    def __init__(self, master, filename=None, xcol=None, ycol=None, plotfilelist=None):
        self.top = master

        self.xcolumn = IntVar()
        self.ycolumn = IntVar()
        self.xaxisVar = StringVar()
        self.xaxisVar.set("defocus")
        self.yaxisVar = StringVar()
        self.yaxisVar.set("micrographs")
        self.filename = ""
        self.D = {}
        self.vx = ()
        self.vy = ()
        self.x0 = 0   # coordinates for mouse rectangle
        self.y0 = 0
        self.x1 = 0
        self.y1 = 0
        self.dragging = 0
        self.eVar1 = StringVar()
        self.emptycolor = getcolor(0)
        self.framecolor = '#aaaaee'
        self.gmax = 0
        self.shapes = ["numbers","circle","diamond","square","triangle"]
        self.symbolcolor = 'white'
        self.activecolor = '#ff0000'
        self.activeshape = 'diamond'   #'plus'
        #self.activegroup = ''
        self.symbolshape = 'circle'
        self.symbolsizes = [15,13,11,9,7,5]
        self.symbolsize = 11
        self.showDefocusSpread = 0

        #   D : keys are 1..n,
        #    each element: ( x, y, groupno, colorstring, mic )

        self.filedata = []   # list of lines from the file [key, micro, defocus, defgrp]

        if filename != None:
            fn, self.ext = os.path.splitext(filename)
            # sets self.D, self.vx, self.vy
            res = self.getdata(filename)
            if res == -1:
                sys.exit()

        keys = self.D.keys()
        keys.sort()

        self.Dkeys = self.D.keys()
        self.Dkeys.sort()
        self.nkeys = len(self.Dkeys)

        self.oldvalue = StringVar()
        self.oldvalue.set("")

        # plotting
        if plotfilelist != None and len(plotfilelist) > 0:
            self.plotfiletemplate = name2template(plotfilelist[0])
        else:
            self.plotfiletemplate = "power/ctf***" + self.ext
        self.plotfile = StringVar()
        self.plotfile.set(self.plotfiletemplate)
        self.gnuplot = ""

        # ------- widgets start here -------

        self.makeMenus(master)

        fg = Frame(master)
        fg.pack(side='top', fill='both', expand=1)
        
        #self.g = Pmw.Blt.Graph(fg, plotbackground="black", takefocus=1) 
        self.g = Pmw.Blt.Graph(fg, takefocus=1)
        self.g.grid_on()
        # use an invisible line to set the plot window size
        self.g.line_create("scatter", xdata=self.vx, ydata=self.vy,
                           linewidth=0, symbol='')
        self.g.legend_configure(hide=1)
        self.g.axis_configure("x", title=self.xaxisVar.get())
        self.g.axis_configure("y", title=self.yaxisVar.get())

        self.g.bind(sequence="<ButtonPress-1>",   func=self.mouseDown)
        self.g.bind(sequence="<ButtonRelease-1>", func=self.mouseUp  )
        self.g.pack(side='top', fill='both', expand=1)

        # put micrograph & defocus entries in bottom of main window
        fe = Frame(fg)
        fe.pack(side='bottom', fill='x')
        self.micEntry = Pmw.EntryField(fe, labelpos='w', label_text='micrograph')
        self.micEntry.component('entry').configure(width=5)
        self.micEntry.pack(side='left', padx=4)
        self.defEntry = Pmw.EntryField(fe, labelpos='w', label_text='defocus')
        self.defEntry.component('entry').configure(width=5)
        self.defEntry.pack(side='right', padx=4)
       
        # bottom frame
        fb = Frame(master, relief='raised', borderwidth=2,
                   background=self.framecolor)
        fb.pack(side='top', fill='x')
        
        # entry for plotfiles
        self.plotEntry = Entry(fb, textvariable=self.plotfile)
        self.plotEntry.pack(side='right', pady=4)
        bpf = Button(fb, text="plot files:", command=self.getplotfiles)
        bpf.pack(side='right', pady=4) 

        self.putsymbols()

       
        # ------- end init -------

    def makeMenus(self,master):
        brelief = 'flat'
        self.mBar = Frame(master, relief='raised', borderwidth=2,
                          background=self.framecolor)        
        self.mBar.pack(side='top',fill = 'x')

        # Make the File menu 
        File_button = Menubutton(self.mBar, text='File', underline=0,
                                 relief=brelief, background=self.framecolor)
        File_button.pack(side='left', padx=5, pady=5)
        File_button.menu = Menu(File_button, tearoff=0)
        File_button['menu'] = File_button.menu

        File_button.menu.add_command(label='Save as', underline=0,
                                     command=self.savefile)
        File_button.menu.add_command(label='Show defocus spread', underline=0,
                                     command=self.getParams)
        
        File_button.menu.add_separator()
        File_button.menu.add_command(label='Quit', underline=0,
                                     command=master.quit)
        # Symbol menu
        Sym_button = Menubutton(self.mBar, text='Symbols', underline=0,
                                 relief=brelief, background=self.framecolor)
        Sym_button.pack(side='left', padx=5, pady=5)
        Sym_button.menu = Menu(Sym_button, tearoff=0)
        Sym_button['menu'] = Sym_button.menu

        # the Symbol --> Shape submenu
        Sym_button.menu.shapes = Menu(Sym_button.menu, tearoff=0)
        self.shapeVar = IntVar()
        if self.symbolshape in self.shapes:
            snum = self.shapes.index(self.symbolshape)
            self.shapeVar.set(snum)
        else:
            self.shapeVar.set(1)
        num = 0
        for shape in self.shapes:
            Sym_button.menu.shapes.add_radiobutton(label=shape, value=num,
                                                   variable=self.shapeVar,
                                                   command=self.getsymbolshape)
            num += 1     
        Sym_button.menu.add_cascade(label='Shapes', menu=Sym_button.menu.shapes)

        # the Symbol --> Size submenu
        Sym_button.menu.sizes = Menu(Sym_button.menu, tearoff=0)
        self.sizeVar = IntVar()
        if self.symbolsize in self.symbolsizes:
            snum = self.symbolsizes.index(self.symbolsize)
            self.sizeVar.set(snum)
        else:
            self.sizeVar.set(1)
        num = 0
        for size in self.symbolsizes:
            Sym_button.menu.sizes.add_radiobutton(label=str(size), value=num,
                                                   variable=self.sizeVar,
                                                   command=self.getsymbolshape)
            num += 1     
        Sym_button.menu.add_cascade(label='Sizes', menu=Sym_button.menu.sizes)
        
        Sym_button.menu.add_command(label='Clear', underline=0,
                                    command=self.cleargroups)

        # put help menu on the right
        Help_button = Menubutton(self.mBar, text='Help', underline=0,
                                 relief=brelief, background=self.framecolor)
        Help_button.pack(side='right', padx=5, pady=5)
        Help_button.menu = Menu(Help_button, tearoff=0)
        Help_button['menu'] = Help_button.menu

        Help_button.menu.add_command(label='Help', underline=0,
                                     command=self.help)

    def cleargroups(self):
        keys = self.D.keys()
        keys.sort()
        # D[k] = (x, y, group, color, mic)
        for k in keys:
            self.D[k][2] = 0
            self.D[k][3] = "white"
        self.putsymbols()

    def getsymbolshape(self):
        self.symbolshape = self.shapes[self.shapeVar.get()]
        self.symbolsize = self.symbolsizes[self.sizeVar.get()]
        self.putsymbols()
        
    def putsymbols(self):
        if self.symbolshape == 'numbers':
            self.puttextsymbols()
            return
        keys = self.D.keys()
        for k in keys:
            mk = str(k)
            currentshape = self.symbolshape,
            if self.g.element_exists(mk):
                currentshape = self.g.element_cget(mk, 'symbol')
                self.g.element_delete(mk)
            if self.g.marker_exists(mk):
                self.g.marker_delete(mk)
            x = (self.D[k][0])
            y = (self.D[k][1])
            group = (self.D[k][2])
            color = (self.D[k][3])
            symbol = self.symbolshape
            if currentshape == self.activeshape:
                symbol = currentshape
            else:
                symbol = self.symbolshape
                
            self.g.line_create(mk, xdata=x, ydata=y,
                               symbol = symbol,
                               pixels = self.symbolsize,
                               fill=color,
                               outline="black")
            
            self.g.element_bind(mk,"<Button-3>",
                               func=lambda event, g=group: self.plotgroup(g))
            self.g.element_bind(mk,"<Enter>",
                               func=lambda event, m=mk: self.entermarker(m))
            self.g.element_bind(mk,"<Leave>",
                               func=lambda event, m=mk: self.leavemarker(m))
            """
            self.g.element_bind(mk,"<Double-Button-1>",
                               func=lambda event, m=mk: self.selectmarker(m))
            """

    def puttextsymbols(self):
        keys = self.D.keys()
        for k in keys:
            mk = str(k)
            #currentshape = self.symbolshape,
            if self.g.marker_exists(mk):
                #currentshape = self.g.marker_cget(mk, 'symbol')
                self.g.marker_delete(mk)
            if self.g.element_exists(mk):
                self.g.element_delete(mk)
            x = (self.D[k][0])
            y = (self.D[k][1])
            group = (self.D[k][2])
            color = (self.D[k][3])
            mic = (self.D[k][4])

            self.g.marker_create("text", name = mk)
            if color != self.emptycolor:
                self.g.marker_configure(mk, 
                                        coords = (x, y),
                                        text = mic,
                                        foreground=color,
                                        shadow = "black")
            else:
                self.g.marker_configure(mk, 
                                        coords = (x, y),
                                        text = mic,
                                        foreground=color,
                                        shadow = "black")
                #print self.g.marker_cget(mk, 'font')

            self.g.marker_bind(mk,"<Button-3>",
                               func=lambda event, g=group: self.plotgroup(g))
            self.g.marker_bind(mk,"<Enter>",
                               func=lambda event, m=mk: self.entermarker(m))
            self.g.marker_bind(mk,"<Leave>",
                               func=lambda event, m=mk: self.leavemarker(m))
        
    def entermarker(self, mk):
        " numbers are markers. all other symbols are elements"
        mic = self.D[int(mk)][4]
        defocus = self.D[int(mk)][0]
        if self.g.element_exists(mk):
            # set the micrograph number in the entry box
            self.micEntry.setentry(mic)
            self.defEntry.setentry(defocus)
            self.g.element_configure(mk, fill=self.activecolor)   
        elif self.g.marker_exists(mk):
            self.micEntry.setentry(mic)
            self.defEntry.setentry(defocus)
            self.g.marker_configure(mk, outline=self.activecolor) #foreground=self.activecolor)
        else:
            print "element %s not found" % str(mic)
            return

    def leavemarker(self, mk):
        color = self.D[int(mk)][3]
        if self.g.element_exists(mk):
            self.g.element_configure(mk, fill=color)
        elif self.g.marker_exists(mk):
            self.g.marker_configure(mk, outline=color)

    def showactivemarkers(self, markers):
        keys = self.D.keys()
        for k in keys:
            mk = str(k)
            color = self.D[k][3]
            
            if mk in markers:   # change properties of the active markers
                if self.g.element_exists(mk):
                    self.g.element_configure(mk, symbol=self.activeshape) 
                elif self.g.marker_exists(mk):
                    self.g.marker_configure(mk, foreground=self.activecolor)

            else:   # reset properties of inactive markers
                if self.g.element_exists(mk):
                    self.g.element_configure(mk, symbol=self.symbolshape)
                elif self.g.marker_exists(mk):
                    self.g.marker_configure(mk, foreground=color)
 

    def assigngroups(self):
        #print self.D
        groups = []
        n = self.nkeys
        for i in range(n):
            k = self.Dkeys[i]
            g = self.D[k][2]
            groups.append(g)
                    
        #print "groupin %s" % str(groups)
        self.gmax = max(groups) + 1

        for i in range(n):
            if groups[i] < 0:
                groups[i] = self.gmax

        #print "groupout %s" % str(groups)
        for i in range(n):
            k = self.Dkeys[i]
            grp = groups[i]
            self.D[k][2] = grp
            self.D[k][3] = getcolor(grp)
            
        self.putsymbols()

    def plotgroup(self,groupno):
        #print "plotting group %d" % groupno
        if groupno == 0:
            return
        self.plotfiletemplate = self.plotfile.get()

        n = self.nkeys
        group = []
        activemarkers = []
        defoci = []
        for i in range(n):
            k = self.Dkeys[i]
            item = self.D[k]
            if item[2] == groupno:
                mic = self.D[k][4]
                fn = makeSpiderFilename(self.plotfiletemplate, mic)
                group.append(fn)
                mk = str(k)
                
                activemarkers.append(mk)
                defocus = self.D[k][0]
                defoci.append(defocus)

        hasplot = 0
        if hasattr(self.gnuplot, 'topwindow') and self.gnuplot.topwindow.winfo_exists():
            hasplot = 1

        if not hasplot:
            self.gnuplotwin = Toplevel(self.top)
            self.gnuplot = pyplot.MyGnuplotPlot(self.gnuplotwin,group)
            self.gnuplot.set_columns("1:5")
            self.gnuplot.selectAll()
            self.gnuplot.plotSelected()
        else:
            self.gnuplot.replace(group)

        self.showactivemarkers(activemarkers)
        if self.showDefocusSpread:
            self.displayDefocusSpread(defoci)

    def getplotfiles(self):
        file = askfilename()
        if file == "":
            return ""
        # shorten the name to below the current working directory
        cwd = os.getcwd() + "/"
        cwd = cwd.replace("/tmp_mnt","")
        file = file.replace(cwd,"")
        
        filetmp = name2template(file)
        self.plotfile.set(filetmp)

    # Mouse functions for selection -----------

    def mouseDrag(self, event):
        (self.x1, self.y1) = self.g.invtransform(event.x, event.y)
             
        self.g.marker_configure("marking rectangle", 
            coords = (self.x0, self.y0, self.x1, self.y0, self.x1, self.y1, self.x0, self.y1, self.x0, self.y0))
    
    def mouseUp(self, event):
        x0 = self.x0
        y0 = self.y0
        x1 = self.x1
        y1 = self.y1
        
        if self.dragging:
            self.g.unbind(sequence="<Motion>")
            self.g.marker_delete("marking rectangle")
            
            if x0 <> x1 and y0 <> y1:

                # make sure the coordinates are sorted
                if x0 > x1: x0, x1 = x1, x0
                if y0 > y1: y0, y1 = y1, y0

                x0 = int(x0)
                y0 = int(y0) + 1
                x1 = int(x1)
                y1 = int(y1)

                group = []

                for y in range(y0,y1+1):
                    if self.D.has_key(y):
                        defocus = int(self.D[y][0])
                        if defocus >= x0 and defocus <= x1:
                            group.append( [y, defocus] )
                            self.D[y][2] = -1  # assign to new group

                #print "%d %d %d %d --> %s" % (x0,y0,x1,y1, str(group))
                self.assigngroups()
                       
    def mouseDown(self, event):
        self.dragging = 0
        
        if self.g.inside(event.x, event.y):
            self.dragging = 1
            (self.x0, self.y0) = self.g.invtransform(event.x, event.y)
            
            self.g.marker_create("line", name="marking rectangle", dashes=(2, 2))
            self.g.bind(sequence="<Motion>",  func=self.mouseDrag)
            
    # File I/O functions -----------

    def getdata(self, filename, xcol=None, ycol=None):
        # expects doc file with columns: MICROGRAPH  DEFOCUS [DEF.GROUP]
        # i.e., if it has group column, use it
        if filename == None: return -1
        if not os.path.exists(filename):
            print "unable to find %s" % filename
            return -1

        n_cols = numberOfColumns(filename)
        if n_cols < 2:
            print "doc files must have at least 2 columns: micrographs, defocus"
            return -1

        K,M,D = readdoc(filename, column=[0,1,2])  # keys,micrographs,defocus

        maxdiff = 800  # maximum range of defocus in a group

        if n_cols == 2:
            # need to sort defocus, generate defocus group numbers
            B = def_sort(M,D,maxdiff)
            
        elif n_cols > 2:
            G = readdoc(filename, column=3)  # group
            
            if not self.isDefGroup(G):
                # 'group' numbers not valid, call def_sort
                B = def_sort(M,D,maxdiff)
            else:
                # make B = [ [mic, def, grp], ..]
                nmic = len(M)
                B =[]
                for i in range(nmic):
                    B.append( [M[i], D[i], G[i]] )
                              
        # convert B lines [mic,def,grp] to filedata lines [def key mic grp]
        n = len(B)
        for i in range(n):
            X = B[i]  # each element is a list
            self.filedata.append( [ X[1], i+1, X[0], X[2] ] )
            
        self.filedata.sort()
        #for i in self.filedata: print i

        # put defocus in x-axis, key in y-axis
        X = []
        Y = []
        for i in range(n):
            entry = self.filedata[i]
            x = entry[0]
            y = entry[1]
            mic = entry[2]
            if len(entry) > 3:
                group = int(entry[3])
            else:
                group = 0
            color = getcolor(group)
            if group > self.gmax:
                self.gmax = group
            if self.D.has_key(y):
                print "Error: Y axis numbers must be unique: %s" % y
            # use the yaxis as the keys. NUMBERS MUST BE UNIQUE
            self.D[y] = [x, y, group, color, mic] 
            X.append(x)
            Y.append(y)

        keys = self.D.keys()
        if len(keys) < 1:
            print "unable to get data from %s" % filename
            return -1

        extra_x = int((max(X) - min(X)) * 0.05)
        extra_y = int((max(Y) - min(Y)) * 0.05)
        
        self.vx = (min(X)-extra_x, max(X)+extra_x)
        self.vy = (min(Y)-extra_y, max(Y)+extra_y)
        return 0

    def isDefGroup(self,G):
        " defocus group numbers must be positive ints <= len(G) "
        n = len(G)
        for i in G:
            if i > 0 and isInt(i) and i <= n:
                pass
            else:
                return 0
        return 1

    def remap(self,g):
        " sorts the values in g; zero is a special case"
        n = len(g)
        val = 1
        oldg = g[0]
        if g[0] == 0:
            newg = [0]
        else:
            newg = [1]
        
        for i in range(1,n):
            if g[i] == 0:
                newg.append(0)
            else:
                if g[i] != oldg and max(newg) != 0:
                    val += 1
                newg.append(val)
            oldg = g[i]
            
        return newg

    def savefile(self):
        hdrstring = " ; /     MICROGRAPH   DEFOCUS    DEF.GROUP\n"
        groups = []
        keys = self.D.keys()
        keys.sort()
        for k in keys:
            groups.append(int(self.D[k][2]))
            
        if min(groups) < 1:
            if not askokcancel('Warning', 'Not all micrographs have been assigned!'):
                return
        outfile = asksavefilename()
        if len(outfile) == 0:
            return

        lines = [hdrstring]

        newgroups = self.remap(groups)
        i = 0
            
        for k in keys:
            defocus = self.D[k][0]
            group   = int(self.D[k][2])
            mic     = int(self.D[k][4])
            n = newgroups[i]
            i += 1
            #print " %3d 3   %6d       %s       %6d %6d" % (k,mic,defocus,group,n)
            line = " %3d 3   %6d       %s       %6d\n" % (k,mic,defocus,n)
            lines.append(line)

        fp = open(outfile, 'w')
        fp.writelines(lines)
        fp.close()

    def openfile(self):
        pass

    def help(self):
        w = Toplevel(self.top)
        s1 = "Drag Left mouse button to group \ntogether points of similar defocus"
        s2 = "Right mouse button: plot all files in a group"
        #s3 = "file://localhost/net/bali/usr1/spider/docs/spire/guitools/defgroup.html"
        #e = StringVar()
        #e.set(s3)
        Label(w, text=s1).pack(side='top', anchor='w', fill='x',padx=2,pady=2)
        Label(w, text=s2).pack(side='top', anchor='w', fill='x',padx=2,pady=2)
        #Entry(w, textvariable=e).pack(side='top', anchor='w', fill='x')
        Button(w, text='Ok', command=w.destroy).pack(side='bottom', pady=5)

    # stuff added for showing the defocus spread -----------------------

    def displayDefocusSpread(self, defoci):
        
        self.computeEnvelope(defoci)
        Env = self.Envelope
        # rescale to graphic
        wd = self.ndefspread_points
        ht = int(wd/2.0)
        emin = 0  # assume squared min(Env)
        emax = max(Env)
        m = ht / (emax - emin)
        b = ht - (m*emax)
        F = []
        for x in Env:
            y = int(m*x + b)
            y = ht - y  # flip it vertically for Canvas display
            F.append(y)

        if not hasattr(self, 'dfspread'):       
            self.dfspread = Canvas(self.g, width=wd, height=ht, background='white')
            self.g.marker_create("window", name="dfspreadWindow")
            #cx,cy = self.g.extents("leftmargin"), self.g.extents("topmargin")
            self.g.marker_configure("dfspreadWindow",
                                    coords = (self.vx[0], self.vy[1]),
                                    window = self.dfspread)
            #self.dfspread.pack(side='top', anchor='nw')
            self.canvasitems = []

        for item in self.canvasitems:
            self.dfspread.delete(item)
            
        for i in range(wd-1):
            x1 = i
            x2 = i+1
            y1 = F[x1]
            y2 = F[x2]
            item = self.dfspread.create_line(x1,y1,x2,y2)
            self.canvasitems.append(item)

    def computeEnvelope(self, defoci):
        " put values in self.Envelope "
        defspread = max(defoci) - min(defoci)
        sum = 0
        for defocus in defoci:
            sum += defocus
        defocus = sum / float(len(defoci))
        
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
        dk    = km1 / float(self.ndefspread_points)
        q1    = (float(self.srcsize.get()) * f2)**2
        gep = float(self.gep.get())
        if gep != 0:
            env   = 1.0/gep**2
        else:
            env = self.infinity
        env1  = env/f2**2
        f     = -pi**2
        ds1   = f1 * float(defspread)
        kappa = ds1 * self.kappa
        dz1   = f1 * float(defocus)

        acr = float(self.acr.get())
        squared = 1  # self.squared.get()

        for i in range(self.ndefspread_points):
            ak = i * dk
            p  = ak**3 - dz1 * ak
            ch = exp(ak**4 * kappa)
            self.Envelope[i] = (exp(f*q1*p**2)*ch)*2*exp(-env1*ak**2)
            #qqt = 2.0*pi*(0.25*ak**4 - 0.5*dz1*ak**2)
            #qqt1 = (1.0-acr)*sin(qqt)-acr*cos(qqt)
            #self.Y[i] = self.E[i] * qqt1
            if squared:
                self.Envelope[i] = (self.Envelope[i])**2
                #self.Y[i] = (self.Y[i])**2
        

    def getParams(self):
        # present a form to the user to get ctf parameters
        self.cs = StringVar()
        self.cs.set('2.0')
        self.kev = StringVar()
        self.kev.set('200')
        self.pixsize = StringVar()
        self.pixsize.set('2.82')
        self.srcsize = StringVar()
        self.srcsize.set('0.0')
        self.acr = StringVar()
        self.acr.set('0.1')
        self.gep = StringVar()
        self.gep.set('1000')
        self.infinity = 1e50
        self.kappa = -pi**2 / (16.0 * log(2.0))
        self.ndefspread_points = 250
        self.Envelope = []
        for i in range(self.ndefspread_points):
            self.Envelope.append(0)
        self.showDefocusSpread = 1

        pwin = Toplevel(self.top)
        fr1 = Frame(pwin)
        cslabel = Label(fr1,text="Cs:")
        csentry = Entry(fr1, textvariable=self.cs, background='white')
        kvlabel = Label(fr1,text="keV:")
        kventry = Entry(fr1, textvariable=self.kev, background='white')
        pixlabel = Label(fr1,text="pixel size:")
        pixentry = Entry(fr1, textvariable=self.pixsize, background='white')
        srclabel = Label(fr1,text="source size:")
        srcentry = Entry(fr1, textvariable=self.srcsize, background='white')
        acrlabel = Label(fr1,text="amplitude contrast ratio:")
        acrentry = Entry(fr1, textvariable=self.acr, background='white')

        cslabel.grid(row=0, column=0, sticky='e')
        csentry.grid(row=0, column=1) #, sticky='e')
        kvlabel.grid(row=1, column=0, sticky='e')
        kventry.grid(row=1, column=1) #, sticky='e')
        pixlabel.grid(row=2, column=0, sticky='e')
        pixentry.grid(row=2, column=1) #, sticky='e')
        srclabel.grid(row=3, column=0, sticky='e')
        srcentry.grid(row=3, column=1) #, sticky='e')
        acrlabel.grid(row=4, column=0, sticky='e')
        acrentry.grid(row=4, column=1) #, sticky='e')
        fr1.pack(side='top')

        fr2 = Frame(pwin)
        Button(fr2, text='Done', command=pwin.destroy).pack()
        fr2.pack(side='top')

    
def printHelp():
    print "Usage: ctfgroup.py -d doc_sort.dat ctfplotfiles*"

# sort the defocus data   
def def_sort(Micro, Defocus, maxdiff=1000):
    " def_sort takes COLUMNS and returns LINES (sorted by defocus)"
    " inputs are lists of micrographs and defocus values."
    " output : a list of lines (sorted by defocus): micrograph, defocus, group"

    n = len(Micro)

    N = []              # create a new list, of (defocus, micrograph) pairs
    for i in range(n):
        N.append( (Defocus[i],Micro[i]) )
    N.sort()           # sorts by defocus, since that's the first element

    # create the output, [micrograph, defocus, def_group]
    prevmic = N[0][1]
    prevdef = N[0][0]
    group = 1    # current def_group
    G = [ (prevmic, prevdef, group) ]

    for i in range(1,n):
        thismic = N[i][1]
        thisdef = N[i][0]
        if (thisdef-prevdef) > maxdiff:
            group = group + 1  # incr. grp no. if difference exceeds maxdiff
            prevdef = thisdef
        G.append( (thismic, thisdef, group) )

    return G

# ----------------------------------------------------------
# input doc file must have micrographs and defocus in 1st, 2nd columns.
# Sorts them with the default 1000 max diff.
# Use -d to indicate defsort file; remainder are plot files.
# (If -d flag not used, it assumes 1st file is def_sort file)
#
# Usage: ctfgroup.py -d def_sort.dat ctf*

if __name__ == '__main__':
    master = Tk()
    args = sys.argv[1:]
    numargs = len(args)
    smallfont = 0
    defsortfile = ""
    plotfiles = []

    if numargs == 0:
        defsortfile = askfilename(title='Select a def_sort doc file')
        if defsortfile == "":
            printHelp()
            sys.exit()
        
    elif numargs > 0:
        if args[0] == "-h" or args[0] == "-help":
            printHelp()
            sys.exit()
            
        while args:
            if args[0] == '-s':
                smallfont = 1
                args = args[1:]
            elif args[0] == '-d':
                defsortfile = args[1]
                args = args[2:]
            else:
                plotfiles.append(args[0])
                args = args[1:]

                
    if defsortfile == "":
        # Use the first plot file as the def_sort file
        defsortfile = plotfiles[0]
        plotfiles = plotfiles[1:]
        if not os.path.exists(defsortfile):
            print "cannot find defsort file %s" % defsortfile
            sys.exit()

    if not smallfont:
        Pmw.initialise(master, fontScheme = 'pmw1')
    master.title("Select defocus groups")

    c = scatter(master, filename=defsortfile, plotfilelist=plotfiles)
    master.mainloop() 
