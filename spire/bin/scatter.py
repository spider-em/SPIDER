#!/usr/bin/env python

# SOURCE:  scatter.py  
# PURPOSE: Scatterplot display
#
# Spider Python Library
# Copyright (C) 2006  Health Research Inc., Menands, NY
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

import                os,string,sys
import                re
from commands         import getoutput

from Tkinter          import  *
import                        Pmw, tkMessageBox
from tkFileDialog     import  askopenfilename, asksaveasfilename
from PIL              import  Image, ImageTk, ImageChops
from Spider           import  Spiderutils

import webbrowser
webpage = "http://www.wadsworth.org/spider_doc/spider/spire/tools-docs/scatter.html"

# Returns "" if user hits 'cancel'
def askfilename():
    " If you hit 'cancel', askopenfilename returns an empty tuple "
    f = askopenfilename()
    if len(f) == 0:
        return ""
    else:
        f = string.replace(f,"/tmp_mnt","")
        return f

shapes = ["square","circle","diamond","triangle","plus","cross","numbers"]

class polygonLines:
    def __init__(self):
        self.coords = []
        self.names = []  # marker names
        self.pccords = []
    def clear(self):
        self.coords = []
        self.names = []
    def pairs2oneList(self, coords=None):
        " ( (1,2), (3,4), (5,6) ) --> (1,2,3,4,5,6) "
        if coords == None:
            coords = self.coords
        n = ()
        for c in coords:
            n = n + c
        return n
    def makecopy(self):
        c = []
        n = []
        for cc in self.coords:
            c.append(cc)
        for nom in self.names:
            n.append(nom)
        return c,n

class Polygon:
    def __init__(self, coords=None):
        if coords == None:
            self.coords = []
        else:
            self.coords = coords
        self.elements = []  # keys for data elements

###############################################################

class scatter:
    """
    D : dictionary of (float) data points,
        dictionary keys are keys from doc file (as strings)
    PD : same, but with integer values to corresponding screen pixels.
    self.pixlist: list of pixel values sorted by x coordinate (x,y,key)
    self.vx,vy: 2 points of min, max data values for initializing graph size

    self.g : main graph object
    self.gpoint : graph with 1-pixel points, used as a mask
    self.gpoly : graph with filled polygon region, used as a mask
    self.impoint, self.impoly : the graphs as Image.images

    ------------------------------------------
    Finding points inside the selected polygon:
    - make 2 copies of the graph, one for data points, the other for
      the filled polygon.
    - use snap to create binary image masks,
    - multiply polygon mask x data point mask,
    - get pixel coords from result, compare to sorted pixelist

    Data point pixels -> generate the image and pixelist whenever a
    new data file is read in. Graph can be discarded.
    Need to generate polygon graph, image as needed, since graph objects
    cannot be read unless they are mapped to the display (i.e. packed).
    
    #  putsymbols: draws symbols in the main graph
    #  mkDataMask: makes graph self.gpoint
    #  putpoints: put data points in self.gpoint (called in mkGraphCopies)
    #  mkPixelist: make sorted list of pixels 

    """
    def __init__(self, master, filename=None, xcol=None, ycol=None):
        self.top = master

        self.xcolumn = IntVar()
        self.ycolumn = IntVar()
        self.xaxisVar = StringVar()
        self.yaxisVar = StringVar()
        self.filename = ""
        self.D = {}
        self.vx = ()
        self.vy = ()

        if filename != None and xcol != None and ycol != None:
            self.D, self.vx, self.vy = self.getdata(filename, xcol, ycol)
        else:
            self.openfile(filename)

        self.state = "None"   # (can be 'Magnify' | 'Polygon' | 'None' )
        self.activestatecolor = "#60E090"
        self.activestatecolor2 = "#70F0A0"

        self.symbolcolor = 'green'
        self.activecolor = 'red'
        self.symbolsize = 9
        self.symbolshape = 'circle'
        self.eVar1 = StringVar()
        self.eVar2 = StringVar()
        self.eVar3 = StringVar()
        
        self.keyvalue = StringVar()  # entry in lower left
        self.keyvalue.set("")
        self.oldvalue = StringVar()
        self.oldvalue.set("")
        self.command = StringVar()
        ###self.command.set("qview.py file***.dat")  July 2018 al
        self.command.set("qview file***.dat")

        # magnifying variables
        self.magnifyVar = IntVar() ; self.magnifyVar.set(0)
        self.polygonVar = IntVar() ; self.polygonVar.set(0)
        self.dragging = 0
        self.x0=0; self.x1=0; self.y0=0; self.y1=0

        # polygon selection variables
        self.poly = polygonLines()
        self.Polygons = {}  # dictionary of polygons
        self.newpoly = 1
        self.uniqname = "1"
        self.selectedPolygon = ""
        self.selectedPolygonColor = "yellow"
        self.PolygonColor = "lightblue"
        self.polyElementsColor = "yellow"
        self.magicnumber= 5  # whooo...but see plotborderwidth,plotpadxy
        self.graph_wd = 0
        self.graph_ht = 0

        # ------- widgets start here -------

        self.mkmenu(master)
        
        fg = Frame(master)
        fg.pack(side='top', fill='both', expand=1)
        
        self.g = Pmw.Blt.Graph(fg, plotbackground="black", takefocus=1) 
        self.g.pack(fill='both', expand=1)
        # use an invisible line to set the plot window size
        self.g.line_create("scatter", xdata=self.vx, ydata=self.vy,
                           linewidth=0, symbol='')
        self.g.legend_configure(hide=1)
        self.g.axis_configure("x", title=self.xaxisVar.get())
        self.g.axis_configure("y", title=self.yaxisVar.get())

        self.mkbottom(master)

        self.plotTheData()

        #self.g.bind("<ButtonPress>", self.transformer)

        # ------- end init -------

    def mkmenu(self,master):
        menuBar = Pmw.MenuBar(master, hull_relief = 'raised',
                              hull_borderwidth = 1)        
        menuBar.pack(fill = 'x')

        # Make the File menu 
        menuBar.addmenu('File', 'helptxt')
        menuBar.addmenuitem('File', 'command', 'open file',
                             label='Open', command=self.openfile)
        menuBar.addmenuitem('File', 'separator')
        menuBar.addmenuitem('File', 'command', 'quit',
                             label='Quit', command=master.quit)
        # Make the Option menu
        menuBar.addmenu('Options', '')
        menuBar.addmenuitem('Options', 'checkbutton', 'turn zoom on/off',
                            label='Magnify', variable=self.magnifyVar,
                            command=self.magnifyOn)
        menuBar.addmenuitem('Options', 'command', 'Reset magnification',
                            label='Reset view', command=self.magReset)
        menuBar.addmenuitem('Options', 'command', 'change symbols',
                            label='Symbols', command=self.symbolsMenu)  
        # Make the Region menu
        menuBar.addmenu('Regions', '')
        menuBar.addmenuitem('Regions', 'checkbutton', 'polygon draw',
                            label='Polygon select', variable=self.polygonVar,
                            command=self.polygonOn)
        menuBar.addmenuitem('Regions', 'command', 'save selected area',
                            label='Save selected polygon', command=self.savePolygon)
        menuBar.addmenuitem('Regions', 'command', 'delete selected area',
                            label='Delete selected polygon', command=self.deletePolygon)

        # Magnify, Region Buttons
        self.magButton = Button(menuBar.component('hull'),
                                text="Magnifier", command=self.magnifyOn)
        self.magButton.pack(side='left',padx=10, pady=10)
        self.magButton.bind('<Double-Button-1>', self.magReset)
        self.polyButton = Button(menuBar.component('hull'),
                                text="Region select", command=self.polygonOn)
        self.polyButton.pack(side='left',padx=10, pady=10)
        
        self.backgroundcolor = self.magButton.cget('background')
        self.activebackgroundcolor = self.magButton.cget('activebackground')

        # Make the Help menu
        menuBar.addmenu('Help', '', side='right')
        menuBar.addmenuitem('Help', 'command', 'view help page',
                            label='Help', command=self.help)
        menuBar.addmenuitem('Help', 'command', 'about this program',
                            label='About', command=self.about)

    def mkbottom(self, master):
        fb = Frame(master)

        # Pack from right to left
        cmdentry = Entry(fb, textvariable=self.command, background='white',
                         width = 40)
        cmdentry.pack(side='right', padx=5, pady=10)
        cmdentry.bind('<Return>', self.callprog_check)
        cmdlabel = Label(fb,text="       Unix command: ")
        cmdlabel.pack(side='right',pady=10,padx=1)
        
        keyentry = Entry(fb, textvariable=self.keyvalue, width=4,
                         background='white')
        keyentry.pack(side='right', padx=5, pady=10)
        keyentry.bind('<Return>', self.callprogram)
        keyentry.bind('<Leave>', self.keyleave)
        keylabel = Label(fb,text="key: ")
        keylabel.pack(side='right',pady=10,padx=1)
        # status area at bottom
        fs = Frame(master)
        Label(fs,text='status:').pack(side='left', padx=4, pady=6)
        self.statusbar = Label(fs,text="", anchor='w', relief='groove')
        self.statusbar.pack(side='right', fill='x', expand=1, pady=6, padx=4)
        
        fs.pack(side='bottom', fill='x') #, expand=1)
        fb.pack(side='bottom', fill='x') #,expand=1)
        
    def statusMessage(self,msg,timeout=5):
        " Cleared after t=timeout seconds : 0 means never clear"
        self.statusbar.configure(text=msg)
        if timeout > 0:
            t = timeout * 1000 # ms
            self.statusbar.after(t, self.statusClear)
            
    def statusClear(self):
        self.statusbar.configure(text="")

    def help(self):
        webbrowser.open(webpage,1)

    def about(self):
        s = "Scatter.py 1.0\n\n" + \
            "A scatterplot tool for reading " +\
            "Spider doc files and executing " +\
            "commands.\n" + \
            "Written in Python, Pmw.Blt\n\n" +\
            "July 28, 2004    Bill Baxter"
        tkMessageBox.showinfo("About scatter.py 1.0", s)

    # ------ plot the data -------
    # NB in Blt lingo, data points marked by circles, triangles, etc
    # are 'elements', while data points represented by numbers (text)
    # are 'markers'. Hence, routines have to check for both.

    def plotTheData(self):
        self.putsymbols()
            
    def putsymbols(self):
        if self.symbolshape == 'numbers':
            self.puttextsymbols()
            return
        keys = self.D.keys()
        for k in keys:
            mk = str(k)
            if self.g.element_exists(mk):
                self.g.element_delete(mk)
            if self.g.marker_exists(mk):
                self.g.marker_delete(mk)
            x = (self.D[k][0])
            y = (self.D[k][1])
                
            self.g.line_create(mk, xdata=x, ydata=y,
                               symbol=self.symbolshape,
                               pixels = self.symbolsize,
                               fill=self.symbolcolor,
                               outline="black")
            self.g.element_bind(mk,"<Enter>",
                               func=lambda event, m=mk: self.entermarker(m))
            self.g.element_bind(mk,"<Leave>",
                               func=lambda event, m=mk: self.leavemarker(m))
            self.g.element_bind(mk,"<Double-Button-1>",
                               func=lambda event, m=mk: self.selectmarker(m))

    def puttextsymbols(self):
        keys = self.D.keys()
        for k in keys:
            mk = str(k)
            if self.g.marker_exists(mk):
                self.g.marker_delete(mk)
            if self.g.element_exists(mk):
                self.g.element_delete(mk)
            x = (self.D[k][0])
            y = (self.D[k][1])

            self.g.marker_create("text", name = mk)
            self.g.marker_configure(mk, 
                                    coords = (x, y),
                                    text = mk,
                                    fill="", # make bgd transparent
                                    foreground=self.symbolcolor)

            self.g.marker_bind(mk,"<Enter>",
                               func=lambda event, m=mk: self.entermarker(m))
            self.g.marker_bind(mk,"<Leave>",
                               func=lambda event, m=mk: self.leavemarker(m))
            self.g.marker_bind(mk,"<Double-Button-1>",
                               func=lambda event, m=mk: self.selectmarker(m))

            
    def selectmarker(self,mk):
        " calls a program when a data point is selected "
        self.keyvalue.set(mk)  # entry box at the bottom
        self.substitute()
        self.callprogram()

    def keyleave(self, event=None):
        self.entermarker(self.keyvalue.get())
            
    def entermarker(self, mk):
        " Numbers are markers. all other symbols are elements"
        if self.g.element_exists(mk):
            self.g.element_configure(mk, fill=self.activecolor)   
        elif self.g.marker_exists(mk):
            self.g.marker_configure(mk, outline=self.activecolor) #foreground=self.activecolor)
        else:
            self.statusMessage("element %s not found" % mk, timeout=0)
            return
        self.leavemarker(self.oldvalue.get())
        self.keyvalue.set(mk)
        self.oldvalue.set(mk)
        self.substitute()     # puts marker number in filename

    def leavemarker(self, mk):
        if self.g.element_exists(mk):
            self.g.element_configure(mk, fill=self.symbolcolor)
        elif self.g.marker_exists(mk):
            self.g.marker_configure(mk, outline=self.symbolcolor)

    def putpoints(self, graph, pixsize=1, symbol='circle', color='white'):
        " puts single pixels in data point graph."
        keys = self.D.keys()
        for k in keys:
            mk = str(k)
            x = (self.D[k][0])
            y = (self.D[k][1])
            graph.line_create(mk, xdata=x, ydata=y,
                              symbol = symbol,
                              pixels = pixsize,
                              outlinewidth = 0,
                              fill = color)

    def mkPixelist(self, graph=None):
        " Creates x sorted list of pixels (x,y,key) in self.pixlist "
        fcoords = []
        keys = self.D.keys()
        for k in keys:
            mk = str(k)
            x = (self.D[k][0])
            y = (self.D[k][1])
            fcoords.append(x)
            fcoords.append(y)

        pcoords = self.getPixelCoords(fcoords, graph=graph)
        pixlist = []
        for k in keys:      # keys list still corresponds to pcoords
            p = pcoords[0]
            x = p[0]
            y = p[1]
            pix = (x,y,k)
            pixlist.append(pix)
            pcoords = pcoords[1:]
        pixlist.sort()
        self.pixlist = pixlist

    # ---------------------------------------------------
    def selectedPolyMarkers(self):
        " colors the elements in the selected polygon "
        polygon = self.Polygons[self.selectedPolygon]
        if hasattr(polygon,'elements') and len(polygon.elements) > 0:
            pass
        else:
            return
        pf = polygon.elements[0]
        keys = self.D.keys()
        
        for key in keys:         
            if self.symbolshape != 'numbers':
                if self.g.element_exists(key):
                    mk = int(key)   # polygon.elements are numbers, keys are strings
                    if mk in polygon.elements:  # why are str() needed for next call?
                        self.g.element_configure(str(key), fill=str(self.polyElementsColor))
                    else:
                        self.g.element_configure(str(key), fill=str(self.symbolcolor))
                else:
                    print "can't find element %s" % key
            else:
                if self.g.marker_exists(key):
                    mk = int(key)   # polygon.elements are numbers, keys are strings
                    if mk in polygon.elements:  
                        self.g.marker_configure(str(key), outline=str(self.polyElementsColor))
                    else:
                        self.g.marker_configure(str(key), outline=str(self.symbolcolor))
                else:
                    print "can't find marker %s" % key         

        
    def symbolsMenu(self):
        w = Toplevel(self.top)
        w.title("Symbols")
        fs = Frame(w)
        
        Label(fs,text='symbols:').grid(row=0,column=0)
        sbox = Pmw.ComboBox(fs, scrolledlist_items=shapes, dropdown=1)
        sbox.selectitem(self.symbolshape)
        sbox.grid(row=0,column=1)
        sbox.bind('<Return>', lambda event, w=w: self.setsymbols(win=w))

        Label(fs,text='size:').grid(row=1,column=0)
        self.eVar1.set(str(self.symbolsize))
        esize = Entry(fs, textvariable=self.eVar1)
        esize.grid(row=1,column=1)
        esize.bind('<Return>', lambda event, w=w: self.setsymbols(win=w))

        Label(fs,text='color:').grid(row=2,column=0)
        self.eVar2.set(self.symbolcolor)
        ecol = Entry(fs, textvariable=self.eVar2)
        ecol.grid(row=2,column=1)
        ecol.bind('<Return>', lambda event, w=w: self.setsymbols(win=w))

        Label(fs,text='active color:').grid(row=3,column=0)
        self.eVar3.set(self.activecolor)
        acol = Entry(fs, textvariable=self.eVar3)
        acol.grid(row=3,column=1)
        acol.bind('<Return>', lambda event, w=w: self.setsymbols(win=w))

        self.newsymbols = (sbox, esize, ecol, acol)
        b = Button(fs,text='ok',command = lambda w=w: self.setsymbols(win=w))
        b.grid(row=4,column=1,sticky='e', padx=4, pady=4)
        fs.pack()       
        
    def setsymbols(self, win=None):
        sbox = self.newsymbols[0]
        self.symbolshape = sbox.get()

        try:
            self.symbolsize = int(self.eVar1.get())
        except:
            pass

        self.symbolcolor = self.eVar2.get()
        self.activecolor = self.eVar3.get()
        win.destroy()
        self.putsymbols()

    def getPixelCoords(self, coords, graph=None):
        " given floating point data coords, returns pairs of coords, "
        " where pixel origin (0,0) is upper left corner of graph area" 
        c = coords[1]
        if type(c) == type( (1,2) ):
            print "coords must be list of numbers, not a list of pairs"
            return
        n = len(coords)
        if n%2 ==1:
            print "error: odd number of coordinates"
            return

        if graph == None:
            graph = self.g
        xedge = graph.extents("leftmargin") + self.magicnumber
        yedge = graph.extents("topmargin") + self.magicnumber

        pixels = []

        for p in range(n/2):
            x = coords[0]
            y = coords[1]
            i,j = graph.transform(x, y)
            pixels.append( (i-xedge, j-yedge) )
            if len(coords) > 1:
                coords = coords[2:]
        return pixels
        

    ###########################################################
    # Magnification vs Polygon selection functions and bindings

    def magnifyOn(self):
        if self.state != 'Magnify':  # turn on magnifier
            self.state = 'Magnify'
            self.magnifyVar.set(1)
            self.polygonVar.set(0)
            self.unbind_poly()
            self.polygonVar.set(0)
            self.bind_mag()
            self.setMagPolyButtons()
            self.statusMessage("Left mouse button selects a region for zooming in", 0)
        else:
            self.state = 'None'
            self.unbind_mag()
            self.setMagPolyButtons()
            self.magnifyVar.set(0)
            self.statusClear()

    def polygonOn(self):
        if self.state != 'Polygon':  # turn on polygon selection
            self.state = 'Polygon'
            self.polygonVar.set(1)
            self.magnifyVar.set(0)
            self.unbind_mag()
            self.bind_poly()
            self.setMagPolyButtons()
            self.polyUnfinish()
            self.statusMessage("Outline a polygon with Left mouse button, finish with Right mouse button",0) 
        else:
            self.state = 'None'
            self.polygonVar.set(0)
            self.unbind_poly()
            self.setMagPolyButtons()
            self.statusClear()
   
    def setMagPolyButtons(self):
        if self.state == 'Magnify':
            self.magButton.configure(background=self.activestatecolor)
            self.magButton.configure(activebackground=self.activestatecolor2)
            self.polyButton.configure(background=self.backgroundcolor)
            self.polyButton.configure(activebackground=self.activebackgroundcolor)
        elif self.state == 'Polygon':
            self.polyButton.configure(background=self.activestatecolor)
            self.polyButton.configure(activebackground=self.activestatecolor2)
            self.magButton.configure(background=self.backgroundcolor)
            self.magButton.configure(activebackground=self.activebackgroundcolor)
        else:  # state == 'None'
            self.polyButton.configure(background=self.backgroundcolor)
            self.polyButton.configure(activebackground=self.activebackgroundcolor)
            self.magButton.configure(background=self.backgroundcolor)
            self.magButton.configure(activebackground=self.activebackgroundcolor)
            
    def bind_mag(self):
        self.g.bind(sequence="<ButtonPress>", func=self.magMouseDown)
        self.g.bind(sequence="<ButtonRelease>", func=self.magMouseUp)
    def unbind_mag(self):
        self.g.unbind(sequence="<ButtonPress>")
        self.g.unbind(sequence="<ButtonRelease>")
            
    def bind_poly(self):
        self.g.bind(sequence="<Button-1>",   func=self.polyMouseDown)
        self.g.bind(sequence="<Button-3>",   func=self.polyFinish)
    def unbind_poly(self):
        self.g.unbind(sequence="<Button-1>")
        self.g.unbind(sequence="<Button-3>")

    
    # Polygon Select Functions ========================================
    
    def polyMouseDown(self, event):
        x, y = self.g.invtransform(event.x, event.y)
        if len(self.poly.coords) > 0:
            xprev, yprev = self.poly.coords[-1]
        else:
            xprev, yprev = x,y
        self.draw_line(xprev,yprev,x,y)
        self.poly.coords.append( (x, y) )
    
    def polyMouseUp(self, event):
        xnext, ynext = self.g.invtransform(event.x,event.y)
        xprev, yprev = self.poly.coords[-1]
        self.poly.coords.append( (xnext, ynext) )
        self.rubberbandLine = self.draw_line(xprev,yprev,xnext,ynext)
        self.poly.names.append(self.rubberbandLine)

    def polyMouseMotion(self, event):
        x, y = self.g.invtransform(event.x, event.y)
	lastx, lasty = self.poly.coords[-1]

	if (lastx != event.x)  and (lasty != event.y) :
            if hasattr(self, 'rubberbandLine'):
                self.g.marker_delete(self.rubberbandLine)
	    self.rubberbandLine = self.draw_line(lastx, lasty, x, y)
	    self.g.update_idletasks()

    def polyUnfinish(self):
        " if user selects another function while partial polygon remains "
        # erase the single lines
        for n in self.poly.names:
            if self.g.marker_exists(n):
                self.g.marker_delete(n)
        self.poly.clear()
        
    def polyFinish(self, event):
        if len(self.poly.coords) == 0:
            return
        # erase the single lines
        for n in self.poly.names:
            if self.g.marker_exists(n):
                self.g.marker_delete(n)
                
        coords = self.poly.pairs2oneList()
                
        if len (self.poly.coords) < 3:
            self.poly.clear()
            return
                
        polyname = self.mknewname(prefix='poly')
        self.g.marker_create("polygon", name=polyname, coords=coords,
                             outline=self.selectedPolygonColor,
                             linewidth=1, fill="")
        
        p = Polygon()  # class instance
        p.coords = coords
        self.Polygons[polyname] = p
        self.selectedPolygon = polyname
        self.polygonSelect(polyname)  # to deselect the others
        self.poly.clear()
        
        self.g.marker_bind(polyname, "<Button-1>",
                           lambda event, n=polyname:self.polygonSelect(n))
        self.g.marker_bind(polyname, "<Double-Button-1>",
                           lambda event, n=polyname:self.polygonSelect(n))

    def draw_line(self,x1,y1,x2,y2):
        name = self.mknewname(prefix="line")
        self.g.marker_create("line", name=name)
        coords = (x1,y1,x2,y2)
        self.g.marker_configure(name, coords=coords, outline='red')
        self.poly.names.append(name)
        return name

    def mknewname(self, prefix=None):
        name = self.uniqname
        p = int(self.uniqname)
        p = p + 1
        self.uniqname = str(p)
        if prefix != None:
            name = prefix + name
        return name

    def polygonSelect(self, name):
        polys = self.Polygons.keys()
        for p in polys:
            if p != name:
                self.g.marker_configure(p,outline=self.PolygonColor)
            else:
                self.g.marker_configure(p,outline=self.selectedPolygonColor)
        self.selectedPolygon = name

    def polybounds(self, coords):
        " Returns 4 tuple of polygon's bounding box"
        c = coords[0]
        xmax = xmin = c[0]
        ymin = ymax = c[1]
        for c in coords:
            x = c[0]
            y = c[1]
            if x > xmax: xmax = x
            if x < xmin: xmin = x
            if y > ymax: ymax = y
            if y < ymin: ymin = y
        return (int(xmin), int(ymin), int(xmax), int(ymax))

    def savePolygon(self):
        #print "saving polygon " + self.selectedPolygon
        if self.selectedPolygon == "":
            print "no polygon selected"
            return
        pobj = self.Polygons[self.selectedPolygon]
        screencoords = self.getPixelCoords(pobj.coords)
        polyname = self.selectedPolygon
       
        # get PIL image of datapoints, polygon
        Pmw.showbusycursor()
        im, imsk = self.mkDataMasks(polyname, pobj.coords)
        
        im = im.convert("1")       # binarize the images
        imsk = imsk.convert("1")

        res = ImageChops.multiply(im,imsk)  

        rp = self.getpixels(res)   # rp = list of pixels found in image
        fk = self.findPixelsInList(rp)  # fk = list of matched elements
        #remove duplicates
        dup = []
        for f in fk:
            if f not in dup:
                dup.append(f)
        fk = dup

        # add the found elements to the Polygon instance
        pobj.elements = fk
        self.Polygons[self.selectedPolygon] = pobj
        self.selectedPolyMarkers()  # colors the selected markers

        self.saveSelectedMarkers()

        del(im, imsk, res)  # does this speed up or slow down next call?
        Pmw.hidebusycursor()

    def saveSelectedMarkers(self, filename=None):
        if self.selectedPolygon == "":
            self.statusMessage("no polygon selected")
            return
            
        if filename == None:
            filename = asksaveasfilename(title="Save points into a doc file")
            if len(filename) == 0:
                return
        print filename
        if type(filename) == type( (1,2) ):
            print 'filename is a tuple!'
            
        poly = self.Polygons[self.selectedPolygon]
        e = poly.elements
        e.sort()

        E = []
        # Make header for doc file and column headings
        E.append(Spiderutils.makeDocfileHeader(filename))
        xaxis = self.xaxisVar.get()
        if xaxis == "": xaxis = 'x axis'
        if len(xaxis) > 11:
            xaxis = xaxis[:11] # just use first 11 chars
        else:
            xaxis = string.center(xaxis,11)
            
        yaxis = self.yaxisVar.get()
        if yaxis == "": yaxis = 'y axis'
        if len(yaxis) > 11:
            yaxis = yaxis[:11] # just use first 11 chars
        else:
            yaxis = string.center(yaxis,11)
        keyhdr = string.center('keys',11)
        E.append(" ; /        %s %s %s\n" %(keyhdr, xaxis, yaxis))
            
        n = len(e)
        for i in range(n): 
            key = i+1
            val = e[i]
            x,y = self.D[val][0], self.D[val][1]
            line =  "%5d %2d %11.1f %11f %11f\n" % (key, 3, float(e[i]), x, y)
            E.append(line)

        try:
            fp = open(filename,'w')
            fp.writelines(E)
            fp.close()
            self.statusMessage('Data written to %s' % filename)
        except:
            self.statusMessage('ERROR: unable to write to %s' % filename, timeout=0)

    def getpixels(self, image):
        " Returns list of pixels above threshold "
        I = list(image.getdata())
        wd,ht = image.size
        if wd == 0 or ht == 0:
            return []
        threshold = 0
        coords = []
        k = 0
        for j in range(ht):
            for i in range(wd):
                p = I[k]
                k = k + 1
                if p > threshold:
                    coords.append( (i,j,p) )
        return coords

    def findPixelsInList(self, plist):
        " plist is list of coords (x,y,intensity) from image.getdata() "
        " compare to sorted pixels in self.pixlist (x,y,key)    "
        if len(plist) == 0:
            return
        fkeys = []
        limit = 2  # i.e. within 1 pixel
        for pix in plist:
            xp, yp = pix[0],pix[1]
            for sp in self.pixlist:
                x, y, key = sp[0], sp[1], sp[2]
                if abs(x-xp) < limit and abs(y-yp) < limit:
                    fkeys.append(key)
        return fkeys

    def Bltgraph2Image(self, graph):
        " Returns a PIL image of graph region only"
        wd = int(graph.cget('width'))
        ht = int(graph.cget('height'))
        
        newname = self.mknewname(prefix='photo')
        photo = PhotoImage(name=newname, master=self.top, width=wd, height=ht)
        #photo = ImageTk.PhotoImage(name=newname, master=self.top, width=wd, height=ht)
        a,b = int(photo.cget('width')), int(photo.cget('height'))
        #print " photoimage size %d %d" % (a,b)
        
        graph.snap(photo)
        c,d = int(photo.cget('width')), int(photo.cget('height'))
        #print " size after snap %d %d" % (c,d)
        #if c != a or d != b:
            #print "ERROR: size mismatch inBltgraph2Image"
            #print " photoimage size %d %d" % (a,b)
            #print " size after graph.snap %d %d" % (c,d)

        x1 = graph.extents("leftmargin") + self.magicnumber
        y1 = graph.extents("topmargin") + self.magicnumber
        x2 = x1 + graph.extents("plotwidth") - 1
        y2 = y1 + graph.extents("plotheight") - 1
        coords = (x1,y1,x2,y2)

        # Image can't read a TkPhotoImage (strange, huh?). So we write it
        # out as a tmp file and read it back with a format Image can grok.       
        tmpfile = 'tmp67egdf.ppm'
        photo.write(tmpfile, from_coords=coords)
        im = Image.open(tmpfile)
        os.remove(tmpfile)

        return im
            
    def mkDataMasks(self, polyname, coords):
        " returns Image.image of datapoints "
        self.statusMessage("Finding points inside selected region...",timeout=0)
        geo = self.top.geometry()
        tmpwin = Toplevel(self.top)
        tmpwin.title("scamdoolie")
        tmpwin.lower()
        tmpwin.geometry(geo)  # try to put it behind the main window
        tmpwin.update_idletasks()
        
        # -- make the graph with data points ---
        g1 = Pmw.Blt.Graph(tmpwin,   #width=g_width, height=g_height,
                           #plotborderwidth = 0,
                           plotbackground="black")
        
        g1.line_create("scatter1", xdata=self.vx, ydata=self.vy,
                       linewidth=0, symbol='')
        g1.legend_configure(hide=1)
        g1.axis_configure("x", title=self.xaxisVar.get())
        g1.axis_configure("y", title=self.yaxisVar.get())
        self.putpoints(g1)
        
        # has to be packed, ow snap returns a 400x 400 image
        g1.pack(expand=1, fill='both')
        g1.update_idletasks()

        # do stuff while the graph is still around
        #print "plotborderwidth: %s" % (str(g1.cget('plotborderwidth')))
        im = self.Bltgraph2Image(g1)
        self.graph_wd = int(g1.cget('width'))
        self.graph_ht = int(g1.cget('height'))
        #print "mkDataMask: graph size = %d %d" % (self.graph_wd, self.graph_ht)

        # ------ create the list of pixel coords for the data ---
        self.mkPixelist(graph=g1)

        g1.pack_forget()
        #########################################################
        # -- make the graph with the filled polygon ---
        g2 = Pmw.Blt.Graph(tmpwin, width=self.graph_wd, height=self.graph_ht,
                           #plotborderwidth = 0,
                           plotbackground="black")
        g2.line_create("scatter2", xdata=self.vx, ydata=self.vy,
                       linewidth=0, symbol='')
        g2.legend_configure(hide=1)
        g2.axis_configure("x", title=self.xaxisVar.get())
        g2.axis_configure("y", title=self.yaxisVar.get())

        g2.marker_create("polygon", name='polyname')
        g2.marker_configure('polyname',
                           coords= coords,
                           fill="white")
        g2.pack(expand=1, fill='both')
        g2.update_idletasks()

        graph_wd = int(g2.cget('width'))
        graph_ht = int(g2.cget('height'))
        #print "mkPolyMask: graph size = %d %d" % (graph_wd, graph_ht)
        
        imsk = self.Bltgraph2Image(g2)
        graph_wd = int(g2.cget('width'))
        graph_ht = int(g2.cget('height'))
        #print "mkPolyMaskafter BLT: graph size = %d %d" % (graph_wd, graph_ht)
        
        tmpwin.destroy() # bye
        self.statusClear()
        return im, imsk    

    def deletePolygon(self):
        pname = self.selectedPolygon
        if self.g.marker_exists(pname):
            self.g.marker_delete(pname)
        del(self.Polygons[pname])
        self.selectedPolygon = ""

    # Magnifier Functions =============================================

    def magMouseDown(self, event):
        self.dragging = 0
        
        if self.g.inside(event.x, event.y):
            self.dragging = 1
            (self.x0, self.y0) = self.g.invtransform(event.x, event.y)
            
            self.g.marker_create("line", name="marking rectangle",dashes=(2, 2),
                                 outline='green')
            
            self.g.bind(sequence="<Motion>",  func=self.magMouseDrag)
    
    def magMouseUp(self, event):
        #global x0, y0, x1, y1
        x0 = self.x0; x1 = self.x1
        y0 = self.y0; y1 = self.y1
        
        if self.dragging:
            self.g.unbind(sequence="<Motion>")
            self.g.marker_delete("marking rectangle")
            
            if x0 <> x1 and y0 <> y1:

                # make sure the coordinates are sorted
                if x0 > x1: x0, x1 = x1, x0
                if y0 > y1: y0, y1 = y1, y0
         
                if event.num == 1:
                   self.zoom(x0, y0, x1, y1) # zoom in
                else:
                   (X0, X1) = self.g.xaxis_limits()
                   k  = (X1-X0)/(x1-x0)
                   x0 = X0 -(x0-X0)*k
                   x1 = X1 +(X1-x1)*k
                   
                   (Y0, Y1) = self.g.yaxis_limits()
                   k  = (Y1-Y0)/(y1-y0)
                   y0 = Y0 -(y0-Y0)*k
                   y1 = Y1 +(Y1-y1)*k
                   
                   self.zoom(x0, y0, x1, y1) # zoom out
        
    def zoom(self, x0, y0, x1, y1):
        self.g.xaxis_configure(min=x0, max=x1)
        self.g.yaxis_configure(min=y0, max=y1)
        
    def magReset(self, event=None):
        xmin = self.vx[0]
        xmax = self.vx[1]
        ymin = self.vy[0]
        ymax = self.vy[1]
        self.g.xaxis_configure(min=xmin, max=xmax)
        self.g.yaxis_configure(min=ymin, max=ymax)
        
    def magMouseDrag(self, event):
        #global x0, y0, x1, y1
        (self.x1, self.y1) = self.g.invtransform(event.x, event.y)
        x0 = self.x0; x1 = self.x1
        y0 = self.y0; y1 = self.y1
             
        self.g.marker_configure("marking rectangle", 
            coords = (x0, y0, x1, y0, x1, y1, x0, y1, x0, y0))

    def transformer(self, event):
        " debugging tool "
        (x,y) = (event.x, event.y)
        print "(%d,%d)" %(x,y),

        # inversetransform clicked pos to get axis-coordinates.
        (x,y) = (self.g.xaxis_invtransform(x), self.g.yaxis_invtransform(y))
        print "--> (%f,%f)" %(x,y),

        # transform back to get original window-coordinates.
        # (not always exact because of round-off errors.)
        (x,y) = (self.g.xaxis_transform(x), self.g.yaxis_transform(y))
        print "--> (%d,%d)" %(x,y)

    ###########################################################
    # Functions for executing commands

    def callprog_check(self, event=None):
        " doesn't substitute if no ***s , bound to cmd entry "
        cmd = self.command.get()
        if string.find(cmd, '*') < 0 and string.find(cmd, '?') < 0:
            cmd = self.substitute()
        self.execprogram(cmd)
    
    def callprogram(self, event=None):
        "executes the text in the entry box self.command"
        cmd = self.substitute()
        #print "executing %s" % (cmd)
        self.execprogram(cmd)

    def execprogram(self, cmd):
        self.statusMessage("starting: %s" % (cmd), 4)
        newpid = os.fork()
        if newpid == 0: # child
            os.system(cmd + " &")    # fork a new process and then go away
            self.top.quit()


    # substitute changes the command in self.command
    def substitute(self, event=None):
        " ( 'out***.dat', '2' ) --> 'out002.dat' writes to self.command"
        " if no stars are found, tries to match numbers before '.'     "
        cmd = self.command.get()
        key = self.keyvalue.get()
        if cmd == "" or key == "":
            return ""

        restq = re.compile('[*]+|[?]+')
        renums = re.compile('[0-9]+\.')  # get digits before dot
        
        a = restq.search(cmd)  # look for stars (or ???'s)
        if not a:
            a = renums.search(cmd)  # no stars, try for digits
            
        if a:
            b = a.span()
            start = b[0]
            end = b[1]
            if cmd[end-1] == '.':
                end = end-1
            length = end - start
            num = string.zfill(key,length)
            cmd = cmd[:start] + num + cmd[end:]

        self.command.set(cmd)
        return cmd
    
       
    ###########################################################
    # Functions for File I/O
    
    def getdata(self, filename=None, xcol=None, ycol=None):
        if filename == None or filename == "":
            filename = self.filename
            if self.filename == None or self.filename == "":
                return

        try:
            #D = spiderutils.readSpiderDocFile(filename, [xcol,ycol])
            D = Spiderutils.readdoc(filename, keys='all')
        except:
            print "Unable to get data from %s" % filename
            return

        keys = D.keys()
        if len(keys) < 1:
            print "%s contains no data" % filename
            return

        f = D[keys[0]]
        xmin = f[0]
        xmax = f[0]
        ymin = f[1]
        ymax = f[1]
        
        for k in keys:
            x = D[k][0]
            y = D[k][1]
            if x < xmin: xmin = x
            elif x > xmax: xmax = x
            if y < ymin: ymin = y
            elif y > ymax: ymax = y

        vector_x = (xmin,xmax)
        vector_y = (ymin,ymax)
        return (D, vector_x, vector_y)

    # ------------------------------
    # display head and tail of file to user
    # user selects columns to be plotted

    def openfile(self, filename=None):
        if filename == None or len(filename) == 0:
            filename = askfilename()
            if len(filename) == 0:
                return
        if not os.path.exists(filename):
            print "cannot find %s" % filename
            return

        self.filename = filename
        w = Toplevel(self.top)
        w.title(os.path.basename(filename))
        
        ft = Frame(w)
        ft.pack(side='top', fill='x', expand=1)
        dxlabelstr = "%s : select columns to plot" % os.path.basename(filename)
        dx = Label(ft, text=dxlabelstr, font="Helvetica 12 bold")
        dx.pack(side='top', padx=10, pady=10)
        dx2 = Label(ft, text="Use '1' for first data column, '0' for key column",
                    font="Helvetica 10")
        dx2.pack(side='top', padx=10, pady=10)

        textbox = Pmw.ScrolledText(w)
        textbox.pack(side='top', fill='both', expand=1)
        width = int(textbox.component('text').cget('width'))
        height = int(textbox.component('text').cget('height'))

        # see if file has more lines than the window
        wccmd = "wc -l %s" % (filename)
        wc = string.strip(getoutput(wccmd))
        if string.find(wc," ") > -1:
            wc = string.split(wc)[0]
        wc = int(wc)

        if wc < 2 * height:  # then use entire file
            lines = getoutput("cat %s" % (filename))
            maxwd = 0
            h = string.split(lines,"\n")
            for line in h:
                if len(line) > maxwd:
                    maxwd = len(line)
            if maxwd > width:
                textbox.component('text').configure(width = maxwd)
            textbox.component('text').insert(END, lines)
            
        else:    # just display the head and tail of the file
            nlines = int(height/3)
            headcmd = "head -%d %s" % (nlines, filename)
            tailcmd = "tail -%d %s" % (nlines, filename)
            hd = getoutput(headcmd)
            tl = getoutput(tailcmd)
            maxwd = 0
            h = string.split(hd,"\n")
            for line in h:
                if len(line) > maxwd:
                    maxwd = len(line)
            if maxwd > width:
                textbox.component('text').configure(width = maxwd)
            
            textbox.component('text').insert(END, hd)
            textbox.component('text').insert(END, '\n    :\n')
            for i in range(nlines-1):
                textbox.component('text').insert(END, '                 :\n')
            textbox.component('text').insert(END, tl)
            
        fr = Frame(w)
        fr.pack(side='bottom', fill='x', expand=1, padx=5, pady=5)
        xce = Pmw.EntryField(fr, labelpos='w',
                             label_text='column to use for X axis:',
                             value=str(self.xcolumn.get()),
                             validate= {'validator':'numeric'})
        xce.component('entry').configure(textvariable=self.xcolumn)
        yce = Pmw.EntryField(fr, labelpos='w',
                             label_text='column to use for Y axis:',
                             value=str(self.ycolumn.get()),
                             validate= {'validator':'numeric'})
        yce.component('entry').configure(textvariable=self.ycolumn)
        xax = Pmw.EntryField(fr, labelpos='w',
                             label_text='Label for X axis:',
                             value=self.xaxisVar.get())
        xax.component('entry').configure(textvariable=self.xaxisVar)
        yax = Pmw.EntryField(fr, labelpos='w',
                             label_text='Label for Y axis:',
                             value=self.yaxisVar.get())
        yax.component('entry').configure(textvariable=self.yaxisVar)
        xce.grid(row=0,column=0, sticky='w', padx=10, pady=4)
        yce.grid(row=1,column=0, sticky='w', padx=10, pady=4)
        xax.grid(row=0,column=1, sticky='e',padx=10, pady=4)
        yax.grid(row=1,column=1, sticky='e',padx=10, pady=4)

        cancel = IntVar()
        cancel.set(0)

        ok = Button(fr, text='Ok',
                    command=lambda x=w,t='ok': self.okdone(x,t))
        ok.grid(row=2, column=0, sticky='sw', padx=5, pady=5)
        quit = Button(fr, text='Cancel',
                      command=lambda x=w, t='cancel': self.okdone(x,t))
        quit.grid(row=2, column=1, sticky='se', padx=5, pady=5)

        w.lift(self.top)
        self.top.lower(w)
        w.update_idletasks()
        self.top.wait_window(w)

    def okdone(self, win, button):
        win.destroy()
        if button == 'ok':
            self.callgetdata()
        else:
            self.top.destroy()
 
    def callgetdata(self):
        xcol = self.xcolumn.get()
        ycol = self.ycolumn.get()
        self.D, self.vx, self.vy = self.getdata(self.filename, xcol=xcol, ycol=ycol)

def printHelp():
    print "Usage: scatter.py docfile [x y]"
    print "where x,y specify which Spider docfile columns to plot."
    print "(1 = first data column, 0 = keys)"

if __name__ == '__main__':
    length = len(sys.argv[1:])

    filename = ""
    
    x = -1
    y = -1

    if length > 0:
        filename = sys.argv[1]
        if filename == '-h' or filename == '-help':
            printHelp()
            sys.exit()
        if not os.path.exists(filename):
            print "cannot find file %s" % filename
            sys.exit()
    if length > 2:
        try:
            x = int(sys.argv[2])
            y = int(sys.argv[3])
        except:
            printHelp()
            sys.exit()

    master = Tk()
    Pmw.initialise(master, fontScheme = 'pmw1')
    master.title("scatter plot")

    if x > -1 and y > -1:
        c = scatter(master, filename, xcol=x, ycol=y)
    else:
        c = scatter(master, filename)
    
    master.mainloop() 
