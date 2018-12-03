#!/usr/bin/env python
#
# SOURCE:  binarytree.py 
#
# MODIFICATIONS:
#    TO DO: save settings
#    TO DO: draw from bottom
#    2015-09-15 -- added menu
#    2015-09-15 -- changes depth on the fly
#    2009-06-04 -- reads existing selection file
#    2009-06-03 -- saves selection file
#    2009-06-02 -- allows skipped nodes
#    print "tree.py, Modified 2015 Sep 15"
#
# Spider Python Library
# Copyright (C) 2006-2018  Health Research Inc., Menands, NY
# Email:    spider@health.ny.gov

import sys
import os
import Pmw 

from   Spider  import Spiderutils
from   Tkinter import *
from   PIL     import Image
from   PIL     import ImageTk

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

def clickOK(parent,win):
    win.destroy()

class BinaryTreeCanvas:
    def __init__(self, master, classTemplate, max_depth=5, savefilename='outfile.dat', 
        margin_width=3, canvasWidth=1600, labelBorder=3) :

        # store passed values
        self.master        = master
        self.classTemplate = classTemplate
        self.max_depth     = max_depth
        self.savefilename  = savefilename
        self.margin_width  = margin_width
        self.canvasWidth   = canvasWidth
        self.labelBorder   = labelBorder

        # initialize
        self.coords_dictionary={}
        self.photo_list=[]
        self.select_dictionary = {}

        # get image dimensions:  ix,iy = im.size  # dimensions
        classname = Spiderutils.template2filename(self.classTemplate,n=1)
        classavg = Image.open(classname)
        self.xdim = classavg.size[0] + 2*self.labelBorder
        self.ydim = classavg.size[1] + 2*self.labelBorder + self.margin_width

        # colors
        self.good_color = 'green'
        self.select_flag = 1
        self.bad_color = 'red'

    def makeMenus(self):
        # ------- create the menu bar -------
        self.mBar = Frame(self.master, relief='raised', borderwidth=1)
        self.balloon = Pmw.Balloon(self.master)

        self.leftframe = Frame(self.mBar) # , relief='raised', borderwidth=1)

        # Make the File menu
        Filebtn = Menubutton(self.leftframe, text='File', relief='flat')
        Filebtn.pack(side=LEFT, padx=5, pady=5)
        Filebtn.menu = Menu(Filebtn, tearoff=0)

        Filebtn.menu.add_command(label='Save selection', underline=0,
                                command=self.saveSelections)
        Filebtn.menu.add_command(label='Read selection', underline=0,
                                command=self.readSelections)
        Filebtn.menu.add_command(label='Close', underline=0,
                                     command=self.master.destroy)
        Filebtn['menu'] = Filebtn.menu

        # Make help menu
        Helpbtn = Menubutton(self.mBar, text='Help', underline=0, relief='flat')
        Helpbtn.pack(side=RIGHT, padx=5, pady=5)
        Helpbtn.menu = Menu(Helpbtn, tearoff=0)
        Helpbtn.menu.add_command(label='Keyboard shortcuts', underline=0, 
                                 command=self.shortcuts)
        Helpbtn['menu'] = Helpbtn.menu

        # Pack menu bar
        self.leftframe.pack(side=LEFT) # , padx=5, pady=5)
        self.mBar.pack(side='top', fill = 'x')

    def drawTree(self) :
        # calculate canvas dimensions
        self.canvasx = 2**(self.max_depth-1)*(self.xdim + self.margin_width)
        self.canvasy = self.max_depth * self.ydim

        if hasattr(self,'tree_canvas'): self.tree_canvas.destroy()
        if hasattr(self,'HscrollBar'):  self.HscrollBar.destroy()
        if hasattr(self,'VscrollBar'):  self.VscrollBar.destroy()

        # define canvas
        self.tree_canvas = Canvas(self.master,height=self.canvasy)
        self.tree_canvas.configure(width=self.canvasx)
        self.tree_canvas.configure(scrollregion=(0,0,self.canvasx,self.canvasy))
        print 'Canvas dimensions:',self.canvasx,self.canvasy

        # scrollbars
        self.HscrollBar = Scrollbar(self.master, command=self.tree_canvas.xview, orient=HORIZONTAL)
        self.tree_canvas.configure(xscrollcommand=self.HscrollBar.set)
        self.HscrollBar.pack(side=BOTTOM, fill=X)

        self.VscrollBar = Scrollbar(self.master, command=self.tree_canvas.yview, orient=VERTICAL)
        self.tree_canvas.configure(yscrollcommand=self.VscrollBar.set)
        self.VscrollBar.pack(side=RIGHT, fill=Y)

        self.class2label = {}  # lookup table to find an image
        
        # loop through depths, from bottom
        for depthCounter in range(self.max_depth):  # range is 0..max_depth-1
            currentDepth = depthCounter + 1
            firstNode = 2**(currentDepth-1)
            lastNode = 2**(currentDepth) - 1
            nodesRow = lastNode-firstNode+1
            bandWidth = self.canvasx / nodesRow
            nodeCoordy = currentDepth * self.ydim

            # loop through nodes
            for rowNode in range(nodesRow):  # range is 0..nodesRow-1
                current_node = rowNode+firstNode  # range is firstNode..lastNode
                nodeCoordx = (rowNode+0.5)*bandWidth

                # write coordinates to dictionary
                self.coords_dictionary[str(current_node)] = [nodeCoordx,nodeCoordy]

                # if not top row, then draw line to parent node
                if currentDepth != 1 :
                    # calculate node# for daughter nodes
                    parentNode = current_node/2  # will convert to INT (rightly so) if odd

                    # get coordinates for daughter nodes
                    parentCoordx = self.coords_dictionary[str(parentNode)][0]
                    parentCoordy = self.coords_dictionary[str(parentNode)][1] - self.ydim + self.margin_width
            
                # check if image exists
                classname = Spiderutils.template2filename(self.classTemplate,n=current_node)
            
                # if file exists, then draw lines and paste image
                if os.path.exists(classname):
                    # draw line to parent
                    if currentDepth != 1 :
                        self.tree_canvas.create_line(nodeCoordx,nodeCoordy-self.ydim, parentCoordx,parentCoordy, tag='parent_line')
                        self.tree_canvas.lower('parent_line')  # otherwise lines are on top of images

                    # prepare image window and store info
                    classavg = Image.open(classname)
                    Tkimage = ImageTk.PhotoImage(classavg.convert2byte(), palette=256)
                    image_label = Label(self.tree_canvas, image=Tkimage, borderwidth=self.labelBorder)
                    image_label.photo = Tkimage
                    image_label.filenum = current_node
                    self.class2label[current_node] = image_label
                    image_label.bind('<Button-1>', lambda event, w=image_label : self.selectClass(w))
                    self.photo_list.append(image_label)

                    # draw image window
                    self.tree_canvas.create_window(nodeCoordx,nodeCoordy, window=image_label, anchor=S)
            
                    last_node = classname  # necessary?

        # Finish drawing
        self.tree_canvas.pack()
        self.master.bind('<Control-t>', self.test)
        self.master.bind('<Control-s>', self.saveSelections)
        self.master.bind('<Control-r>', self.readSelections)
        self.master.bind('<Control-w>', self.closeWindow)
        self.master.bind('<plus>', self.levelUp)
        self.master.bind('<minus>', self.levelDown)
 
    def shortcuts(self):
        sc = Toplevel(self.master)
        sc.title("Shortcuts")
        rownum = 0

        scmain = Frame(sc, borderwidth=2, relief=RIDGE)
        Label(scmain, text='Main window:').grid(row=rownum, sticky=W)

        rownum += 1
        Label(scmain, text='+'                 ).grid(row=rownum, column=0, sticky=W)
        Label(scmain, text='Add row'           ).grid(row=rownum, column=1, sticky=W)

        rownum += 1
        Label(scmain, text='-'                 ).grid(row=rownum, column=0, sticky=W)
        Label(scmain, text='Remove row'        ).grid(row=rownum, column=1, sticky=W)

        rownum += 1
        Label(scmain, text='Control-s'         ).grid(row=rownum, column=0, sticky=W)
        Label(scmain, text='Save selections'   ).grid(row=rownum, column=1, sticky=W)

        rownum += 1
        Label(scmain, text='Control-r'         ).grid(row=rownum, column=0, sticky=W)
        Label(scmain, text='Read selections'   ).grid(row=rownum, column=1, sticky=W)

        rownum += 1
        Label(scmain, text='Control-w'         ).grid(row=rownum, column=0, sticky=W)
        Label(scmain, text='Close window'      ).grid(row=rownum, column=1, sticky=W)

        scmain.pack(padx=5, pady=5, expand=1)

        okframe = Frame(sc)
        Button(okframe, text="OK", command=sc.destroy).pack()
        sc.bind('<Return>', lambda p=self.master, w=sc: clickOK(p,w)) 
        okframe.pack()

    def closeWindow(self, event=None):
        self.master.destroy()

    def levelDown(self, event=None) :
        self.max_depth -= 1
        print 'levelDown, max_depth:', self.max_depth
        self.drawTree()

    def levelUp(self, event=None) :
        self.max_depth += 1
        print 'levelup, max_depth:', self.max_depth
        self.drawTree()

    def test(self, event=None) :
        print "select_dictionary length:", len(self.select_dictionary)
        print "select_dictionary:", self.select_dictionary
        #print "master:", hasattr(self,'master')
        #print "tree_canvas:", hasattr(self,'tree_canvas')

    def selectClass(self, widget) :
        # get current color
        currentColor = widget.cget('background')
        
        # if current color isn't good color, then select
        if currentColor != self.good_color :
            widget.configure(background=self.good_color)
            self.select_dictionary[widget.filenum] = self.select_flag
        else : # deselect
            widget.configure(background=self.bad_color)
            del self.select_dictionary[widget.filenum]

    def saveSelections(self, event=None) :
        headers = ['class_num']
        spiderDictionary = {}
        key = 0
        selectKeys = self.select_dictionary.keys()
        selectKeys.sort()

        for counter in selectKeys:
            key = key + 1
            spiderDictionary[key] = [counter, self.select_dictionary[counter]]

        # if non-empty, backup prior versions
        if len(self.select_dictionary) > 0 : backup(self.savefilename)

        if Spiderutils.writeSpiderDocFile(self.savefilename, spiderDictionary, headers=headers, append=1):
            print 'Wrote', key, 'keys to %s' % os.path.basename(self.savefilename)
        else:
            print "Unable to write to", self.savefilename

    def readSelections(self, event=None) :
        if os.path.exists(self.savefilename):
            goodclassdoc = Spiderutils.readSpiderDocFile(self.savefilename)
            goodclasskeys = goodclassdoc.keys()
            found_counter = 0

            for key in goodclasskeys :
                filenumber = int(goodclassdoc[key][0])
                classnum   = int(goodclassdoc[key][1])

                # map particle to label and goodclass_label.configure
                if self.class2label.has_key(filenumber) :
                    goodclass_label = self.class2label[filenumber]
                    goodclass_label.configure(background=self.good_color)
                    self.select_dictionary[filenumber] = self.select_flag
                    found_counter = found_counter + 1
                else :
                    print "WARNING readSelections: image", key, filenumber, "not found"
            
            print found_counter, "keys found"

        else:
            print outfile, 'does not exist'

if __name__ == "__main__":

    argCounter = 0

    argCounter = argCounter + 1
    if sys.argv[argCounter:] :
        file_example = sys.argv[argCounter:]
        nodeTemplate = Spiderutils.name2template(file_example[0])
    #    print 'file template:', nodeTemplate
    else :
        print
        print "syntax: tree.py node_img001.ext {selectfile.ext max_depth margin_width canvas_width}"
        print
        sys.exit()

    argCounter = argCounter + 1
    if sys.argv[argCounter:] :
        max_depth = int(sys.argv[argCounter])
    else :
	max_depth = 5
    print 'max_depth:', max_depth

    argCounter = argCounter + 1
    if sys.argv[argCounter:] :
        savefilename = sys.argv[argCounter]
    else :
    #    print "template:", nodeTemplate
        extension = os.path.splitext(nodeTemplate)[1]
        savefilename = 'goodclasses' + extension
    print 'savefilename:', savefilename

    argCounter = argCounter + 1
    if sys.argv[argCounter:] :
        margin_width = int(sys.argv[argCounter])
	print 'margin_width:', margin_width

    argCounter = argCounter + 1
    if sys.argv[argCounter:] :
        canvasWidth = int(sys.argv[argCounter])
	print 'canvas_width:', canvasWidth

    argCounter = argCounter + 1
    if sys.argv[argCounter:] :
        labelBorder = int(sys.argv[argCounter])
	print 'image_border:', labelBorder

    root = Tk()
    root.title('tree.py')

    tree = BinaryTreeCanvas(root, nodeTemplate, max_depth=max_depth, savefilename=savefilename)
    tree.makeMenus()
    tree.drawTree()

    root.mainloop()
