#!/usr/bin/env python

# Spider Python Library
# Copyright (C) 2006  Health Research Inc.
#
# HEALTH RESEARCH INCORPORATED (HRI),
# ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455
#
# Email:  spider@wadsworth.org
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

import sys
from Spider import Spiderutils

if sys.argv[1:] :
    filelist = sys.argv[1:]
    NodeTemplate = Spiderutils.name2template(filelist[0])
#    print 'file template:', NodeTemplate
else : 
    print
    print "syntax: tree.py node_img001.ext {max_depth margin_width canvas_width}"
    print
    sys.exit()

if sys.argv[2:] :
    maxDepth = int(sys.argv[2])
else : 
    maxDepth = 6
print 'max_depth:', maxDepth

last_node = 2**(maxDepth) - 1
last_file = Spiderutils.template2filename(NodeTemplate,n=last_node)
print 'last_file:', last_file

if sys.argv[3:] :
    marginWidth = int(sys.argv[3])
else : 
    marginWidth = 2
print 'margin_width:', marginWidth

if sys.argv[4:] :
    canvasWidth = int(sys.argv[4])
else : 
    canvasWidth = 1600
print 'canvas_width:', canvasWidth

from Tkinter import *
import Image
import ImageTk

root = Tk()
root.title('tree.py')
Coords={}
Photolist=[]

# get image dimensions:  ix,iy = im.size  # dimensions
Nodename = Spiderutils.template2filename(NodeTemplate,n=1)
Node = Image.open(Nodename)
xdim,ydim = Node.size

# calculate canvas dimensions
canvasx = 2**(maxDepth-1)*(xdim+marginWidth)
canvasy = maxDepth*ydim
Tree = Canvas(root,height=canvasy)
Tree.configure(width=canvasx)
Tree.configure(scrollregion=(0,0,canvasx,canvasy))
print 'Canvas dimensions:',canvasx,canvasy

# loop through depths, from bottom
for depthCounter in range(maxDepth):  # range is 0..maxDepth-1
    currentDepth = maxDepth-depthCounter  # range is maxDepth..1 by -1
    firstNode = 2**(currentDepth-1)
    lastNode = 2**(currentDepth) - 1
    nodesRow = lastNode-firstNode+1
    bandWidth = canvasx/nodesRow
    nodeCoordy = currentDepth*ydim
    
    # loop through nodes
    for rowNode in range(nodesRow):  # range is 0..nodesRow-1
        currentNode = rowNode+firstNode  # range is firstNode..lastNode
        nodeCoordx = (rowNode+0.5)*bandWidth

        # write coordinates to dictionary
        Coords[str(currentNode)] = [nodeCoordx,nodeCoordy]
        
        # if not bottom row, then draw lines to daughter nodes
        if currentDepth != maxDepth:
            # calculate node# for daughter nodes
            daughterNode1 = 2*currentNode
            daughterNode2 = 2*currentNode + 1
            
            # get coordinates for daughter nodes
            daughter1Coordx = Coords[str(daughterNode1)][0]
            daughter2Coordx = Coords[str(daughterNode2)][0]
            daughterCoordy = Coords[str(daughterNode1)][1] - ydim
            
            # draw lines from parent to daughters
            Tree.create_line(nodeCoordx,nodeCoordy-ydim, daughter1Coordx,daughterCoordy)
            Tree.create_line(nodeCoordx,nodeCoordy-ydim, daughter2Coordx,daughterCoordy)
            
        # draw in node
        Nodename = Spiderutils.template2filename(NodeTemplate,n=currentNode)
        Node = Image.open(Nodename)
        Tkimage = ImageTk.PhotoImage(Node.convert2byte(), palette=256)  # , master=root)
        Photolist.append(Tkimage)
        Tree.create_image(nodeCoordx,nodeCoordy, image=Tkimage, anchor=S)  # CENTER)
        last_node = Nodename

# Finish drawing
scrollBar = Scrollbar(root, command=Tree.xview, orient=HORIZONTAL)
Tree.configure(xscrollcommand=scrollBar.set)
Tree.pack()
Tree.postscript(file="test")
scrollBar.pack(side=BOTTOM, fill=X)
root.mainloop()
