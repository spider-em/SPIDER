#!/usr/bin/env python
#
# A simple SPIDER image viewer.
# Image (Python Imaging Library) used to convert SPIDER image,
# Tkinter used for display.

from Tkinter import *
import Image, ImageTk
import sys

if not sys.argv[1:]:
    print "Usage: python simpleview.py filename"
    sys.exit()
filename = sys.argv[1]

# open the file with the Python Imaging library (PIL)
im = Image.open(filename)

root = Tk()  # initialize Tkinter

# convert the PIL image into a Tk-compatible image
tkimage = ImageTk.PhotoImage(im.convert2byte(), palette=256)

# Paste the Tk image into a Tkinter Label
Label(root, image=tkimage).pack()

root.mainloop()  # start the GUI
