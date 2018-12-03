#!/usr/bin/env python

from   Tkinter import *
from   PIL     import Image
from   PIL     import ImageTk

import os, sys

class Viewer:
    def __init__(self, master, filename):
        self.top = master
        self.filename = filename
        self.index = 0
        #display first image
        if not os.path.exists(filename):
            print "Unable to find %s" % filename
            self.top.quit()

        labeltext = os.path.basename(filename)

        im = Image.open(filename)
        self.im = im
        if im.format == "SPIDER":
            if im.istack != 0:
                self.nimages = im.nimages
                labeltext += " : stack file with %d images" % im.nimages
            bim = im.convert2byte()
        self.size = im.size

        self.title = Label(text=labeltext)
        self.title.pack()
                
        self.tkimage = ImageTk.PhotoImage(bim, palette=256)
        self.tmp = ""
        
        self.lbl = Label(master, image=self.tkimage)
        self.lbl.pack(side='top')

        # the button frame
        fr = Frame(master)
        fr.pack(side='top', expand=1, fill='x')
        back = Button(fr, text="back", command=self.backframe)
        back.grid(row=0, column=0, sticky="w", padx=4, pady=4)

        ilabel = Label(fr, text="image number:")
        ilabel.grid(row=0, column=1, sticky="e", pady=4)

        self.evar = IntVar()
        self.evar.set(1)
        entry = Entry(fr, textvariable=self.evar, width=6)
        entry.grid(row=0, column=2, sticky="w", pady=4)
        entry.bind('<Return>', self.getimgnum)
        
        next = Button(fr, text="next", command=self.nextframe)
        next.grid(row=0, column=3, sticky="e", padx=4, pady=4)

    def backframe(self):
        index = self.im.tell()
        index = index - 1  # back up one frame
        if index < 0:
            index = self.nimages-1
        self.im.seek(index)
        self.evar.set(index+1)
        self.toframe()

    def nextframe(self):
        index = self.im.tell()
        index = index + 1  # forward one frame
        if index >= self.nimages:
            index = 0
        self.im.seek(index)
        self.evar.set(index+1)
        self.toframe()

    def toframe(self):
        # this line needs to be in a separate function (?)
        self.tkimage.paste(self.im.convert2byte())

    def getimgnum(self, event=None):
        index = self.evar.get() - 1
        if index < 0 or index >= self.nimages:
            index = self.im.tell()
            self.evar.set(index+1)
        self.im.seek(index)
        self.toframe()

# --------------------------------------------------------------------
if __name__ == "__main__":

    if not sys.argv[1:]:
        print "Usage: viewstack.py stackfile"
        sys.exit()
    filename = sys.argv[1]

    root = Tk()
    app = Viewer(root, filename)
    root.mainloop()
