#!/usr/bin/env python
import os, glob
from Spider.Spiderutils import *
import Tkinter
import Image, ImageTk
"""
This module writes text to a file, which is then run in spider as a batch
file. The results of processing are then displayed in a window. It looks
for a set of files matching the template 'sar***', then takes the average
of these files with the spider AS R command. 
"""
def view(filename):
    " display an image "
    im = Image.open(filename)
    root = Tkinter.Tk()
    tkimage = ImageTk.PhotoImage(im.convert2byte(), palette=256)
    Tkinter.Label(root, image=tkimage).pack()
    root.mainloop()
    
def runBatchString(bstring, ext, delete=1):
    # write string to a batchfile and run it in spider
    batchfile = 'tmp123456.bat'
    fp = open(batchfile, 'w')
    fp.write(bstring)
    fp.close()

    # run spider
    spider = findSpider()
    if len(spider) > 0: print runSpider(spider, batchfile, ext)
    else: print "unable to find spider"

    if delete: os.remove(batchfile)

# template for a batch file that calls SPIDER's averaging command, AS R
asr = """
AS R
_input
_filenums
A
_avg
_var
EN D
"""

# ---------------------------------------------------------------------
if __name__ == '__main__':

    template = "sar***"

    # create a string of the file numbers, e.g., "1-3,5,8"
    filelist = glob.glob(template) # get files that match template
    filenums = []         
    for filename in filelist:
        filenums.append(filenumber(filename))
    filestring = numberlist2string(filenums)
    print "file numbers: %s" % filestring
    file, ext = os.path.splitext(filelist[0])  # get data extension

    avg = 'avg001'
    var = 'var001'
    # replace the placeholders with filename strings
    asr = asr.replace("_input", template)
    asr = asr.replace("_filenums", filestring)
    asr = asr.replace("_avg", avg)
    asr = asr.replace("_var", var)
    print asr
    
    runBatchString(asr, ext)  # runBatchString defined above
    
    avg_file = avg + ext
    if os.path.exists(avg_file) and isSpiderImage(avg_file):
        print "Sending %s to a window..." % avg_file
        view(avg_file)
