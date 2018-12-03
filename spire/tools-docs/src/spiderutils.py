import string, sys
import re, os, struct, time
from string import split, find

###################################################################
#
# Checking file types
#
#     istextfile()       boolean
#     isSpiderDocFile()  boolean
#     isSpiderImage()    boolean
#     isSpiderBin()      returns "image","volume","Fourier"  or 0

noNumbers = re.compile("[^\d^\s^\.^\+^\-Ee]")
text_characters = "".join(map(chr, range(32, 127)) + list("\n\r\t\b"))
_null_trans = string.maketrans("", "")

# Test for text vs binary should call istextfile(filename)
# returns 1 for text, 0 for binary
def istextfile(filename, blocksize = 512):
    if os.path.isdir(filename):
        return 0
    return istext(open(filename).read(blocksize))

def istext(s):
    if "\0" in s:
        return 0

    if not s:  # Empty files are considered text
        return 1

    # Get the non-text characters (maps a character to itself then
    # use the 'remove' option to get rid of the text characters.)
    t = s.translate(_null_trans, text_characters)

    # If more than 30% non-text characters, then
    # this is considered a binary file
    if len(t)/len(s) > 0.30:
        return 0
    return 1

# quits as soon as it gets a good data line
def isSpiderDocfile(file):
    try:
        fp = open(file, 'r')
    except:
        print 'unable to open %s' % (file)
        return 0

    comments = 0
    isDoc = 0
    blank = 0
    while 1:
        s = fp.readline()
        if s == "":  # only EOF should return blank
            break

        if len(s) > 2 and s[0] == " " and s[1] == ';':   # Spider comment line
            continue

        if noNumbers.match(s):  # if find any nondigits, +, _ etc
            isDoc = 0
            break

        ss = split(s)
        # test for new format: nums divided by blanks, 1st value is an int,
        try:
            i = int(ss[0])
            # and there are N data columns, where N = s[1]
            n = int(ss[1])
            if len(ss[2:]) == n:
                isDoc = 1
                break         # then it's new (SPIDER 11.0 Feb 2004)
        except:
            pass
        
        # see if it's the older fixed column format
        if len(s) < 13:
            isDoc = 0
            break
        try:
            key = int(s[0:6])   # 1st 6 chars are key
            n = int(s[6])       # 7th char is N
            f = float(s[7:13])   # see if there's 1 good data value
            isDoc = 1
            break
        except:
            isDoc = 0
            break
    fp.close()
    return isDoc

# -------------------- routines for binary files 

def isInt(f):
    try:
        i = int(f)
        if f-i == 0: return 1
        else:        return 0
    except:
        return 0

iforms = [1,3,-11,-12,-21,-22]

# returns iform, if t is a valid Spider header,
# otherwise returns 0
def isSpiderHeader(t):
    h = (99,) + t   # add 1 value so can use spider header index start=1
    # header values 1,2,5,12,13,22,23 should be integers
    if not isInt(h[1]): return 0
    if not isInt(h[2]): return 0
    if not isInt(h[5]): return 0
    if not isInt(h[12]): return 0
    if not isInt(h[13]): return 0
    if not isInt(h[22]): return 0
    if not isInt(h[23]): return 0
    # check iform
    iform = int(h[5])
    if not iform in iforms: return 0
    # check other header values
    labrec = int(h[13])   # no. records in file header
    labbyt = int(h[22])   # total no. of bytes in header
    lenbyt = int(h[23])   # record length in bytes
    #print "labrec = %d, labbyt = %d, lenbyt = %d" % (labrec,labbyt,lenbyt)
    if labbyt != (labrec * lenbyt): return 0
    # looks like a valid header
    return iform

# returns "image","volume","Fourier"  or 0
def isSpiderBin(filename):
    if not os.path.exists(filename):
        return 0
    minsize = 92  # 23 * 4 bytes
    if os.path.getsize(filename) < minsize:
        return 0
    try:
        fp = open(filename,'rb')
        f = fp.read(minsize)   # read 23 * 4 bytes
        fp.close()
    except:
        return 0
    bigendian = 1
    t = struct.unpack('>23f',f)    # try big-endian first
    iform = isSpiderHeader(t)
    if iform == 0:
        bigendian = 0
        t = struct.unpack('<23f',f)  # little-endian
        iform = isSpiderHeader(t)
    if iform == 0:
        return 0
    elif iform == 1:
        return "image"
    elif iform == 3:
        return "volume"
    elif iform in [-11,-12,-21,-22]:
        return "Fourier"

def isSpiderImage(file):
    if isSpiderBin(file) == "image":
        return 1
    else:
        return 0
    
# -----------------------------------------------------
def getfilenumber(filename):
    " returns string with leading zeroes "
    filename = os.path.basename(filename)
    fname,ext = os.path.splitext(filename)
    number = ""
    done = 0
    while not done:
        c = fname[-1]
        try:
            int(c)
            number = c + number
            fname = fname[:-1]
        except:
            done = 1
    return number

def fileReadLines(filename):
    try:
        fp = open(filename,'r')
        B = fp.readlines()
        fp.close()
        return B
    except IOError, e:
        print 'Unable to open file \n' + filename, e
        return None

# returns current time as tuple of 3 strings: (date, time, ID)
# e.g., ('16-OCT-03', '13:08:16', '031016130816')
def nowisthetime():
    tt = time.localtime(time.time())
    # localtime return format: (2003, 10, 16, 12, 48, 30, 3, 289, 1)
    t = string.split(time.asctime(tt))
    # asctime return format: 'Thu Oct 16 12:50:17 2003'
    mo = string.upper(t[1])
    day = t[2]
    timestr = t[3]
    yr = t[4]
    datestr = "%s-%s-%s" % (day, mo, yr)

    yr = yr[-2:]
    # this is just to get the month as a number
    d = map(str,tt)   # stringify all numbers in the tuple
    mon = d[1]
    if len(mon) < 2: mon = '0' + mon
    (h,m,s) = string.split(timestr,':')
    idstr = "%s%s%s%s%s%s" % (yr,mon,day,h,m,s)

    return (datestr, timestr, idstr)


###################################################################
#
# Reading, writing Spider document files
#
"""
Has to handle old format FORTRAN FORMAT: (I6,I1,9G12.3),
and new format: 1PGE12.5 (puts space b/w columns)

Input: filename, [list of columns] (default = read all columns)
Returns a dictionary: D[key] = [list of requested column values]

Column 0 refers to keys, column 1 to 1st SPIDER data column

TO DO: use Numpy arrays
"""
def readSpiderDocFile(filename, col_list=None):
    format = "unknown"
    if type(col_list) == type(1) and col_list != 0:
        col_list = list(col_list)
    if col_list == None: col_list = ["read_all"]

    F = fileReadLines(filename)
    if F == None: return None

    # get the first line with data
    for line in F:
        if ';' in line:
            continue
        s = split(line)
        break

    # test for new format divided by blanks
    # if 1st value is an int,
    try:
        i = int(s[0])
        # and there are N data columns, where N = s[1]
        n = int(s[1])
        if len(s[2:]) == n:
            format = "new"    # then it's new (SPIDER 11.0 Feb 2004)
        if col_list[0] == "read_all":
            col_list = range(1,n+1)
    except:
        pass
    # see if it's old format
    if format == "unknown":
        try:
            key = int(line[0:6])   # 1st 6 chars are key
            n = int(line[6])       # 7th char is N
            f = float(line[7:13])   # see if there's 1 good data value
            format = "old"
            if col_list[0] == "read_all":
                col_list = range(1,n+1)
        except:
            pass

    if format == "unknown":   # don't know what to do with it
        return None

    # ------- get the data ------
    D = {}
    if format == "new":
        for line in F:
            if ';' in line:
                continue
            s = split(line)
            key = int(s[0])
            values = []
            for column in col_list:
                i = column+1   # python index = spider column+1
                values.append(float(s[i]))
            D[key] = values

    elif format == "old":
        for line in F:
            if ';' in line:
                continue
            key = int(line[0:6])
            #n = line[6]
            values = []
            for col in col_list:
                start = 7 + (col-1)*12
                end = start+12
                values.append(float(line[start:end]))
            D[key] = values

    if len(D) == 0:
        return None
    else:
        return D

def makeDocfileHeader(filename):
    filename = os.path.basename(filename)
    ext = os.path.splitext(filename)
    ext = ext[1:]
    date,time,idstr = nowisthetime()
    h = " ;%s/%s   %s AT %s   %s\n" % (ext,ext,date,time,filename)
    return h

# -----------------------------------------------------
# data must be in dictionary form: D[key] = [list of column values]
def writeSpiderDocFile(filename, data, headers=None, append=0):
    try:
        if append == 0:
            fp = open(filename, 'w')
        else:
            fp = open(filename, 'w+')
        fname = os.path.basename(filename)
        ext = os.path.splitext(filename)[1]
        ext = ext[1:] # remove dot
        date,time,idstr = nowisthetime()
        h = " ;%s/%s   %s AT %s   %s\n" % (ext,ext,date,time,fname)
        fp.write(h)
        
        if headers != None and len(headers) > 0:
            h = " ; /    "
            for header in headers:
                hlen = len(header)
                if hlen > 11:
                    header = header[:11]
                elif hlen < 11:
                    pad = 11 - hlen
                    p = int(pad/2)
                    hd = ' '*p + header + ' '*p
                    if pad > p+p:
                        hd += ' '
                    header = hd
                h += header
            fp.write(h+"\n")

        # write data
        keys = data.keys()
        keys.sort()
        for key in keys:
            values = data[key]
            n = len(values)
            h = "%d %2d " % (key,n)
            for value in values:
                h += "%12.5f " % (float(value))
            fp.write(h+"\n")
        fp.close()
    except:
       print "unable to create %s" % filename
       return 0
    return 1

###################################################################
#
# File number routines

re_nums = re.compile('\d+\D')  # ints followed by one non-int char

def filenumber(file):
    " given 'mic021.dat' --> returns 21 "
    if len(file) == 0: return

    m = re_nums.findall(file)
    if len(m) == 1:
        n = m[0][:-1] # remove trailing char
        return int(n)
    elif len(m) > 1:  # if more than one, use n nearest the dot
        for item in m:
            if string.find(item,'.') > -1:
                n = item[:-1]
                return int(n)
            
def name2template(file):
    " given 'mic021.dat' --> returns mic***.dat "
    if len(file) == 0: return

    m = re_nums.findall(file)  # m = ['021.']
    if len(m) < 1:
        return
    elif len(m) > 1:  # if more than one, use n nearest the dot
        for item in m:
            if string.find(item,'.') > -1:
                numstr = item
    elif len(m) == 1:
        numstr = m[0]
        
    stars = "*" * (len(numstr) - 1)
    stars += "."
    a = string.find(file,numstr)  # find index
    tmp = string.replace(file,numstr,stars,1)
    return tmp

# template should have ***'s, numfile should be a numbered filename
# Input: supply either number (N) or numfile
def template2filename(template, numfile=None, n=None):
    " given (pic***.dat, doc003.dat) --> returns pic003.dat "
    if numfile != None:
        n = filenumber(numfile)
    elif n != None:
        n = int(n)
    else:
        print "template2filename: n=None, numfile=None"
        return
    nstars = string.count(template,"*")
    numstr = string.zfill(str(n), nstars)
    sts = "*" * nstars
    filename = string.replace(template,sts,numstr)
    return filename

###################################################################
#
if __name__ == '__main__':

    nargs = len(sys.argv[1:])
	
    if nargs == 0:
	    print "Usage: readSpiderDocFile.py filename '[columnlist]'"
	    sys.exit()
    elif nargs == 1:
	    filename = sys.argv[1]
	    D = readSpiderDocFile(filename)
    elif nargs > 1:
	    filename = sys.argv[1]
	    try:
                col_list = eval(sys.argv[2])
	        D = readSpiderDocFile(filename,col_list)
	    except:
                D = readSpiderDocFile(filename)
    if D != None:
        keys = D.keys()
        keys.sort()
        for key in keys:
            print "%d: %s" % (key, str(D[key]))

    writeSpiderDocFile('test.out',D,['A VERY LONG HEADER', 'shorty', '__ELEVEN___'])

        
