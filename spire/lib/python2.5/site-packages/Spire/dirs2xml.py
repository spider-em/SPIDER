
import os, re

from Spider.Spiderutils import fileReadLines

def stripComment(line, strip=1):
    " removes Spider comment, and end whitespace from remaining text "
    n = line.find(";")
    if n > -1:
        line = line[:n].rstrip()
        
    if strip:
        line = line.strip()
    return line
# Batch file routines

re_hdr = re.compile('END +BATCH +HEADER')  
re_reg = re.compile('[xX][0-9][0-9] *=')   # "x11 =" patterns
re_bat1 = re.compile(':: *SPIDER +BATCH +FILE *::') # lower case, ignores too many spaces
re_bat2 = re.compile(':: *SPIDER +PROCEDURE +FILE *::')

# isSpiderBatchfile
# only looks for ":: spider batch file ::" or "[x11,x12]"  in first line .
# No other parsing
def isSpiderBatchfile(file):
    """ only checks first few lines of text.
        only ":: spider batch file ::"
          or ; --- End batch header ---
          or "FR G/L" pattern followed by [symbol]
        are checked
    " or ?symbol? pattern
    """
    try:
        fp = open(file, 'r')
    except:
        print 'unable to open %s' % (file)
        return 0

    max = 40
    for i in range(max):
        line = fp.readline()
        if not line:
            fp.close()
            return 0
        line = line.strip()
        line = line.upper()
        cmd = ""
        if len(line) > 3:
            cmd = line[0:4]
            
        if len(line) == 0:
            continue
        elif re_bat1.search(line) or re_bat2.search(line):
            fp.close()
            return 1
        elif re_hdr.search(line):
            fp.close()
            return 1
        # coment must come after header
        elif line[0] == ";":
            continue
        elif re_reg.match(line):
            fp.close()
            return 1
        # FR must be followed by symbol
        elif cmd == "FR G" or cmd == "FR L":
            fp.close()
            return 1
            #line = fp.readline()
            #line = line.strip()
            #if re_sym.match(line):
        # too much work to eliminate a text file
        #else:
            #fp.close()
            #return 0

re_proc_hdr  = re.compile('\[([xX][0-9][0-9],?)+\]')  # "[x12,x13,..]" patterns
re_proc_hdr2 = re.compile('\(([xX][0-9][0-9],?)+\)')  # "(x12,x13,..)" patterns
re_proc_hdr3 = re.compile('\( *(( *, *)?\[[\._\w-]+\])+ *\)')  # ([a],[b]) patterns
re_sym = re.compile('\?[\w ]+\?')   # "?symbol?" patterns
def isSpiderProcedurefile(file):
    """
    Checks 2 conditions:
    if top has either [x12] register or FR ?XXX? pattern,
    if there's a RE statement.
    
    Only checking the first few lines of text doesn't work. procs can
    start with arbitrary Spider ops.
    [x11,x12] (or (x11,x12)) pattern 
    or "FR" pattern followed by ?symbol?
    or comments or blank lines are allowed
    """
    B = fileReadLines(file)
    if B == None:
        return 0

    has_top = 0
    has_RE = 0

    # check top for registers or FR
    maxlin = min(22, len(B))
    for i in range(maxlin):
        line = B[i]
        line = stripComment(line)
        if len(line) == 0:
            continue
        # note use of "re.match" (beginning of line) vs "re.search"
        elif re_proc_hdr2.match(line) or re_proc_hdr.match(line) or re_proc_hdr3.match(line):
            has_top = 1
            break
        # FR must be followed by "?symbol?"
        elif line[0:2].upper() == "FR":
            next = B[i+1]
            next = next.strip()
            #if re_sym.match(line):  # doesn't work for ? (input) ?
            if next[0] == '?' and next[1:].find('?') > 0:
                has_top = 1
                break

    if not has_top:
        return 0

    B.reverse()
    for line in B:
        s = stripComment(line).upper()
        if s == 'RE':
            has_RE = 1
            break

    if has_top and has_RE:
        return 1
    else:
        return 0

def isBatchFile(f):
    if isSpiderBatchfile(f) or isSpiderProcedurefile(f):
        return True
    else:
        return False

def getDirFiles(dirname):
    os.chdir(dirname)
    dlist = []
    blist = []
    flist = os.listdir(os.getcwd())

    for f in flist:
        if os.path.isdir(f):
            dlist.append(f)
        elif isBatchFile(f) and f[-1] != '~':
            blist.append(f)

    d = [dirname, blist]
    for subdir in dlist:
        d[1].append([subdir])   # for now just support 1 level of recursion
    #    d.append(getDirFiles(subdir))
    os.chdir('..')
    return d


" dir is of form [dirname [ f, f, [dir], [dir],..]] "

projectDir = "/usr1/bbaxter/hcc"

os.chdir(projectDir)

dlist = []
blist = []
flist = os.listdir(os.getcwd())

for f in flist:
    if os.path.isdir(f):
        dlist.append(f)
    elif isBatchFile(f) and f[-1] != '~':
        blist.append(f)

projdirs = ['.', blist]

for subdir in dlist:
    projdirs[1].append(getDirFiles(subdir))

print projdirs


