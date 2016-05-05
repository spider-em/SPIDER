ss = ["INT([diam]/(2.0*[pixsize]))", "[radius] + [alignsh]", "INT([winsize]/2) - 1",
      "[maxrad] - [alignsh] - 1", "INT(([win-frac]*[winsize])/2.0 )",
      "1.0 / (2.0 * [pixsize])", "12.398 / SQR([kev] * (1022.0 + [kev]))"]
      

def removeSymbols(s):
    n = ""
    insym = 0
    for ch in s:
        if ch == '[':
            insym = 1
        elif ch == ']':
            insym = 0
            continue

        if not insym and ch != ' ':
            n += ch
    return n

for d in ss:
    print "%s : %s" % (d, removeSymbols(d))
                       
    
