#
# Spire: SPIDER Reconstruction Engine
# $Id: spiderLocalMssql.py,v 1.1 2012/07/03 14:21:59 leith Exp $
#
# functions for connecting to Microsoft SQL Server
#
# This file should contain site-specific code for establishing
# a database connection.
#
"""
Database API:

To set up communication between Spire and the external database, you must
create a SpiderDatabase class object. This object must have the methods
   isExtDatabaseAlive
   sendQuery
   getProjectInfo

from spiderClasses import SpiderDatabase

class spiderDatabase(SpiderDatabase):
    " inherits from spiderClasses.SpiderDatabase"
    def __init__(self):
        self.databaseName = "e.g., the MySQL database"
        self.id = GB.P.ID
        self.testquery = "an example SQL query"
        self.projectquery = "a query to get project info"
        self.connectProg = ""
   
    def isExtDatabaseAlive(self):
        if (some test): return 0
        else:  return 1

    def sendQuery(self, query):
        " send SQL query to database "
        pass

    def getProjectInfo(self):
        pass


isExtDatabaseAlive returns 1 if connection to database can be established.
==================
def isExtDatabaseAlive():
    if sometest():
        return 1
    else:
        return 0

getProjectInfo put data from database into SpiderProject object
==============
def getProjectInfo(project = SpiderProject.object):
    The project ID should be set before the call.
    This function is responsible for the SQL query that
    retrieves the required information.

"""
import os
import string,sys
from commands import getoutput

import GB
import spiderUtils
from spiderClasses import SpiderDatabase
import dbupload

"""
Information specific to Wadsworth

This implementation relies on 2 java programs, Connect and Upload,
which are tailored for connecting to our Ms SQLServer database.
"""

# set the CLASSPATH environmental variable.
# This only happens once, when this module is first imported
"""
classpath = "/net/sylt/usr3/yuchen/msSQLjdbc/lib:" \
            "/msbase.jar:" \
            "/net/sylt/usr3/yuchen/msSQLjdbc/lib/msutil.jar:" \
            "/net/sylt/usr3/yuchen/msSQLjdbc/lib/mssqlserver.jar"

cp = os.environ['CLASSPATH']
if classpath not in cp:
    cp = cp + ":.:" + classpath
os.environ['CLASSPATH'] = cp
"""

Fields = ['Title',
          'File_Extension',
          'Working_Directory',
          'Pixel_Size',
          'Voltage',
          'Name',
          'Short_Description',
          'Date',
          'First_Name',
          'Last_Name',
          'Specimen_Name',
          'Specimen_Description']

# extract the first set of chars that corresp to float number.
# e.g., "2.82 (scanned at 1.8)" returns 2.82
# input: string
# Returns: float, or 'None' if none found
def extractFloat(s):
    try:
        i = float(s)
        return i
    except:
        intchars = string.digits + "-."
        n = ""
        for c in s:
            if c in intchars:   # add char to the good string
                n += c
            elif len(n) == 0:  # non-int, but haven't found any yet
                continue
            else:
                break
        try:
            i = float(n)
            return i
        except:
            return None

# needs information about a specific project
testquery = 'select Title, File_Extension, ' \
        'Working_Directory, ' \
        'Pixel_Size,Voltage, Name, Short_Description, Date, First_Name, ' \
        'Last_Name,Specimen_Name, Specimen_Description from Projects, ' \
        'People, EM where Author_ID=People_ID and Microscope=EM_ID and ' \
        'Project_ID='

#######################################################################
# spiderDatabase class - must inherit from spiderClasses.SpiderDatabase

class spiderDatabase(SpiderDatabase):
    " inherits from spiderClasses.SpiderDatabase"
    def __init__(self):
        self.databaseName = "the Project Archive"
        self.id = GB.P.ID
        self.testquery = testquery + self.id
        self.projectquery = self.testquery
        self.connectProg = "java Connect"
        self.uploadProg = "java Upload"
        # Usage (with quotes): java Upload "filename" "projID" "column"
   
    def isExtDatabaseAlive(self):
        test = self.connectProg + ' "select Title from Projects where Project_ID=23"'
        x = getoutput(test)
        if string.find(x, "Exception") > -1:
            GB.errstream(x)
            return 0
        else:
            return 1

    def sendQuery(self, query):
        " send SQL query to database "
        cmd = self.connectProg + ' "%s"' % query
        x = getoutput(cmd)
        x = string.replace(x,"\t", "\n")
        x = string.replace(x,"\r", "\n")
        return(x)

    # operates on the global project (GB.P). assumes project ID has been set.
    def getProjectInfo(self):
        sp = GB.P
        id = sp.ID

        GB.outstream("Trying to establish a connection to the database for Project " + id + "...")

        # construct the command
        cmd = self.connectProg + ' "select '  # start of double quote
        for item in Fields:
            cmd += item + ", "

        cmd = cmd[:-2]  # get rid of last comma (and space)
        cmd += " from Projects, People, EM where Author_ID=People_ID and Microscope=EM_ID"
        cmd += ' and Project_ID=' + id + '"'  # end of double quote
        #print cmd

        x = getoutput(cmd)
        if string.find(x,"Exception") > -1:
            GB.errstream(" ********* DATABASE ACCESS ERROR ******* ")
            GB.errstream(x)
            return
        s = string.split(x,'\t')
        if len(s) < 6:
            GB.errstream(" ********* database info error ******* ")
            for item in s: GB.errstream(item)
            return

        #for item in s: print item

        # get the required fields (obviously this could be more efficient)
        title = s[0]
        if title != 'null' and title != "" and title != None:
            sp.title = title
            
        ext = s[1]
        if ext != 'null' and ext != "" and ext != None:
            sp.dataext = ext
            
        workdir = s[2]
        if workdir != 'null' and workdir != "" and workdir != None:
            sp.projdir = workdir
            projdir = string.split(workdir,"/")
            if projdir[0] == "" and projdir[1] == 'net':
                sp.host = projdir[2]
            elif projdir[0] == 'net':
                sp.host = projdir[1]

        pixelsize = s[3]
        if pixelsize != 'null' and pixelsize != "" and pixelsize != None:
            k = extractFloat(pixelsize)
            if k != None:
                sp.pixelsize = k
            else:
                sp.pixelsize = None

        kv = s[4]
        if kv != 'null' and kv != "" and kv != None:
            k = extractFloat(kv)
            if k != None:
                sp.kv = k
            else:
                sp.kv = None
                 
        microscope = s[5]
        """ choices in Archive Database drop-down menu:
        Philips Tecnai F20
        Philips Tecnai F30
        Philips EM 420
        HVEM
        IVEM
        Zeiss910
        """
        if microscope != 'null' and microscope != "" and microscope != None:
            if string.find(microscope, "F20") > -1:
                sp.Cs = 2.0
            elif string.find(microscope, "F30") > -1:
                sp.Cs = 2.26

        # Put the remainder in an attribute called dbinfo
        Fields[Fields.index('Name')] = 'Microscope'  # replace 'name' w/ 'microscope'
        sp.dbinfo = {'fields': Fields }
        k = 0
        for field in Fields:
            sp.dbinfo[field] = s[k]
            k += 1

    def upload(self, uploadobj, outputfunction):
        output = outputfunction

        for item in uploadobj:
            table = item[0]
            column = item[1]
            dtype = item[3]
            value = item[4]
            if dtype == 'text_file' or dtype == 'binary':
                filename = value
            
            if dtype == 'string':
                sqlstr = "UPDATE %s set %s = %s where Project_ID = %s"\
                         ";" % (table, column, value, GB.P.ID)
                res = self.sendQuery(sqlstr)
                if res != "": output(res)
                else: output("sent %s = %s" % (column, value))
                
            elif dtype == 'text_file' or dtype == 'binary':
                if not os.path.exists(filename):
                    output("unable to find %s" % filename)
                    continue
                
                cmd = '%s "%s" "%s" "%s"' % (self.uploadProg, filename, GB.P.ID, column)
                print cmd
                x = getoutput(cmd)
                if string.find(x,"Exception") > -1:
                    GB.errstream(" ********* DATABASE ACCESS ERROR ******* ")
                    GB.errstream(x)
                    output(x)
                    return
                else:
                    output("sent %s = %s" % (column, value))

##########################################################################
# upload section : transfer data to the external database
#
# specify data to be loaded. Each item is a list:
#
# [table, column, prompt, type, default]
#
#  table : database table to insert into
#  column: column_name to insert into
#  prompt: a more user-friendly descriptive label
#  type : must be one of ("string", "text_file", "binary")
#         string is a value, such as "15.6A",
#         text_file is a text file to be inserted in its entirety,
#         binary is a binary file (image, volume) to be loaded.
#  default: default file for text & binary, leave "" for string
table = 'Projects'
upload_data = [
    [table, 'Resolution',    'Resolution',    "string",    ""],
    [table, 'order_file',    'order',         "text_file", "Particles/order_picked"],
    [table, 'order_defocus', 'order defocus', "text_file", "Power_Spectra/order_defgrps"],
    [table, 'order_select',  'order select',  "text_file", "Refinement/order_select"],
    [table, 'ctfs',          'ctfs',          'text_file', "Power_Spectra/def_avg"],
    [table, 'Volume',        'Volume',        'binary',    'Refinement/final/bpr**'],
    [table, 'Final_Resolution_Curve', 'Resolution curve', 'binary', 'Refinement/plot01.gif']
    ]
    
