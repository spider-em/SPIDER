#!/usr/bin/env python

"""
Example program that connects to a MySQL database. First it makes
sure the database exists (if not, it is created). Then it creates a
table in the database, and puts some data in the table.
*** If the table already exists, this program will try to delete it! ***

Usage: testmysql.py [ -u user -p password -d database -t table]
Usage: testmysql.py [ table database user password ]

All arguments are optional; they must be in the above order. The defaults are:
table = "projects"
database = "test_spire"
user = ""   i.e., anonymous login must be enabled for mysql
password = ""

You must first: 
 - have MySQL installed, 
 - have the MySQL server started, 
 - have the Python mySQL interface, MySQLdb, installed
   (download from sourceforge.net/projects/mysql-python),
 - the MySQL database should have anonymous login allowed,
   or the username must already exist.

Using MySQLdb in Python:
   After importing MySQLdb, create a connection object:
    import MySQLdb
    conn=MySQLdb.connect(db=_DATABASE, user=_USER, passwd=_PASSWORD)

   Then you can perform queries and fetch the results:
    c = conn.cursor()
    query = "SELECT spam, eggs, sausage FROM breakfast_table"
    c.execute("%s" % (query,) )
    results = c.fetchall()

   See the MySQLdb documentation.
"""

import MySQLdb
import os, sys

_VERBOSE = 1
_DATABASE = "test_spire"
_TABLE = "projects"
_USER = ""
_PASSWORD = ""

class TestData:
    def __init__(self, database=_DATABASE, table=_TABLE):
        self.database = database
        self.tablename = table

        self.table  =  "(project_id    VARCHAR(12)  NULL,\
                    pixelsize     VARCHAR(12)  NULL,\
                    kev           VARCHAR(12)  NULL,\
                    title         VARCHAR(64)  NULL,\
                    project_file  VARCHAR(24)  NULL,\
                    data_ext      VARCHAR(3)   NULL,\
                    host          VARCHAR(24)  NULL,\
                    projdir       VARCHAR(64)  NULL)" 


        self.data  = [('project_id', '149'),
                 ('pixelsize', '4.78'),
                 ('kev', '100'),
                 ('title', "Project for testing MySQL connection"),
                 ('project_file', 'proj149'),
                 ('data_ext', 'xyz'),
                 ('host', os.uname()[1]),
                 ('projdir', os.path.join(os.getcwd(),'xyz'))
                 ]


def execute_query(db, qstring):
	c=db.cursor()
	c.execute("%s" % (qstring,) )
	q = c.fetchall()
	c.close()
	return q

   
def hasDatabase(db, database):
    ret = execute_query(db, "SHOW DATABASES")
    # returns list of tuples:   (('mysql',), ('sample',), ('test',))
    dstr = "('%s',) in ret" % database
    return eval(dstr)

def createDatabase(db, database, verbose=_VERBOSE):
    " create database, if it doesn't exist "
    if not hasDatabase(db, database):
        if verbose:
            print "creating database '%s'..." % database
        qstr = "CREATE DATABASE %s" % database
        execute_query(db, qstr)
    else:
        if verbose:
            print "database '%s' is present." % database

# Returns a list of tables 
def get_tables(db):	
    T = execute_query(db, "SHOW TABLES")
    Tables = []
    for item in T:
            Tables.append(item[0])   # list of tables
    return Tables
           
# create table. If it already exists, overwrite it.

def createTable(db, tablename, table, verbose=_VERBOSE):
    " create a table in database db. If it already exists, overwrite it "
    # well, check first
    tabs = get_tables(db)
    if tablename in tabs:
        print "Database %s already has a table named %s" % (td.database, tablename)
        res = raw_input("Delete %s? (y/n): " % tablename)
        res = res.strip()
        if len(res) > 0:
            if res[0] == 'n' or res[0] == 'N':
                return 0

    execute_query(db, "DROP TABLE IF EXISTS %s" % (tablename,) )
    # table structure
    qstr = "CREATE TABLE %s " % tablename
    qstr += table

    if verbose:
        print "creating table '%s' :" % tablename
    execute_query(db, qstr)
    # show the created table
    if verbose:
        show_columns(db, tablename)

    return 1


def show_columns(db, table):
    hline  = "+--------------+-------------+------+-----+---------+-------+"
    header = "| Field        | Type        | Null | Key | Default | Extra |"
    blank  = "|              |             |      |     |         |       |"

    print hline
    print header
    print hline

    ret = execute_query(db, "show columns from %s" % table)
    b = blank[1:-1].split('|')  # get a list of blank fields

    # each item in d gets put in a blank
    for item in ret:
        print insert_value(b, item)

    print hline
    
def insert_value(blank, values):
    # blank is string of blanks: ['  ', '  ', '   ']
    # values is tuple of values:
    line = ""
    k = 0
    for value in values:
        spaces = blank[k]
        if value == None:
            value = "NULL"
        ns = len(spaces)
        nv = len(value)
        trail = ns - (nv+1)
        new = "| " + value + trail*" "
        line += new
        k += 1
    line += "|"
    return line

# Load data into table. ----------------------------------
valueColumns = ['project_id', 'pixelsize', 'kev', 'resolution']
textFileColumns = ['res_plot', 'res_curve']
binaryFileColumns = ['res_image']
textColumns = ['title', 'project_file', 'data_ext', 'host']

def data2list(data):
    # convert 'project_id  149\n  pixelsize   4.78\n ...'
    # to      [['project_id', '149'], ['pixelsize', '4.78']]
    d = data.split("\n")
    B = []
    for item in d:
        p = item.split()
        if len(p) > 0:
            B.append(p)
    return B

def load_table(db, table, data):
    " Put some test data into the database table"
    #B = data2list(data)  # data can be read from a text file
    B = data
    # B is a list of [column_name, value] pairs
    
    # INSERT INTO table (column1, column2,..) VALUES (val1, val2, ..)
    curs=db.cursor()
    columns = "("
    values = "("
    projid = None
    #textFiles = {}  # dictionary of textfiles
    
    for item in B:
        col = item[0]
        val = item[1]
        if col in valueColumns:   
            columns += '%s,' % col
            values += '"%s",' % val
            if col == 'project_id':
                projid = val
                
        elif col in textColumns:   
            columns += '%s,' % col
            values += '"%s",' % val
            if col == 'project_id':
                projid = val
                
        elif col in textFileColumns:
            filename = val
            if not os.path.exists(filename):
                print filename + " not found"
                continue
            fp = open(filename)
            F = fp.read()  # read in entire file as a string
            fp.close()
            #textFiles[col] = F # put in dictionary, accessed by column name
            columns += '%s,' % col
            values += '"%s",' % MySQLdb.string_literal(F)
            
        elif col in binaryFileColumns:
            filename = val
            if not os.path.exists(filename):
                print filename + " not found"
                continue
            fp = open(filename, 'rb')
            F = fp.read()  # read in entire file as a string
            fp.close()
            #textFiles[col] = F # put in dictionary, accessed by column name
            columns += '%s,' % col
            values += '"%s",' % MySQLdb.string_literal(F)

    # insert the values into the table
    columns = columns[:-1] + ")"
    values = values[:-1] + ")"  
    sqlstr = 'INSERT INTO %s %s VALUES %s' % (table, columns, values)
    #print sqlstr
    curs.execute("%s" % (sqlstr,) )
   
    curs.close()
    # show data in table
    showdata(db, table)
    
def showdata(db, table):
    hline  = "+------------+-----------+------+--------------+----------+"
    header = "| project_id | pixelsize | kev  | project_file | data_ext |"
    blank  = "|            |           |      |              |          |"

    print hline
    print header
    print hline

    ret = execute_query(db, "select project_id, pixelsize, kev, project_file, data_ext from %s" % table)
    b = blank[1:-1].split('|')

    # each item in d gets put in a blank
    for item in ret:
        print insert_value(b, item)
    print hline

def processargs(arglist):
    " -d database -uuser -> {'-d': 'database', '-u': 'user'} "
    D = {}
    while 1:
        if len(arglist) < 1:
            return D
        
        a = arglist[0]
        if a[0] == '-':
            if len(a) == 2:
                if len(arglist) < 2:
                    return D
                D[a] = arglist[1]
                arglist = arglist[2:]
            elif len(a) > 2:
                # there might not be a space between flag and arg (-uuser)
                flag = a[0:2]
                val = a[2:]
                D[flag] = val
                arglist = arglist[1:]
        else:
            # there's no flag
            arglist = arglist[1:]

# --------------------------------------------------------	
if __name__ == '__main__':

    args = sys.argv[1:]

    if '-h' in args or '-help' in args:
        print "Usage: testmysql.py [-u user -p password -d database -t table]"
        sys.exit()

    if len(args) > 0:
        # default values for _USER, etc. are set above
        # dictionary has {flag : value} pairs
        argdict = processargs(sys.argv[1:])
        if '-u' in argdict: _USER = argdict['-u']
        if '-p' in argdict: _PASSWORD = argdict['-p']
        if '-d' in argdict: _DATABASE = argdict['-d']
        if '-t' in argdict: _TABLE = argdict['-t']

    # create an instance of the test data to put in  the new table
    td = TestData(database=_DATABASE, table=_TABLE)
        
    # create a connection object
    db=MySQLdb.connect(user=_USER, passwd=_PASSWORD)

    createDatabase(db, td.database)
    execute_query(db, "USE " + td.database)
    created = createTable(db, td.tablename, td.table)

    if created:
        print "\nSome data loaded into table '%s':" %  td.tablename
        load_table(db, td.tablename, td.data)
    
    db.close()
