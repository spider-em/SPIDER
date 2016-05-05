"""
This module demonstrates how to provide database connectivity for Spire.

1) there must be a class called spiderDatabase, which inherits from
   spiderClasses.SpiderDatabase. You need to import spiderClasses to
   do this. Also the Spire module GB must be imported.

2) The spiderDatabase class must have the following:
   Attributes:
         databaseName
         id
         projectquery
   Methods:
         __init__()
         isExtDatabaseAlive()
         sendQuery()
         getProjectInfo()
         upload()

3) This file must also import the Python package for the particular
   database you are using. E.g., for MySQL you need MySQLdb, for ODBC
   use the mxODBC package. See www.python.org/topics/database/modules.html
   for a list of supported databases

4) This file must be copied to spiderLocalDB.py in the
   Python_path/site-packages/Spire directory.

"""

# projectquery is a string with an SQL query for testing the database connection.
# 'Project_Title', et al. are columns in a table called Projects
projectquery = 'SELECT Project_Title, File_Extension, Working_Directory, ' \
               'Pixel_Size, Electron_Voltage, Spherical_Aberration ' \
               'FROM Projects WHERE Project_ID=112233'

import GB      # module that holds global variables for a Spire project.
import spiderClasses

# All instances of 'PythonDatabaseLib' below should be replaced by the
# Python database package you are using (e.g., MySQLdb for MySQL).
import PythonDatabaseLib 

class spiderDatabase(spiderClasses.SpiderDatabase):

    def __init__(self):
        " sets the four required attributes "
        self.databaseName = "the Laboratory Project Database"
        self.id = GB.P.ID  # GB.P holds project global variables for Spire
        self.projectquery = projectquery

    def isExtDatabaseAlive(self):
        " Whatever your database library does to check for a connection."
        " Returns 0 or 1 "
        # replace next line with actual library function
        if PythonDatabaseLib.TestConnection(): 
            return 1
        else:
            return 0

    def sendQuery(self,query):
        " Whatever your database library does to send a SQL query."
        " Returns query results. "
        # replace next line with actual library function
        query_results = PythonDatabaseLib.SendQuery(query)
        return query_results

    def getProjectInfo(self):
        """
        This function fetches project information used by Spire. It must
        return a dictionary with the keys:
        'title', 'dataext', 'projdir', 'pixelsize', 'kv', 'Cs'
        where
        title = title of the project
        dataext = SPIDER 3-letter data extension
        projdir = working directory of the project
        pixelsize = pixel size (A)
        kv = EM accelerating voltage
        Cs = spherical aberration
        """
        output = self.sendQuery(self.projectquery)
        D = {}
        D['title'] = output[0]
        D['dataext'] = output[1]
        D['projdir'] = output[2]
        D['pixelsize'] = output[3]
        D['kv'] = output[4]
        D['Cs'] = output[5]
        return D

    def upload(self, upload_object):
        """
        Given an uploadable object, establish a connection and insert data
        into the appropriate place in the external database.
        """
        # replace next line with actual library function
        PythonDatabaseLib.InsertData(upload_object)
