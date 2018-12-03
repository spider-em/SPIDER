"""
This module demonstrates how to provide database connectivity for Spire.

1) there must be a class called spiderDatabase, which inherits from
   spiderClasses.SpiderDatabase. You need to import spiderClasses to
   do this. Also the Spire module GB must be imported.

2) The spiderDatabase class must have the following:
   Attributes:
         databaseName
         id
         testquery
         projectquery
   Methods:
         isExtDatabaseAlive()
         sendQuery()
         getProjectInfo()
         upload()

3) This file must also import the Python package for the particular
   database you are using. E.g., for MySQL you need MySQLdb, for ODBC
   use the mxODBC package. See www.python.org/topics/database/modules.html
   for a list of supported databases

4) This file must be copied to spiderLocalDB.py in the
   Python path/site-packages/Spire directory.

"""

# testquery is a string with an SQL query for testing the database connection.
# 'Project_Title', et al. are columns in a table called PROJECTS
testquery = 'SELECT Project_Title, File_Extension, Working_Directory, ' \
            'Pixel_Size, Electron_Voltage, Spherical_Aberration ' \
            'FROM PROJECTS WHERE Project_ID=112233'

import GB      # module that holds global variables for a Spire project.
import spiderClasses

# All instances of 'PythonDatabaseLib' below should be replaced by the
# Python database package you are using (e.g., MySQLdb for MySQL).
import PythonDatabaseLib 

class spiderDatabase(spiderClasses.SpiderDatabase):

    def processArgs(self):
        " sets the four required attributes "
        self.databaseName = "the laboratory project database"
        self.id = GB.P.ID  # GB.P holds project global variables for Spire
        self.testquery = testquery
        self.projectquery = self.testquery

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
        This function uses the results of the object's sendQuery function
        to obtain the following items for use by Spire:
        - title of the project
        - data extension
        - working directory of the project
        - pixel size (A)
        - electron voltage
        - spherical aberration
        These items are then copied into the relevant GB.P attributes.
        """
        output = self.sendQuery(self.projectquery)
        GB.P.title = output[0]
        GB.P.dataext = output[1]
        GB.P.projdir = output[2]
        GB.P.pixelsize = output[3]
        GB.P.kv = output[4]
        GB.P.Cs = output[5]

    def upload(self, upload_object):
        """
        Given an uploadable object, establish a connection and insert data
        into the appropriate place in the external database.
        """
        # replace next line with actual library function
        PythonDatabaseLib.InsertData(upload_object)
