# -*- coding: utf-8 -*-
"""
This module provides a set of tools for reading NETL's material property
databse (*.accdb file) or a convenient json file.
"""

try:
    import pypyodbc
except ImportError:
    pypyodbc = False

import json
import os
import datetime

class ParticleDatabase(object):
    """
    This class provides a set of methods for reading NETL's material property
    database either in the MS Access (*.accdb) or the convenient json format.
    """
    def __init__(self, databasepath=None):
        self.databasepath = databasepath
        self.database = {}
        self.__idindex = {}
        self.materials = None

        if self.databasepath:
            self.readdatabase()

    def save2json(self, path='./Material_Properties.json'):
        """
        Save the mmaterial databse to a *.json file.
        """

        dthandler = lambda obj: obj.isoformat() if isinstance(obj, datetime.datetime) else json.JSONEncoder().default(obj)

        json.dump(self.database, open(path, 'w'), default=dthandler)

    def readdatabase(self, databasepath=None):
        """
        Read an MS Access (*.accdb) database file or json (*.json) file.
        """
        if databasepath:
            self.databasepath = databasepath

        # Check extension
        name, ext = os.path.splitext(self.databasepath)

        if ext == '.accdb' or ext == '.mdb':
            if not pypyodbc:
                raise IOError('pypyodbc library not avaliable.')

            self.__connect2database()
            self.__buildindex()
            self.__addsize()
            self.__adddensity()
            self.__addumf()

            self.connection.close()

        elif ext == '.json':
            self.__readjson()

        else:
            raise IOError('File extension not recognized.')

        self.materials = self.database.keys()

    def __readjson(self):
        self.database = json.load(open(self.databasepath,'r'))

    def __connect2database(self):
        connection_string = 'Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=' + self.databasepath
        self.connection = pypyodbc.connect(connection_string)

        self.cursor = self.connection.cursor()

    def __buildindex(self):

        self.cursor.execute('''SELECT * FROM Mat_Index;''')

        cols = [d[0] for d in self.cursor.description]

        idindx = cols.index('id')
        matindx = cols.index('material')
        geldindx = cols.index('geldart')
        noteindx = cols.index('notes')

        for row in self.cursor.fetchall():
            name = ': '.join(['{0:03d}'.format(row[idindx]), row[matindx]])
            self.database[name] = {'id':       row[idindx],
                                   'material': row[matindx],
                                   'geldart':  row[geldindx],
                                   'notes':    row[noteindx],
                                   }

            self.__idindex[row[idindx]] = name

    def __addsize(self):
        self.cursor.execute('''SELECT * FROM Size;''')
        cols = [d[0] for d in self.cursor.description]

        colindx = {}
        for i, col in enumerate(cols):
            colindx[col] = i

        cols.pop(cols.index('id'))

        for row in self.cursor.fetchall():
            for col in cols:
                self.database[self.__idindex[row[colindx['id']]]][col] = row[colindx[col]]

    def __addumf(self):
        self.cursor.execute('''SELECT * FROM Fluidization;''')
        cols = [d[0] for d in self.cursor.description]

        colindx = {}
        for i, col in enumerate(cols):
            colindx[col] = i

        cols.pop(cols.index('id'))

        for row in self.cursor.fetchall():
            for col in cols:
                self.database[self.__idindex[row[colindx['id']]]][col] = row[colindx[col]]

    def __adddensity(self):
        self.cursor.execute('''SELECT * FROM Densities;''')
        cols = [d[0] for d in self.cursor.description]

        colindx = {}
        for i, col in enumerate(cols):
            colindx[col] = i

        cols.pop(cols.index('id'))

        for row in self.cursor.fetchall():
            for col in cols:
                self.database[self.__idindex[row[colindx['id']]]][col] = row[colindx[col]]

    def __getitem__(self, key):
        return self.database[key]

    def keys(self):
        return self.database.keys()


if __name__ == '__main__':

    # Read MS Access database
    databasepath = r'K:\Common\PSDF\Carbon Capture\C2U\Experiments\Bed Materials\Material_Properties.accdb'
    partdatabase = ParticleDatabase(databasepath)

    # Save it to a json file
    partdatabase.save2json('./Material_Properties.json')

    # Read json file
    databasepath = './Material_Properties.json'
    partdatabase = ParticleDatabase(databasepath)

    # Print material names
    print partdatabase.materials

    # Print material data
    print partdatabase['079: Sorbent BN v2']
