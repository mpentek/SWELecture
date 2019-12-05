#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Other utilities

Author: philipp.bucher@tum.de, a.winterstein@tum.de, mate.pentek@tum.de, anoop.kodakkal@tum.de
      
Note: ...

Created on:  16.01.2018
Last update: 13.11.2018
'''
#===============================================================================

import sys
import time
import math


class FileWriter:
    def __init__(self, FileName, DataNames, OpenMode="w"):
        if type(DataNames) is not list:
            raise Exception("The result column names have to be passed as list!")
        
        self.file = open(FileName, OpenMode)

        self.num_results = len(DataNames)

        for name in DataNames:
            self.file.write(str(name) + "\t")
        self.file.write("\n")
        self.file.flush()

    def WriteToFile(self, Results):
        if type(Results) is not list:
            raise Exception("The results  have to be passed as list!")
        if self.num_results != len(Results):
            raise Exception("Wrong number of results passed")
        
        for result in Results:
            self.file.write(str(result) + "\t")
        self.file.write("\n")
        self.file.flush()

    def CloseFile(self):
        self._close_file()
        
    def __del__(self): # in case the user forgets to close the file
        self._close_file()

    def _close_file(self):
        try:
            if not self.file.closed:
                self.file.close()
                print("Result file was closed")
        except AttributeError:
            print("File Object did not exist")
        

def TimeRoundValue(delta_time):
    return abs(int(math.log10(delta_time))) + 2
