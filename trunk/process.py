"""
$process.py

PULSAR Project:

 Copyright (C) 2006-2007 Jean-Paul Amoureux, Christian Fernandez
 JPA - Unite de Catalyse et Chimie du Solide, Lille, France.
 CF  - Laboratoire Catalyse et Spectrochimie, Caen, France.

LICENSE:

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License (GPL)
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 To read the license please visit http://www.gnu.org/copyleft/gpl.html

PURPOSE OF THIS FILE:

 Class for the Pulsar Script Processing

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""

import wx
import os
import sys
import string
import thread
from config import PULSARPATH
from functions import *

from debug import *

#------------------------------------------------------------------------------
def runScript(name):
#------------------------------------------------------------------------------

    [dirname,filename]=os.path.split(name)
    [scriptname,extension]=os.path.splitext(filename)
    DEBUG_MSG("process RunScript: "+filename)
    # Read file
    #-----------
    WRITE_STRING(' Reading script "'+scriptname+'"\n')              
    try:
        f=open(name,"r")
        text = f.read()
    except IOError, (errno, strerror):
        WRITE_STRING("I/O error(%s): %s" % (errno, strerror))
    f.close()

    # delete old spectra *.spe
    #-------------------------
    spectrumpath=dirname+"/"+scriptname+".spe"
    if os.path.exists(spectrumpath):
        os.remove(spectrumpath)
        WRITE_STRING("Old spectrum file: "+spectrumpath+" has been removed")
    
    # Define some default variables
    #------------------------------
    os.chdir(PULSARPATH)
    imports  = 'from functions import *\n'    

    # Perform before run commands
    #----------------------------

    header  = '_script="' + string.join(dirname.split('\\'),'/') + '/' + scriptname + '"'+ '\n'
    # to avoid a bug with filename all \ are replaced with /.
    header += "\nWRITE_STRING('SCRIPT PROCESS: Started...')\n"

    # Add the header to the script
    #-----------------------------
    text = header + text
    
    # Perform after run commands
    #---------------------------
    # text += "\nsim1.write_spectra(_script)\n"
    text +="\ntry:"
    text +="\n   os.rename('scratch.tmp',_script+'.spe')"
    text +="\nexcept:"
    text +="\n   WRITE_STRING('ERROR: '+_script+'.spe cannot be created.')" 
    text +="\n   WRITE_STRING('   May be you forgot the instruction sim.write_spectra() at the end of your file')"
    text +="\n   cancel=True"
    text +="\nWRITE_STRING('SCRIPT PROCESS: End')\n"
    text = text.replace('\r', '\n')
    text = text.replace('\n\n', '\n')    

    # create a temporary copy of the executed file
    #---------------------------------------------
    WRITE_STRING(' Creating temporary script file"'+scriptname+'.tmp"\n')              
    try:
        f=open("workspace\\"+scriptname+".tmp","wb") 
        f.write(text)
    except IOError, (errno, strerror):
        WRITE_STRING("I/O error(%s): %s" % (errno, strerror))
        f.close()
        return False

    f.close()

    max = 2000   #TODO: try to find a correct estimate of the running time from the parameters readed in the script
    dlg = wx.ProgressDialog("Script execution",
                           "Please wait...",
                           maximum = max,
                           parent= None,
                           style = wx.PD_CAN_ABORT
                            | wx.PD_APP_MODAL
                            | wx.PD_ELAPSED_TIME
                            #| wx.PD_ESTIMATED_TIME
                            #| wx.PD_REMAINING_TIME
                            )

    count = 0
    keepGoing = True
    iserror=False

    # Wait the script finish its execution
    #--------------------------------------
    thrd=PulsarProcess(text)
    thrd.Start()

    try:
        while keepGoing:
            count += 1
            if count>=max : count=max-1999
            wx.MilliSleep(1)
            if thrd.IsRunning():
                keepGoing=dlg.Update(count)
            else:
                keepGoing=False
            
    # handle runtime errors (TODO: Complete this)
    #--------------------------------------------
    except PulsarError, e:
        if e.value==1000:
            WRITE_STRING("***** User asked for process termination *****\n")
            iserror=True
        if e.value==10:
            iserror=True

    except:
        WRITE_STRING("*** Unknown error ***\n")
        iserror=True
        
    # Remove dialog box
    #------------------
    dlg.Destroy()

    # delete the temporary file
    #--------------------------    
    #os.remove(scriptname+".tmp")

    DEBUG_MSG("RunScript Ended")

    return not iserror   # return true if there is no error
      
#===================
class PulsarProcess:
#===================
    """
    Run script
    """
    
    def __init__(self, script):
        self.script=script
        
    def Start(self):
        self.running = True
        thread.start_new_thread(self.Run, ())

    def Stop(self):
        self.running = False

    def IsRunning(self):
        return self.running

    def Run(self):
        try:
            code=compile(self.script,'<script>','exec')
            iserror=False
            exec(code)
            self.running=False
            
        except SyntaxError, (message,(filename, lineno, offset, line)):
            #atttributes filename, lineno, offset and text for easier access to the details.
            #str() of the exception instance returns only the message. 
            WRITE_STRING("*** ERROR *** "+message+" script")
            WRITE_STRING("\tthe error is in line "+str(lineno)+" or above!")
            WRITE_STRING("\t\t  "+line.replace('\n',''))
            WRITE_STRING("\t\t"+offset*" "+"^")
            iserror=True
            self.running=False   

        except NameError, details:
            #atttributes filename, lineno, offset and text for easier access to the details.
            #str() of the exception instance returns only the message. 
            WRITE_STRING("*** ERROR *** "+str(details))
            iserror=True
            self.running=False   

        except:
            WRITE_STRING("*** ERROR *** non-handled error (please report this to the developpers)")
            iserror=True
            self.running=False   
            raise          

        
        
