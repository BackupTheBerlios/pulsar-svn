"""
$debug.py

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

 Code to handle DEBUGGING and MESSAGE information
 (adapted from some mathplotlib source code)

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""
__revision__ = "1"
__revision_date__="2007.02.02" 


import sys
import wx
from f95pulsar import parameters
from string import rstrip, lstrip, strip

DEBUG = parameters.debug
verbose = parameters.verbose

if DEBUG :
        import traceback
        import inspect
        
#------------------------------------------------------------------------------
def pydebug_msg():
#------------------------------------------------------------------------------    
# print debug-"message" from the F95PULSAR module
# Callback function called by Fortran
    try:
        verifie=inspect.isclass() # check if inspect is loaded
    except:
        import inspect
        
    string=rstrip(parameters.message.tostring())
    if parameters.debug:
        wx.LogMessage("Debug>  %s" % (string))
    
#------------------------------------------------------------------------------
def pywrite_string():
#------------------------------------------------------------------------------    
# print "message" from the F95PULSAR module
# Callback function called by Fortran
    string=rstrip(parameters.message.tostring())
    WRITE_STRING(string)

#------------------------------------------------------------------------------
def DEBUG_MSG(string):
#------------------------------------------------------------------------------
# Display a message if debug=.true.
    if DEBUG:
        wx.LogMessage("Debug>  %s" % (string))
         
#------------------------------------------------------------------------------
def debug_on_error(type, value, tb):
#------------------------------------------------------------------------------
    traceback.print_exception(type, value, tb)
    wx.LogMessage(" ") 

if DEBUG:
    sys.excepthook = debug_on_error

#------------------------------------------------------------------------------
def WRITE_VALUE(string, value, format):
#------------------------------------------------------------------------------
    wx.LogMessage("%s "+format % (string, value))

#------------------------------------------------------------------------------
def WRITE_STRING(string):
#------------------------------------------------------------------------------
    wx.LogMessage("%s " % string)
