"""
$main.py

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
  
 Welcome screen and starting procedure

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""



import wx
import os
import sys
from images import opj
from editor import PulsarNotebookFrame


#=====================================
class MySplashScreen(wx.SplashScreen):
#=====================================
    """ Introduction screen of Pulsar """
    def __init__(self, filename=None, parent=None):
        self.filename = filename
        self.parent = parent
        try:
            #it may be a standalone application
            #in this case, the PYTHONPATH is changed
            fil=opj(os.path.join(os.path.dirname(sys.executable), 'splash.bmp'))
            test=open(fil, "rb").read()
        except:
            #if not then load it as usual
            fil=opj("images\splash.bmp")
	bmp = wx.Bitmap(fil)
        wx.SplashScreen.__init__(self,
                                 bmp,
                                 wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT,
                                 3000, None, -1)
        wx.InitAllImageHandlers()        
        # create a notebook
        self.frame = PulsarNotebookFrame(filename=self.filename)
        self.Bind(wx.EVT_CLOSE, self.onclose)
        self.futurecall = wx.FutureCall(2000, self.showmain)

    def onclose(self, evt): 
        """Make sure the default handler runs too
        so this window gets destroyed"""
        evt.Skip()
        self.Hide()
        
        # if the timer is still running then go ahead and show the
        # main frame now
        if self.futurecall.IsRunning():
            self.futurecall.Stop()
            self.showmain()

    def showmain(self):
        """ Show the main frame """
        self.frame.Show()
        if self.parent is not None:
            self.parent.SetTopWindow(self.frame)
        if self.futurecall.IsRunning():
            self.Raise()
        
#==========================
class PulsarApp(wx.App):
#==========================    
    """Pulsar main application."""
    #---------------------------------
    def __init__(self, filename=None):
        self.filename = filename
        wx.App.__init__(self, redirect=False)

    #----------------    
    def OnInit(self):
        splash = MySplashScreen(self.filename, parent=self)
        splash.Show()
        return True

#---------------------------------------------------------------------------   
def main(argv):
    if len(argv) < 2:
        #Execute the interactive window application
        #------------------------------------------
    	app = PulsarApp()
    	app.MainLoop()
    
    else: 
        # Else load the file then execute the application with this name
        #----------------------------------------------------------------
        name = argv[1]
        [dirname,filename] = os.path.split(name)
        [scriptname,extension] = os.path.splitext(filename)

        if extension != '.pul':
            print "\n*** ERROR: only file with the *.pul extension can be executed by Pulsar.py ***"
            raise SystemExit
        app = Pulsar(filename=name)
        app.MainLoop()        
