"""
$editor.py

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

 Class for the Pulsar Script Editor

TODO: Make a lot of cleaning of this file (that was mostly adapted from PyCrust)

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""



import wx

if wx.Platform == '__WXMSW__':
    import wx.lib.iewin    as  iewin  #ActiveX IEtml

import wx.html
    
import sys, os
import frame as frame
import version
import editwindow as editwindow
import dispatcher as dispatcher
from buffer import Buffer
from plot import *

DEBUG = False
MOREDEBUG = False
SERIOUSDEBUG = False
PULSARPATH,PULSARNAME=os.path.split(sys.argv[0])

#==============================
class EditorFrame(frame.Frame):
    """Frame containing one editor."""

    def __init__(self, parent=None, id=-1, title='pyPULSAR '+version.__PROJECT_VERSION__+" rev."+version.__PROJECT_REVISION__,
                 pos=wx.DefaultPosition, size=(800, 600), 
                 style=wx.DEFAULT_FRAME_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE,
                 filename=None):
        """Create EditorFrame instance."""
        if DEBUG : print "EditorFrame"

        frame.Frame.__init__(self, parent, id, title, pos, size, style)
        self.buffers = {}
        self.buffer = None  # Current buffer.
        self.editor = None
        self._defaultText = title + ' - A Solid-state NMR simulation package'
        self._statusText = self._defaultText
        self.SetStatusText(self._statusText)
        wx.EVT_IDLE(self, self.OnIdle)
        self._setup()
        if filename:
            self.bufferCreate(filename)

    def _setup(self):
        """Setup prior to first buffer creation.
        Useful for subclasses."""
        if DEBUG : print "EditorFrame","_setup" 
        pass

    def setEditor(self, editor):
        """setEditor"""
        if DEBUG : print "EditorFrame","setEditor"
         
        self.editor = editor
        self.buffer = self.editor.buffer
        self.buffers[self.buffer.id] = self.buffer

##    def OnAbout(self, event):
##        """Display an About window."""
##        if DEBUG : print "EditorFrame","onAbout"
##         
##        title = 'About pyPulsar'
##        text = title='pyPULSAR '+version.__PROJECT_VERSION__+" rev."+version.__PROJECT_REVISION__+'\n\n'+ \
##             'A solid State NMR Simulation program\n'+ \
##             'by C.Fernandez and J.P. Amoureux.\n\n\n'+ \
##             '(www-lcs.ensicaen.fr)\n\n'+ \
##             '(For more information: Look at the documentation) '
##        
##        dialog = wx.MessageDialog(self, text, title,
##                                  wx.OK | wx.ICON_INFORMATION)
##        dialog.ShowModal()
##        dialog.Destroy()

    def OnClose(self, event):
        """Event handler for closing."""
        if DEBUG : print "EditorFrame","onClose"
 
        for buffer in self.buffers.values():
            self.buffer = buffer
            if buffer.hasChanged():
                cancel = self.bufferSuggestSave()
                if cancel and event.CanVeto():
                    event.Veto()
                    return
        self.Destroy()

    def OnIdle(self, event):
        """Event handler for idle time."""
        if MOREDEBUG : print "EditorFrame","OnIdle",'#'
         
        self._updateStatus()
        if hasattr(self, 'notebook'):
            self._updateTabText()
        self._updateTitle()
        event.Skip()

    def _updateStatus(self):
        """Show current status information."""
        if MOREDEBUG : print "EditorFrame","_updateStatus",'#'
        if self.editor and hasattr(self.editor, 'getStatus'):
            status = self.editor.getStatus()
            text = 'File: %s  |  Line: %d  |  Column: %d' % status
        else:
            text = self._defaultText
        if text != self._statusText:
            self.SetStatusText(text)
            self._statusText = text

    def _updateTabText(self,modified=None):
        """Show current buffer information on notebook tab."""
        if MOREDEBUG : print "EditorFrame","_updateTabText",'#'
        suffix = ' *'
        size=2
        notebook = self.notebook
        selection = notebook.GetSelection()
        if selection == -1:
            return
        text = notebook.GetPageText(selection)

        if modified is not None:
            modified=os.path.split(modified)[1]
            notebook.SetPageText(selection, modified)
            
        window = notebook.GetPage(selection)
        if window.editor and window.editor.buffer.hasChanged():
            if text.endswith(suffix):
                pass
            else:
                notebook.SetPageText(selection, text + suffix)
        else:
            if text.endswith(suffix):
                notebook.SetPageText(selection, text[:(len(text)-size)])

    def _updateTitle(self):
        """Show current title information."""
        if DEBUG : print "EditorFrame","_updateTitle" 
        title = self.GetTitle()
        if self.bufferHasChanged():
            if title.startswith('* '):
                pass
            else:
                self.SetTitle('* ' + title)
        else:
            if title.startswith('* '):
                self.SetTitle(title[2:])
        
    def hasBuffer(self):
        """Return True if there is a current buffer."""
        if DEBUG : print "EditorFrame","hasBuffer",'#'
         
        if self.buffer:
            return True
        else:
            return False

    def bufferClose(self):
        """Close buffer."""
        if DEBUG : print "EditorFrame","bufferClose"
        if self.bufferHasChanged():
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        self.bufferDestroy()
        cancel = False
        return cancel

    def bufferCreate(self, filename=None):
        """Create new buffer."""
        if DEBUG : print "EditorFrame","bufferCreate"
         
        self.bufferDestroy()
        buffer = Buffer()
        self.panel = panel = wx.Panel(parent=self, id=-1)
        wx.EVT_ERASE_BACKGROUND(panel, lambda x: x)        
        editor = Editor(parent=panel)
        panel.editor = editor
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(editor.window, 1, wx.EXPAND)
        panel.SetSizer(sizer)
        panel.SetAutoLayout(True)
        sizer.Layout()
        buffer.addEditor(editor)
        buffer.open(filename)
        self.setEditor(editor)
        self.editor.setFocus()
        self.SendSizeEvent()      

    def bufferDestroy(self):
        """Destroy the current buffer."""
        if DEBUG : print "EditorFrame","bufferDestroy" 
        if self.buffer:
            for editor in self.buffer.editors.values():
                editor.destroy()
            self.editor = None
            del self.buffers[self.buffer.id]
            self.buffer = None
            self.panel.Destroy()

    def bufferHasChanged(self):
        """Return True if buffer has changed since last save."""
        if MOREDEBUG : print "EditorFrame","bufferHasChanged","#" 
        if self.buffer:
            return self.buffer.hasChanged()
        else:
            return False

    def bufferNew(self):
        """Create new buffer."""
        if DEBUG : print "EditorFrame","bufferNew" 
        if self.bufferHasChanged():
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        self.bufferCreate()
        cancel = False
        return cancel

    def bufferOpen(self):
        """Open file in buffer."""
        if DEBUG : print "EditorFrame","bufferOpen" 
        if self.bufferHasChanged():
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        filedir = ''
        if self.buffer and self.buffer.doc.filedir:
            filedir = self.buffer.doc.filedir
        result = openSingle(directory=filedir)
        if result.path:
            self.bufferCreate(result.path)
        cancel = False
        return cancel

    def bufferSave(self):
        """Save buffer to its file."""
        if DEBUG : print "EditorFrame","bufferSave"
        if self.buffer.doc.filepath:
            self.buffer.save()
            cancel = False
        else:
            cancel = self.bufferSaveAs()
        return cancel

    def bufferSaveAs(self):
        """Save buffer to a new filename."""
        if DEBUG : print "EditorFrame","bufferSaveAs" 
        if self.bufferHasChanged() and self.buffer.doc.filepath:
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        filedir = ''
        if self.buffer and self.buffer.doc.filedir:
            filedir = self.buffer.doc.filedir
        result = saveSingle(directory=filedir)
        if result.path:
            self.buffer.saveAs(result.path)
            #suppress *.spe
            [scriptname,extension] = os.path.splitext(result.path)
            f=scriptname+".spe"
            if os.path.exists(f):
                print f+" removed"
                os.remove(f)
            else:
                print f+" not found"
            self.bufferClose()
            self.bufferCreate(result.path)
            cancel = False
        else:
            cancel = True
        return cancel

    def bufferSuggestSave(self):
        """Suggest saving changes.  Return True if user selected Cancel."""
        if DEBUG : print "EditorFrame","bufferSuggestSave" 
        result = messageDialog(parent=None,
                               message='%s has changed.\n'
                                       'Would you like to save it first'
                                       '?' % self.buffer.name,
                               title='Save current file?')
        if result.positive:
            cancel = self.bufferSave()
        else:
            cancel = result.text == 'Cancel'
        return cancel


#======================================
class EditorNotebookFrame(EditorFrame):
    """Frame containing one or more editors in a notebook."""

    def __init__(self, parent=None, id=-1, title='pyPULSAR ' \
		    +version.__PROJECT_VERSION__+" rev."+version.__PROJECT_REVISION__,
                 pos=wx.DefaultPosition, size=(800, 600), 
                 style=wx.DEFAULT_FRAME_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE,
                 filename=None):
        """Create EditorNotebookFrame instance."""
        if DEBUG : print "EditorNotebookFrame"

        self.notebook = None
        EditorFrame.__init__(self, parent, id, title, pos,
                             size, style, filename)
        if self.notebook:
            dispatcher.connect(receiver=self._editorChange,
                               signal='EditorChange', sender=self.notebook)
        

    def _editorChange(self, editor):
        """Editor change signal receiver."""
        if DEBUG : print "EditorNotebookFrame","_editorChange" 
        self.setEditor(editor)

    def _updateTitle(self):
        if MOREDEBUG : print "EditorNotebookFrame","_updateTitle","#"
        title = self.GetTitle()
        if self.bufferHasChanged():
            if title.startswith('* '):
                pass
            else:
                self.SetTitle('* ' + title)
        else:
            if title.startswith('* '):
                self.SetTitle(title[2:])
        
    def bufferCreate(self, filename=None):
        """Create new buffer."""
        if DEBUG : print "EditorNotebookFrame","bufferCreate" 
        buffer = Buffer()
        panel = wx.Panel(parent=self.notebook, id=-1)
        wx.EVT_ERASE_BACKGROUND(panel, lambda x: x)        
        editor = Editor(parent=panel)
        panel.editor = ed



        itor
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(editor.window, 1, wx.EXPAND)
        panel.SetSizer(sizer)
        panel.SetAutoLayout(True)
        sizer.Layout()
        buffer.addEditor(editor)
        buffer.open(filename)
        self.setEditor(editor)
        self.notebook.AddPage(page=panel, text=self.buffer.name, select=True)
        self.editor.setFocus()

    def bufferDestroy(self):
        """Destroy the current buffer."""
        if DEBUG : print "EditorNotebookFrame","bufferDestroy" 
        selection = self.notebook.GetSelection()
        if DEBUG: print "Destroy Selection:", selection
        if selection > 0:  # Don't destroy the PyCrust tab.
            if self.buffer:
                del self.buffers[self.buffer.id]
                self.buffer = None  # Do this before DeletePage().
            self.notebook.DeletePage(selection)

    def bufferNew(self):
        """Create new buffer."""
        if DEBUG : print "EditorNotebookFrame","bufferNew" 
        self.bufferCreate()
        cancel = False
        return cancel

    def bufferOpen(self):
        """Open file in buffer."""
        if DEBUG : print "EditorNotebookFrame","bufferOpen" 
        filedir = ''
        if self.buffer and self.buffer.doc.filedir:
            filedir = self.buffer.doc.filedir
        result = openMultiple(directory=filedir)
        for path in result.paths:
            self.bufferCreate(path)
        cancel = False
        return cancel

#=================================
class EditorNotebook(wx.Notebook):
    """A notebook containing a page for each editor."""

    def __init__(self, parent):
        """Create EditorNotebook instance."""
        if DEBUG : print "EditorNotebook"
        wx.Notebook.__init__(self, parent, id=-1, style=wx.CLIP_CHILDREN)
        wx.EVT_NOTEBOOK_PAGE_CHANGING(self, self.GetId(),
                                      self.OnPageChanging)
        wx.EVT_NOTEBOOK_PAGE_CHANGED(self, self.GetId(),
                                     self.OnPageChanged)
        wx.EVT_IDLE(self, self.OnIdle)
        
    def OnIdle(self, event):
        """Event handler for idle time."""
        if MOREDEBUG : print "EditorNotebook","OnIdle","#"  
        self._updateTabText()
        event.Skip()

    def _updateTabText(self):
        """Show current buffer display name on all but first tab."""
        if MOREDEBUG : print "EditorNotebook","_updateTabText",'#'  
        size = 2
        changed = ' *'
        unchanged = '  '
        selection = self.GetSelection()
        if selection < 1:
            return
        text = self.GetPageText(selection)
        window = self.GetPage(selection)
        if not window.editor:
            return
        if text.endswith(changed) or text.endswith(unchanged):
            name = text[:-size]
        else:
            name = text
        if name != window.editor.buffer.name:
            text = window.editor.buffer.name
        if window.editor.buffer.hasChanged():
            if text.endswith(changed):
                text = None
            elif text.endswith(unchanged):
                text = text[:-size] + changed
            else:
                text += changed
        else:
            if text.endswith(changed):
                text = text[:-size] + unchanged
            elif text.endswith(unchanged):
                text = None
            else:
                text += unchanged
        if text is not None:
            self.SetPageText(selection, text)
            self.Refresh()  # Needed on Win98.

    def OnPageChanging(self, event):
        """Page changing event handler."""
        if DEBUG : print "EditorNotebook","OnPageChanging"  
        event.Skip()

    def OnPageChanged(self, event):
        """Page changed event handler."""
        if DEBUG : print "EditorNotebook","OnPageChanged"  
        new = event.GetSelection()
        window = self.GetPage(new)
        dispatcher.send(signal='EditorChange', sender=self,
                        editor=window.editor)
        window.SetFocus()
        event.Skip()

#==============================================
class PulsarNotebookFrame(EditorNotebookFrame):
    """Frame containing a notebook containing EditorShellNotebooks."""

    def __init__(self, parent=None, id=-1, title='pyPULSAR ' \
		    +version.__PROJECT_VERSION__+" rev."+version.__PROJECT_REVISION__,
                 pos=wx.DefaultPosition, size=(750, 550), 
                 style=wx.DEFAULT_FRAME_STYLE,
                 filename=None, singlefile=False):
        
        """Create EditorShellNotebookFrame instance. *** """
        if DEBUG : print "PulsarNotebookFrame"
        self._singlefile = singlefile
        
        EditorNotebookFrame.__init__(self, parent, id, title, pos,
                                     size, style, filename)
       
    def _setup(self):
        """Setup prior to first buffer creation.

        Called automatically by base class during init. *** """
        if DEBUG : print "PulsarNotebookFrame._setup"  
        if not self._singlefile:
            self.notebook = EditorNotebook(parent=self)

##    def OnAbout(self, event):
##        """Display an About window."""
##        if DEBUG : print "PulsarNotebookFrame","OnAbout"  
##        title = 'About pyPulsar'
##        text ='pyPULSAR '+version.__PROJECT_VERSION__+" rev."+version.__PROJECT_REVISION__+'\n\n'+ \
##             '<b>A solid state NMR Simulation program</b>\n'+ \
##             'by C.Fernandez and J.P. Amoureux.\n\n\n'+ \
##             '(http:\\www-lcs.ensicaen.fr\pyPulsar)\n\n'+ \
##             'Documentation: (Documentation) '
##        dialog = wx.MessageDialog(self, text, title,
##                                  wx.OK | wx.ICON_INFORMATION)
##        dialog.ShowModal()
##        dialog.Destroy()

    def OnAbout(self, event):
        dlg = PulsarAbout(self)
        dlg.ShowModal()
        dlg.Destroy()

    def OnDoc(self, event):
        dlg = PulsarDoc(self)
        dlg.Show(True)
        
    def bufferCreate(self, filename=None):
        """Create new buffer."""
        if DEBUG : print "PulsarNotebookFrame","bufferCreate" 
        if self._singlefile:
            self.bufferDestroy()
            notebook = EditorShellNotebook(parent=self,
                                           filename=filename)
            self.notebook = notebook
        else:
            notebook = EditorShellNotebook(parent=self.notebook,
                                           filename=filename)
        self.setEditor(notebook.editor)
        if not self._singlefile:
            self.notebook.AddPage(page=notebook, text=self.buffer.name,
                                  select=True)
        self.editor.setFocus()

    def bufferDestroy(self):
        """Destroy the current buffer."""
        if DEBUG : print "PulsarNotebookFrame","bufferDestroy"  
        if self.buffer:
            self.editor = None
            del self.buffers[self.buffer.id]
            self.buffer = None  # Do this before DeletePage().
        if self._singlefile:
            self.notebook.Destroy()
            self.notebook = None
        else:
            selection = self.notebook.GetSelection()
            if DEBUG: print "Destroy Selection:", selection
            self.notebook.DeletePage(selection)

    def bufferNew(self):
        """Create new buffer."""
        if DEBUG : print "PulsarNotebookFrame","bufferNew"  
        if self._singlefile and self.bufferHasChanged():
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        self.bufferCreate()
        cancel = False
        return cancel

    def bufferOpen(self):
        """Open file in buffer."""
        if DEBUG : print "PulsarNotebookFrame","bufferOpen"  
        if self._singlefile and self.bufferHasChanged():
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        filedir = PULSARPATH+"/workspace"
        if self.buffer and self.buffer.doc.filedir:
            filedir = self.buffer.doc.filedir
        if self._singlefile:
            result = openSingle(directory=filedir)
            if result.path:
                self.bufferCreate(result.path)
        else:
            result = openMultiple(directory=filedir)
            for path in result.paths:
                self.bufferCreate(path)
        cancel = False
        return cancel

    def simulationCompute(self):
	"""
	Run the simulation and show the resulting spectrum
	"""
	if DEBUG : print "PulsarNotebookFrame","simulationCompute" 
        nb=self.notebook.GetCurrentPage()
        nb.SetSelection(1)
        self.editor.log.ClearAll()
        # Make this log active
        wx.Log_SetActiveTarget(MyLog(self.editor.log, 0))
        
        #Run the simulation
        #------------------
	if self.buffer.simulationCompute():
	    self.SetStatusText('Computing done!')
	    nb.SetSelection(1)

            #Show the resulting spectra
            #----------------------
	    try:
                if  self.editor.plot.plot_data(mode="reset") :
                    nb.SetSelection(2)
                    self.editor.plot.Show()
                    self.SetStatusText('Plotting done!')
            except:
                #plot_data returns False
                wx.LogMessage("There was an error when trying of plotting the spectra")

	else:
            #simulationCompute returns False
	    self.SetStatusText('Error: unable to perform computing')


#============================================================================
class EditorShellNotebook(wx.Notebook):
    """A notebook containing an editor page and a shell page."""

    #-----------------
    def __init__(self, parent, filename=None):
        """Create EditorShellNotebook instance."""
        if DEBUG : print "EditorShellNotebook"

        self.parent=parent
        wx.Notebook.__init__(self, parent, id=-1)
        editorparent = editorpanel = wx.Panel(self, -1)
        logparent = logpanel = wx.Panel(self, -1)
        plotparent = plotpanel = wx.Panel(self, -1)
        self.buffer = Buffer()
        
        self.editor = Editor(parent=editorparent)
        self.buffer.addEditor(self.editor)
        self.buffer.open(filename)
        self.log= ViewLog(parent=logparent)
        self.editor.log = self.log
        
        self.plot = PlotFigure(parent=plotparent,filename=filename)
        self.editor.plot = self.plot
        
        # Set the log target to be this textctrl using our own class
        wx.Log_SetActiveTarget(MyLog(self.log, 0))
        
        # for serious debugging
        if SERIOUSDEBUG:
            wx.Log_SetActiveTarget(wx.LogStderr())
            wx.Log_SetTraceMask(wx.TraceMessages)

        self.AddPage(page=editorpanel, text='Editor', select=True)
        self.AddPage(page=logpanel, text='Log')
        self.AddPage(page=plotpanel, text='Plot')

        # Setup sizers

        editorsizer = wx.BoxSizer(wx.VERTICAL)
        editorsizer.Add(self.editor.window, 1, wx.EXPAND)
        editorpanel.SetSizer(editorsizer)
        editorpanel.SetAutoLayout(True)

        logsizer = wx.BoxSizer(wx.VERTICAL)
        logsizer.Add(self.log.window, 1, wx.EXPAND)
        logpanel.SetSizer(logsizer)
        logpanel.SetAutoLayout(True)

        plotsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        plotsizer1.Add(self.plot.canvas, 1, wx.EXPAND)
        plotsizer1.Add(self.plot.choice, 0, wx.EXPAND)
        plotsizer = wx.BoxSizer(wx.VERTICAL)
        plotsizer.Add(plotsizer1, 1, wx.EXPAND)
        plotsizer.Add(self.plot.toolbar, 0, wx.GROW)
        plotpanel.SetSizer(plotsizer)

        self.plot.Fit()

        # Plot the spectra read from filename
        if  self.plot.plot_data() :
            self.SetSelection(2)
            self.plot.Show()
        
        wx.EVT_NOTEBOOK_PAGE_CHANGED(self, self.GetId(), self.OnPageChanged)
        
    #-----------------
    def OnPageChanged(self, event):
        """Page changed event handler."""
        if DEBUG : print "EditorShellNotebook","OnPageChanged"  
        selection = event.GetSelection()
        if selection == 0:
            self.editor.setFocus()
##        if selection == 1:
##            self.log.setFocus()
##        if selection == 2:
##            self.plot.setFocus()
        event.Skip()

    #-----------------
    def SetFocus(self):
        if DEBUG : print "EditorShellNotebook","SetFocus"  
        wx.Notebook.SetFocus(self)
        selection = self.GetSelection()
        if selection == 0:
            self.editor.setFocus()
##        if selection == 1:
##            self.log.setFocus()
##        if selection == 2:
##            self.plot.setFocus()


#============================================================================
class ViewLog:
#============================================================================
    """ A class to handle log from Pulsar """

    #----------------------
    def __init__(self, parent, id=-1, size=wx.DefaultSize):
        """Create ViewLog instance."""
        if DEBUG : print "Viewlog"

        self.window = wx.TextCtrl(parent, id, '',style=wx.TE_MULTILINE|wx.TE_RICH)
        self.window.SetInsertionPoint(0)
        points = self.window.GetFont().GetPointSize()
        f = wx.Font(points, wx.MODERN, wx.NORMAL, wx.NORMAL)
        self.window.SetDefaultStyle(wx.TextAttr("BLACK", wx.NullColour, f))
        # self.window.SetStyle(1, 1000, wx.TextAttr("BLUE", wx.NullColour, f))

        wx.EVT_CHAR(self.window, self.EvtChar)
        wx.EVT_WINDOW_DESTROY(self.window, self.OnWindowDestroy)

    def OnSetFocus(self, evt):
        if DEBUG : print "Viewlog","onSetFocus"  
        evt.Skip()

    #------------------
    def setFocus(self):
        """Set the input focus to the editor window."""
        if DEBUG : print "Viewlog","setFocus"
        
    def OnKillFocus(self, evt):
       if DEBUG : print "Viewlog","OnKillFocus"  
       evt.Skip()

    def OnWindowDestroy(self, evt):
# TODO: May be a saving of the log file would be great.
        if DEBUG : print "Viewlog","OnWindowDestroy"  
        evt.Skip()

    def EvtChar(self, event):
        if DEBUG : print "Viewlog","EvtChar"  
        event.Skip()

    def AppendText(self,text):
        if DEBUG : print "Viewlog","AppendText"  
        self.window.AppendText(text)
        self.window.SetInsertionPoint(self.window.GetLastPosition())
        
    #-----------------
    def ClearAll(self):
        if DEBUG : print "Viewlog","ClearAll"  
        self.window.Clear()

#===============================================================================
class MyLog(wx.PyLog):
    """
     A custom PyLog class                                                          
    """
    
    def __init__(self, textCtrl, logTime=0):
        if DEBUG : print "Mylog"

        wx.PyLog.__init__(self)
        self.tc = textCtrl
        self.logTime = logTime

    def DoLogString(self, message, timeStamp):
        if DEBUG : print "Mylog","DoLogString"  
        if self.logTime:
            message = time.strftime("%X", time.localtime(timeStamp)) + \
                      ": " + message
        if self.tc:
            self.tc.AppendText(message + '\n')

        
#============================================================================
class Editor:
    """Editor having an EditWindow."""

    #-----------------
    def __init__(self, parent, id=-1, pos=wx.DefaultPosition,
                 size=wx.DefaultSize,
                 style=wx.CLIP_CHILDREN | wx.SUNKEN_BORDER):
        """Create Editor instance."""
        if DEBUG : print "Editor"

        self.window = EditWindow(self, parent, id, pos, size, style)
        self.id = self.window.GetId()
        self.buffer = None
        # Assign handlers for keyboard events.
        wx.EVT_CHAR(self.window, self.OnChar)
        wx.EVT_KEY_DOWN(self.window, self.OnKeyDown)
         
    #-----------------
    def _setBuffer(self, buffer, text):
        """Set the editor to a buffer.  Private callback called by buffer."""
        if DEBUG : print "Editor","_setBuffer"  
        self.buffer = buffer
        self.clearAll()
        self.setText(text)
        self.emptyUndoBuffer()
        self.setSavePoint()

    #-----------------
    def destroy(self): 
        """Destroy all editor objects."""
        if DEBUG : print "Editor","destroy" 
        self.window.Destroy()

    #-----------------
    def clearAll(self):
        if DEBUG : print "Editor","clearAll"  
        self.window.ClearAll()

    #-----------------
    def emptyUndoBuffer(self):
        if DEBUG : print "Editor","emptyUndoBuffer"  
        self.window.EmptyUndoBuffer()

    #-----------------
    def getStatus(self):
        """Return (filepath, line, column) status tuple."""
        if MOREDEBUG : print "Editor","getStatus","###"  
        if self.window:
            pos = self.window.GetCurrentPos()
            line = self.window.LineFromPosition(pos) + 1
            col = self.window.GetColumn(pos)
            if self.buffer:
                name = self.buffer.doc.filepath or self.buffer.name
            else:
                name = ''
            status = (name, line, col)
            return status
        else:
            return ('', 0, 0)

    #-----------------
    def getText(self):
        """Return contents of editor."""
        if DEBUG : print "Editor","getText"  
        return self.window.GetText()

    #--------------------
    def hasChanged(self):
        """Return True if contents have changed."""
        if MOREDEBUG : print "Editor","hasChanged","###"  
        return self.window.GetModify()

    #------------------
    def setFocus(self):
        """Set the input focus to the editor window."""
        if DEBUG : print "Editor","setFocus"  
        self.window.SetFocus()

    #----------------------
    def setSavePoint(self):
        if DEBUG : print "Editor","setSavePoint"  
        self.window.SetSavePoint()

    #-----------------------
    def setText(self, text):
        """Set contents of editor."""
        if DEBUG : print "Editor","setText"  
        self.window.SetText(text)

    #-----------------------
    def OnChar(self, event):
        """Keypress event handler.
        
        Only receives an event if OnKeyDown calls event.Skip() for the
        corresponding event."""
        if DEBUG : print "Editor","OnChar" 
        event.Skip()
            
    #--------------------------
    def OnKeyDown(self, event):
        """Key down event handler."""
        if DEBUG : print "Editor","OnKeyDown" 
        key = event.KeyCode()
        controlDown = event.ControlDown()
        altDown = event.AltDown()
        shiftDown = event.ShiftDown()
        # Let Ctrl-Alt-* get handled normally.
        if controlDown and altDown:
            event.Skip()
        # Increase font size.
        elif controlDown and key in (ord(']'),):
            dispatcher.send(signal='FontIncrease')
        # Decrease font size.
        elif controlDown and key in (ord('['),):
            dispatcher.send(signal='FontDecrease')
        # Default font size.
        elif controlDown and key in (ord('='),):
            dispatcher.send(signal='FontDefault')
        else:
            event.Skip()

#============================================================================
class EditWindow(editwindow.EditWindow):
    """EditWindow based on StyledTextCtrl."""

    #----------------------------------------------------------------
    def __init__(self, editor, parent, id=-1, pos=wx.DefaultPosition,
                 size=wx.DefaultSize,
                 style=wx.CLIP_CHILDREN | wx.SUNKEN_BORDER):
        """Create EditWindow instance."""
        if DEBUG : print "EditWindow"
        editwindow.EditWindow.__init__(self, parent, id, pos, size, style)
        self.editor = editor


#============================================================================
class DialogResults:
    """DialogResults class."""

    #----------------------------
    def __init__(self, returned):
        """Create wrapper for results returned by dialog."""
        if DEBUG : print "DialogResults"
        self.returned = returned
        self.positive = returned in (wx.ID_OK, wx.ID_YES)
        self.text = self._asString()
        
    #------------------
    def __repr__(self):
        if DEBUG : print "DialogResults","__repr__" 
        return str(self.__dict__)

    #-------------------
    def _asString(self):
        if DEBUG : print "DialogResults","_asString" 
        returned = self.returned
        if returned == wx.ID_OK:
            return "Ok"
        elif returned == wx.ID_CANCEL:
            return "Cancel"
        elif returned == wx.ID_YES:
            return "Yes"
        elif returned == wx.ID_NO:
            return "No"


#----------------------------------------------------------------------------
def fileDialog(parent=None, title='Open', directory=PULSARPATH, filename='',
               wildcard='All Files (*.*)|*.*',
               style=wx.OPEN | wx.MULTIPLE):
    """File dialog wrapper function."""
    if DEBUG : print "fileDialog ",directory 
    dialog = wx.FileDialog(parent, title, directory, filename,
                           wildcard, style)
    result = DialogResults(dialog.ShowModal())
    if result.positive:
        result.paths = dialog.GetPaths()
    else:
        result.paths = []
    dialog.Destroy()
    if DEBUG: print "file dialog result ", result
    return result


#----------------------------------------------------------------------------
def openSingle(parent=None, title='Open', directory=PULSARPATH, filename='',
               wildcard='PULSAR files (*.pul)|*.pul', style=wx.OPEN):
    """File dialog wrapper function."""
    if DEBUG : print "openSingle" 
    dialog = wx.FileDialog(parent, title, directory, filename,
                           wildcard, style)
    result = DialogResults(dialog.ShowModal())
    if result.positive:
        result.path = dialog.GetPath()
    else:
        result.path = None
    dialog.Destroy()
    return result

#----------------------------------------------------------------------------
def openMultiple(parent=None, title='Open', directory=PULSARPATH, filename='',
                 wildcard='PULSAR files (*.pul)|*.pul',
                 style=wx.OPEN | wx.MULTIPLE):
    """File dialog wrapper function."""
    if DEBUG : "openMultiple" 
    return fileDialog(parent, title, directory, filename, wildcard, style)


#----------------------------------------------------------------------------
def saveSingle(parent=None, title='Save', directory=PULSARPATH, filename='',
               wildcard='PULSAR files (*.pul)|*.pul',
               style=wx.SAVE | wx.HIDE_READONLY | wx.OVERWRITE_PROMPT):
    """File dialog wrapper function."""
    if DEBUG : print "saveSingle" 
    dialog = wx.FileDialog(parent, title, directory, filename,
                           wildcard, style)
    result = DialogResults(dialog.ShowModal())
    if result.positive:
        result.path = dialog.GetPath()
    else:
        result.path = None
    dialog.Destroy()
    return result


#----------------------------------------------------------------------------
def directory(parent=None, message='Choose a directory', path=PULSARPATH, style=0,
              pos=wx.DefaultPosition, size=wx.DefaultSize):
    """Dir dialog wrapper function."""
    if DEBUG : print "directory" 
    dialog = wx.DirDialog(parent, message, path, style, pos, size)
    result = DialogResults(dialog.ShowModal())
    if result.positive:
        result.path = dialog.GetPath()
    else:
        result.path = None
    dialog.Destroy()
    return result


#----------------------------------------------------------------------------
def messageDialog(parent=None, message='', title='Message box',
                  style=wx.YES_NO | wx.CANCEL | wx.CENTRE | wx.ICON_QUESTION,
                  pos=wx.DefaultPosition):
    """Message dialog wrapper function."""
    if DEBUG : print "messageDialog" 
    dialog = wx.MessageDialog(parent, message, title, style, pos)
    result = DialogResults(dialog.ShowModal())
    dialog.Destroy()
    return result

#----------------------------------------------------------------------

class PulsarAbout(wx.Dialog):
    """ An about box that uses an HTML window """

    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, -1, 'About pyPulsar',
                          size=(600, 400) )
        
        html = wx.html.HtmlWindow(self, -1)
        html.LoadPage("README.html")
        button = wx.Button(self, wx.ID_OK, "Ok")
       
        # constraints for the html window
        lc = wx.LayoutConstraints()
        lc.top.SameAs(self, wx.Top, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.bottom.SameAs(button, wx.Top, 5)
        lc.right.SameAs(self, wx.Right, 5)
        html.SetConstraints(lc)

        # constraints for the button
        lc = wx.LayoutConstraints()
        lc.bottom.SameAs(self, wx.Bottom, 5)
        lc.centreX.SameAs(self, wx.CentreX)
        lc.width.AsIs()
        lc.height.AsIs()
        button.SetConstraints(lc)

        self.SetAutoLayout(True)
        self.Layout()
        self.CentreOnParent(wx.BOTH)

class PulsarDoc(wx.Dialog):
    """ An about box that uses an IEHtml window """

    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, -1, 'Loading ONLINE documentation...be patient',
                          size=(999, 640) )
        
        self.html = iewin.IEHtmlWindow(self, -1) #, style = wx.NO_FULL_REPAINT_ON_RESIZE)
        self.html.LoadUrl("http://www-lcs.ensicaen.fr/pyPulsar/index.php/Documentation")
        
        self.Bind(iewin.EVT_DocumentComplete, self.OnDocumentComplete, self.html)
        
        # constraints for the html window
        lc = wx.LayoutConstraints()
        lc.top.SameAs(self, wx.Top, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.bottom.SameAs(self, wx.Bottom, 5)
        lc.right.SameAs(self, wx.Right, 5)
        self.html.SetConstraints(lc)

        self.SetAutoLayout(True)
        self.Layout()
        self.CentreOnParent(wx.BOTH)

    def OnDocumentComplete(self, evt):
        if evt.URL!="about:blank":
            self.SetTitle( "ONLINE Documentation:"+evt.URL)
            
