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

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""
__revision__ = "1"
__revision_date__="2007.02.02" 



import wx
import sys, os
import frame as frame
import version
import editwindow as editwindow
import dispatcher as dispatcher
from buffer import Buffer
from plot import *


DEBUG = False
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
        if DEBUG : print self.__class__

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
         
        pass

    def setEditor(self, editor):
         
        self.editor = editor
        self.buffer = self.editor.buffer
        self.buffers[self.buffer.id] = self.buffer

    def OnAbout(self, event):
        """Display an About window."""
         
        title = 'About pyPulsar'
        text = title='pyPULSAR '+version.__PROJECT_VERSION__+" rev."+version.__PROJECT_REVISION__+'\n\n'+ \
             'A solid State NMR Simulation program\n'+ \
             'by C.Fernandez and J.P. Amoureux.\n\n\n'+ \
             '(www-lcs.ensicaen.fr)\n\n'+ \
             '(For more information: Look at the documentation) '
        
        dialog = wx.MessageDialog(self, text, title,
                                  wx.OK | wx.ICON_INFORMATION)
        dialog.ShowModal()
        dialog.Destroy()

    def OnClose(self, event):
        """Event handler for closing."""
         
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
         
        self._updateStatus()
        if hasattr(self, 'notebook'):
            self._updateTabText()
        self._updateTitle()
        event.Skip()

    def _updateStatus(self):
        """Show current status information."""
         
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
        suffix = ' *'
        size=2
        notebook = self.notebook
        selection = notebook.GetSelection()
        if selection == -1:
            return
        text = notebook.GetPageText(selection)

##CF ADDED Modified argument to handle SaveAs change of file
        if modified is not None:
            modified=os.path.split(modified)[1]
            notebook.SetPageText(selection, modified)
            
        window = notebook.GetPage(selection)
##        if not (text.endswith("  ") or text.endswith(suffix)):
##                 notebook.SetPageText(selection, text + "  ")
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
         
        if self.buffer:
            return True
        else:
            return False

    def bufferClose(self):
        """Close buffer."""
         
        if self.bufferHasChanged():
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        self.bufferDestroy()
        cancel = False
        return cancel

    def bufferCreate(self, filename=None):
        """Create new buffer."""
         
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
         
        if self.buffer:
            for editor in self.buffer.editors.values():
                editor.destroy()
            self.editor = None
            del self.buffers[self.buffer.id]
            self.buffer = None
            self.panel.Destroy()


    def bufferHasChanged(self):
        """Return True if buffer has changed since last save."""
         
        if self.buffer:
            return self.buffer.hasChanged()
        else:
            return False

    def bufferNew(self):
        """Create new buffer."""
         
        if self.bufferHasChanged():
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        self.bufferCreate()
        cancel = False
        return cancel

    def bufferOpen(self):
        """Open file in buffer."""
        print "editor frame buffer open" 
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

##     def bufferPrint(self):
##         """Print buffer."""
##         pass

##     def bufferRevert(self):
##         """Revert buffer to version of file on disk."""
##         pass

    def bufferSave(self):
        """Save buffer to its file."""
         
        if self.buffer.doc.filepath:
            self.buffer.save()
            cancel = False
        else:
            cancel = self.bufferSaveAs()
        return cancel

    def bufferSaveAs(self):
        """Save buffer to a new filename."""
         
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
            #print filedir, self.buffer.doc.filedir,self.buffer.doc.filepath
            self._updateTabText(modified=self.buffer.doc.filepath)
            cancel = False
        else:
            cancel = True
        return cancel

    def bufferSuggestSave(self):
        """Suggest saving changes.  Return True if user selected Cancel."""
         
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
        if DEBUG : print self.__class__

        self.notebook = None
        EditorFrame.__init__(self, parent, id, title, pos,
                             size, style, filename)
        if self.notebook:
            dispatcher.connect(receiver=self._editorChange,
                               signal='EditorChange', sender=self.notebook)
        

    def _editorChange(self, editor):
        """Editor change signal receiver."""
         
        self.setEditor(editor)

##    def OnAbout(self, event):
##        """Display an About window."""
##         
##        title = 'About PyAlaMode'
##        text = 'Another fine, flaky program.'
##        dialog = wx.MessageDialog(self, text, title,
##                                  wx.OK | wx.ICON_INFORMATION)
##        dialog.ShowModal()
##        dialog.Destroy()

    def _updateTitle(self):
##        """Show current title information."""
##         
##        pass
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
         
        buffer = Buffer()
        panel = wx.Panel(parent=self.notebook, id=-1)
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
        self.notebook.AddPage(page=panel, text=self.buffer.name, select=True)
        self.editor.setFocus()

    def bufferDestroy(self):
        """Destroy the current buffer."""
         
        selection = self.notebook.GetSelection()
##         print "Destroy Selection:", selection
        if selection > 0:  # Don't destroy the PyCrust tab.
            if self.buffer:
                del self.buffers[self.buffer.id]
                self.buffer = None  # Do this before DeletePage().
            self.notebook.DeletePage(selection)

    def bufferNew(self):
        """Create new buffer."""
         
        self.bufferCreate()
        cancel = False
        return cancel

    def bufferOpen(self):
        """Open file in buffer."""
         
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
        if DEBUG : print self.__class__
        wx.Notebook.__init__(self, parent, id=-1, style=wx.CLIP_CHILDREN)
        wx.EVT_NOTEBOOK_PAGE_CHANGING(self, self.GetId(),
                                      self.OnPageChanging)
        wx.EVT_NOTEBOOK_PAGE_CHANGED(self, self.GetId(),
                                     self.OnPageChanged)
        wx.EVT_IDLE(self, self.OnIdle)
        
    def OnIdle(self, event):
        """Event handler for idle time."""
         
        self._updateTabText()
        event.Skip()

    def _updateTabText(self):
        """Show current buffer display name on all but first tab."""
         
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
         
        event.Skip()

    def OnPageChanged(self, event):
        """Page changed event handler."""
         
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
        
        """Create EditorShellNotebookFrame instance."""
        if DEBUG : print self.__class__
        self._singlefile = singlefile
        
        EditorNotebookFrame.__init__(self, parent, id, title, pos,
                                     size, style, filename)
       
    def _setup(self):
        """Setup prior to first buffer creation.

        Called automatically by base class during init."""
         
        if not self._singlefile:
            self.notebook = EditorNotebook(parent=self)

    def OnAbout(self, event):
        """Display an About window."""
         
        title = 'About pyPulsar'
        text ='pyPULSAR '+version.__PROJECT_VERSION__+" rev."+version.__PROJECT_REVISION__+'\n\n'+ \
             'A solid state NMR Simulation program\n'+ \
             'by C.Fernandez and J.P. Amoureux.\n\n\n'+ \
             '(www-lcs.ensicaen.fr)\n\n'+ \
             '(For more information: Look at the documentation) '
        dialog = wx.MessageDialog(self, text, title,
                                  wx.OK | wx.ICON_INFORMATION)
        dialog.ShowModal()
        dialog.Destroy()

    def bufferCreate(self, filename=None):
        """Create new buffer."""

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
         
        if self.buffer:
            self.editor = None
            del self.buffers[self.buffer.id]
            self.buffer = None  # Do this before DeletePage().
        if self._singlefile:
            self.notebook.Destroy()
            self.notebook = None
        else:
            selection = self.notebook.GetSelection()
##             print "Destroy Selection:", selection
            self.notebook.DeletePage(selection)

    def bufferNew(self):
        """Create new buffer."""
         
        if self._singlefile and self.bufferHasChanged():
            cancel = self.bufferSuggestSave()
            if cancel:
                return cancel
        self.bufferCreate()
        cancel = False
        return cancel

    def bufferOpen(self):
        """Open file in buffer."""
         
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

 
    #---------------------------
    def simulationCompute(self):
	"""
	Run the simulation and show the resulting spectrum
	"""
        nb=self.notebook.GetCurrentPage()
        nb.SetSelection(1)
        self.editor.log.ClearAll()

        #Run the simulation
        #------------------
	if self.buffer.simulationCompute():
	    self.SetStatusText('Computing done!')
	    nb.SetSelection(1)

            #Show the resulting spectra
            #----------------------
            if  self.editor.plot.plot_data(mode="reset") :
                nb.SetSelection(2)
                self.editor.plot.Show()
                self.SetStatusText('Plotting done!')
            else:
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
        if DEBUG : print self.__class__

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
        """Create Editor instance."""
        if DEBUG : print self.__class__

        self.window = wx.TextCtrl(parent, id, '',style=wx.TE_MULTILINE)
        self.window.SetInsertionPoint(0)
        wx.EVT_CHAR(self.window, self.EvtChar)
        wx.EVT_WINDOW_DESTROY(self.window, self.OnWindowDestroy)

    def OnSetFocus(self, evt):
         
        evt.Skip()

    #------------------
    def setFocus(self):
        """Set the input focus to the editor window."""
#        self.window.SetFocus()
         
#        wx.LogMessage('focus')

    def OnKillFocus(self, evt):
#        wx.LogMessage("OnKillFocus")
        
       evt.Skip()

    def OnWindowDestroy(self, evt):
#        wx.LogMessage("OnWindowDestroy")
# TODO: MAy be a saving of the log file would be great.
         
        evt.Skip()

    def EvtChar(self, event):
#        wx.LogMessage('EvtChar: %d' % event.GetKeyCode())
         
        event.Skip()

    def AppendText(self,text):
         
        self.window.AppendText(text)

    #-----------------
    def ClearAll(self):
         
        self.window.Clear()

#===============================================================================
class MyLog(wx.PyLog):
    """
     A custom PyLog class                                                          
    """
    
    def __init__(self, textCtrl, logTime=0):
        if DEBUG : print self.__class__

        wx.PyLog.__init__(self)
        self.tc = textCtrl
        self.logTime = logTime

    def DoLogString(self, message, timeStamp):
         
        if self.logTime:
##            message = time.strftime("%X", time.localtime(timeStamp)) + \
##                      ": " + message
            pass
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
        if DEBUG : print self.__class__

        self.window = EditWindow(self, parent, id, pos, size, style)
        self.id = self.window.GetId()
        self.buffer = None
        # Assign handlers for keyboard events.
        wx.EVT_CHAR(self.window, self.OnChar)
        wx.EVT_KEY_DOWN(self.window, self.OnKeyDown)
         
    #-----------------
    def _setBuffer(self, buffer, text):
        """Set the editor to a buffer.  Private callback called by buffer."""
         
        self.buffer = buffer
        self.clearAll()
        self.setText(text)
        self.emptyUndoBuffer()
        self.setSavePoint()

    #-----------------
    def destroy(self):
         
        """Destroy all editor objects."""
        self.window.Destroy()

    #-----------------
    def clearAll(self):
         
        self.window.ClearAll()

    #-----------------
    def emptyUndoBuffer(self):
         
        self.window.EmptyUndoBuffer()

    #-----------------
    def getStatus(self):
        """Return (filepath, line, column) status tuple."""
         
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
         
        return self.window.GetText()

    #--------------------
    def hasChanged(self):
        """Return True if contents have changed."""
         
        return self.window.GetModify()

    #------------------
    def setFocus(self):
        """Set the input focus to the editor window."""
         
        self.window.SetFocus()

    #----------------------
    def setSavePoint(self):
         
        self.window.SetSavePoint()

    #-----------------------
    def setText(self, text):
        """Set contents of editor."""
         
        self.window.SetText(text)

    #-----------------------
    def OnChar(self, event):
        """Keypress event handler.
        
        Only receives an event if OnKeyDown calls event.Skip() for the
        corresponding event."""
        event.Skip()
            
    #--------------------------
    def OnKeyDown(self, event):
        """Key down event handler."""

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
        if DEBUG : print self.__class__
        editwindow.EditWindow.__init__(self, parent, id, pos, size, style)
        self.editor = editor


#============================================================================
class DialogResults:
    """DialogResults class."""

    #----------------------------
    def __init__(self, returned):
        """Create wrapper for results returned by dialog."""
        if DEBUG : print self.__class__
        self.returned = returned
        self.positive = returned in (wx.ID_OK, wx.ID_YES)
        self.text = self._asString()
        
    #------------------
    def __repr__(self):
        return str(self.__dict__)

    #-------------------
    def _asString(self):
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
    dialog = wx.FileDialog(parent, title, directory, filename,
                           wildcard, style)
    result = DialogResults(dialog.ShowModal())
    if result.positive:
        result.paths = dialog.GetPaths()
    else:
        result.paths = []
    dialog.Destroy()
    return result


#----------------------------------------------------------------------------
def openSingle(parent=None, title='Open', directory=PULSARPATH, filename='',
               wildcard='PULSAR files (*.pul)|*.pul', style=wx.OPEN):
    """File dialog wrapper function."""
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
    return fileDialog(parent, title, directory, filename, wildcard, style)


#----------------------------------------------------------------------------
def saveSingle(parent=None, title='Save', directory=PULSARPATH, filename='',
               wildcard='PULSAR files (*.pul)|*.pul',
               style=wx.SAVE | wx.HIDE_READONLY | wx.OVERWRITE_PROMPT):
    """File dialog wrapper function."""
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
    dialog = wx.MessageDialog(parent, message, title, style, pos)
    result = DialogResults(dialog.ShowModal())
    dialog.Destroy()
    return result



