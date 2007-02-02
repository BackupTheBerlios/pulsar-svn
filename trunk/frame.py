"""
$frame.py

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

 Frames and Menus

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""
__revision__ = "1"
__revision_date__="2007.02.02"

import wx
import images as images

ID_NEW = wx.ID_NEW
ID_OPEN = wx.ID_OPEN
ID_REVERT = wx.ID_REVERT
ID_CLOSE = wx.ID_CLOSE
ID_SAVE = wx.ID_SAVE
ID_SAVEAS = wx.ID_SAVEAS
ID_PRINT = wx.ID_PRINT
ID_EXIT = wx.ID_EXIT
ID_UNDO = wx.ID_UNDO
ID_REDO = wx.ID_REDO
ID_CUT = wx.ID_CUT
ID_COPY = wx.ID_COPY
ID_PASTE = wx.ID_PASTE
ID_CLEAR = wx.ID_CLEAR
ID_SELECTALL = wx.ID_SELECTALL
ID_ABOUT = wx.ID_ABOUT
ID_COPY_PLUS = wx.NewId()
ID_COMPUTE = wx.NewId()
ID_PASTE_PLUS = wx.NewId()


class Frame(wx.Frame):
    """Frame with standard menu items."""

    def __init__(self, parent=None, id=-1, title='Editor',
                 pos=wx.DefaultPosition, size=wx.DefaultSize, 
                 style=wx.DEFAULT_FRAME_STYLE):
        """Create a Frame instance."""
	wx.Frame.__init__(self, parent, id, title, pos, size, style)
        self.CreateStatusBar()
        self.SetStatusText('Frame')
        self.SetIcon(images.getPulsarIcon())
        self.__createMenus()
        wx.EVT_CLOSE(self, self.OnClose)
        self.Centre(wx.BOTH)
        
    def OnClose(self, event):
        """Event handler for closing."""
        self.Destroy()

    def __createMenus(self):
        m = self.fileMenu = wx.Menu()
        m.Append(ID_NEW, '&New \tCtrl+N',
                 'New file')
        m.Append(ID_OPEN, '&Open... \tCtrl+O',
                 'Open file')
        m.AppendSeparator()
        m.Append(ID_REVERT, '&Revert \tCtrl+R',
                 'Revert to last saved version')
        m.Append(ID_CLOSE, '&Close \tCtrl+W',
                 'Close file')
        m.AppendSeparator()
        m.Append(ID_SAVE, '&Save... \tCtrl+S',
                 'Save file')
        m.Append(ID_SAVEAS, 'Save &As \tShift+Ctrl+S',
                 'Save file with new name')
        m.AppendSeparator()
        m.Append(ID_PRINT, '&Print... \tCtrl+P',
                 'Print file')
        m.AppendSeparator()
        m.Append(ID_EXIT, 'E&xit', 'Exit Program')

        m = self.editMenu = wx.Menu()
        m.Append(ID_UNDO, '&Undo \tCtrl+Z',
                 'Undo the last action')
        m.Append(ID_REDO, '&Redo \tCtrl+Y',
                 'Redo the last undone action')
        m.AppendSeparator()
        m.Append(ID_CUT, 'Cu&t \tCtrl+X',
                 'Cut the selection')
        m.Append(ID_COPY, '&Copy \tCtrl+C',
                 'Copy the selection')
        m.Append(ID_COPY_PLUS, 'Cop&y Plus \tShift+Ctrl+C',
                 'Copy the selection - retaining prompts')
        m.Append(ID_PASTE, '&Paste \tCtrl+V', 'Paste from clipboard')
        m.Append(ID_PASTE_PLUS, 'Past&e Plus \tShift+Ctrl+V',
                 'Paste and run commands')
        m.AppendSeparator()
        m.Append(ID_CLEAR, 'Cle&ar',
                 'Delete the selection')
        m.Append(ID_SELECTALL, 'Select A&ll \tCtrl+A',
                 'Select all text')
        m = self.computeMenu = wx.Menu()
        m.Append(ID_COMPUTE, '&Compute \tShift+Ctrl+C',
                 'Compute using the current script - If modified, the current script must be saved before the menu becomes active')
         
        m = self.helpMenu = wx.Menu()
        m.AppendSeparator()
        m.Append(ID_ABOUT, '&About...', 'About this program')

        b = self.menuBar = wx.MenuBar()
        b.Append(self.fileMenu, '&File')
        b.Append(self.editMenu, '&Edit')
        b.Append(self.computeMenu, '&Simulation')
        b.Append(self.helpMenu, '&Help')
        self.SetMenuBar(b)

        wx.EVT_MENU(self, ID_NEW, self.OnFileNew)
        wx.EVT_MENU(self, ID_OPEN, self.OnFileOpen)
        wx.EVT_MENU(self, ID_REVERT, self.OnFileRevert)
        wx.EVT_MENU(self, ID_CLOSE, self.OnFileClose)
        wx.EVT_MENU(self, ID_SAVE, self.OnFileSave)
        wx.EVT_MENU(self, ID_SAVEAS, self.OnFileSaveAs)
        wx.EVT_MENU(self, ID_COMPUTE, self.OnSimulationCompute)
        wx.EVT_MENU(self, ID_PRINT, self.OnFilePrint)
        wx.EVT_MENU(self, ID_EXIT, self.OnExit)
        wx.EVT_MENU(self, ID_UNDO, self.OnUndo)
        wx.EVT_MENU(self, ID_REDO, self.OnRedo)
        wx.EVT_MENU(self, ID_CUT, self.OnCut)
        wx.EVT_MENU(self, ID_COPY, self.OnCopy)
        wx.EVT_MENU(self, ID_COPY_PLUS, self.OnCopyPlus)
        wx.EVT_MENU(self, ID_PASTE, self.OnPaste)
        wx.EVT_MENU(self, ID_PASTE_PLUS, self.OnPastePlus)
        wx.EVT_MENU(self, ID_CLEAR, self.OnClear)
        wx.EVT_MENU(self, ID_SELECTALL, self.OnSelectAll)
        wx.EVT_MENU(self, ID_ABOUT, self.OnAbout)
        
        wx.EVT_UPDATE_UI(self, ID_NEW, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_OPEN, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_REVERT, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_CLOSE, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_SAVE, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_SAVEAS, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_COMPUTE, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_PRINT, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_UNDO, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_REDO, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_CUT, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_COPY, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_COPY_PLUS, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_PASTE, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_PASTE_PLUS, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_CLEAR, self.OnUpdateMenu)
        wx.EVT_UPDATE_UI(self, ID_SELECTALL, self.OnUpdateMenu)
 
    def OnFileNew(self, event):
        self.bufferNew()

    def OnFileOpen(self, event):
        self.bufferOpen()

    def OnFileRevert(self, event):
        self.bufferRevert()

    def OnFileClose(self, event):
        self.bufferClose()

    def OnFileSave(self, event):
        self.bufferSave()

    def OnFileSaveAs(self, event):
        self.bufferSaveAs()

    def OnSimulationCompute(self, event):
        #use to run the PULSAR program
        self.simulationCompute()

    def OnFilePrint(self, event):
        self.bufferPrint()

    def OnExit(self, event):
        self.Close(False)

    def OnUndo(self, event):
        win = wx.Window_FindFocus()
        win.Undo()

    def OnRedo(self, event):
        win = wx.Window_FindFocus()
        win.Redo()

    def OnCut(self, event):
        win = wx.Window_FindFocus()
        win.Cut()

    def OnCopy(self, event):
        win = wx.Window_FindFocus()
        win.Copy()

    def OnCopyPlus(self, event):
        win = wx.Window_FindFocus()
        win.CopyWithPrompts()

    def OnPaste(self, event):
        win = wx.Window_FindFocus()
        win.Paste()

    def OnPastePlus(self, event):
        win = wx.Window_FindFocus()
        win.PasteAndRun()

    def OnClear(self, event):
        win = wx.Window_FindFocus()
        win.Clear()

    def OnSelectAll(self, event):
        win = wx.Window_FindFocus()
        win.SelectAll()

    def OnAbout(self, event):
        """Display an About window."""
        title = 'About'
        text = 'Your message here.'
        dialog = wx.MessageDialog(self, text, title,
                                  wx.OK | wx.ICON_INFORMATION)
        dialog.ShowModal()
        dialog.Destroy()

    def OnUpdateMenu(self, event):
        """Update menu items based on current status and context."""
        win = wx.Window_FindFocus()
        id = event.GetId()
        event.Enable(True)
        try:
            if id == ID_NEW:
                event.Enable(hasattr(self, 'bufferNew'))
            elif id == ID_OPEN:
                event.Enable(hasattr(self, 'bufferOpen'))
            elif id == ID_REVERT:
                event.Enable(hasattr(self, 'bufferRevert')
                             and self.hasBuffer())
            elif id == ID_CLOSE:
                event.Enable(hasattr(self, 'bufferClose')
                             and self.hasBuffer())
            elif id == ID_SAVE:
                event.Enable(hasattr(self, 'bufferSave')
                             and self.bufferHasChanged())
            elif id == ID_SAVEAS:
                event.Enable(hasattr(self, 'bufferSaveAs')
                             and self.hasBuffer())
            elif id == ID_COMPUTE:
                event.Enable(hasattr(self, 'simulationCompute')
                             and self.hasBuffer()
                             and not self.bufferHasChanged())
            elif id == ID_PRINT:
                event.Enable(hasattr(self, 'bufferPrint')
                             and self.hasBuffer())
            elif id == ID_UNDO:
                event.Enable(win.CanUndo())
            elif id == ID_REDO:
                event.Enable(win.CanRedo())
            elif id == ID_CUT:
                event.Enable(win.CanCut())
            elif id == ID_COPY:
                event.Enable(win.CanCopy())
            elif id == ID_COPY_PLUS:
                event.Enable(win.CanCopy() and hasattr(win, 'CopyWithPrompts'))
            elif id == ID_PASTE:
                event.Enable(win.CanPaste())
            elif id == ID_PASTE_PLUS:
                event.Enable(win.CanPaste() and hasattr(win, 'PasteAndRun'))
            elif id == ID_CLEAR:
                event.Enable(win.CanCut())
            elif id == ID_SELECTALL:
                event.Enable(hasattr(win, 'SelectAll'))
            else:
                event.Enable(False)
        except AttributeError:
            # This menu option is not supported in the current context.
            event.Enable(False)
