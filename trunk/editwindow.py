"""
$editwindow.py

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

 Edit Windows

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""

import wx
from wx import stc

import keyword
import os
import sys
import time

import dispatcher  as dispatcher


if 'wxMSW' in wx.PlatformInfo:
    FACES = { 'times'     : 'Times New Roman',
              'mono'      : 'Courier New',
              'helv'      : 'Arial',
              'lucida'    : 'Lucida Console',
              'other'     : 'Comic Sans MS',
              'size'      : 10,
              'lnsize'    : 8,
              'backcol'   : '#FFFFFF',
              'calltipbg' : '#FFFFB8',
              'calltipfg' : '#404040',
            }

elif 'wxGTK' in wx.PlatformInfo and 'gtk2' in wx.PlatformInfo:
    FACES = { 'times'     : 'Serif',
              'mono'      : 'Monospace',
              'helv'      : 'Sans',
              'other'     : 'new century schoolbook',
              'size'      : 10,
              'lnsize'    : 9,
              'backcol'   : '#FFFFFF',
              'calltipbg' : '#FFFFB8',
              'calltipfg' : '#404040',
            }

elif 'wxMac' in wx.PlatformInfo:
    FACES = { 'times'     : 'Lucida Grande',
              'mono'      : 'Courier New',
              'helv'      : 'Geneva',
              'other'     : 'new century schoolbook',
              'size'      : 13,
              'lnsize'    : 10,
              'backcol'   : '#FFFFFF',
              'calltipbg' : '#FFFFB8',
              'calltipfg' : '#404040',
            }

else: # GTK1, etc.
    FACES = { 'times'     : 'Times',
              'mono'      : 'Courier',
              'helv'      : 'Helvetica',
              'other'     : 'new century schoolbook',
              'size'      : 12,
              'lnsize'    : 10,
              'backcol'   : '#FFFFFF',
              'calltipbg' : '#FFFFB8',
              'calltipfg' : '#404040',
            }


class EditWindow(stc.StyledTextCtrl):
    """EditWindow based on StyledTextCtrl."""

    def __init__(self, parent, id=-1, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.CLIP_CHILDREN | wx.SUNKEN_BORDER):
        """Create EditWindow instance."""
        stc.StyledTextCtrl.__init__(self, parent, id, pos, size, style)
        self.__config()
        stc.EVT_STC_UPDATEUI(self, id, self.OnUpdateUI)
        dispatcher.connect(receiver=self._fontsizer, signal='FontIncrease')
        dispatcher.connect(receiver=self._fontsizer, signal='FontDecrease')
        dispatcher.connect(receiver=self._fontsizer, signal='FontDefault')

    def _fontsizer(self, signal):
        """Receiver for Font* signals."""
        size = self.GetZoom()
        if signal == 'FontIncrease':
            size += 1
        elif signal == 'FontDecrease':
            size -= 1
        elif signal == 'FontDefault':
            size = 0
        self.SetZoom(size)

    def __config(self):
        """Configure shell based on user preferences."""
        self.SetMarginType(1, stc.STC_MARGIN_NUMBER)
        self.SetMarginWidth(1, 40)

        self.SetLexer(stc.STC_LEX_PYTHON)
        self.SetKeyWords(0, ' '.join(keyword.kwlist))

##CF
        # Enable folding
        self.SetProperty("fold", "1" ) 

        # Highlight tab/space mixing (shouldn't be any)
        #self.SetProperty("tab.timmy.whinge.level", "1")

        # Set left and right margins
        self.SetMargins(2,2)

        # Set up the numbers in the margin for margin #1
        self.SetMarginType(1, wx.stc.STC_MARGIN_NUMBER)
        # Reasonable value for, say, 4-5 digits using a mono font (40 pix)
        self.SetMarginWidth(1, 40)

        # Indentation and tab stuff
        self.SetIndent(4)               # Proscribed indent size for wx
        self.SetIndentationGuides(True) # Show indent guides
        self.SetBackSpaceUnIndents(True)# Backspace unindents rather than delete 1 space
        self.SetTabIndents(True)        # Tab key indents
        self.SetTabWidth(4)             # Proscribed tab size for wx
        self.SetUseTabs(False)          # Use spaces rather than tabs, or
                                        # TabTimmy will complain!    
        # White space
        self.SetViewWhiteSpace(False)   # Don't view white space

        # EOL: Since we are loading/saving ourselves, and the
        # strings will always have \n's in them, set the STC to
        # edit them that way.            
        self.SetEOLMode(wx.stc.STC_EOL_LF)
        self.SetViewEOL(False)
        
	# No right-edge mode indicator
        self.SetEdgeMode(stc.STC_EDGE_NONE)

        # Setup a margin to hold fold markers
        self.SetMarginType(2, stc.STC_MARGIN_SYMBOL)
        self.SetMarginMask(2, stc.STC_MASK_FOLDERS)
        self.SetMarginSensitive(2, True)
        self.SetMarginWidth(2, 12)

        # and now set up the fold markers
        self.MarkerDefine(stc.STC_MARKNUM_FOLDEREND,     stc.STC_MARK_BOXPLUSCONNECTED,  "white", "black")
        self.MarkerDefine(stc.STC_MARKNUM_FOLDEROPENMID, stc.STC_MARK_BOXMINUSCONNECTED, "white", "black")
        self.MarkerDefine(stc.STC_MARKNUM_FOLDERMIDTAIL, stc.STC_MARK_TCORNER,  "white", "black")
        self.MarkerDefine(stc.STC_MARKNUM_FOLDERTAIL,    stc.STC_MARK_LCORNER,  "white", "black")
        self.MarkerDefine(stc.STC_MARKNUM_FOLDERSUB,     stc.STC_MARK_VLINE,    "white", "black")
        self.MarkerDefine(stc.STC_MARKNUM_FOLDER,        stc.STC_MARK_BOXPLUS,  "white", "black")
        self.MarkerDefine(stc.STC_MARKNUM_FOLDEROPEN,    stc.STC_MARK_BOXMINUS, "white", "black")

##/CF
        self.setStyles(FACES)
        self.SetViewWhiteSpace(False)
        self.SetTabWidth(4)
        self.SetUseTabs(False)
        self.SetWrapMode(False)
        try:
            self.SetEndAtLastLine(False)
        except AttributeError:
            pass

    def setStyles(self, faces):
        """Configure font size, typeface and color for lexer."""

        # Default style
        self.StyleSetSpec(stc.STC_STYLE_DEFAULT,
                          "face:%(mono)s,size:%(size)d,back:%(backcol)s" % \
                          faces)

        self.StyleClearAll()

        # Built in styles
##        self.StyleSetSpec(stc.STC_STYLE_LINENUMBER,
##                          "back:#C0C0C0,face:%(mono)s,size:%(lnsize)d" % faces)
##        self.StyleSetSpec(stc.STC_STYLE_CONTROLCHAR,
##                          "face:%(mono)s" % faces)
##        self.StyleSetSpec(stc.STC_STYLE_BRACELIGHT,
##                          "fore:#0000FF,back:#FFFF88")
##        self.StyleSetSpec(stc.STC_STYLE_BRACEBAD,
##                          "fore:#FF0000,back:#FFFF88")

        # Python styles
##        self.StyleSetSpec(stc.STC_P_DEFAULT,
##                          "face:%(mono)s" % faces)
##        self.StyleSetSpec(stc.STC_P_COMMENTLINE,
##                          "fore:#007F00,face:%(mono)s" % faces)
##        self.StyleSetSpec(stc.STC_P_NUMBER,
##                          "")
##        self.StyleSetSpec(stc.STC_P_STRING,
##                          "fore:#7F007F,face:%(mono)s" % faces)
##        self.StyleSetSpec(stc.STC_P_CHARACTER,
##                          "fore:#7F007F,face:%(mono)s" % faces)
##        self.StyleSetSpec(stc.STC_P_WORD,
##                          "fore:#00007F,bold")
##        self.StyleSetSpec(stc.STC_P_TRIPLE,
##                          "fore:#7F0000")
##        self.StyleSetSpec(stc.STC_P_TRIPLEDOUBLE,
##                          "fore:#000033,back:#FFFFE8")
##        self.StyleSetSpec(stc.STC_P_CLASSNAME,
##                          "fore:#0000FF,bold")
##        self.StyleSetSpec(stc.STC_P_DEFNAME,
##                          "fore:#007F7F,bold")
##        self.StyleSetSpec(stc.STC_P_OPERATOR,
##                          "")
##        self.StyleSetSpec(stc.STC_P_IDENTIFIER,
##                          "")
##        self.StyleSetSpec(stc.STC_P_COMMENTBLOCK,
##                          "fore:#7F7F7F")
##        self.StyleSetSpec(stc.STC_P_STRINGEOL,
##                          "fore:#000000,face:%(mono)s,back:#E0C0E0,eolfilled" % faces)

       # Caret color
        self.SetCaretForeground("BLUE")
        # Selection background
        self.SetSelBackground(1, '#66CCFF')

        self.SetSelBackground(True, wx.SystemSettings_GetColour(wx.SYS_COLOUR_HIGHLIGHT))
        self.SetSelForeground(True, wx.SystemSettings_GetColour(wx.SYS_COLOUR_HIGHLIGHTTEXT))

                # Line numbers in margin
        self.StyleSetSpec(wx.stc.STC_STYLE_LINENUMBER,'fore:#000000,back:#99A9C2')    
        # Highlighted brace
        self.StyleSetSpec(wx.stc.STC_STYLE_BRACELIGHT,'fore:#00009D,back:#FFFF00')
        # Unmatched brace
        self.StyleSetSpec(wx.stc.STC_STYLE_BRACEBAD,'fore:#00009D,back:#FF0000')
        # Indentation guide
        self.StyleSetSpec(wx.stc.STC_STYLE_INDENTGUIDE, "fore:#CDCDCD")

        # Python styles
        self.StyleSetSpec(wx.stc.STC_P_DEFAULT, 'fore:#000000')
        # Comments
        self.StyleSetSpec(wx.stc.STC_P_COMMENTLINE,  'fore:#008000,back:#F0FFF0')
        self.StyleSetSpec(wx.stc.STC_P_COMMENTBLOCK, 'fore:#008000,back:#F0FFF0')
        # Numbers
        self.StyleSetSpec(wx.stc.STC_P_NUMBER, 'fore:#008080')
        # Strings and characters
        self.StyleSetSpec(wx.stc.STC_P_STRING, 'fore:#800080')
        self.StyleSetSpec(wx.stc.STC_P_CHARACTER, 'fore:#800080')
        # Keywords
        self.StyleSetSpec(wx.stc.STC_P_WORD, 'fore:#000080,bold')
        # Triple quotes
        self.StyleSetSpec(wx.stc.STC_P_TRIPLE, 'fore:#800080,back:#FFFFEA')
        self.StyleSetSpec(wx.stc.STC_P_TRIPLEDOUBLE, 'fore:#800080,back:#FFFFEA')
        # Class names
        self.StyleSetSpec(wx.stc.STC_P_CLASSNAME, 'fore:#0000FF,bold')
        # Function names
        self.StyleSetSpec(wx.stc.STC_P_DEFNAME, 'fore:#008080,bold')
        # Operators
        self.StyleSetSpec(wx.stc.STC_P_OPERATOR, 'fore:#800000,bold')
        # Identifiers. I leave this as not bold because everything seems
        # to be an identifier if it doesn't match the above criterae
        self.StyleSetSpec(wx.stc.STC_P_IDENTIFIER, 'fore:#000000')



    def OnUpdateUI(self, event):
        """Check for matching braces."""
        # If the auto-complete window is up let it do its thing.
        if self.AutoCompActive() or self.CallTipActive():
            return
        braceAtCaret = -1
        braceOpposite = -1
        charBefore = None
        caretPos = self.GetCurrentPos()
        if caretPos > 0:
            charBefore = self.GetCharAt(caretPos - 1)
            styleBefore = self.GetStyleAt(caretPos - 1)

        # Check before.
        if charBefore and chr(charBefore) in '[]{}()' \
        and styleBefore == stc.STC_P_OPERATOR:
            braceAtCaret = caretPos - 1

        # Check after.
        if braceAtCaret < 0:
            charAfter = self.GetCharAt(caretPos)
            styleAfter = self.GetStyleAt(caretPos)
            if charAfter and chr(charAfter) in '[]{}()' \
            and styleAfter == stc.STC_P_OPERATOR:
                braceAtCaret = caretPos

        if braceAtCaret >= 0:
            braceOpposite = self.BraceMatch(braceAtCaret)

        if braceAtCaret != -1  and braceOpposite == -1:
            self.BraceBadLight(braceAtCaret)
        else:
            self.BraceHighlight(braceAtCaret, braceOpposite)

    def CanCopy(self):
        """Return True if text is selected and can be copied."""
        return self.GetSelectionStart() != self.GetSelectionEnd()

    def CanCut(self):
        """Return True if text is selected and can be cut."""
        return self.CanCopy() and self.CanEdit()

    def CanEdit(self):
        """Return True if editing should succeed."""
        return not self.GetReadOnly()

    def CanPaste(self):
        """Return True if pasting should succeed."""
        return stc.StyledTextCtrl.CanPaste(self) and self.CanEdit()
