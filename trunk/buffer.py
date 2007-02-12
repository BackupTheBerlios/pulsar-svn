"""
$buffer.py

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

import os
import sys
import document as document
import process
from f95pulsar import reset

DEBUG=False
MOREDEBUG=False

class Buffer:
    """Buffer class."""
    id = 0

    def __init__(self, filename=None):
        """Create a Buffer instance."""
        Buffer.id += 1
        if DEBUG : print "Buffer",Buffer.id
        self.id = Buffer.id
        self.name = ''
        self.editors = {}
        self.editor = None
        self.syspath = sys.path[:]
        self.confirmed= True
        while True:
            try:
                self.syspath.remove('')
            except ValueError:
                break
        while True:
            try:
                self.syspath.remove('.')
            except ValueError:
                break
        self.open(filename)

    def addEditor(self, editor):
        """Add an editor."""
        if DEBUG : print "Buffer", "addEditor"
        self.editor = editor
        self.editors[editor.id] = editor

    def hasChanged(self):
        """Return True if text in editor has changed since last save."""
        if MOREDEBUG : print "Buffer", "hasChanged","###"
        if self.editor:
            return self.editor.hasChanged()
        else:
            return False

    def new(self, filepath):
        """New empty buffer."""
        if DEBUG : print "Buffer", "new"
        if not filepath:
            return
        if os.path.exists(filepath):
            self.confirmed = self.overwriteConfirm(filepath)
        else:
            self.confirmed = True

    def open(self, filename):
        """Open file into buffer."""
        if DEBUG : print "Buffer", "open"
        self.doc = document.Document(filename)
        self.name = self.doc.filename or ('Untitled:' + str(self.id))
        if self.doc.filedir and self.doc.filedir not in self.syspath:
            # To create the proper context for updateNamespace.
            self.syspath.insert(0, self.doc.filedir)
        if self.doc.filepath and os.path.exists(self.doc.filepath):
            self.confirmed = True
        if self.editor:
            text = self.doc.read()
            self.editor._setBuffer(buffer=self, text=text)
        if DEBUG: print "filename: ",filename
        
    def overwriteConfirm(filepath):
        """Confirm overwriting an existing file."""
        if DEBUG : print "Buffer", "overwriteconfirm"
        return False

    def save(self):
        """Save buffer."""
        if DEBUG : print "Buffer", "save"
        filepath = self.doc.filepath
        if not filepath:
            return  # XXX Get filename
        if not os.path.exists(filepath):
            self.confirmed = True
        if not self.confirmed:
            self.confirmed = self.overwriteConfirm(filepath)
        if self.confirmed:
            self.doc.write(self.editor.getText())
            if self.editor:
                self.editor.setSavePoint()

    def saveAs(self, filename):
        """Save buffer."""
        if DEBUG : print "Buffer", "saveAs"        
        self.doc = document.Document(filename)
        self.name = self.doc.filename
        self.save()

    def simulationCompute(self):
        """Run Simulation.
        Return True if updated, False if there was an error."""
        if DEBUG : print "Buffer", "SimulationCompute"     
        syspath = sys.path
        sys.path = self.syspath
        text = self.editor.getText()
        return process.runScript(self.doc.filepath)
    
