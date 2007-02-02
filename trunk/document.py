"""
$document.py

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

 Document Class

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""
__revision__ = "1"
__revision_date__="2007.02.02" 


import os


class Document:
    """Document class."""

    def __init__(self, filename=None):
        """Create a Document instance."""
        self.filename = filename
        self.filepath = None
        self.filedir = None
        self.filebase = None
        self.fileext = None
        if self.filename:
            self.filepath = os.path.realpath(self.filename)
            self.filedir, self.filename = os.path.split(self.filepath)
            self.filebase, self.fileext = os.path.splitext(self.filename)

    def read(self):
        """Return contents of file."""
        if self.filepath and os.path.exists(self.filepath):
            f = file(self.filepath, 'rb')
            try:
                return f.read()
            finally:
                f.close()
        else:
            return ''

    def write(self, text):
        """Write text to file."""
        try:
            f = file(self.filepath, 'wb')
            f.write(text)
        finally:
            if f:
                f.close()
