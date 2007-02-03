"""
$images.py

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

 Support for icons.

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""


import wx
import os
import cStringIO
import sys

def getPulsarIcon():
    icon = wx.EmptyIcon()
    icon.CopyFromBitmap(getPulsarBitmap())
    return icon

def getPulsarBitmap():
    return wx.BitmapFromImage(getPulsarImage())

def getPulsarImage():
    stream = cStringIO.StringIO(getPulsarData())
    return wx.ImageFromStream(stream)

def getPulsarData():
    try:
        #it may be a standalone application
        #in this case, the PYTHONPATH is changed
        fil=opj(os.path.join(os.path.dirname(sys.executable), 'pulsar.png'))
        return open(fil, "rb").read()
    except:
        #if not then load it as usual
        fil="images\pulsar.png"
        return open(fil, "rb").read()

def opj(path):
    """Convert paths to the platform-specific separator"""
    str = apply(os.path.join, tuple(path.split('/')))
    # HACK: on Linux, a leading / gets lost...
    if path.startswith('/'):
        str = '/' + str
    return str

