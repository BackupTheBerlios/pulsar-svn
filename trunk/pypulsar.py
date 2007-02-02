"""
$pypulsar.py

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

 Starting python code

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""
__revision__ = "21"
__revision_date__="2007.02.02" 

from main import *
import sys
import pytz.zoneinfo.UTC  # NEEDED BY PYINSTALLER (I don't know why?)

#---------------------------------------------------------------------------   
if __name__ == "__main__":
    import psyco
    try:
        psyco.full()
        print "psyco started."
    except ImportError:
        print "psyco error."
        pass
    main(sys.argv)
