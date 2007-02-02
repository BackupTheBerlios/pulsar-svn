::$Build_pyPulsar.bat

::PULSAR Project:

:: Copyright (C) 2006-2007 Jean-Paul Amoureux, Christian Fernandez
:: JPA - Unite de Catalyse et Chimie du Solide, Lille, France.
:: CF  - Laboratoire Catalyse et Spectrochimie, Caen, France.

::LICENSE:

:: This program is free software; you can redistribute it and/or
:: modify it under the terms of the GNU General Public License (GPL)
:: as published by the Free Software Foundation; either version 2
:: of the License, or (at your option) any later version.
:: This program is distributed in the hope that it will be useful,
:: but WITHOUT ANY WARRANTY; without even the implied warranty of
:: MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
:: GNU General Public License for more details.
:: To read the license please visit http://www.gnu.org/copyleft/gpl.html

::PURPOSE OF THIS FILE:

:: BUILD EXE FILE WITH PYINSTALLER (PYI)
::****************************************************************************************************
:: I was unable to make a standalone exe using Py2EXE (because of MAtplotLib and Scipy usage)
:: For PYinstaller it is also a problem, but a work around this problem, showed me that it is
:: only necessary to copy the missing import folder of Scipy and the mpl-data directory of Matplotlib
:: The later has been renamed to "matplotlib" to work properly.
::**************************************************************************************************** 

@echo off

set version=0.1.rc1

:: debug=0  (debug false)
:: debug =1 (debug true)

set debug=0

:: CHANGE THESE LINES ACCORDING TO THE INSTALLATION DIRECTORIES
::-------------------------------------------------------------
set PYIDIR=D:\program_files\Pulsar\pyinstaller
set SRCDIR=D:\program_files\Pulsar\trunk
set DESTDIR=%SRCDIR%\pypulsar-%version%

:: Building Command
::-----------------
set BUILD=%PYIDIR%\build.py

:: Build...
::---------
python %BUILD% %SRCDIR%\pypulsar.spec

:: Copy some missing files (needed to make the executable working with MatplotLib and Scipy
::------------------------------------------------------------------------------------------
xcopy %SRCDIR%\Add_to_dist\*.* "%DESTDIR%\"  /S /D /Y

:: Execute the application
::------------------------
%DESTDIR%\pypulsar.exe
