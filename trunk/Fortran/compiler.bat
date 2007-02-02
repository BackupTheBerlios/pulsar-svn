@echo off
::$compiler.bat

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

:: COMPILATION OF THE PROGRAM PULSAR AS A PYTHON MODULE

:: ****************************************************
:: * This file start the compilation of f95pulsar.pyd *
:: * (an extension module for python)                 *
:: * It creates Fortran callable functions            *
:: * from the pulsar.f95 program compiled with g95    *
:: * or lahey                                         *
:: ****************************************************

cls

:: Compiler= lahey  - use lahey lf95 compiler
:: Compiler= g95    - use gnu g95 compiler

	set Compiler=lahey

:: Debug?

	set debug=0
	if %debug%==1	set debugoptions= --debug

:: delete the existing file

	if exist f95pulsar.pyd 			del f95pulsar.pyd 
        if exist *.mod                          del *.mod
        if exist *.*~                           del *.*~
        if exist *.bak                          del *.bak

:: compile

	f2py.py --verbose --fcompiler=%Compiler% %debugoptions% -m f95pulsar -c f95pulsar.f95

:: test result
 
	if not exist "f95pulsar.pyd"		goto:eof
 
:: copy the resulting file
 
	copy f95pulsar.pyd ..\f95pulsar.pyd 
        if exist *.mod                          del *.mod

:: execute
	cd ..
	Run_pyPulsar.cmd
	cd fortran

:eof
