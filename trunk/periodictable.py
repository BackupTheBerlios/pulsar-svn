"""
$periodictable.py

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

 List of isotopes

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""
__revision__ = "1"
__revision_date__="2007.02.02" 


#TODO: put more precision on the data
# Dictionnary of nucleus for nmr
isotopes={'1H':{'name':'Proton 1H','spin':0.5,'larmor':-100,'abundance':99.98,'absolute':1.00,'relative':1.00,'quadrupole':0},
'2D':{'name':'Deuteron 2D','spin':1,'larmor':-15.351,'abundance':0.015,'absolute':1.45E-6,'relative':9.65E-3,'quadrupole':2.73E-3},
'3He':{'name':'Helium 3He','spin':0.5,'larmor':76.178,'abundance':1.3E-4,'absolute':5.75E-7,'relative':0.44,'quadrupole':0},
'6Li':{'name':'Lithium 6Li','spin':1,'larmor':-14.716,'abundance':7.42,'absolute':6.31E-4,'relative':8.50E-3,'quadrupole':-8E-4},
'7Li':{'name':'Lithium 7Li','spin':1.5,'larmor':-38.863,'abundance':92.58,'absolute':0.27,'relative':0.29,'quadrupole':-0.045},
'9Be':{'name':'Beryllium 9Be','spin':1.5,'larmor':14.053,'abundance':100,'absolute':1.39E-2,'relative':1.39E-2,'quadrupole':0.052},
'10B':{'name':'Boron 10B','spin':3,'larmor':-10.746,'abundance':19.58,'absolute':3.90E-3,'relative':1.99E-2,'quadrupole':0.074},
'11B':{'name':'Boron 11B','spin':1.5,'larmor':-32.084,'abundance':80.42,'absolute':0.13,'relative':0.17,'quadrupole':0.0355},
'13C':{'name':'Carbon 13C','spin':0.5,'larmor':-25.144,'abundance':1.108,'absolute':1.76E-4,'relative':1.59E-2,'quadrupole':0},
'14N':{'name':'Nitrogen 14N','spin':1,'larmor':-7.224,'abundance':99.63,'absolute':1.01E-3,'relative':1.01E-3,'quadrupole':0.016},
'15N':{'name':'Nitrogen 15N','spin':0.5,'larmor':10.133,'abundance':0.37,'absolute':3.85E-6,'relative':1.04E-3,'quadrupole':0},
'17O':{'name':'Oxygen 17O','spin':2.5,'larmor':13.557,'abundance':0.037,'absolute':1.08E-5,'relative':2.91E-2,'quadrupole':-0.026},
'19F':{'name':'Fluorine 19F','spin':0.5,'larmor':-94.077,'abundance':100,'absolute':0.83,'relative':0.83,'quadrupole':0},
'21Ne':{'name':'Neon 21Ne','spin':1.5,'larmor':7.894,'abundance':0.257,'absolute':6.43E-6,'relative':2.50E-3,'quadrupole':0.09},
'23Na':{'name':'Sodium 23Na','spin':1.5,'larmor':-26.451,'abundance':100,'absolute':9.25E-2,'relative':9.25E-2,'quadrupole':0.12},
'25Mg':{'name':'Magnesium 25Mg','spin':2.5,'larmor':6.1195,'abundance':10.13,'absolute':2.71E-4,'relative':2.67E-3,'quadrupole':0.22},
'27Al':{'name':'Aluminum 27Al','spin':2.5,'larmor':-26.057,'abundance':100,'absolute':0.21,'relative':0.21,'quadrupole':0.149},
'29Si':{'name':'Silicon 29Si','spin':0.5,'larmor':19.865,'abundance':4.7,'absolute':3.69E-4,'relative':7.84E-3,'quadrupole':0},
'31P':{'name':'Phosphorus 31P','spin':0.5,'larmor':-40.481,'abundance':100,'absolute':6.63E-2,'relative':6.63E-2,'quadrupole':0},
'33S':{'name':'Sulfur 33S','spin':1.5,'larmor':-7.67,'abundance':0.76,'absolute':1.72E-5,'relative':2.26E-3,'quadrupole':-0.055},
'35Cl':{'name':'Chlorine 35Cl','spin':1.5,'larmor':-9.798,'abundance':75.53,'absolute':3.55E-3,'relative':4.70E-3,'quadrupole':-0.08},
'37Cl':{'name':'Chlorine 37Cl','spin':1.5,'larmor':-8.156,'abundance':24.47,'absolute':6.63E-4,'relative':2.71E-3,'quadrupole':-0.0632},
'39K':{'name':'Potassium 39K','spin':1.5,'larmor':-4.667,'abundance':93.1,'absolute':4.73E-4,'relative':5.08E-4,'quadrupole':0.055},
'41K':{'name':'Potassium 41K','spin':1.5,'larmor':-2.561,'abundance':6.88,'absolute':5.78E-6,'relative':8.40E-5,'quadrupole':0.067},
'43Ca':{'name':'Calcium 43Ca','spin':3.5,'larmor':6.728,'abundance':0.145,'absolute':9.28E-6,'relative':6.40E-3,'quadrupole':-0.05},
'45Sc':{'name':'Scandium 45Sc','spin':3.5,'larmor':-24.29,'abundance':100,'absolute':0.30,'relative':0.30,'quadrupole':-0.22},
'47Ti':{'name':'Titanium 47Ti','spin':2.5,'larmor':-5.637,'abundance':7.28,'absolute':1.52E-4,'relative':2.09E-3,'quadrupole':0.29},
'49Ti':{'name':'Titanium 49Ti','spin':3.5,'larmor':-5.638,'abundance':5.51,'absolute':2.07E-4,'relative':3.76E-3,'quadrupole':0.24},
'50V':{'name':'Vanadium 50V','spin':6,'larmor':-9.97,'abundance':0.24,'absolute':1.33E-4,'relative':5.55E-2,'quadrupole':0.21},
'51V':{'name':'Vanadium 51V','spin':3.5,'larmor':-26.289,'abundance':99.76,'absolute':0.38,'relative':0.38,'quadrupole':-0.052},
'53Cr':{'name':'Chromium 53Cr','spin':1.5,'larmor':5.652,'abundance':9.55,'absolute':8.62E-3,'relative':9.03E-4,'quadrupole':0.030},
'55Mn':{'name':'Manganese 55Mn','spin':2.5,'larmor':-24.664,'abundance':100,'absolute':0.18,'relative':0.18,'quadrupole':0.55},
'57Fe':{'name':'Iron 57Fe','spin':0.5,'larmor':-3.231,'abundance':2.19,'absolute':7.38E-7,'relative':3.37E-5,'quadrupole':0},
'59Co':{'name':'Cobalt 59Co','spin':3.5,'larmor':-23.614,'abundance':100,'absolute':0.28,'relative':0.28,'quadrupole':0.40},
'61Ni':{'name':'Nickel 61Ni','spin':1.5,'larmor':8.936,'abundance':1.19,'absolute':4.25E-5,'relative':3.57E-3,'quadrupole':0.16},
'63Cu':{'name':'Copper 63Cu','spin':1.5,'larmor':-26.505,'abundance':69.09,'absolute':6.43E-2,'relative':9.31E-2,'quadrupole':-0.211},
'65Cu':{'name':'Copper 65Cu','spin':1.5,'larmor':-28.394,'abundance':30.91,'absolute':3.52E-2,'relative':0.11,'quadrupole':-0.195},
'67Zn':{'name':'Zinc 67Zn','spin':2.5,'larmor':-6.254,'abundance':4.11,'absolute':1.17E-4,'relative':2.85E-3,'quadrupole':0.15},
'69Ga':{'name':'Gallium 69Ga','spin':1.5,'larmor':-24.003,'abundance':60.4,'absolute':4.17E-2,'relative':6.91E-2,'quadrupole':0.178},
'71Ga':{'name':'Gallium 71Ga','spin':1.5,'larmor':-30.495,'abundance':39.6,'absolute':5.62E-2,'relative':0.14,'quadrupole':0.112},
'73Ge':{'name':'Germanium 73Ge','spin':4.5,'larmor':3.488,'abundance':7.76,'absolute':1.08E-4,'relative':1.4E-3,'quadrupole':-0.2},
'75As':{'name':'Arsenic 75As','spin':1.5,'larmor':-17.126,'abundance':100,'absolute':2.51E-2,'relative':2.51E-2,'quadrupole':0.3},
'77Se':{'name':'Selenium 77Se','spin':0.5,'larmor':-19.067,'abundance':7.58,'absolute':5.25E-4,'relative':6.93E-3,'quadrupole':0},
'79Br':{'name':'Bromine 79Br','spin':1.5,'larmor':-25.053,'abundance':50.54,'absolute':3.97E-2,'relative':7.86E-2,'quadrupole':0.33},
'81Br':{'name':'Bromine 81Br','spin':1.5,'larmor':-27.006,'abundance':49.46,'absolute':4.87E-2,'relative':9.85E-2,'quadrupole':0.28},
'83Kr':{'name':'Krypton 83Kr','spin':4.5,'larmor':3.847,'abundance':11.55,'absolute':2.17E-4,'relative':1.88E-3,'quadrupole':0.15},
'85Rb':{'name':'Rubidium 85Rb','spin':2.5,'larmor':-9.655,'abundance':72.15,'absolute':7.57E-3,'relative':1.05E-2,'quadrupole':0.25},
'87Rb':{'name':'Rubidium 87Rb','spin':1.5,'larmor':-32.721,'abundance':27.85,'absolute':4.87E-2,'relative':0.17,'quadrupole':0.12},
'87Sr':{'name':'Strontium 87Sr','spin':4.5,'larmor':4.333,'abundance':7.02,'absolute':1.88E-4,'relative':2.69E-3,'quadrupole':0.36},
'89Y':{'name':'Yttrium 89Y','spin':0.5,'larmor':4.899,'abundance':100,'absolute':1.18E-4,'relative':1.18E-4,'quadrupole':0},
'91Zr':{'name':'Zirconium 91Zr','spin':2.5,'larmor':9.33,'abundance':11.23,'absolute':1.06E-3,'relative':9.48E-3,'quadrupole':-0.21},
'93Nb':{'name':'Niobium 93Nb','spin':4.5,'larmor':-24.442,'abundance':100,'absolute':0.48,'relative':0.48,'quadrupole':-0.2},
'95Mo':{'name':'Molybdenum 95Mo','spin':2.5,'larmor':-6.514,'abundance':15.72,'absolute':5.07E-4,'relative':3.23E-3,'quadrupole':0.12},
'97Mo':{'name':'Molybdenum 97Mo','spin':2.5,'larmor':6.652,'abundance':9.46,'absolute':3.24E-4,'relative':3.43E-3,'quadrupole':1.1},
'99Ru':{'name':'Ruthenium 99Ru','spin':2.5,'larmor':4.605,'abundance':12.72,'absolute':2.48E-5,'relative':1.95E-4,'quadrupole':-0.19},
'101Ru':{'name':'Ruthenium 101Ru','spin':2.5,'larmor':5.161,'abundance':17.07,'absolute':2.40E-4,'relative':1.41E-3,'quadrupole':0.076},
'103Rh':{'name':'Rhodium 103Rh','spin':0.5,'larmor':3.147,'abundance':100,'absolute':3.11E-5,'relative':3.11E-5,'quadrupole':0},
'105Pd':{'name':'Palladium 105Pd','spin':2.5,'larmor':4.576,'abundance':22.23,'absolute':2.49E-4,'relative':1.12E-3,'quadrupole':-0.8},
'107Ag':{'name':'Silver 107Ag','spin':0.5,'larmor':4.046,'abundance':51.82,'absolute':3.43E-5,'relative':6.62E-5,'quadrupole':0},
'109Ag':{'name':'Silver 109Ag','spin':0.5,'larmor':4.652,'abundance':48.18,'absolute':4.86E-5,'relative':1.01E-4,'quadrupole':0},
'111Cd':{'name':'Cadmium 111Cd','spin':0.5,'larmor':21.205,'abundance':12.75,'absolute':1.21E-3,'relative':9.54E-3,'quadrupole':0},
'113Cd':{'name':'Cadmium 113Cd','spin':0.5,'larmor':22.182,'abundance':12.26,'absolute':1.33E-3,'relative':1.09E-2,'quadrupole':0},
'113In':{'name':'Indium 113In','spin':4.5,'larmor':-21.866,'abundance':4.28,'absolute':1.47E-2,'relative':0.34,'quadrupole':1.14},
'115In':{'name':'Indium 115In','spin':4.5,'larmor':-21.914,'abundance':95.72,'absolute':0.33,'relative':0.34,'quadrupole':0.83},
'117Sn':{'name':'Tin 117Sn','spin':0.5,'larmor':35.625,'abundance':7.61,'absolute':3.44E-3,'relative':4.52E-2,'quadrupole':0},
'119Sn':{'name':'Tin 119Sn','spin':0.5,'larmor':37.272,'abundance':8.58,'absolute':4.44E-3,'relative':5.18E-2,'quadrupole':0},
'121Sb':{'name':'Antimony 121Sb','spin':2.5,'larmor':-23.93,'abundance':57.25,'absolute':9.16E-2,'relative':0.16,'quadrupole':-0.53},
'123Sb':{'name':'Antimony 123Sb','spin':3.5,'larmor':-12.959,'abundance':42.75,'absolute':1.95E-2,'relative':4.57E-2,'quadrupole':-0.68},
'123Te':{'name':'Tellurium 123Te','spin':0.5,'larmor':26.207,'abundance':0.87,'absolute':1.56E-4,'relative':1.80E-2,'quadrupole':0},
'125Te':{'name':'Tellurium 125Te','spin':0.5,'larmor':31.596,'abundance':6.99,'absolute':2.20E-3,'relative':3.15E-2,'quadrupole':0},
'127I':{'name':'Iodine  127I','spin':2.5,'larmor':-20.007,'abundance':100,'absolute':9.34E-2,'relative':9.34E-2,'quadrupole':-0.79},
'129Xe':{'name':'Xenon 129Xe','spin':0.5,'larmor':27.66,'abundance':26.44,'absolute':5.60E-3,'relative':2.12E-2,'quadrupole':0},
'131Xe':{'name':'Xenon 131Xe','spin':1.5,'larmor':-8.199,'abundance':21.18,'absolute':5.84E-4,'relative':2.76E-3,'quadrupole':-0.12},
'133Cs':{'name':'Cesium 133Cs','spin':3.5,'larmor':-13.117,'abundance':100,'absolute':4.74E-2,'relative':4.74E-2,'quadrupole':-0.0030},
'135Ba':{'name':'Barium 135Ba','spin':1.5,'larmor':-9.934,'abundance':6.59,'absolute':3.22E-4,'relative':4.90E-3,'quadrupole':0.18},
'137Ba':{'name':'Barium 137Ba','spin':1.5,'larmor':-11.113,'abundance':11.32,'absolute':7.76E-4,'relative':6.86E-3,'quadrupole':0.28},
'138La':{'name':'Lanthanum 138La','spin':5,'larmor':-13.193,'abundance':0.089,'absolute':8.18E-5,'relative':9.19E-2,'quadrupole':-0.47},
'139La':{'name':'Lanthanum 139La' ,'spin':3.5,'larmor':-14.126,'abundance':99.91,'absolute':5.91E-2,'relative':5.92E-2,'quadrupole':0.21},
'177Hf':{'name':'Hafnium 177Hf','spin':3.5,'larmor':-3.12,'abundance':18.5,'absolute':1.18E-4,'relative':6.38E-4,'quadrupole':4.5},
'179Hf':{'name':'Hafnium 179Hf','spin':4.5,'larmor':1.869,'abundance':13.75,'absolute':2.97E-5,'relative':2.16E-4,'quadrupole':5.1},
'181Ta':{'name':'Tantalum 181Ta','spin':3.5,'larmor':-11.97,'abundance':99.98,'absolute':3.60E-2,'relative':3.60E-2,'quadrupole':3.0},
'183W':{'name':'Tungsten 183W','spin':0.5,'larmor':-4.161,'abundance':14.4,'absolute':1.03E-5,'relative':7.20E-4,'quadrupole':0},
'185Re':{'name':'Rhenium 185Re','spin':2.5,'larmor':-22.513,'abundance':37.07,'absolute':4.93E-2,'relative':0.13,'quadrupole':2.8},
'187Re':{'name':'Rhenium 187Re','spin':2.5,'larmor':-22.744,'abundance':62.93,'absolute':8.62E-2,'relative':0.13,'quadrupole':2.6},
'187Os':{'name':'Osmium 187Os','spin':0.5,'larmor':-2.303,'abundance':1.64,'absolute':2.00E-7,'relative':1.22E-5,'quadrupole':0},
'189Os':{'name':'Osmium 189Os','spin':1.5,'larmor':-7.758,'abundance':16.1,'absolute':3.76E-4,'relative':2.34E-3,'quadrupole':0.8},
'191Ir':{'name':'Iridium 191Ir','spin':1.5,'larmor':-1.718,'abundance':37.3,'absolute':9.43E-6,'relative':2.53E-5,'quadrupole':1.5},
'193Ir':{'name':'Iridium 193Ir','spin':1.5,'larmor':-1.871,'abundance':62.7,'absolute':2.05E-5,'relative':3.27E-5,'quadrupole':1.4},
'195Pt':{'name':'Platinum 195Pt','spin':0.5,'larmor':-21.499,'abundance':33.8,'absolute':3.36E-3,'relative':9.94E-3,'quadrupole':0},
'197Au':{'name':'Gold 197Au','spin':1.5,'larmor':-1.712,'abundance':100,'absolute':2.51E-5,'relative':2.51E-5,'quadrupole':0.58},
'199Hg':{'name':'Mercury 199Hg','spin':0.5,'larmor':-17.827,'abundance':16.84,'absolute':9.54E-4,'relative':5.67E-3,'quadrupole':0},
'201Hg':{'name':'Mercury 201Hg','spin':1.5,'larmor':6.599,'abundance':13.22,'absolute':1.90E-4,'relative':1.44E-3,'quadrupole':0.5},
'203Tl':{'name':'Thallium 203Tl','spin':0.5,'larmor':-57.149,'abundance':29.5,'absolute':5.51E-2,'relative':0.18,'quadrupole':0},
'205Tl':{'name':'Thallium 205Tl','spin':0.5,'larmor':-57.708,'abundance':70.5,'absolute':0.13,'relative':0.19,'quadrupole':0},
'207Pb':{'name':'Lead 207Pb','spin':0.5,'larmor':-20.921,'abundance':22.6,'absolute':2.07E-3,'relative':9.16E-3,'quadrupole':0},
'209Bi':{'name':'Bismuth 209 Bi','spin':4.5,'larmor':-16.069,'abundance':100,'absolute':0.13,'relative':0.13,'quadrupole':-0.4}}
