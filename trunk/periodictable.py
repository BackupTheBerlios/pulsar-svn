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


#TODO: put more precision on the data
# Dictionnary of nucleus for nmr
isotopes={
'1H':{ 'abundance':99.98, 'name':'Hydrogen 1H', 'relative':1.0, 'mass':1, 'quadrupole':0, 'larmor':-100, 'spin':0.5, 'absolute':1.0, },
'2D':{ 'abundance':0.015, 'name':'Deuterium 2D', 'relative':0.00965, 'mass':2, 'quadrupole':0.00273, 'larmor':-15.351, 'spin':1, 'absolute':1.45e-006, },
'3He':{ 'abundance':0.00013, 'name':'Helium 3He', 'relative':0.44, 'mass':3, 'quadrupole':0, 'larmor':76.178, 'spin':0.5, 'absolute':5.75e-007, },
'6Li':{ 'abundance':7.42, 'name':'Lithium 6Li', 'relative':0.0085, 'mass':6, 'quadrupole':-0.0008, 'larmor':-14.716, 'spin':1, 'absolute':0.000631, },
'7Li':{ 'abundance':92.58, 'name':'Lithium 7Li', 'relative':0.29, 'mass':7, 'quadrupole':-0.045, 'larmor':-38.863, 'spin':1.5, 'absolute':0.27, },
'9Be':{ 'abundance':100, 'name':'Beryllium 9Be', 'relative':0.0139, 'mass':9, 'quadrupole':0.052, 'larmor':14.053, 'spin':1.5, 'absolute':0.0139, },
'10B':{ 'abundance':19.58, 'name':'Boron 10B', 'relative':0.0199, 'mass':10, 'quadrupole':0.074, 'larmor':-10.746, 'spin':3, 'absolute':0.0039, },
'11B':{ 'abundance':80.42, 'name':'Boron 11B', 'relative':0.17, 'mass':11, 'quadrupole':0.0355, 'larmor':-32.084, 'spin':1.5, 'absolute':0.13, },
'13C':{ 'abundance':1.108, 'name':'Carbon 13C', 'relative':0.0159, 'mass':13, 'quadrupole':0, 'larmor':-25.144, 'spin':0.5, 'absolute':0.000176, },
'14N':{ 'abundance':99.63, 'name':'Nitrogen 14N', 'relative':0.00101, 'mass':14, 'quadrupole':0.016, 'larmor':-7.224, 'spin':1, 'absolute':0.00101, },
'15N':{ 'abundance':0.37, 'name':'Nitrogen 15N', 'relative':0.00104, 'mass':15, 'quadrupole':0, 'larmor':10.133, 'spin':0.5, 'absolute':3.85e-006, },
'17O':{ 'abundance':0.037, 'name':'Oxygen 17O', 'relative':0.0291, 'mass':17, 'quadrupole':-0.026, 'larmor':13.557, 'spin':2.5, 'absolute':1.08e-005, },
'19F':{ 'abundance':100, 'name':'Fluorine 19F', 'relative':0.83, 'mass':19, 'quadrupole':0, 'larmor':-94.077, 'spin':0.5, 'absolute':0.83, },
'21Ne':{ 'abundance':0.257, 'name':'Neon 21Ne', 'relative':0.0025, 'mass':21, 'quadrupole':0.09, 'larmor':7.894, 'spin':1.5, 'absolute':6.43e-006, },
'23Na':{ 'abundance':100, 'name':'Sodium 23Na', 'relative':0.0925, 'mass':23, 'quadrupole':0.12, 'larmor':-26.451, 'spin':1.5, 'absolute':0.0925, },
'25Mg':{ 'abundance':10.13, 'name':'Magnesium 25Mg', 'relative':0.00267, 'mass':25, 'quadrupole':0.22, 'larmor':6.1195, 'spin':2.5, 'absolute':0.000271, },
'27Al':{ 'abundance':100, 'name':'Aluminum 27Al', 'relative':0.21, 'mass':27, 'quadrupole':0.149, 'larmor':-26.057, 'spin':2.5, 'absolute':0.21, },
'29Si':{ 'abundance':4.7, 'name':'Silicon 29Si', 'relative':0.00784, 'mass':29, 'quadrupole':0, 'larmor':19.865, 'spin':0.5, 'absolute':0.000369, },
'31P':{ 'abundance':100, 'name':'Phosphorus 31P', 'relative':0.0663, 'mass':31, 'quadrupole':0, 'larmor':-40.481, 'spin':0.5, 'absolute':0.0663, },
'33S':{ 'abundance':0.76, 'name':'Sulfur 33S', 'relative':0.00226, 'mass':33, 'quadrupole':-0.055, 'larmor':-7.67, 'spin':1.5, 'absolute':1.72e-005, },
'35Cl':{ 'abundance':75.53, 'name':'Chlorine 35Cl', 'relative':0.0047, 'mass':35, 'quadrupole':-0.08, 'larmor':-9.798, 'spin':1.5, 'absolute':0.00355, },
'37Cl':{ 'abundance':24.47, 'name':'Chlorine 37Cl', 'relative':0.00271, 'mass':37, 'quadrupole':-0.0632, 'larmor':-8.156, 'spin':1.5, 'absolute':0.000663, },
'39K':{ 'abundance':93.1, 'name':'Potassium 39K', 'relative':0.000508, 'mass':39, 'quadrupole':0.055, 'larmor':-4.667, 'spin':1.5, 'absolute':0.000473, },
'41K':{ 'abundance':6.88, 'name':'Potassium 41K', 'relative':8.4e-005, 'mass':41, 'quadrupole':0.067, 'larmor':-2.561, 'spin':1.5, 'absolute':5.78e-006, },
'43Ca':{ 'abundance':0.145, 'name':'Calcium 43Ca', 'relative':0.0064, 'mass':43, 'quadrupole':-0.05, 'larmor':6.728, 'spin':3.5, 'absolute':9.28e-006, },
'45Sc':{ 'abundance':100, 'name':'Scandium 45Sc', 'relative':0.3, 'mass':45, 'quadrupole':-0.22, 'larmor':-24.29, 'spin':3.5, 'absolute':0.3, },
'47Ti':{ 'abundance':7.28, 'name':'Titanium 47Ti', 'relative':0.00209, 'mass':47, 'quadrupole':0.29, 'larmor':-5.637, 'spin':2.5, 'absolute':0.000152, },
'49Ti':{ 'abundance':5.51, 'name':'Titanium 49Ti', 'relative':0.00376, 'mass':49, 'quadrupole':0.24, 'larmor':-5.638, 'spin':3.5, 'absolute':0.000207, },
'50V':{ 'abundance':0.24, 'name':'Vanadium 50V', 'relative':0.0555, 'mass':50, 'quadrupole':0.21, 'larmor':-9.97, 'spin':6, 'absolute':0.000133, },
'51V':{ 'abundance':99.76, 'name':'Vanadium 51V', 'relative':0.38, 'mass':51, 'quadrupole':-0.052, 'larmor':-26.289, 'spin':3.5, 'absolute':0.38, },
'53Cr':{ 'abundance':9.55, 'name':'Chromium 53Cr', 'relative':0.000903, 'mass':53, 'quadrupole':0.03, 'larmor':5.652, 'spin':1.5, 'absolute':0.00862, },
'55Mn':{ 'abundance':100, 'name':'Manganese 55Mn', 'relative':0.18, 'mass':55, 'quadrupole':0.55, 'larmor':-24.664, 'spin':2.5, 'absolute':0.18, },
'57Fe':{ 'abundance':2.19, 'name':'Iron 57Fe', 'relative':3.37e-005, 'mass':57, 'quadrupole':0, 'larmor':-3.231, 'spin':0.5, 'absolute':7.38e-007, },
'59Co':{ 'abundance':100, 'name':'Cobalt 59Co', 'relative':0.28, 'mass':59, 'quadrupole':0.4, 'larmor':-23.614, 'spin':3.5, 'absolute':0.28, },
'61Ni':{ 'abundance':1.19, 'name':'Nickel 61Ni', 'relative':0.00357, 'mass':61, 'quadrupole':0.16, 'larmor':8.936, 'spin':1.5, 'absolute':4.25e-005, },
'63Cu':{ 'abundance':69.09, 'name':'Copper 63Cu', 'relative':0.0931, 'mass':63, 'quadrupole':-0.211, 'larmor':-26.505, 'spin':1.5, 'absolute':0.0643, },
'65Cu':{ 'abundance':30.91, 'name':'Copper 65Cu', 'relative':0.11, 'mass':65, 'quadrupole':-0.195, 'larmor':-28.394, 'spin':1.5, 'absolute':0.0352, },
'67Zn':{ 'abundance':4.11, 'name':'Zinc 67Zn', 'relative':0.00285, 'mass':67, 'quadrupole':0.15, 'larmor':-6.254, 'spin':2.5, 'absolute':0.000117, },
'69Ga':{ 'abundance':60.4, 'name':'Gallium 69Ga', 'relative':0.0691, 'mass':69, 'quadrupole':0.178, 'larmor':-24.003, 'spin':1.5, 'absolute':0.0417, },
'71Ga':{ 'abundance':39.6, 'name':'Gallium 71Ga', 'relative':0.14, 'mass':71, 'quadrupole':0.112, 'larmor':-30.495, 'spin':1.5, 'absolute':0.0562, },
'73Ge':{ 'abundance':7.76, 'name':'Germanium 73Ge', 'relative':0.0014, 'mass':73, 'quadrupole':-0.2, 'larmor':3.488, 'spin':4.5, 'absolute':0.000108, },
'75As':{ 'abundance':100, 'name':'Arsenic 75As', 'relative':0.0251, 'mass':75, 'quadrupole':0.3, 'larmor':-17.126, 'spin':1.5, 'absolute':0.0251, },
'77Se':{ 'abundance':7.58, 'name':'Selenium 77Se', 'relative':0.00693, 'mass':77, 'quadrupole':0, 'larmor':-19.067, 'spin':0.5, 'absolute':0.000525, },
'79Br':{ 'abundance':50.54, 'name':'Bromine 79Br', 'relative':0.0786, 'mass':79, 'quadrupole':0.33, 'larmor':-25.053, 'spin':1.5, 'absolute':0.0397, },
'81Br':{ 'abundance':49.46, 'name':'Bromine 81Br', 'relative':0.0985, 'mass':81, 'quadrupole':0.28, 'larmor':-27.006, 'spin':1.5, 'absolute':0.0487, },
'83Kr':{ 'abundance':11.55, 'name':'Krypton 83Kr', 'relative':0.00188, 'mass':83, 'quadrupole':0.15, 'larmor':3.847, 'spin':4.5, 'absolute':0.000217, },
'85Rb':{ 'abundance':72.15, 'name':'Rubidium 85Rb', 'relative':0.0105, 'mass':85, 'quadrupole':0.25, 'larmor':-9.655, 'spin':2.5, 'absolute':0.00757, },
'87Rb':{ 'abundance':27.85, 'name':'Rubidium 87Rb', 'relative':0.17, 'mass':87, 'quadrupole':0.12, 'larmor':-32.721, 'spin':1.5, 'absolute':0.0487, },
'87Sr':{ 'abundance':7.02, 'name':'Strontium 87Sr', 'relative':0.00269, 'mass':87, 'quadrupole':0.36, 'larmor':4.333, 'spin':4.5, 'absolute':0.000188, },
'89Y':{ 'abundance':100, 'name':'Yttrium 89Y', 'relative':0.000118, 'mass':89, 'quadrupole':0, 'larmor':4.899, 'spin':0.5, 'absolute':0.000118, },
'91Zr':{ 'abundance':11.23, 'name':'Zirconium 91Zr', 'relative':0.00948, 'mass':91, 'quadrupole':-0.21, 'larmor':9.33, 'spin':2.5, 'absolute':0.00106, },
'93Nb':{ 'abundance':100, 'name':'Niobium 93Nb', 'relative':0.48, 'mass':93, 'quadrupole':-0.2, 'larmor':-24.442, 'spin':4.5, 'absolute':0.48, },
'95Mo':{ 'abundance':15.72, 'name':'Molybdenum 95Mo', 'relative':0.00323, 'mass':95, 'quadrupole':0.12, 'larmor':-6.514, 'spin':2.5, 'absolute':0.000507, },
'97Mo':{ 'abundance':9.46, 'name':'Molybdenum 97Mo', 'relative':0.00343, 'mass':97, 'quadrupole':1.1, 'larmor':6.652, 'spin':2.5, 'absolute':0.000324, },
'99Ru':{ 'abundance':12.72, 'name':'Ruthenium 99Ru', 'relative':0.000195, 'mass':99, 'quadrupole':-0.19, 'larmor':4.605, 'spin':2.5, 'absolute':2.48e-005, },
'101Ru':{ 'abundance':17.07, 'name':'Ruthenium 101Ru', 'relative':0.00141, 'mass':101, 'quadrupole':0.076, 'larmor':5.161, 'spin':2.5, 'absolute':0.00024, },
'103Rh':{ 'abundance':100, 'name':'Rhodium 103Rh', 'relative':3.11e-005, 'mass':103, 'quadrupole':0, 'larmor':3.147, 'spin':0.5, 'absolute':3.11e-005, },
'105Pd':{ 'abundance':22.23, 'name':'Palladium 105Pd', 'relative':0.00112, 'mass':105, 'quadrupole':-0.8, 'larmor':4.576, 'spin':2.5, 'absolute':0.000249, },
'107Ag':{ 'abundance':51.82, 'name':'Silver 107Ag', 'relative':6.62e-005, 'mass':107, 'quadrupole':0, 'larmor':4.046, 'spin':0.5, 'absolute':3.43e-005, },
'109Ag':{ 'abundance':48.18, 'name':'Silver 109Ag', 'relative':0.000101, 'mass':109, 'quadrupole':0, 'larmor':4.652, 'spin':0.5, 'absolute':4.86e-005, },
'111Cd':{ 'abundance':12.75, 'name':'Cadmium 111Cd', 'relative':0.00954, 'mass':111, 'quadrupole':0, 'larmor':21.205, 'spin':0.5, 'absolute':0.00121, },
'113Cd':{ 'abundance':12.26, 'name':'Cadmium 113Cd', 'relative':0.0109, 'mass':113, 'quadrupole':0, 'larmor':22.182, 'spin':0.5, 'absolute':0.00133, },
'113In':{ 'abundance':4.28, 'name':'Indium 113In', 'relative':0.34, 'mass':113, 'quadrupole':1.14, 'larmor':-21.866, 'spin':4.5, 'absolute':0.0147, },
'115In':{ 'abundance':95.72, 'name':'Indium 115In', 'relative':0.34, 'mass':115, 'quadrupole':0.83, 'larmor':-21.914, 'spin':4.5, 'absolute':0.33, },
'117Sn':{ 'abundance':7.61, 'name':'Tin 117Sn', 'relative':0.0452, 'mass':117, 'quadrupole':0, 'larmor':35.625, 'spin':0.5, 'absolute':0.00344, },
'119Sn':{ 'abundance':8.58, 'name':'Tin 119Sn', 'relative':0.0518, 'mass':119, 'quadrupole':0, 'larmor':37.272, 'spin':0.5, 'absolute':0.00444, },
'121Sb':{ 'abundance':57.25, 'name':'Antimony 121Sb', 'relative':0.16, 'mass':121, 'quadrupole':-0.53, 'larmor':-23.93, 'spin':2.5, 'absolute':0.0916, },
'123Sb':{ 'abundance':42.75, 'name':'Antimony 123Sb', 'relative':0.0457, 'mass':123, 'quadrupole':-0.68, 'larmor':-12.959, 'spin':3.5, 'absolute':0.0195, },
'123Te':{ 'abundance':0.87, 'name':'Tellurium 123Te', 'relative':0.018, 'mass':123, 'quadrupole':0, 'larmor':26.207, 'spin':0.5, 'absolute':0.000156, },
'125Te':{ 'abundance':6.99, 'name':'Tellurium 125Te', 'relative':0.0315, 'mass':125, 'quadrupole':0, 'larmor':31.596, 'spin':0.5, 'absolute':0.0022, },
'127I':{ 'abundance':100, 'name':'Iodine 127I', 'relative':0.0934, 'mass':127, 'quadrupole':-0.79, 'larmor':-20.007, 'spin':2.5, 'absolute':0.0934, },
'129Xe':{ 'abundance':26.44, 'name':'Xenon 129Xe', 'relative':0.0212, 'mass':129, 'quadrupole':0, 'larmor':27.66, 'spin':0.5, 'absolute':0.0056, },
'131Xe':{ 'abundance':21.18, 'name':'Xenon 131Xe', 'relative':0.00276, 'mass':131, 'quadrupole':-0.12, 'larmor':-8.199, 'spin':1.5, 'absolute':0.000584, },
'133Cs':{ 'abundance':100, 'name':'Cesium 133Cs', 'relative':0.0474, 'mass':133, 'quadrupole':-0.003, 'larmor':-13.117, 'spin':3.5, 'absolute':0.0474, },
'135Ba':{ 'abundance':6.59, 'name':'Barium 135Ba', 'relative':0.0049, 'mass':135, 'quadrupole':0.18, 'larmor':-9.934, 'spin':1.5, 'absolute':0.000322, },
'137Ba':{ 'abundance':11.32, 'name':'Barium 137Ba', 'relative':0.00686, 'mass':137, 'quadrupole':0.28, 'larmor':-11.113, 'spin':1.5, 'absolute':0.000776, },
'138La':{ 'abundance':0.089, 'name':'Lanthanum 138La', 'relative':0.0919, 'mass':138, 'quadrupole':-0.47, 'larmor':-13.193, 'spin':5, 'absolute':8.18e-005, },
'139La':{ 'abundance':99.91, 'name':'Lanthanum 139La', 'relative':0.0592, 'mass':139, 'quadrupole':0.21, 'larmor':-14.126, 'spin':3.5, 'absolute':0.0591, },
'177Hf':{ 'abundance':18.5, 'name':'Hafnium 177Hf', 'relative':0.000638, 'mass':177, 'quadrupole':4.5, 'larmor':-3.12, 'spin':3.5, 'absolute':0.000118, },
'179Hf':{ 'abundance':13.75, 'name':'Hafnium 179Hf', 'relative':0.000216, 'mass':179, 'quadrupole':5.1, 'larmor':1.869, 'spin':4.5, 'absolute':2.97e-005, },
'181Ta':{ 'abundance':99.98, 'name':'Tantalum 181Ta', 'relative':0.036, 'mass':181, 'quadrupole':3.0, 'larmor':-11.97, 'spin':3.5, 'absolute':0.036, },
'183W':{ 'abundance':14.4, 'name':'Tungsten 183W', 'relative':0.00072, 'mass':183, 'quadrupole':0, 'larmor':-4.161, 'spin':0.5, 'absolute':1.03e-005, },
'185Re':{ 'abundance':37.07, 'name':'Rhenium 185Re', 'relative':0.13, 'mass':185, 'quadrupole':2.8, 'larmor':-22.513, 'spin':2.5, 'absolute':0.0493, },
'187Os':{ 'abundance':1.64, 'name':'Osmium 187Os', 'relative':1.22e-005, 'mass':187, 'quadrupole':0, 'larmor':-2.303, 'spin':0.5, 'absolute':2e-007, },
'187Re':{ 'abundance':62.93, 'name':'Rhenium 187Re', 'relative':0.13, 'mass':187, 'quadrupole':2.6, 'larmor':-22.744, 'spin':2.5, 'absolute':0.0862, },
'189Os':{ 'abundance':16.1, 'name':'Osmium 189Os', 'relative':0.00234, 'mass':189, 'quadrupole':0.8, 'larmor':-7.758, 'spin':1.5, 'absolute':0.000376, },
'191Ir':{ 'abundance':37.3, 'name':'Iridium 191Ir', 'relative':2.53e-005, 'mass':191, 'quadrupole':1.5, 'larmor':-1.718, 'spin':1.5, 'absolute':9.43e-006, },
'193Ir':{ 'abundance':62.7, 'name':'Iridium 193Ir', 'relative':3.27e-005, 'mass':193, 'quadrupole':1.4, 'larmor':-1.871, 'spin':1.5, 'absolute':2.05e-005, },
'195Pt':{ 'abundance':33.8, 'name':'Platinum 195Pt', 'relative':0.00994, 'mass':195, 'quadrupole':0, 'larmor':-21.499, 'spin':0.5, 'absolute':0.00336, },
'197Au':{ 'abundance':100, 'name':'Gold 197Au', 'relative':2.51e-005, 'mass':197, 'quadrupole':0.58, 'larmor':-1.712, 'spin':1.5, 'absolute':2.51e-005, },
'199Hg':{ 'abundance':16.84, 'name':'Mercury 199Hg', 'relative':0.00567, 'mass':199, 'quadrupole':0, 'larmor':-17.827, 'spin':0.5, 'absolute':0.000954, },
'201Hg':{ 'abundance':13.22, 'name':'Mercury 201Hg', 'relative':0.00144, 'mass':201, 'quadrupole':0.5, 'larmor':6.599, 'spin':1.5, 'absolute':0.00019, },
'203Tl':{ 'abundance':29.5, 'name':'Thallium 203Tl', 'relative':0.18, 'mass':203, 'quadrupole':0, 'larmor':-57.149, 'spin':0.5, 'absolute':0.0551, },
'205Tl':{ 'abundance':70.5, 'name':'Thallium 205Tl', 'relative':0.19, 'mass':205, 'quadrupole':0, 'larmor':-57.708, 'spin':0.5, 'absolute':0.13, },
'207Pb':{ 'abundance':22.6, 'name':'Lead 207Pb', 'relative':0.00916, 'mass':207, 'quadrupole':0, 'larmor':-20.921, 'spin':0.5, 'absolute':0.00207, },
'209Bi':{ 'abundance':100, 'name':'Bismuth 209Bi', 'relative':0.13, 'mass':209, 'quadrupole':-0.4, 'larmor':-16.069, 'spin':4.5, 'absolute':0.13, },
}
