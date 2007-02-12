"""
$functions.py

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

 This file contains some of the routines used by pulsar 
 for the input of the starting parameters

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""

#-------------------------------------------------------------------------------
#  Imports:
#-------------------------------------------------------------------------------
# Major package imports.
import wx
import os 
import sys 
import re 
import types 
import string
import math
from Numeric import *
from editor import *
from f95pulsar import *
from periodictable import *
from scipy.fftpack import ifft,fft

import os.path, StringIO

from debug import *

# module constants
#------------------
#logicals
TRUE = true = YES = Yes = yes = Y = y = True
FALSE = false = NO = No = no = N = n = False

#units
HZ = Hz = hz = 1.0
KHZ = kHz = khz = 1.0e3
MHZ = MHz = mhz = 1.0e6
USEC = usec = 1.0
MSEC = msec = 1.0e3
SEC = sec = 1.0e6
TESLA = Tesla = tesla = 1
ANGSTROM = ANGS = angs = angstrom = 1           
NM = nm = 10   #nanometer
K = k = 1024   #kpoint

#miscellaneous
MAS = mas = 54.735610
STATIC = static = 0
GAMMA_H = gamma_H = 26.7510e7/(2.0*pi)
AUTO = auto = -999
FULL = full = -998
ALL = all = 0
CENTRAL = Central = central = 1
SATELLITE = Satellite = satellite = 2

#==============================================================================
class Simulation:
#==============================================================================
    """
    A class to handle the parameters relative to a single pyPulsar simulation
    """
    
    current_spectrum = []
    collection_spectra = []
    collection_data = []
    collection_label=[]
    sample=None
    sequence=None

    def __del__(self):
        WRITE_STRING("DESTROY this  simulation instance")
###        reset()
##        self.current_spectrum = []
##        self.collection_spectra = []
##        self.collection_data = []
##        self.collection_label=[]
        
    #----------------------------------------------------------------------------------------
    def __init__(self, protonfrequency = "400 MHz", verbose=False, debug=False):
    #----------------------------------------------------------------------------------------
        """
        Class constructor -additionally set default parameters for the simulation
        """
        #WRITE_STRING("INIT a simulation instance")
        #it is mandatory to reset f95pulsar when doing a new simulation
        #that means zero some arrays that can have different size in
        #different simulation

        #Reset most of the already allocated parameters
        self.verbose=verbose or debug
        self.debug=debug
        self.reset_pulsar()
        self.reset_nucleus()
        self.reset_spectra()
        
        self.current_spectrum = []
        self.collection_spectra = []
        self.collection_data = []
        self.collection_label=[]
       
        # Verbose?
        self.verbose=verbose or debug
        self.debug=debug
        
        #Fortan flags
        parameters.verbose=self.verbose
        parameters.debug=self.debug
        
        # initialise some basic parameters for the simulation
        self.set_protonfrequency(protonfrequency)            
        self.spinningspeed=0.
        self.spinningangle=MAS
        self.qfactor=0.001
        self.accuracy=8
        self.rfstep=5
        self.npts=4*K
        self.sw=10*kHz
        self.aliasing=no
        self.nsb=0
        self.nall=0
        self.flag_sum=False
        
        # load all new parameters in Pulsar
        # the most important here is the protonfrequency
        self.validate()
        
        # initialise a "void" sample
        #   zero nuclei inside
        self.sample=Sample()

        # initialise a "void" pulse sequence
        self.sequence=PulseSequence()
               
    #----------------------
    def reset_pulsar(self):
    #----------------------
        """
        Deallocate most of the important arrays in the f95Pulsar module
        """
        DEBUG_MSG("RESET PULSAR")
        reset()      #CALL TO A F95PULSAR reset function
        
    #-----------------------
    def reset_spectra(self):
    #-----------------------
        """
        Delete all spectra in a collection to zero (initialisation)  
        """
        DEBUG_MSG("RESET SPECTRA")
        self.collection_spectra = []
        self.list_lb=[]
        self.list_gb=[]
        self.list_concentration=[]
        
    #-----------------------
    def reset_pulse(self):
    #-----------------------
        """
        Deallocate pulse arrays (initialisation) 
        """
        DEBUG_MSG("RESET PULSE SEQUENCE")
        parameters.pulse = []
        parameters.delay = []
        parameters.coher = []
        

    #-----------------------    
    def reset_nucleus(self):
    #-----------------------
        """
        Deallocate the parameter arrays concerning the nuclei
        """
        DEBUG_MSG("RESET NUCLEUS")
        parameters.nucleus=[]
        
    #-------------------------------------------------------
    def execute_pulsar(self,TEST=False):
    #-------------------------------------------------------
        """
        Run the simulation using the f95Pulsar module
        """
        DEBUG_MSG("EXECUTE PULSAR")
        # When running the full program we need to reset all arrays to None
        # and load all the new parameters
        self.validate()

        #reset_share()
        #reset_operators()
        control_parameters(pydebug_msg)
        diagonalize.dealloc()
        
        # No nucleus?
        if not self.sample.nuclei:
            return
        
        # Settings in case of spinning sideband aliasing
        if self.aliasing:
            DEBUG_MSG("ALIASING=YES")
            self.set_nsb("FULL")
            parameters.nsb = self.nsb
            parameters.sw = self.spinningspeed*float(self.nsb*2)
            parameters.npts = int(self.npts*(parameters.sw/self.sw))         

        # Display the current parameters if needed
        self.display()
        
        #Execute the simulation
        if not TEST : compute(pydebug_msg,pywrite_string)
           
    #-----------------
    def display(self):
    #-----------------
        if verbose:
            WRITE_STRING("\nSPECTROMETER:")
            WRITE_STRING(  "-------------")
            WRITE_STRING("\tSpectrometer field: " + strclean(parameters.spectrometerfield,5)+" Tesla")
            WRITE_STRING("\tProton frequency: " + strclean(parameters.protonfrequency/mhz,5)+" MHz")
            WRITE_STRING("\nPROBEHEAD:")
            WRITE_STRING("------------")
            WRITE_STRING("\tspinning speed: " + str(parameters.spinningspeed/khz)+" kHz")
            WRITE_STRING("\tspinning angle: " + strclean(parameters.spinningangle,8)+" deg.")
            WRITE_STRING("\tquality factor: " + str(parameters.qfactor))
            WRITE_STRING("\nSIMULATION PARAMETERS:")
            WRITE_STRING("----------------------")
            WRITE_STRING("\tnumber of points: " + str(parameters.npts))
            WRITE_STRING("\tspectral width: " + str(parameters.sw/khz)+" kHz")  
            WRITE_STRING("\tnumber of spinning sidebands: " + str(parameters.nsb*2))
            WRITE_STRING("\tpowder integration accuracy: " + str(parameters.accuracy))
            WRITE_STRING("\trf integration step: " + str(parameters.rfstep))
            if (parameters.nall == 0): strg = "-- detect all transitions"
            if (parameters.nall == 1): strg = "-- detect central transition only"
            if (parameters.nall == 2): strg = "-- detect satellite transitions only"    
            WRITE_STRING("\tDetection: "+ strg)
            WRITE_STRING("\nOTHER INFORMATIONS:")
            WRITE_STRING("-------------------")
            WRITE_STRING("\tAliasing of the spinning sidebands: "+str(self.aliasing))

    #------------------    
    def validate(self):
    #------------------    
        """
        set parameters in the fortran module
        """
        DEBUG_MSG("VALIDATE SIMULATION")

        parameters.spectrometerfield = self.field
        parameters.protonfrequency = self.protonfrequency
        parameters.spinningangle = self.spinningangle
        parameters.spinningspeed = self.spinningspeed
        parameters.qfactor = self.qfactor
        parameters.rfstep = self.rfstep
        parameters.npts = self.npts
        parameters.nsb = self.nsb
        parameters.accuracy = self.accuracy
        parameters.sw = self.sw
        parameters.nall = self.nall            
                        
    #--------------------------    
    def set_field(self, field):
    #--------------------------
        """
        Change the spectrometer field
        and calculate the corresponding proton frequency
        """
        DEBUG_MSG("SET_FIELD")
        try:
            self.field = abs(evaluate(field))
        except:
            self.field=0
        if self.field == 0 :
            WRITE_STRING("\n\terror in set_field: reset to default 9.395*tesla\n")
            self.field = 9.395*tesla
        self.protonfrequency = self.field*gamma_H

    #----------------------------------------
    def set_frequency(self, protonfrequency):
    #----------------------------------------
        """ similar to set_protonfrequency
            kept for backward compatibility """
        self.set_protonfrequency(protonfrequency)
        
    #----------------------------------------
    def set_protonfrequency(self, protonfrequency):
    #----------------------------------------
        """
        Change the proton frequency
        and calculate the corresponding spectrometer field   
        """
        DEBUG_MSG("SET_PROTONFREQUENCY")
        try:
            self.protonfrequency = abs(evaluate(protonfrequency))
        except:
            WRITE_STRING("\n\terror in set_proton_frequency: reset to default 400*mhz\n")
            self.protonfrequency = 400*mhz
        self.field = self.protonfrequency/gamma_H
        
    #------------------------------------------
    def set_spinningspeed(self, spinningspeed):
    #------------------------------------------
        """
        Change the spinning speed
        """
        DEBUG_MSG("SET_SPINNINGSPEED")
        self.spinningspeed = abs(evaluate(spinningspeed))
        if self.spinningspeed<1*hz and self.nsb > 0:
            WRITE_STRING("\n\tWARNING: the spinning speed being 0, nsb has also been set to 0\n")
            self.nsb = 0

    #------------------------------------------
    def set_spinningangle(self, spinningangle):
    #------------------------------------------
        """
        Change the spinning angle  
        """
        DEBUG_MSG("SET_SPINNINGANGLE")
        self.spinningangle = abs(evaluate(spinningangle))
        
    #----------------------    
    def set_nsb(self, nsb):
    #----------------------    
        """
        Change the number of spinning sideband - can be AUTO
        """
        DEBUG_MSG("SET_NSB")
        self.nsb = evaluate(nsb)
        if self.nsb == -999 and self.spinningspeed > 0.1: 
            # calculate number of spinning sidebands for satellite transition
            try:
                s = parameters.nucleus[0,1]
                wQ=3.*parameters.quadrupole[0,1]/(2.*s*(2*s-1))
                width=min(self.sw*1.1,wQ*(2*s-1.)*1.25)
                self.nsb = int(width/self.spinningspeed/2.)
            except:
                WRITE_STRING("\n\t*** ERROR *** the 'AUTO' flag in 'set_nsb' can be used only after a nucleus S has been defined")
                WRITE_STRING("\t (Additionnaly, the sw (using set_sw) and the spinningspeed (using set_spinningspeed) must ")
                WRITE_STRING("\t  also be defined prior to this calculation)\n")
        elif self.nsb == -998 and self.spinningspeed > 0.1: 
            # calculate the maximum number of spinning sidebands for satellite transition
            try:
                s = parameters.nucleus[0,1]
                wQ=3.*parameters.quadrupole[0,1]/(2.*s*(2*s-1))
                width=wQ*(2*s-1.)*1.25
                self.nsb = int(width/self.spinningspeed/2.)
                
            except:
                WRITE_STRING("error with FULL calculation")
        else:
            if self.spinningspeed<0.1:
                self.nsb = 0
            else:
                self.nsb = evaluate(nsb)
                
        self.nsb=abs(self.nsb)        

    #------------------------            
    def set_npts(self, npts):
    #------------------------
        """
        Change the number of points
        """
        DEBUG_MSG("SET_NPTS")
        self.npts = abs(evaluate(npts))
        self.npts = int(2**round(log(self.npts+1)/log(2)))
        self.npts = max(min(self.npts,parameters.nptsmax),128)

    #---------------------------------
    def set_sw(self, sw, aliasing=no):
    #---------------------------------
        """
        Change the spectral width
        """
        DEBUG_MSG("SET_SW")
        self.sw = abs(abs(evaluate(sw)))
        self.sw = max(self.sw, 100*Hz)
        self.aliasing=evaluate(aliasing)
        
    #--------------------------------
    def set_accuracy(self, accuracy):
    #--------------------------------
        """
        Change the accuracy
        """
        DEBUG_MSG("SET_ACCURACY")
        self.accuracy = abs(evaluate(accuracy))
        self.accuracy = max(min(self.accuracy,parameters.accuracymax),1)

    #----------------------------
    def set_rfstep(self, rfstep):
    #----------------------------
        """
        Change the rf integration step
        """
        DEBUG_MSG("SET_RFSTEP")
        self.rfstep = abs(evaluate(rfstep))
        self.rfstep = max(self.rfstep, .1)

    #----------------------------
    def set_detect(self, detect):
    #----------------------------
        """ 
        Set detection parameters
        args:
            detect choice: "all", "central", "satellite"
        """
        DEBUG_MSG("SET_DETECT")
        self.nall = max(min(evaluate(detect),2),0)
        try:
            spins = parameters.nucleus[0,1]
            if spins<1 : self.nall = 0
        except:
            WRITE_STRING("\n\t!!! WARNING: 'set_detect' can be used only after a nucleus S has been defined\n")
    
    #------------------------------
    def set_qfactor(self, qfactor):
    #------------------------------
        """
        Change the probe quality factor  
        """
        DEBUG_MSG("SET_QFACTOR")
        self.qfactor = max(evaluate(qfactor),0.0001)

    #-----------------------------
    def store_spectrum(self,data="undefined",label=""):
    #--------------------------------------------------
        """
        add a spectrum to the collection of spectra  
        """
        DEBUG_MSG("Store Spectrum")
        spec = array(parameters.spec)
        s = self.collection_spectra
        self.collection_spectra.append(spec)
        if data != "undefined":
            self.collection_data.append(data)      #identifier
            DEBUG_MSG("Stored : "+str(data))
        if label !="":
            self.collection_label=label
            DEBUG_MSG("Stored : "+label)
        else:
            self.collection_label="variable index"            
        self.current_spectrum = array(spec)
        self.list_lb.append(self.lb)
        self.list_gb.append(self.gb)
        self.list_concentration.append(self.concentration)
        
        WRITE_STRING("\nA spectrum has been added to collection spectra ("+str(len(self.collection_spectra))+")") 

    #------------------------------------------
    def sum_spectra(self, do_sum=True):
    #------------------------------------------
        """
        set a flag for future summation of spectra  
        """
        DEBUG_MSG("SUM_SPECTRA "+str(do_sum))
        self.flag_sum=do_sum
                    
    #-------------------------------------------
    def write_spectra(self):
    #-------------------------------------------
        """
        Write spectra into a file in various format
        args: filename - name of the output file (default:'pulsar.spe'
        """
        DEBUG_MSG("WRITE_SPECTRA ")

        # TODO : implement various format
        n = self.npts
        sw= self.sw
        vl = parameters.nucleus[0,2]
        np1 = len(self.collection_spectra)
        sw1 = 0
        vl1 = 0
        lines = []
        lines.append("PULSAR"+'\n')
        lines.append("NP="+string.strip(str(n))+'\n')
        lines.append("SW="+string.strip(str(sw))+'\n')
        lines.append("VL="+string.strip(str(vl))+'\n')
        if np1>1:
            lines.append("NI="+string.strip(str(np1))+'\n')
            lines.append("SW1="+string.strip(str(sw1))+'\n')
            lines.append("VL1="+string.strip(str(vl1))+'\n')
            lines.append("SUM="+string.strip(str(self.flag_sum))+'\n')
        for id in range(np1):
            lines.append("LB="+string.strip(str(self.list_lb[id]))+"\n")
            lines.append("GB="+string.strip(str(self.list_gb[id]))+"\n")
            lines.append("CONCENTRATION="+string.strip(str(self.list_concentration[id]))+"\n")
            
        lines.append("TYPE=SPE"+'\n')
        if len(self.collection_data)==np1:
            print self.collection_label
            lines.append("YDATA "+self.collection_label+'\n')
            for id in range(np1):
                line = string.strip(repr(self.collection_data[id]))+'\n'
                lines.append(line)
        lines.append("XDATA"+'\n')
        for id in range(np1):
            spec=[]
            spec=self.collection_spectra[id]
            
            #case of an aliasing
            if self.aliasing:
                i=0
                specs=[]
                np=parameters.npts
                for j in range(n):
                    sp=0.
                    for i in range(np/n):
                        sp=sp+spec[j+i*n]
                    specs.append(sp)
                spec=[0]*n
                spec[1:n/2-1]=specs[n/2:n-1]
                spec[n/2:n-1]=specs[1:n/2-1]
                
            #Make a correction of the spec by creating actual Abs and Disp component
            fid=[]
            try:
                fid = fft(array(spec,'D'),None, -1,0)
                for i in range(n/2,n):
                    fid[i]=complex(0.,0.)  #dissymetrisation
                   
                spec= ifft(fid,None, -1,1)

            except:
                WRITE_STRING("ERROR in fft")

            for i in range(len(spec)):
                line = string.strip(repr(spec[i].imag))+" "+string.strip(repr(spec[i].real))+'\n'
                lines.append(line)
        lines.append("END"+'\n')
        try:
            os.remove("scratch.tmp")
        except:
            pass
        f = open("scratch.tmp","wb")
        f.writelines(lines)
        f.close()
        if (verbose) :
            WRITE_STRING("\nSaving of spectra data done.")


    #-----------------------------------
    def add_nucleus(self, nucleus=None):
    #-----------------------------------
        """
        add a nucleus to the current simulation
        """
        DEBUG_MSG("ADD NUCLEUS")
        if not nucleus:
            WRITE_STRING("\nError in add_nucleus - check the argument of the function")
        # we add this nucleus to the sample
        self.sample.add_nucleus(nucleus)
        
    #-----------------------------
    def set_channel(self,*channel):
    #------------------------------
        """
        define the channels
        """
        DEBUG_MSG("SET_CHANNEL : "+channel[0])
        # the channel correspond to a Nucleus selection
        try:
            self.S_channel=channel[0]
            if len(channel)>1:
                self.I_channel=channel[1]
                DEBUG_MSG("SET_CHANNEL : "+channel[1])
        except:
            WRITE_STRING("\nError in set_channel")
            
    #-----------------------------------------------------
    def select_nucleus(self,nucleus=None,index=None):
    #-----------------------------------------------------
        """
        define the observed nucleus used for the calculation (use only its index)
        arg:
            nucleus: String
            index: integer - index start at 0
        """
        DEBUG_MSG("SELECT_NUCLEUS")

        #First we have to handle a possible errors arising when no observed channels have been defined
        try:
            essai=self.S_channel
        except:
            print  self.sample.observed
            self.S_channel=self.sample.observed
        
        if (nucleus is None) and (index is None):
            # select_nucleus()
            #   we assumed that the observing channel in this case
            #   corresponds to the first set nucleus on the observed nucleus  
            nucleus=self.S_channel
            index=0
            
        elif type(nucleus)is not types.StringType:
            # select_nucleus(1)
            #    The nucleus name was not given - probably only the index
            index=nucleus
            nucleus=self.S_channel
            
        elif type(index) is types.IntType :
            if (nucleus is None):
                nucleus=self.S_channel
                
        else:
            nucleus=self.S_channel
            index=0

        WRITE_STRING("\n\n*Nucleus '"+ nucleus +"' with index : "+str(index)+" has been selected\n" +
                         "*********************************************************************")

        try:
            if type(nucleus) is types.StringType:
                #in the case where the nucleus is not the same as the observed channel, we have to correct
                #the chemical shift to take into account the frequency difference
                delta = isotopes[nucleus]["larmor"] - isotopes[self.S_channel]["larmor"]
                delta=delta*mhz   
                abundance=isotopes[nucleus]["abundance"]
                
            #Validation
            self.sample.validate(nucleus,index,delta, abundance)

            #Read the lb, gb and concentration from the current nucleus
            self.lb=self.sample.nuclei[nucleus][index].lb
            self.gb=self.sample.nuclei[nucleus][index].gb
            self.concentration=self.sample.nuclei[nucleus][index].concentration
        
        except KeyError:
            WRITE_STRING("\n\n*** ERROR *** '"+nucleus+ \
                  "' was not found in the nucleus list - Check the name and/or the mass number")
        except:
            WRITE_STRING("\nError in select_nucleus")
            
    #-------------------------------------------------------------
    def get_numberofnuclei(self,isotope):
    #-------------------------------------------------------------
        """
        Get the number of nuclei of a given type
        Parameters
        * isotope: required isotope name 
        Usage
            ...
            get_numberofnuclei("27Al")   # determine the number of 27Al nuclei
            ...
        """
        DEBUG_MSG("GET NUMBER OF NUCLEI")
        try:
            nb= len(self.sample.nuclei[isotope])
            if (verbose) : WRITE_STRING("\nNumber of nuclei:"+str(nb))
        except:
            WRITE_STRING("\nError in get_numberofnuclei")
        return nb
    
    #------------------------
    def set_idealpulse(self):
    #------------------------
        """
        The pulse is considered as an ideal pulse
        """
        DEBUG_MSG("SET_IDEAL")
        self.sequence.idealpulse()
        
    #------------------------------------------------------------------------------------
    def set_pulse(self,length=0,powerS=0,offsetS=0,phaseS=0,powerI=0,offsetI=0,phaseI=0):
    #------------------------------------------------------------------------------------
        """
         SYNTAX: set_pulse(self,length=0,powerS=0,offsetS=0,phaseS=0,powerI=0,offsetI=0,phaseI=0)
         'Set pulse parameters'
         optional args:
            length - pulse duration (default is 0)
            powerS, powerI - rf power respectively on spin S (1) and I (2)
            offsetS, offsetI - offset respectively on spin S and I
            phaseS, phaseI - pulse phase instruction resp. on S and I
        """
        DEBUG_MSG("SET_PULSE")
        self.sequence.pulse(length,powerS,offsetS,phaseS,powerI,offsetI,phaseI)

    #-------------------------------------------
    def set_delay(self,length=0,decouple=false):
    #-------------------------------------------
        """
         SYNTAX : set_delay(self,length,decouple) 
         'Set delay parameters'
         args:
            length - delay duration
            decouple - decoupling flag (default=false)
        """
        DEBUG_MSG("SET_DELAY")
        self.sequence.delay(length,decouple)

    #-----------------------------------------
    def select_coherence(self,dp=999,nbphases=1):
    #-----------------------------------------
        """
         SYNTAX : select_coherence(self,dp=999,nbphases=1) 
         'select coherence transfer and mode of selection'
         args:
            nucl - index of the nucleus (SpinS, SpinI)
            dp - coherence transfer jump (999 : no selection)
            nbphases - number of phases (1 by default : no phase cycling)
        """
        DEBUG_MSG("SELECT_COHERENCE")
        self.sequence.coherence(dp,nbphases)
       
#==============================================================================
class Sample:
#==============================================================================
    """
    A class to handle the parameters relative to a set of nucleus for one sample
    """
    #------------------    
    def __init__(self):
    #------------------    
        """
        initialize the class with a zero number of nucleus
        """
        DEBUG_MSG("Sample INIT")
        self.nuclei = {}
        self.observed = None    #no observed nucleus defined (bydefault it will be the first one define)
        
    #---------------------------------    
    def add_nucleus(self, newnucleus):
    #---------------------------------   
        """
        add a nucleus to the sample
        """
        DEBUG_MSG("Sample ADD_NUCLEUS")
        key=newnucleus.isotope # the name of the isotope become the key
        if not self.nuclei:
            # if there is no nuclei already defined
            self.observed=key   # in case we need to know what was the first nuclei entered in this dictionary
            DEBUG_MSG("observed default"+self.observed)

        if not self.nuclei.has_key(key):
            self.nuclei[key] = []  
        self.nuclei[key].append(newnucleus)
        DEBUG_MSG(" -- add_nucleus "+key)
        
    #---------------------------------------------------    
    def validate(self,key,index=0,delta=0, abundance=1):
    #---------------------------------------------------    
        """ validation of one nucleus """
        if not self.nuclei: return 
        if not self.nuclei.has_key(key): return
        #extract the nuclei with the key and the index indicated 
        nucleus=self.nuclei[key][index]                  
        nucleus.validate(delta, abundance)
        
    #----------------------    
    def validate_all(self):
    #----------------------    
        """ validation of all nuclei in the sample """
        DEBUG_MSG("Sample VALIDATE ALL")
        if self.nuclei:

            nuclei=[]
            for key,nuclei in self.nuclei.items():
                for nucleus in nuclei:
                    nucleus.validate()
                 
#==============================================================================
class Nucleus:
#==============================================================================
    """
    A class to define the nucleus parameters
    """
    
    #---------------------------------    
    def __init__(self,isotope = "1H"):
    #---------------------------------    
        """
        set default parameters
          - default isotope = '1H'
        """
        DEBUG_MSG("Nucleus INIT")
        self.isotope = isotope
        self.nucleus=None
        try:

            self.nucleus = isotopes[isotope]      
            self.spin = self.nucleus["spin"]
            self.larmor = abs(self.nucleus["larmor"]) 
            self.name = self.nucleus["name"]
            self.quadmoment = self.nucleus["quadrupole"]
            self.abundance = self.nucleus["abundance"]
            self.relative = self.nucleus["relative"]
            self.absolute = self.nucleus["absolute"]
            
            self.set_chemicalshift()
            self.set_quadrupole()
            self.set_relaxation()
            self.lb=0.01
            self.gb=0.
            self.concentration=1.
            
        except KeyError:
            WRITE_STRING("*** ERROR *** '"+isotope+ \
                  "' was not found in the nucleus list - Check the name and/or the mass number")
        except:
            WRITE_STRING("*** ERROR *** Error in nucleus INIT")
      
    #-----------------------------------    
    def validate(self,delta, abundance):
    #-----------------------------------     
        """
        set parameters in the fortran module
        """
        DEBUG_MSG("Nucleus VALIDATE")
        if not self.nucleus:
            return

        if parameters.nucleus == None:
            # if the 'parameters.nucleus' doesn't exist
            # create the temporary 'ar' array with a  single line
            #----------------------------------------------------
            self.index = 1
            ar = [[self.index,self.spin,self.larmor*parameters.protonfrequency/100.,abundance]]   
                            #conversion made to actual larmor frequency
            
        else:  
            # the 'parameters.nucleus' array exist already
            # then copy the existing 'parameters.nucleus' array into the temporary 'ar' array
            #---------------------------------------------------------------------------------
            nofnuclei  =  len(parameters.nucleus)
            ar = [[0,0,0,0] for i in range(nofnuclei+1)]
            i = 0
            while i < nofnuclei :
                ar[i] = [i+1,parameters.nucleus[i,1],parameters.nucleus[i,2],parameters.nucleus[i,3]]  #copy of the previous elements
                i = i+1
            # and add a new line 
            #-------------------
            self.index = nofnuclei+1
            ar[nofnuclei] = [self.index,self.spin,self.larmor*parameters.protonfrequency/100.,abundance]      
                            #conversion made to actual larmor frequency
            
        # allocate or reallocate and initialize the 'parameters.nucleus' array with the 'ar' array
        #-----------------------------------------------------------------------------------------
        parameters.nucleus = ar
        
        nofnuclei  =  len(parameters.nucleus)

        # chemical shift
        #---------------
        # conversion ppm to hz (depends on the larmor frequency)
        # add the shift due to eventual difference with the radio frequency
        iso = self.iso*parameters.nucleus[self.index-1,2]/mhz + delta
        csa = self.csa*parameters.nucleus[self.index-1,2]/mhz
        if parameters.chemicalshift == None:
            # if the 'parameters.chemicalshift' array doen't exist, create it with one line 
            # (which in principe correspond to index 1 as this function is called immediately
            # after the creation of the each nuclei)
            ar = [[self.index,iso,csa,self.etacsa,self.alphacsa,self.betacsa,self.gammacsa]]
        else:
            # the size of the 'parameters.shift' has to be the same as the size of 'parameters.nucleus'
            # copy of the previous elements 
            ar = [[0,0,0,0,0,0,0] for i in range(nofnuclei)]   
            i = 0
            size = len(parameters.chemicalshift) 
            while i < size :
                ar[i] = [int(parameters.chemicalshift[i,0]),parameters.chemicalshift[i,1], parameters.chemicalshift[i,2], parameters.chemicalshift[i,3],parameters.chemicalshift[i,4],parameters.chemicalshift[i,5],parameters.chemicalshift[i,6]]  
                i = i+1
            # replace the line of the given index with the new one 
            ar[self.index-1] = [self.index,iso,csa,self.etacsa,self.alphacsa,self.betacsa,self.gammacsa]       
        parameters.chemicalshift = ar                                          #allocate or reallocate and initialize    

        if parameters.quadrupole==None:
            # if the 'parameters.quadrupole' array doesn't exist, create it with one line 
            # (which in principe correspond to index 1)        
            ar = [[self.index,self.cq,self.etaq,self.alphaq,self.betaq,self.gammaq]]   
        else:
            #copy of the previous elements
            ar = [[0,0,0,0,0,0] for i in range(nofnuclei)]
            i = 0
            size = len(parameters.quadrupole) 
            while i < size :
                ar[i] = [int(parameters.quadrupole[i,0]),parameters.quadrupole[i,1], parameters.quadrupole[i,2], parameters.quadrupole[i,3],parameters.quadrupole[i,4],parameters.quadrupole[i,5]]  
                i = i+1
            ar[self.index-1] = [self.index,self.cq,self.etaq,self.alphaq,self.betaq,self.gammaq]       
        parameters.quadrupole = ar                                          #allocate or reallocate and initialize       
        # T2
        if parameters.t2 == None:
            # if the 'parameters.t2' array doen't exist, create it with one line 
            # (which in principe correspond to index 1)
            ar = [[self.index,self.t2]]   
        else:
            # copy of the previous elements
            size = len(parameters.t2)
            ar = [[0,0] for i in range(nofnuclei)]
            i = 0
            while i < size :
                ar[i] = [int(parameters.t2[i,0]),parameters.t2[i,1]]  
                i = i+1
            # replace the line of the given index with the new one 
            ar[self.index-1] = [self.index,self.t2]                 # the numbering start at zero (this is why I use index-1)   
        parameters.t2 = ar                                          #allocate or reallocate and initialize        
              
        if verbose:
            #Display the current parameters
            #-------------------------------
            WRITE_STRING("\nNUCLEUS: "+self.name)
            WRITE_STRING(  "------------")
            if self.index == 1 : WRITE_STRING("\t-This is the observed nucleus")
            if self.index == 2 : WRITE_STRING("\t-This is the excited but non observed nucleus I \n\t-only ONE excited nucleus in this version-")  
            if self.index > 2 :  WRITE_STRING("\t-additionnal NON-excited nucleus")  
            i = self.index-1
            spin=parameters.nucleus[i,1]
            if (spin*2 == int(spin)*2):
                WRITE_STRING("\tspin: " + str(int(spin)))
            else:
                WRITE_STRING("\tspin: " + str(int(spin*2))+"/2")   
            WRITE_STRING("\tlarmor frequency: " + str(parameters.nucleus[i,2]/mhz)+" MHz")
            if self.index<3:
                WRITE_STRING("\tisotropic shift: " + str(parameters.chemicalshift[i,1]) + " Hz")
                if parameters.chemicalshift[i,2]!=0 :
                    WRITE_STRING("\tcsa anisotropy: " + str(parameters.chemicalshift[i,2])+" Hz")
                if parameters.chemicalshift[i,2]!=0 :
                    WRITE_STRING("\teta csa: " + str(parameters.chemicalshift[i,3]))
                if parameters.chemicalshift[i,2]!=0 :
                    WRITE_STRING("\teuler csa [alpha,beta,gamma]: ["+ str(parameters.chemicalshift[i,4])+","+str(parameters.chemicalshift[i,5])+","+str(parameters.chemicalshift[i,6])+"]")     
                if parameters.quadrupole[i,1]!=0 :
                    WRITE_STRING("\tQuadrupolar coupling : " + str(parameters.quadrupole[i,1]/mhz)+" MHz")
                if parameters.quadrupole[i,1]!=0 :
                    WRITE_STRING("\teta Quadrupole: " + str(parameters.chemicalshift[i,2]))
                if parameters.quadrupole[i,1]!=0 :
                    WRITE_STRING("\teuler Quadrupole [alpha,beta,gamma]: ["+ str(parameters.quadrupole[i,3])+ \
                                 ","+str(parameters.quadrupole[i,4])+","+str(parameters.quadrupole[i,5])+"]")     
                WRITE_STRING("\tT2 relaxation: "+str(parameters.t2[i,1]/sec)+" sec")
                WRITE_STRING("\tLB: "+str(self.lb)+" Hz")
                WRITE_STRING("\tGB: "+str(self.gb)+" Hz")
                WRITE_STRING("\tconcentration: "+str(self.concentration))
          
    # --------------------------------------------------------------------------------
    def set_chemicalshift(self, iso=0, csa = 0, eta = 0,alpha = 0,beta = 0,gamma = 0):
    #---------------------------------------------------------------------------------
        """
        Set the chemical shift of a given nucleus
        args:
            index - index of the nucleus (defined previously by set_nucleus)
            iso - isotropic chemical shift
            csa - chemical shift anisotropy
            eta - assymetry of the csa
            alpha, beta, gamma - euler angles describing the orientation 
                                 with respect to the molecular frame
        """
        DEBUG_MSG("Nucleus SET CHEMICAL SHIFT")
        self.iso=evaluate(iso,isppm = True)
        self.csa=evaluate(csa,isppm = True)
        self.etacsa=evaluate(eta)
        self.alphacsa=evaluate(alpha)
        self.betacsa=evaluate(beta)
        self.gammacsa=evaluate(gamma)
        
    # set_relaxation---------------------
    def set_relaxation(self,t2 = 1.0e+30):
    #------------------------------------
        """
        Set the T2 of a given nucleus
        args:
           index - index of the nucleus (defined previously by set_nucleus)
           t2 - T2 relaxation time 
        """
        DEBUG_MSG("Nucleus SET RELAXATION")
        self.t2 = evaluate(t2)
            
    # set_quadrupole ----------------------------------------------------
    def set_quadrupole(self,cq = 0,eta = 0,alpha = 0,beta = 0,gamma = 0):
    #--------------------------------------------------------------------
        """
        Set the quadrupolar parameters of a given nucleus
        args:
            index - index of the nucleus (defined previously by set_nucleus)
            cq - quadrupolar coupling constant in Hz
            eta - assymetry of the quadrupolar interaction tensor
            alpha, beta, gamma - euler angles describing the orientation 
                                 with respect to the molecular frame
        """
        DEBUG_MSG("Nucleus SET QUADRUPOLE")
        if not self.nucleus:
            return
        
        self.cq = evaluate(cq)
        self.etaq = evaluate(eta)
        self.alphaq = evaluate(alpha)
        self.betaq = evaluate(beta)
        self.gammaq = evaluate(gamma)     
        if self.spin<1: 
            #if the spin is non quadrupolar
            self.cq = 0.0  
            self.etaq = 0.0
            self.alphaq = 0.0
            self.betaq = 0.0
            self.gammaq = 0.0

    # set_lb-------------------
    def set_lb(self,lb = 0.01):
    #--------------------------
        """
        Set the LB of a given nucleus
        args:
           lb - broadening (lorenzian)
        """
        DEBUG_MSG("Nucleus SET LB")
        self.lb = evaluate(lb)

    # set_gb------------------
    def set_gb(self,gb = 0.0):
    #-------------------------
        """
        Set the GB of a given nucleus
        args:
           gb - broadening (gaussian)
        """
        DEBUG_MSG("Nucleus SET GB")
        self.gb = evaluate(gb)

    # set_concentration-------
    def set_concentration(self,concentration = 1.0):
    #-------------------------
        """
        Set the CONCENTRATION of a given nucleus
        """
        DEBUG_MSG("Nucleus SET CONCENTRATION")
        self.concentration = evaluate(concentration)
        
#==============================================================================    
class PulseSequence:
#==============================================================================    
    """ A class to handle the pulse sequence definition """

    #------------------
    def __init__(self):
    #------------------
        DEBUG_MSG("Pulse INIT")
        pass

    #--------------------
    def idealpulse(self):
    #--------------------
        DEBUG_MSG("Pulse IDEALPULSE")
        parameters.idealpulse= True
        self.addpulse()

    #------------------------
    def pulse(self,length=0,powerS=0,offsetS=0,phaseS=0,powerI=0,offsetI=0,phaseI=0):
    #------------------------
        DEBUG_MSG("pulse PULSE")
        parameters.idealpulse= False
        self.addpulse(length,powerS,offsetS,phaseS,powerI,offsetI,phaseI)

    # ----------------------------------------------------------------------------------
    def addpulse(self,length=0,powerS=0,offsetS=0,phaseS=0,powerI=0,offsetI=0,phaseI=0):
    #-----------------------------------------------------------------------------------
        DEBUG_MSG("pulse ADDPULSE")

        length=evaluate(length)
        powerS=evaluate(powerS)
        powerI=evaluate(powerI)
        offsetS=evaluate(offsetS)
        offsetI=evaluate(offsetI)
        phaseS=evaluate(phaseS)
        phaseI=evaluate(phaseI)     

        if parameters.pulse==None:
            index=1
            ar=[[index,length,powerS,offsetS,phaseS,powerI,offsetI,phaseI]]
            if (verbose) :
                WRITE_STRING("\nPULSE SEQUENCE:")
                WRITE_STRING(  "---------------")
        else:

            size=len(parameters.pulse) 
            ar=[[0,0,0,0,0,0,0,0] for i in range(size+1)]
            i=0
            while i < size :
                ar[i]=[int(parameters.pulse[i,0]), parameters.pulse[i,1], parameters.pulse[i,2], parameters.pulse[i,3], parameters.pulse[i,4], parameters.pulse[i,5], parameters.pulse[i,6], parameters.pulse[i,7]]  
                #copy of the previous elements
                i=i+1
            index=size+1
            ar[size]=[index,length,powerS,offsetS,phaseS,powerI,offsetI,phaseI]       
        parameters.pulse=ar                                          #allocate or reallocate and initialize       
        size=len(parameters.pulse)
        i=size-1
        strg="Pulse"+" P"+str(int(parameters.pulse[i,0]))+" "
        strg= strg + "["
        strg= strg + str(parameters.pulse[i,1])+" usec"  
        if parameters.pulse[i,2]!=0 : 
            strg= strg + "("+str(parameters.pulse[i,2]/khz)+" kHz," \
                       +str(parameters.pulse[i,3])+" Hz," \
                       +str(parameters.pulse[i,4])+"):S"
        if parameters.pulse[i,5]!=0 : 
            strg= strg + "("+str(parameters.pulse[i,5]/khz)+" kHz," \
                        +str(parameters.pulse[i,6])+" Hz," \
                        +str(parameters.pulse[i,7])+"):I"
        strg= strg + "]"
        if (verbose) : WRITE_STRING(strg)
        # put automatically a delay=0 and select_coherence()-> no selection by default
        self.delay()
        self.coherence() 
         
    # -----------------------------------------------
    def delay(self,length=0,decoupling=false):
    #------------------------------------------------
        DEBUG_MSG("pulse DELAY")

        if parameters.pulse==None:
            WRITE_STRING("WARNING: A sequence cannot start with a delay (no effect). It is ignored")
            return
        # the index will be that of the preceeding pulse
        nofpulses=len(parameters.pulse)

        length=evaluate(length)
        decouple=0
        if decoupling : decouple=1

        if parameters.delay==None:
            ar=[[1,length,decouple]]   
        else:
            ar=[[0,0,0] for i in range(nofpulses)]
            #copy of the previous elements        
            i=0
            size=len(parameters.delay)
            while i < size :
                ar[i]=[int(parameters.delay[i,0]), parameters.delay[i,1], parameters.delay[i,2]]  
                i=i+1
            ar[nofpulses-1]=[nofpulses,length,decouple]       
        parameters.delay=ar  #allocate or reallocate and initialize       
        if (parameters.delay[nofpulses-1,1]>0):   
            strg= "Delay "+"D"+str(int(parameters.delay[nofpulses-1,0]))
            strg= strg + "[ "+str(parameters.delay[nofpulses-1,1])+" usec ]  "
            if (parameters.delay[nofpulses-1,2]==0): 
                strg=strg + "-- without decoupling"
            else:
                strg=strg +"-- with decoupling"
            if (verbose) : WRITE_STRING(strg)
        
    # coherence ----------------------------
    def coherence(self,dp=999,nbphases=1):
    #---------------------------------------
        DEBUG_MSG("pulse COHERENCE")
        dp=evaluate(dp)
        nbphases=evaluate(nbphases)
        if parameters.pulse==None:
            WRITE_STRING("WARNING: A sequence cannot start with a coherence selection (no effect). It is ignored")
            return
        if dp==999:
            nbphases==1    #no selection at all in this case                                             
        # the index will be that of the preceeding pulse
        nofpulses=len(parameters.pulse) 
        if parameters.coher==None:
            ar=[[1,dp,nbphases]]   
        else:
            ar=[[0,0,0,0] for i in range(nofpulses)]
            #copy of the previous elements        
            i=0
            size=len(parameters.coher)
            while i < size :
                ar[i]=[int(parameters.coher[i,0]), parameters.coher[i,1], parameters.coher[i,2]]  
                i=i+1
            ar[nofpulses-1]=[nofpulses,dp,nbphases]       
        parameters.coher=ar          # allocate or reallocate and initialize       
        if dp!=999:
            if (verbose) : WRITE_STRING( "coherence selection : Number of phases N="+str(nbphases)+", Required "+chr(127)+"p="+str(dp))  
        if nbphases==1:
            if (verbose) : WRITE_STRING( "no coherence selection")

    # check coherence transfers        
    def check_ctp(self):
    #---------------
        #global nofpulses, ns, ctpall
        
        if (verbose) : WRITE_STRING( "\nLIST OF POSSIBLE COHERENCE TRANSFER PATHWAYS (CTP):")

        self.nofpulses=len(parameters.pulse)
         
        if (self.nofpulses==1): 
            WRITE_STRING( "single pulse, no CTP selection possible")
            return
             
        self.ns=2*parameters.nucleus[0,1]  #maximum order
        pend=0
        ctp=range(self.nofpulses+1)
        self.ctpall=[]
        self.searchp(1,pend,ctp)
        
        if len(self.ctpall)==0 :
            WRITE_STRING("ERROR: ")
            WRITE_STRING( "No valid coherence transfer has been found:")
            WRITE_STRING( "  Check your coherence transfer settings!")
            WRITE_STRING( "   e.g., sum("+chr(127)+"p) must be " )
            WRITE_STRING( "         -1 for the observed spin S and 0 for spin I" )
            exit()

        # store the result
        nctp=size(self.ctpall,0)
        self.ctpall.sort()
        if nctp!=pow(2*self.ns+1,self.nofpulses-1):
            if (verbose) : WRITE_STRING( "number of CTPs: " + nctp)
            self.print_ctp()
            parameters.ctp=self.ctpall
        else:
            if (verbose) : WRITE_STRING(" all allowed") 
            #TODO: check some other situation where it is not needed to
            #      make the calculation of all CTPs 

    # print selected coherence transfers        
    def print_ctp(self):
    #---------------
        #global ctpall, ctpdict
        i=1
        self.ctpall.sort()
        #create a dictionnary (to reference esaly the various coherence)
        self.ctpdict={}
        for c in self.ctpall:
            self.ctpdict[i]=c
            strg = str(int(i)) + ":  {"  
            for j in range(size(self.ctpall,1)-1):
                item=c[j]
                if abs(item)<=999:
                    strg= strg + str(int(item)).rjust(3)+" ->"
                else:
                    strg= strg + "("+str(int(item-10000)).rjust(3)+","+str(-int(item-10000)).ljust(3)+") ->"
            strg=strg + " -1 }"
            if (verbose) : WRITE_STRING(strg)
            i=i+1
          
    # compress coherence transfers----------------------------------------------------        
    def compress_ctp(self):
    #------------------------------------------------------------------------------
        #global nofpulses, ns, ctpall
        
        if (verbose) : WRITE_STRING( "\nREDUCING NUMBER OF CTP :")

        l=len(self.ctpall)
        if l==0 :
            WRITE_STRING( "ERROR: ")
            WRITE_STRING( "No valid coherence transfer has been found:")
            exit()

        #look for pathway where two symmetrical coherences are found 
        for c1 in self.ctpall:
            #print "check ",c1  
            for c2 in self.ctpall:
                if c2!=c1:
                    #print "compare with ",c2     
                    for i in range(1,size(self.ctpall,1)-1):  
                        #look each items except the first one (always 0) and the last (-1).  
                        store=c2
                        store[i]=-store[i]
                        if store==c1:
                            #this two pathway are identical except for one element which is symetric
                            #we can reduce it to only one (we add 10000 to says that two symetric element
                            #must be retained)
                            self.ctpall.remove(c2)
                            self.ctpall.remove(c1)
                            store[i]=abs(store[i])+10000
                            self.ctpall.append(store)
                            #WRITE_STRING( store )
                        break
        if len(self.ctpall)== l :
            if (verbose) : WRITE_STRING( "No reduction possible." )   
        else:
            parameters.ctp=self.ctpall
            self.print_ctp()
        
    # searchp----------------------------------------------------------------------
    def searchp(self,nc,pend,ctpi):
    #------------------------------------------------------------------------------
        """
        recursive function
        used by check_ctp
        """  
        #global nofpulses, ns, ctpall
        ctpi[nc-1]=pend
        if nc==(self.nofpulses+1):  
            if pend==-1 and nc==(self.nofpulses+1):
                c=[]
                for i in ctpi:
                    c.append(i)
                self.ctpall.append(c)
            return
        else:  
            dp=parameters.coher[nc-1,1]
            nbphases=parameters.coher[nc-1,2]
            if dp==999 :   #all jumps allowed
                dp=1 
                nbphases=1
            for n in range(-int(3*self.ns),int(3*self.ns)):
                p=pend+dp+nbphases*n
                if p<=self.ns and p>=-self.ns :
                    self.searchp(nc+1, p,ctpi)
             
    # remove_ctp-------------------------------------------------------------------
    def remove_ctp(self,index):
    #------------------------------------------------------------------------------
        """ remove_ctp : to remove one or more ctp in the list
            arg: index is a scalar
                 or a list in the format (index1, index2...)"""
        #global ctpdict, ctpall
        if (verbose) : WRITE_STRING( "\nREMOVING CTP"+str(index))
        # here the dictionnary is useful
        #(because we cannot know easily the order of the list element)
        try:
            for i in index:    
                del self.ctpdict[i]
        except:
            try:
                del self.ctpdict[index]  #if list is a single number     
            except:
                WRITE_STRING( "ERROR in the remove_ctp command")
                WRITE_STRING( self.remove_ctp.__doc__)
        self.ctpall=ctpdict.values()
        parameters.ctp=self.ctpall
        self.print_ctp()

    # set ctp----------------------------------------------------------------------
    def set_ctp(self,list):
    #------------------------------------------------------------------------------
        """ create one CTP(s) from a list of value """
        #global nofpulses, ns, ctpall
        if (verbose) : WRITE_STRING( "\nADDING A COHERENCE TRANSFER PATHWAY: ")
        if len(list)!= self.nofpulses+1:
            WRITE_STRING( "ERROR: len of the list is not correct.")
            WRITE_STRING( "Must be number of pulse+1")
            exit()
        if list[0]!=0 and list[self.nofpulses]!=-1:   
            WRITE_STRING( "ERROR: List must start with 0 end ending with -1")
            exit()
           
        self.ctpall.append(list)
        parameters.ctp=self.ctpall
        self.print_ctp()


# DIPOLAR COUPLINGS
#====================
      
# dipole_from_distance-------------------
def dipole_from_distance(distance,f1,f2):
#----------------------------------------
    """
    'Compute dipolar coupling from distances' 
    args:
        distance - distances between nucleus
        f1, f2 - larmor frequency of the two nuclei 
    """
    DEBUG_MSG("dipole from distance")
    field = parameters.spectrometerfield
    h = 6.6260693e-34
    mu0 = 4*pi*1.0e-7
    f1 = f1/field  #conversion en gamma
    f2 = f2/field
    dipolar = abs(f1*f2*mu0*h/(distance*1.0e-10)**3/4/pi)
    return dipolar
                 

# indirect j coupling --------------------------------------------------------------
def set_indirect(index,second,j = 0,delta = 0,eta = 0,alpha = 0,beta = 0,gamma = 0):
#-----------------------------------------------------------------------------------
    """
    set_indirect -- 
    'Set the indirect J coupling between a pair of nucleus'
    args:
        index - index of the first nucleus (defined previously by set_nucleus)
                Required to be one.
        second - index of the second nucleus (defined previously by set_nucleus)
        j - J coupling value
        delta - ansotropy of J
        alpha, beta, gamma - euler angles with respect to the molecular frame
    """
    DEBUG_MSG("Set Indirect")
    if parameters.nucleus == None: 
        WRITE_STRING("\n*set_indirect*\nERROR: No nucleus were defined yet. Cannot set indirect j coupling'!\n")
        WRITE_STRING(set_indirect.__doc__)
        exit()
    else:
        if index != 1: 
            WRITE_STRING("\n*set_indirect*\nWARNING: the index of the first nucleus is required to be 1(for the moment)!")
            WRITE_STRING("It has been automatically changed.")
            index = 1
        if second == 1:   
            WRITE_STRING("\n*set_indirect*\nWARNING: the index of the second nucleus is required to be different of the first index!")
            WRITE_STRING("It has been automatically changed to 'index+1'.")
            second = index + 1
        if second > len(parameters.nucleus): 
            WRITE_STRING("\n*set_indirect*\nERROR: the second index does not correspond to an existing nucleus!\n")
            WRITE_STRING(set_indirect.__doc__)
            exit()
    j = evaluate(j)
    delta = evaluate(delta)
    eta = evaluate(eta)
    alpha = evaluate(alpha)
    beta = evaluate(beta) 
    gamma = evaluate(gamma)  
    if parameters.indirect==None:
        if (verbose) : WRITE_STRING("\nindirect J coupling:")
        ar=[[index,second,j,delta,eta,alpha,beta,gamma]]   
    else:
        size=len(parameters.indirect) 
        ar=[[0,0,0,0,0,0,0,0] for i in range(size+1)]
        i=0
        while i < size :
            ar[i]=[int(parameters.indirect[i,0]),int(parameters.indirect[i,1]),parameters.indirect[i,2], parameters.indirect[i,3], parameters.indirect[i,4],parameters.indirect[i,5],parameters.indirect[i,6],parameters.indirect[i,7]]  #copy of the previous elements
            i=i+1
        ar[size]=[index,second,j,delta,eta,alpha,beta,gamma]       
    parameters.indirect=ar                                          #allocate or reallocate and initialize       
    size=len(parameters.indirect)
    i=size-1
    if (verbose) :
        WRITE_STRING("\tIndirect J coupling ("+str(index)+","+str(second)+"): "+str(parameters.indirect[i,2])+" Hz")
        if parameters.indirect[i,3]!=0 : WRITE_STRING("\tdelta J ("+str(index)+","+str(second)+"): "+str(parameters.indirect[i,3])+" Hz")
        if parameters.indirect[i,3]!=0 : WRITE_STRING("\teta J ("+str(index)+","+str(second)+"): "+str(parameters.indirect[i,4]))
        if parameters.indirect[i,3]!=0 : WRITE_STRING("\teuler J ("+str(index)+","+str(second)+"): ["+ str(parameters.indirect[i,5])+","+str(parameters.indirect[i,6])+","+str(parameters.indirect[i,7])+"]")     

# direct dipolar coupling -----------------------------------------
def set_dipole(index,second,dip=0,distance=0,alpha=0,beta=0,gamma=0):
#------------------------------------------------------------------
    """
    -- set_dipole -- 
    'Set the direct dipole coupling between a pair of nucleus'
    args:
        index - index of the first nucleus (defined previously by set_nucleus)
                Required to be one.
        second - index of the second nucleus (defined previously by set_nucleus)
        dip - dipolar coupling value
        distances - distance between nucleus in angstrom (used only if d=0)
        alpha, beta, gamma - euler angles with respect to the molecular frame
    """
    DEBUG_MSG("Set Dipole")
    if parameters.nucleus==None: 
        WRITE_STRING("\n*set_dipole*\nERROR: No nucleus were defined yet. Cannot set dipole coupling'!\n")
        WRITE_STRING(set_dipole.__doc__)
        exit()
    else:
        if index!=1: 
            WRITE_STRING("\n*set_dipole*\nWARNING: the index of the first nucleus is required to be 1(for the moment)!")
            WRITE_STRING("It has been automatically changed.")
            index=1
        if second==1:   
            WRITE_STRING("\n*set_dipole*\nWARNING: the index of the second nucleus is required to be different of the first index!")
            WRITE_STRING("It has been automatically changed to 'index+1'.")
            second=index+1
        if second>len(parameters.nucleus): 
            WRITE_STRING("\n*set_dipole*\nERROR: the second index does not correspond to an existing nucleus!\n")
            WRITE_STRING(set_dipole.__doc__)
            exit()
    dip=evaluate(dip)
    distance=evaluate(distance)
    if (distance!=0) and (dip==0) : dip=dipole_from_distance(distance,parameters.nucleus[index-1,2],parameters.nucleus[second-1,2])
    alpha=evaluate(alpha)
    beta=evaluate(beta) 
    gamma=evaluate(gamma)
    if parameters.dipole==None:
        if (verbose) : WRITE_STRING("\ndipolar coupling:")
        ar=[[index,second,dip,distance,alpha,beta,gamma]]   
    else:
        size=len(parameters.dipole) 
        ar=[[0,0,0,0,0,0,0] for i in range(size+1)]
        #copy of the previous elements
        i=0
        while i < size :
            ar[i]=[int(parameters.dipole[i,0]),int(parameters.dipole[i,1]),parameters.dipole[i,2], parameters.dipole[i,3], parameters.dipole[i,4],parameters.dipole[i,5],parameters.dipole[i,6]] 
            i=i+1
        ar[size]=[index,second,dip,distance,alpha,beta,gamma]       
    parameters.dipole=ar                                          #allocate or reallocate and initialize       
     
    size=len(parameters.dipole)
    i=size-1
    if (verbose) : WRITE_STRING("\tdipolar coupling ("+str(index)+","+str(second)+"): "+str(parameters.dipole[i,2])+" Hz")
    if (verbose) : WRITE_STRING("\tpolar angles [alpha,beta] ("+str(index)+","+str(second)+"): ["+ str(parameters.dipole[i,4])+","+str(parameters.dipole[i,5])+","+str(parameters.dipole[i,6])+"]")     

#  MODULE FUNCTIONS
#==================

def evaluate(st,isppm=False):
#----------------------------
    """
    Allow the evaluation of strings containing units e.g., '10 kHz'
    Return False if there is an error
    """
    res=False
    if type(st) is types.StringType:
        st=string.upper(st)
        p = re.compile('\s+')
        m = p.split(string.strip(st))
        try:
            res = eval(m[0])
        except NameError:
            WRITE_STRING("*** ERROR *** during the evaluation of '"+st+"'  >>   "+m[0]+" is not defined")
            
        except SyntaxError:
            WRITE_STRING("*** ERROR *** Syntax error in '"+st+"'")   
        if len(m) > 1:
            if isppm :
##                print st
##                print m
##                print string.strip(m[1])
                if string.strip(m[1]) != "PPM":
                    WRITE_STRING("ERROR of unit - 'ppm' required")     
            else:
                if string.strip(m[1]) != "PPM":
                    try:
                        res=res*eval(m[1])
                    except NameError:
                        WRITE_STRING("*** ERROR *** during the evaluation of '"+st+"'  >>  "+m[1]+" is not defined")                       
                    except SyntaxError:
                        WRITE_STRING("*** ERROR *** Syntax error in '"+st+"'")        
                else:
                    WRITE_STRING("ERROR of unit - 'ppm' cannot be used here")
                    exit()
    else:
        res = st
    return res      

# strclean---------
def strclean(st,n):
#------------------    
    return str(st)[0:n]
   
# showmatrice-----
def showmatrix(d):
#-----------------
    nr = len(d)
    nc = len(d[0])
    ii = 2
    strg=''
    while ii < nr:
        jj = 2
        while jj < nc:
            strg=strg + " "+str((int(d[ii,jj]*100.))/100.)
            jj = jj + 1
        str=str+"\n"  
        ii = ii + 1
    WRITE_STRING(strg)    

#------------------------------------------------------------------------------
def nowhitespace(st):
#------------------------------------------------------------------------------
    """suppress whitespace in a string"""
    #TODO Can be nicer with regular expression
    st=list(str(st))
    ch=''
    for x in st: ch=ch+string.strip(str(x))
    return string.strip(str(ch))

#------------------------------------------------------------------------------
def TestMe():
#------------------------------------------------------------------------------
    if (verbose) : WRITE_STRING("Here we are!")
        
#------------------------------------------------------------------------------
def HelpMe():
#------------------------------------------------------------------------------
    """ help for the python commands in pulsar"""

    WRITE_STRING( set_nucleus.__doc__)
    WRITE_STRING( set_chemical_shift.__doc__)
    WRITE_STRING( set_quadrupole.__doc__)
    WRITE_STRING( set_dipole.__doc__)
    WRITE_STRING( set_indirect.__doc__)
    WRITE_STRING( set_relaxation.__doc__)
    WRITE_STRING( set_spectrometer.__doc__)
    WRITE_STRING( set_probehead.__doc__)
    WRITE_STRING( init_simulation.__doc__)

#----------------
def Delete_Log():
#----------------
    """Delete the pyPulsar.log file"""
    os.remove("pyPulsar.log")
    




