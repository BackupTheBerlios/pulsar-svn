# -*- coding: cp1252 -*-
# os.system("pulsar.exe test.in")

import os
import string
import sys
import getopt

# Function spawn : Execute the external program in a subprocess and wait
# ----------------------------------------------------------------------
# (Hints from : http://effbot.org/librarybook/os.htm)
# ----------------------------------------------------------------------
def spawn(program, *args):
    if os.name in ("nt", "dos"):
        exefile = ".exe"
    else:
        exefile = ""
    try:
        # check if the os module provides a shortcut
        return os.spawnvp(program, (program,) + args)
    except AttributeError:
        pass
    try:
        spawnv = os.spawnv
    except AttributeError:
        # assume it's unix
        pid = os.fork()
        if not pid:
            os.execvp(program, (program,) + args)
        return os.wait()[0]
    else:
        # got spawnv but no spawnp: go look for an executable
            file = os.path.join(program) + exefile
            print file
            try:
                return spawnv(os.P_WAIT, file, (file,) + args)
            except os.error:
                pass
# end Function spawn
# -----------------------------------------------------------------------

#-------- Main ----------------------------------------------------------
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
             raise Usage(msg)

        # Build the fortran input file
        #-----------------------------

        #init units
        hz=1
        khz=1000
        mhz=1000000
        us=1
        ms=1000
        s=1000000

        #init parameters

        #GENERAL
        #-------
        VERBOSE=1      # 1 VERBOSE OUTPUT                                 
        ACCURACY=5     # 2 ACCURACY ON EULER ANGLES (1-8)                                 
        RF=5           # 3 RF INTEGRATION STEP : VAS ( )                                  
        TD=512         # 4 SPECTRUM POINTS NUMBER  (2**n<16385)                         
        SW=100*khz     # 5 TOTAL SPECTRAL WIDTH (SW:Default Hz)                               
        LB=10*hz       # 6 LORENTZIAN  LINEWIDTH BROADENING (F.W.H.M) (Hz)               
        GB=0           # 7 GAUSSIAN    LINEWIDTH BROADENING (F.W.H.M) (Hz)                
        ANGLE=54.736   # 8 SPINNER  ANGLE / Bo ( )                                       
        SPEED=0        # 9 SPINNER FREQUENCY (Hz)                                     
        QF=0.001       #10 PROBE QUALITY Q-FACTOR (F.W.H.M=Wo/Q)                          
           
        #S OBSERVED NUCLEUS
        #------------------
        VLS=400*mhz    #12 S LARMOR FREQENCY
        SPINS=0.5      #13 S SPIN VALUE <6                                                
        T2S=900000000  #14 T2 temps de relaxation homogène (micro-secondes)       
        CQS=0          #15 QUADRUPOLAR CONSTANT : Cq=eeqQ/h (KHz)1, 3,   12            
        ETAQS=0        #16 QUADRUPOLAR ASYMMETRY  (FROM 0 TO 1)                           
        ISOS=0         #17 CHEMICAL SHIFT (P.P.M) % MIDDLE OF THE S WINDOW 19.883       
        CSAS=0         #18 CSA: TOTAL WIDTH OF THE SPECTRUM(P.P.M)                        
        ETACS=0        #19 CSA: ASYMMETRY  (FROM 0 TO 1)                                  
        ACS=0          #20                  alpha ( ) / quadrupolar PAS                   
        BCS=0          #21                  beta  ( ) / quadrupolar PAS                   
        GCS=0          #22                  gamma ( ) / quadrupolar PAS 
                      
        #NUCLEUS HETEROGENEOUSLY COUPLED TO S (NON EXCITED)
        #--------------------------------------------------
        NDIS=0         #23 NUMBER OF NUCLEI DIPOLARLY-COUPLED TO THIS NUCLEUS <4          
        I1=0           #24 1 SPIN VALUE OF THIS DIPOLARLY-COUPLED NUCLEUS                 
        J1=0           #25 1 J COUPLING (Hz)                                              
        D1=0           #26 1 INDIVIDUAL INHOMOGENEOUS DIPOLAR INTERACTION (Hz)            
        B1=0           #27                   1 POLAR ANGLE (°) / quadrupolar OZ           
        A1=0           #28                   1 POLAR ANGLE (°) / quadrupolar OXZ          
        I2=0           #29 2 SPIN VALUE OF THIS DIPOLARLY-COUPLED NUCLEUS                 
        J2=0           #30 2 J COUPLING                                                   
        D2=0           #31 2 INDIVIDUAL INHOMOGENEOUS DIPOLAR INTERACTION (Hz)            
        B2=0           #32                   2 POLAR ANGLE (°) / quadrupolar OZ           
        A2=0           #33                   2 POLAR ANGLE (°) / quadrupolar OXZ          
        I3=0           #34 3 SPIN VALUE OF THIS DIPOLARLY-COUPLED NUCLEUS                 
        J3=0           #35 3 J COUPLING (Hz)                                              
        D3=0           #36 3 INDIVIDUAL INHOMOGENEOUS DIPOLAR INTERACTION (Hz)            
        B3=0           #37                   3 POLAR ANGLE (°) / quadrupolar OZ           
        A3=0           #38                   3 POLAR ANGLE (°) / quadrupolar OXZ          
        NSB=0          #39 HALF-NUMBER OF SIDEBANDS <21                                   
        NALL=0         #40 ALL TRANSITIONS :1  / ELSE :0 transition (-1/2,1/2)            
       
        # I EXCITED NUCLEUS,NOT OBSERVED
        # ------------------------------
        VLI=0          #42 LARMOR FREQUENCY (MHz)                                       
        SPINI=0        #43 SPIN VALUE OF THIS DIPOLARLY-COUPLED NUCLEUS <6                
        T2I=900000000  #44 T2 temps de relaxation homogène (micro-secondes)       
        JCoupling=0    #45 J COUPLING (Hz)                                               
        Dipole=0       #46 I-S INDIVIDUAL INHOMOGENEOUS DIPOLAR INTERACTION (Hz)        
        Beta=0         #47        I-S POLAR ANGLE (°) / S quadrupolar OZ                  
        Alpha=0        #48        I-S POLAR ANGLE (°) / S quadrupolar OXZ                 
        CQI=0          #49 QUADRUPOLAR CONSTANT : Cq=eeqQ/h (KHz)                         
        ETAQI=0        #50 QUADRUPOLAR ASYMMETRY  (FROM 0 TO 1)                           
        AQI=0          #51                  alpha ( ) / S quadrupolar PAS                 
        BQI=0          #52                  beta  ( ) / S quadrupolar PAS                 
        GQI=0          #53                  gamma ( ) / S quadrupolar PAS                 
        ISOI=0         #54 CHEMICAL SHIFT (P.P.M) % MIDDLE OF THE WINDOW                  
        CSAI=0         #55 CSA: TOTAL WIDTH OF THE SPECTRUM(P.P.M)                      
        ETACI=0        #56 CSA: ASYMMETRY  (FROM 0 TO 1)                                  
        ACI=0          #57                  alpha ( ) / S quadrupolar PAS                 
        BCI=0          #58                  beta  ( ) / S quadrupolar PAS                 
        GCI=0          #59                  gamma ( ) / S quadrupolar PAS                 

        # COHERENCE.SELECTION
        #--------------------
        NPHASING=0     #60
        NBOUCLE=0      #61
        IFASING=1      #62
        NPHASES= (    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
        LEVELS = (0,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)

        # SEQUENCE
        #-----------
        NCYCL=2        #63
        PULSE={}
        PULSE[1]= (   5500.,    0.0,   0.000,  1,  1, 10,       0.,    0.0,   0.000,  15.152,  4377.273, 0)
        PULSE[2]= (   5500.,    0.0,   0.000,  2, -1, 10,       0.,    0.0,   0.000,  30.303,  4384.849, 0)

        #read the input file and execute commands
        f = open("test.in", "r")             #Open the file
        exec(f.read())                       #Execute
        f.close                              #Close
        
        #write parameters to the 'pulsar.dat' file
        f = open("pulsar.dat", "w")     #open the file for writing
        f.write('VERBOSE '+str(VERBOSE)+'\n')      # 1 
        f.write('ACCURACY '+str(ACCURACY)+'\n')     # 2                            
        f.write('RF '+str(RF)+'\n')           # 3                                   
        f.write('TD '+str(TD)+'\n')           # 4                          
        f.write('SW '+str(SW)+'\n')           # 5                                
        f.write('LB '+str(LB)+'\n')           # 6
        f.write('GB '+str(GB)+'\n')           # 7                 
        f.write('ANGLE '+str(ANGLE)+'\n')     # 8                                        
        f.write('SPEED '+str(SPEED)+'\n')     # 9                                      
        f.write('QF '+str(QF)+'\n')           #10                           
        f.write('-  0\n')
        f.write('VLS '+str(VLS)+'\n')          #12 
        f.write('SPINS '+str(SPINS)+'\n')        #13                                                 
        f.write('T2S '+str(T2S)+'\n')          #14        
        f.write('CQS '+str(CQS)+'\n')          #15             
        f.write('ETAQS '+str(ETAQS)+'\n')        #16                            
        f.write('ISOS '+str(ISOS)+'\n')         #17        
        f.write('CSAS '+str(CSAS)+'\n')         #18                         
        f.write('ETACS '+str(ETACS)+'\n')        #19                                   
        f.write('ACS '+str(ACS)+'\n')          #20                                     
        f.write('BCS '+str(BCS)+'\n')          #21                                     
        f.write('GCS '+str(GCS)+'\n')          #22                   
        f.write('NDIS '+str(NDIS)+'\n')         #23           
        f.write('I1 '+str(I1)+'\n')           #24                  
        f.write('J1 '+str(J1)+'\n')           #25                                              
        f.write('D1 '+str(D1)+'\n')           #26             
        f.write('B1 '+str(B1)+'\n')           #27                              
        f.write('A1 '+str(A1)+'\n')           #28                             
        f.write('I2 '+str(I2)+'\n')           #29                  
        f.write('J2 '+str(J2)+'\n')           #30                                                   
        f.write('D2 '+str(D2)+'\n')           #31             
        f.write('B2 '+str(B2)+'\n')           #32                             
        f.write('A2 '+str(A2)+'\n')           #33                            
        f.write('I3 '+str(I3)+'\n')           #34                 
        f.write('J3 '+str(J3)+'\n')           #35                                               
        f.write('D3 '+str(D3)+'\n')           #36             
        f.write('B3 '+str(B3)+'\n')           #37                             
        f.write('A3 '+str(A3)+'\n')           #38                             
        f.write('NSB '+str(NSB)+'\n')          #39                                  
        f.write('NALL '+str(NALL)+'\n')         #40             
        f.write('- '+'0\n')                        
        f.write('VLI '+str(VLI)+'\n')          #42                                      
        f.write('SPINI '+str(SPINI)+'\n')        #43                
        f.write('T2I '+str(T2I)+'\n')          #44        
        f.write('JCoupling '+str(JCoupling)+'\n')    #45                                                
        f.write('Dipole '+str(Dipole)+'\n')       #46         
        f.write('Beta '+str(Beta)+'\n')         #47                         
        f.write('Alpha '+str(Alpha)+'\n')        #48                         
        f.write('CQI '+str(CQI)+'\n')          #49                          
        f.write('ETAQI '+str(ETAQI)+'\n')        #50                           
        f.write('AQI '+str(AQI)+'\n')          #51                                  
        f.write('BQI '+str(BQI)+'\n')          #52                                  
        f.write('GQI '+str(GQI)+'\n')          #53                                   
        f.write('ISOI '+str(ISOI)+'\n')         #54                   
        f.write('CSAI '+str(CSAI)+'\n')         #55                       
        f.write('ETACI '+str(ETACI)+'\n')        #56                                  
        f.write('ACI '+str(ACI)+'\n')          #57                                   
        f.write('BCI '+str(BCI)+'\n')          #58                                   
        f.write('GCI '+str(GCI)+'\n')          #59                                  
        f.write('NPHASING '+str(NPHASING)+'\n')     
        f.write('NBOUCLE '+str(NBOUCLE)+'\n')      
        f.write('IFASING '+str(IFASING)+'\n')
        for i in range(len(NPHASES)):                 
            f.write(str(NPHASES[i])+' ')
        f.write('\n')   
        for i in range(len(LEVELS)):                 
            f.write(str(LEVELS[i])+' ')
        f.write('\n')
        f.write('NCYCL '+str(NCYCL)+'\n')
        for j in range(NCYCL):
            for i in range(len(PULSE[j+1])):                 
                f.write(str(PULSE[j+1][i])+' ')
            f.write('\n')

        f.close()

        # spawn the external Pulsar fortan program
        # ----------------------------------------
        if spawn("pulsar", "pulsar.dat")!= 0:
            print 'ERROR'


        
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":

    argv = sys.argv
    main()
