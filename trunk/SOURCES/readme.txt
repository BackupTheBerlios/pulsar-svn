PULSAR VERSION 0.1
ORIGINAL PROGRAM OF J.P.AMOUREUX
------------------------------------------------------------------------
Modified by J.P. Amoureux and C. Fernandez for a free release under 
GPL licence.
------------------------------------------------------------------------

CALCUL DE LA MATRICE DENSITE TOTALE (S*I) ET DETERMINATION DU SIGNAL SUIVANT S APRES UNE SEQUENCE DE PREPARATION (pulses et attentes) POUR 2 SPINS ISOLES I ET S. CECI EN     
ABSENCE DE TOUT MOUVEMENT DE REORIENTATION MOLECULAIRE.     
PULSES CAN BE SENT IN THE I AND S CHANNELS.                     
HOWEVER,EACH OF THEM DOES NOT INTERACT ON THE OTHER SPIN.  
L'OPERATEUR D'EVOLUTION EST INTEGRE STEP PAR STEP PENDANT LES PULSES: METHODE LEAP-FROG.                              
LORSQU'UN PULSE DURE PLUS QU'UNE PERIODE DU ROTOR (spin-lock ou CP), L'OPERATEUR D'EVOLUTION EST CALCULE SUR UNE SEULE ROTATION.                                     
                                                            
TROIS TYPES D'UTILISATIONS SONT POSSIBLES:         
Elles peuvent etre utilisées séparément ou simultanément.   
Dans les cas (2,3) un seul spectre est calculé               
      (1) le cyclage de phase réellement utilisé sur le spectro 
          --> tous les spectres correspondant sont calculés.    
      (2) celui correspondant simplement à des pulses,          
              sans sélection dans la matrice densité.          
      (3) celui correspondant aux niveaux de quanta sélectionnés
              directement dans la matrice densité. 
      Lorsqu'un seul niveau est sélectionné, il doit l'etre 
                              dans la 1° colonne.                  
                                                                    
LA SEQUENCE DE PREPARATION EST ANALYSEE TEMPORELLEMENT ET CELLE D'OBSERVATION FREQUENTIELLEMENT. En effet, on suppose que les procesus irréversibles (e.g. dipolaire homo) sont négligeables; ce qui signifie qu'il n'y a pas de transfer entre cohérences. Pendant l'observation, ces processus ne sont pris en compte que par un élargissement Lorentzo/Gaussien.                                          
                                                                     
POWDER SAMPLE IS : STATIC (TETA=0) OR SPINNING : V.A.S       
                                                                    
INTERACTIONS (with possible different P.A.Ss) ARE :         
           SPIN S : CSA, QUADRUPOLAR, J and INHOMOGENEOUS-DIPOLAR   
                    with other spins excited (I->S) or not (K->S)  
                                        I spin is perfectly decoupled during acquisition
                    I spin can be perfectly decoupled or not         
                               during delays: 1 or 0 in last column.     
           SPIN I : CSA and QUADRUPOLAR INTERACTIONS                 
                    S->I : J and INHOMOGENEOUS DIPOLAR.             
                    I spin is perfectly decoupled during acquisition
                    I spin can be perfectly decoupled or not        
                               during delays: 1 or 0 in last column.

SPINi  : valeurs des                        les tableaux permettent 
SPINs  : spins etudies                      des spins de valeur <5  
                                                                     
NCYCL  = Nombre de cycles<10000. One cycle = un pulse + une attente.
         Pour chaque valeur de ALPHA, BETA et GAMMA,le programme    
         effectue une boucle sur la valeur de NCYCL.                
                                                                     
ISPEED = 0 : VAS ;  = 1 : STATIC                                     
                                                                     
ITOUR   * PULSES LONGS EN VAS : NOMBRE -1 DE TOURS ENTIERS=0,1,2,3..
        * STATIC OU PULSES COURTS EN VAS                  =-1       
                                                                    
NTP     * PULSES LONGS EN VAS : NOMBRE DE PAS (DTP),POUR OBTENIR UN  
                         TOUR EXACTEMENT + LE "PETIT" DERNIER(DTP1) 
        * STATIC OU PULSES COURTS EN VAS:NOMBRE DE PAS CORRESPONDANT
                                               A LA LONGEUR DU PULSE 
                                                                     
DTP     * PULSES LONGS EN VAS : PAS EN mmS (*2.PI.E-6) POUR OBTENIR 
                                UN NOMBRE ENTIER DE PAS SUR UN TOUR 
        * STATIC OU PULSES COURTS EN VAS : PAS EN mmS (*2.PI.E-6)   
                                                                     
DTP1    * PULSES LONGS EN VAS : le dernier tour est décomposé en    
NTPint    NTPint pas de durée DTP, plus un dernier petit pas de duré
          DTP1 (<DTP) en mms(*2.PI.E-6). Il peut arriver que DTP1=0.   
        * STATIC OU PULSES COURTS EN VAS : DTP1=NTPint=0                

---------------------------------------------------------------------
             LES FREQUENCES SONT EN Hz                              
              LES TEMPS SONT MULTIPLIES PAR 2*pi                     
              LES ANGLES SONT EN RADIANS                             
*********************************************************************

