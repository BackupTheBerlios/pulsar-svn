!     Last change:  CF   15 Feb 2006    1:02 pm
! File: $Id: pulsar.for, v 0.1
! ----------------------------------------------------------------------
! PULSAR Project
!        Copyright (C) 2006 Jean-Paul Amoureux, Christian Fernandez
!        JPA - Unité de Catalyse et Chimie du Solide, Lille, France.
!        CF  - Laboratoire Catalyse et Spectrochimie, Caen, France.
!
! ----------------------------------------------------------------------
! LICENSE
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License (GPL)
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! To read the license please visit http://www.gnu.org/copyleft/gpl.html
! ----------------------------------------------------------------------
! Original Author of file: Jean-Paul Amoureux
! Purpose of file:
!    This is the main program file
! ----------------------------------------------------------------------
! LAST MODIFICATIONS:
!   20060214 CF - Creation of this file from the original Pulsar program
!   20060215 CF - Splitting of the file into MODULES
!
! TODO:
!   . Transform array in allocatable arrays
!   . Write the following header in english
!   . Transform in free form fortran syntax
!
! ----------------------------------------------------------------------
! List of subroutines in the various files
!
! LL.f95
!-------
!   LLcp
!   LLredor
!   LL3spinsRedor
!   LLtedor
!   LL
!
!PULSAR.f95
!----------
!   PULSAR
!   Main
!
!SUBROUTINES.f95
!---------------
!   TRI
!   DIPOLAR
!   SELE
!   PRINT
!   PREPARATION
!   AMcalc
!   AM
!   ORIENTATION
!   AMdif
!   EDGEA
!   EDGEB
!   SPECTREA
!   SPECTREB
!   DIAG
!   HTRIDI
!   TQL2
!   HTRIBK
!   WIGNER
!   FFT
!   FFT1
!   FAST
!   FR4TR
!   FORD1
!   FORD2
!   WIDTH
!   TENT
!   OPERATOR
!   ZERO
!----------------------------------------------------------------------------

PROGRAM  PULSAR
!==============
CHARACTER*256 cl

        CALL getcl(cl)                ! Get command line (input filename)
        OPEN(10,file=cl,status="old") ! Open this file for input
        CALL Main(cl)                 ! Call the main program
        CLOSE(10)                     ! Close the input file

END PROGRAM


SUBROUTINE Main()
!================
! Main subroutine
!================

      USE common_module

      DIMENSION NISPs(10),TL(61,30),Inorm1(16),Inorm2(16)
      DIMENSION FIDX(32800),FIDY(32800),SPIX(32800),SPIY(32800)           
      DIMENSION P(99999),Idec(9999),PHASEi(9999),Iref(9999)
      DIMENSION RorTs(11,11),RoiTs(11,11),RorTi(11,11),RoiTi(11,11)      
      DIMENSION JJ(20),Nphase(20),FASE(20),Level(0:20)         
      DATA Inorm1/2154752,933774,589820,430020,338069,278417,236611,205703, &
                  181927,163068,147749,135059,124378,115254,107382,100517/
      DATA Inorm2/6122925,3121593,2088467,1568265,1255355,1046463,897112, &
                  785115,697918,628170,571061,523524,483263,448748,418821,392639/

!      OPEN(20,file="res.txt",status="replace")
!      OPEN(22,FILE="spectre.txt",STATUS="replace")

      READ(10,1000)  ((TL(I,J),J=1,30),P(I),I=1,60)        
      READ(10,1001) IPRINT,Nphasing,Nboucle,Ifasing,NCYCL  
      READ(10,1002) (Nphase(I),I=1,20)
      READ(10,1007) (Level(J),J=0,20)    
      READ(10,1000)  (TL(61,J),J=1,30),P(61)
      READ(10,1008) (NCfoo,(P(51+J+10*NC),J=1,3),Iref(NC),IQ1(NC,1),IQ2(NC,1),(P(51+J+10*NC),J=5,9),Idec(NC),NC=1,NCYCL)

!     -----------------------------------------------------------
      IF(IPRINT.EQ.1) THEN
!      OPEN(12,file="debug.txt",status="replace")
      END IF

      ISPEED=0
      PI=3.1415926536
      PI2=0.01745329252
      SQ6=SQRT(6.)
!              DONNEES SPECTROMETRE
      NT=5*INT(P(2)+0.0001)
      powder1=FLOAT(Inorm1(NT/5))/10000000000.
      powder2=FLOAT(Inorm2(NT/5))/10000000000.
      NU=NT
      NPTS=INT(P(4)+0.0001)
      FST=-P(5)/2.
      DLB=P(6)
      DGB=P(7)
      WLs=P(12)*1000000.
      IF(P(42).LT.0.4.OR.P(43).LT.0.4) THEN
                p(42)=0.00001 
                P(43)=0.00001
                P(44)=100000000.
                DO 20 I=45,59
   20           P(I)=0.     
!            on élimine le noyau I qui n'a alors pas d'interaction avec celui S
      END IF
      WLi=P(42)*1000000.
 
      IF(P(8).LE.0.0001.OR.P(9).LE.0.0001) THEN                             
             ISPEED=1
             NU=1              
             P(8)=0.
             P(9)=0.
             P(39)=0.
      END IF
      DU=2.*PI/FLOAT(NU)              

      CE=COS(P(8)*PI2)
      SE=SIN(P(8)*PI2)
      DO 30 I=-4,4
      DO 30 J=-4,4
   30 DE(I,J)=0.                                                 
      CALL WIGNER(CE,SE,DE)

      DO 40 KS=1,2
      DO 40 L=0,4
      DO 40 I=L-2,2
   40 DK(KS,L,I)=DE(I,KS)*DE(L-I,-KS)

      WROT=P(9)
      SPINs=P(13)
      SPINi=P(43)
      IF(P(3).LT.0.1) P(3)=0.1
      IF(P(3).GT.360.)P(3)=360.
      IF(ISPEED.EQ.1) RFstep=100000000.
      IF(ISPEED.EQ.0) RFstep=P(3)/WROT/0.00036
      QQ=4.*P(10)*P(10)/WLs/WLs               
!              DONNEES NOYAU ET LECTURE
      Ns=INT(2*SPINs+0.0001)+1
      Ni=INT(2*SPINi+0.0001)+1
      Nsi=Ns*Ni
      IF(MOD(SPINs,1.).LT.0.001) P(40)=1.
      NALL=INT(P(40)+0.0001)
              
      CALL TRI(P,Idec)
!     -------------------------------définition de tous les opérateurs----------------------------- 
      CALL OPERATOR
!     --------------calcul des phases des pulses permettant la sélection d'un chemin de cohérence--
      CALL SELE(P,Level,Iref,Nphase,Ifasing,PHASEi,IDNfase,DNscan)
!*************************BOUCLE SUR LE TEMPS T1 LORS D'UNE EXPERIENCE A 2 DIMENSIONS *************
      DO 995 IUY=10,10,1
      Pfit=0.5*FLOAT(IUY)

      CALL PRINT(TL,P,Nphase,Level,IPRINT,Ifasing,IDNfase)       
      T2s=P(14)*0.000002*PI
      T2i=P(44)*0.000002*PI
      
      CALL PREPARATION(T2s,T2i,P,PI,PI2,RFstep)  
!****************BOUCLE POSSIBLE SUR LES ESPECES S, ICI UNE SEULE ESPECE***************************
      DO 700 KLM=1,1
      CALL ZERO()   

      DJ(4)=P(45)
      Wdip(4)=P(46) 
      A2Ds=-P(47)*PI2
      A3Ds=-P(48)*PI2
      A23Ds(4)=ABS(A2Ds)+ABS(A3Ds)
      IF(A23Ds(4).GE.0.01) CALL ORIENTATION(0.,A2Ds,A3Ds,0.,4,FFDs,GGDs)
 
      COEFPs=0.
      COEFPi=0.
      CQs=P(15)
      CQi=P(49)
      IF(SPINs.GT.0.6) COEFPs=-500.*SQ6*CQs/(4.*SPINs*(2.*SPINs-1.))
      IF(SPINi.GT.0.6) COEFPi=-500.*SQ6*CQi/(4.*SPINi*(2.*SPINi-1.))
      COEFSs=6.*COEFPs**2/WLs
      COEFSi=6.*COEFPi**2/WLi
      ETAQs=P(16)
      ETAQi=P(50)
      DO 310 K=1,Ns-1
      DMs=SPINs+1.-FLOAT(K)
      U1s(K)=-0.50*COEFSs*(24.*DMs*(DMs-1.)-4.*SPINs*(SPINs+1.)+9.)
  310 U2s(K)=-0.25*COEFSs*(12.*DMs*(DMs-1.)-4.*SPINs*(SPINs+1.)+6.)
      A1qi=P(51)*PI2
      A2qi=P(52)*PI2 
      A3qi=P(53)*PI2
      A123qi=ABS(A1qi)+ABS(A2qi)+ABS(A3qi)
      IF(A123qi.GE.0.01) CALL ORIENTATION(A1qi,A2qi,A3qi,ETAQi,1,FFQi,GGQi)
               
      DISOs=P(17)*WLs/1000000.
      ETACs=P(19)
      DELTACs=P(18)*WLs/(3.+ETACs)/1000000.
      A1cs=P(20)*PI2
      A2cs=P(21)*PI2
      A3cs=P(22)*PI2
      A123cs=ABS(A1cs)+ABS(A2cs)+ABS(A3cs)
      IF(A123cs.GE.0.01) CALL ORIENTATION(A1cs,A2cs,A3cs,ETACs,1, FFCs,GGCs)
                     
      DISOi=P(54)*WLi/1000000.
      ETACi=P(56)
      DELTACi=P(55)*WLi/(3.+ETACi)/1000000.
      A1ci=P(57)*PI2
      A2ci=P(58)*PI2
      A3ci=P(59)*PI2
      A123ci=ABS(A1ci)+ABS(A2ci)+ABS(A3ci)
      IF(A123ci.GE.0.01) CALL ORIENTATION(A1ci,A2ci,A3ci,ETACi,1,FFCi,GGCi)
      NSB=INT(P(39))
!-------------------BOUCLE SUR LE DIPOLAIRE NON MANIPULE--(K-S)--------------
      CALL DIPOLAR(3,IN,P,NISPs,COEFdip)
               DO 600 Ii1=-NISPs(1),NISPs(1),2                           
               WdipIs(1)=Wdip(1)*FLOAT(Ii1)/2.                           
               DO 600 Ii2=-NISPs(2),NISPs(2),2                           
               WdipIs(2)=Wdip(2)*FLOAT(Ii2)/2.                           
               DO 600 Ii3=-NISPs(3),NISPs(3),2                           
               WdipIs(3)=Wdip(3)*FLOAT(Ii3)/2.                             
      DISOst=DISOs+(DJ(1)*FLOAT(Ii1)+DJ(2)*FLOAT(Ii2)+DJ(3)*FLOAT(Ii3))/2.
      A123T=A123cs*ABS(DELTACs)+A123ci*ABS(DELTACi)+A23Ds(1)*ABS(Ii1)+ &
            A23Ds(2)*ABS(Ii2)+A23Ds(3)*ABS(Ii3)+A23Ds(4)+A123qi*1000.*ABS(CQi)
!-------------------BOUCLE SUR LES PHASES DES PULSES, SI IL Y A CYCLAGE DE PHASE
      DO 410 IIZ=0,40
      DO 410 JJZ=0,80
      DO 410 KZ=1,11
      DO 410 NSBs=-20,20
      QRT(IIZ,JJZ,KZ,NSBs)=0.
  410 QIT(IIZ,JJZ,KZ,NSBs)=0.
  
      DO 500 J1=0,Nphase(1)-1
      DO 500 J2=0,Nphase(2)-1
      DO 500 J3=0,Nphase(3)-1
      DO 500 J4=0,Nphase(4)-1
      DO 500 J5=0,Nphase(5)-1
      DO 500 J6=0,Nphase(6)-1
      DO 500 J7=0,Nphase(7)-1
      DO 500 J8=0,Nphase(8)-1 
      DO 500 J9=0,Nphase(9)-1
      DO 500 J10=0,Nphase(10)-1
      DO 500 J11=0,Nphase(11)-1
      DO 500 J12=0,Nphase(12)-1
      DO 500 J13=0,Nphase(13)-1 
      DO 500 J14=0,Nphase(14)-1
      DO 500 J15=0,Nphase(15)-1
      DO 500 J16=0,Nphase(16)-1
      DO 500 J17=0,Nphase(17)-1
      DO 500 J18=0,Nphase(18)-1
      DO 500 J19=0,Nphase(19)-1
      DO 500 J20=0,Nphase(20)-1 
      JJ(1)=J1
      JJ(2)=J2
      JJ(3)=J3
      JJ(4)=J4
      JJ(5)=J5
      JJ(6)=J6
      JJ(7)=J7 
      JJ(8)=J8
      JJ(9)=J9
      JJ(10)=J10
      JJ(11)=J11
      JJ(12)=J12
      JJ(13)=J13  
      JJ(14)=J14
      JJ(15)=J15
      JJ(16)=J16
      JJ(17)=J17
      JJ(18)=J18
      JJ(19)=J19
      JJ(20)=J20  
      PHIR=0.
      IF(Ifasing.NE.0) THEN
           DO 450 If=1,IDNfase
           FASE(If)=360.*FLOAT(JJ(If))/FLOAT(Nphase(If))
  450      PHIR=PHIR+FASE(If)*(Level(If-1)-Level(If)) 
!                  on calcule la phase de chaque paquet de pulses ainsi que celle du receiver
           DO 460 NC=1,NCYCL
           PPP=FASE(Iref(NC))-PHIR+PHASEi(NC)       
  460      PULSEs(NC,3)=AMOD(PPP,360.)*PI2 
!                  on calcule la phase de chaque pulse shiftée de celle du receiver (0 dans PULSAR) 
      END IF

      IF(IPRINT.EQ.1) THEN
       WRITE(6,1000)  (TL(61,J),J=1,30),P(61)
       WRITE(6,1008) (NC,(P(51+J+10*NC),J=1,2),PULSEs(NC,3)/PI2,Iref(NC),  &
                      IQ1(NC,1),IQ2(NC,1),(P(51+J+10*NC),J=5,9),Idec(NC),NC=1,NCYCL)
      END IF

      IF(A123T.LT.0.1) THEN
             DNORM=FLOAT(NT)*COEFdip/FLOAT(NU)**2/FLOAT(Ns)/DNscan  
             CALL EDGEA()            
      ELSE        
             DNORM=FLOAT(NT)*COEFdip/FLOAT(NU)**2/FLOAT(Ns)/DNscan/4.
             CALL EDGEB()
      END IF     

  500 CONTINUE

      IF(A123T.LT.0.1) THEN
             FNEXP=FLOAT(NU)*DNscan*powder1/COEFdip             
             CALL SPECTREA()
      ELSE        
             FNEXP=FLOAT(NU)*DNscan*powder2/COEFdip
             CALL SPECTREB()
      END IF     

  600 CONTINUE
  700 CONTINUE
!----------------------------------------------------------------------      
      DO 720 I=1,Ns
      DO 720 J=1,Ns
      RorTs(I,J)=0.
      RoiTs(I,J)=0.
      DO 720 KK=0,Ni-1
      RorTs(I,J)=RorTs(I,J)+RorTt(I+KK*Ns,J+KK*Ns)/FLOAT(Ni)/FNEXP
  720 RoiTs(I,J)=RoiTs(I,J)+RoiTt(I+KK*Ns,J+KK*Ns)/FLOAT(Ni)/FNEXP
           
      DO 730 I=1,Ni
      DO 730 J=1,Ni
      RorTi(I,J)=0. 
      RoiTi(I,J)=0.
      DO 730 KK=1,Ns
      RorTi(I,J)=RorTi(I,J)+RorTt(KK+(I-1)*Ns,KK+(J-1)*Ns)*WLs/WLi/FLOAT(Ns)/FNEXP
  730 RoiTi(I,J)=RoiTi(I,J)+RoiTt(KK+(I-1)*Ns,KK+(J-1)*Ns)*WLs/WLi/FLOAT(Ns)/FNEXP

                                                                  
!----------------------!
!     Broadenings      !
!----------------------!

      CALL WIDTH(NPTS,DLB,DGB,FST,SPX,FIDX,FIDY)
      CALL FFT(2*NPTS,FIDX,SPX,SPIX)
      CALL WIDTH(NPTS,DLB,DGB,FST,SPY,FIDX,FIDY)
      CALL FFT(2*NPTS,FIDX,SPY,SPIY)

      DO 900 IT=1,NPTS
      SPY(IT)=SPY(IT)-SPIX(IT)
      SPX(IT)=SPX(IT)+SPIY(IT)
      SPIY(IT)=0.
  900 SPIY(2*NPTS+1-IT)=0.

      NN2=NPTS/2
      DO 910 I=1,NN2
      SPIY(I)=SPY(NN2+I)
  910 SPIY(2*NPTS+1-I)=SPY(NN2+1-I)

      CALL FFT(2*NPTS,SPIY,FIDY,FIDX)

!-----------------!
!     OUTPUT      !
!-----------------!
! Writing spectra
      OPEN(2,file="test.out",status="replace")
      WRITE (2,*) "<spectre>"
      DO 920 I=1,NPTS
  920 WRITE (2,*) FST+(I-1)*(-2*FST/NPTS),1000000*SPX(I),1000000*SPY(I)
      WRITE (2,*) "</spectre>"

! writing density matrices
      WRITE (2,*) "<MDs>"
      WRITE(2,*) SPINs
      DO 960 I=1,Ns
  960 WRITE(2,5001) (RorTs(I,J),RoiTs(I,J),J=1,Ns)
      WRITE (2,*) "</MDs>"
                 
      IF(SPINi.LT.0.1) GOTO 990
      WRITE (2,*) "<MDi>"
      WRITE(2,*) SPINi
      DO 970 I=1,Ni                               
  970 WRITE(2,5001) (RorTi(I,J),RoiTi(I,J),J=1,Ni)
      WRITE (2,*) "</MDi>"

      WRITE (2,*) "<MDis>"
      WRITE(2,*) SPINs,SPINi
      DO 980 I=1,Nsi                               
  980 WRITE(2,5001) (RorTt(I,J)/FNEXP,RoiTt(I,J)/FNEXP,J=1,Nsi)
      WRITE (2,*) "</MDis>"
  990 CONTINUE 
  995 CONTINUE 

 1000 FORMAT(30A2,F13.3)
 1001 FORMAT(13X,I1,15X,I1,11X,I2,16X,I1,9X,I4)
 1002 FORMAT(13X,20(1X,I2)) 
 1007 FORMAT(10X,21(1X,I2))
 1008 FORMAT(I4,F9.0,F7.1,F8.3,3I3,F9.0,F7.1,2F8.3,F10.3,I2)
 5001 FORMAT(100(2X,F7.1,1X,F7.1))
 5002 FORMAT (F16.8,F16.8,F16.8)

      close(2)


      close(5)
!      close(12)
!      close(20)
!      close(22)
!      STOP
       END
!






