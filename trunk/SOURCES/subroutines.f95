!     Last change:  CF   16 Feb 2006    8:49 am
! File: $Id: subroutines.f95, v 0.1
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
!Last modifications:
!  CF   16 Feb 2006    8:46 am :  some cleaning (not finish of course)
!------------------------------------------------------------------


!***********************************************************************
SUBROUTINE TRI(P,Idec)
!***********************************************************************
      USE common_module
      DIMENSION P(*),Idec(*)
      IF(SPINs.LE.0.6) P(15)=0.
      IF(ABS(P(15)).LE.0.001) THEN
           P(16)=0.
           P(20)=0.
           P(21)=0.
           P(22)=0.
      END IF
      IF(ABS(P(18)).LE.0.001) THEN
           P(19)=0.
           P(20)=0.
           P(21)=0.
           P(22)=0.
      END IF  
      IF((ABS(P(16))+ABS(P(19))+ABS(P(21))).LT.0.001) THEN
           P(20)=0.
           P(22)=0.
      END IF
      Nmin=5*INT(P(23))+1 
      DO 10 INN=Nmin,15
   10 P(23+INN)=0.      
      DO 20 INs=1,3
   20 IF(ABS(P(22+5*INs)).LT.0.001) P(23+5*INs)=0.     
                
      IF(P(43).LE.0.6) P(49)=0.
      IF(ABS(P(49)).LE.0.001) THEN
           P(50)=0.
           P(51)=0. 
           P(52)=0.
           P(53)=0.
      END IF         
      IF(ABS(P(55)).LT.0.001) THEN
           P(56)=0.
           P(57)=0.
           P(58)=0.
           P(59)=0.
      END IF

      DO 30 NC=1,NCYCL
      IF(Idec(NC).EQ.1) Decouple(NC)=0.
   30 IF(Idec(NC).EQ.0) Decouple(NC)=1.

      RETURN
      END
!

!***********************************************************************
      SUBROUTINE DIPOLAR(NATsmax,IN,P,NISPs,COEFdip)
!***********************************************************************       
      USE common_module
      DIMENSION P(*),NISPs(*)
      COEFdip=1.
      NATsd=INT(P(23))
      DO 10 INs=1,NATsmax
      NISPs(INs)=INT(2.*P(19+5*INs)+0.00001)
      COEFdip=COEFdip/(NISPs(INs)+1.)
      DJ(INs)=P(20+5*INs)
      Wdip(INs)=P(21+5*INs)
      A2Ds=-P(22+5*INs)*0.01745329252
      A3Ds=-P(23+5*INs)*0.01745329252
      A23Ds(INs)=ABS(A2Ds)+ABS(A3Ds)
   10 IF(A23Ds(INs).GE.0.01) CALL ORIENTATION(0.,A2Ds,A3Ds,0.,INs,FFDs,GGDs)
      RETURN                                                                   
      END                                                                      
!
!***********************************************************************
      SUBROUTINE SELE(P,Level,Iref,Nphase,Ifasing,PHASEi,IDNfase,DNscan)
!***********************************************************************       
      USE common_module
      DIMENSION P(*),Level(0:20),Iref(*),Nphase(*),PHASEi(*)
      Level(0)=0
      IF(Nboucle.GT.NCYCL.OR.Ifasing.EQ.1) Nboucle=0
      IF(Nphasing.EQ.0) THEN
             Nboucle=0
             DO 50 I=1,9999
             DO 50 J=1,2
                   IQ1(I,J)=10
   50              IQ2(I,J)=10               
      ELSE              
             IQ2t=0
             DO 60 I=1,9999
                   IQ1(I,2)=10
                   IQ2(I,2)=10
                   IF(I.GE.NCYCL) IQ1(I,1)=10   
                   IF(I.GE.NCYCL) IQ2(I,1)=10
                   IF(I.GT.NCYCL) P(59+10*I)=0.
   60              IF(IABS(IQ2(I,1)).LT.Ns) IQ2t=IQ2t+1
             IF(IQ2t.EQ.0) Nboucle=0
      END IF

!     cyclage de phase: cohérence des données et calcul du nombre total de scans expérimentaux 
!     (DNscan) et du nombre de paquets de pulses qui sont phasés (IDNfase<NCYCL).
      DNscan=1.
      IDNfase=0
      DO 70 I=1,20
      IF(Ifasing.EQ.0.OR.I.GE.NCYCL.OR.Nphase(I).LT.1) Nphase(I)=1
      IF(I.GE.NCYCL) Level(I)=-1
      IF(Nphase(I).GT.1) IDNfase=I
   70 DNscan=DNscan*Nphase(I)  

      IF(Nboucle.GT.0) THEN
      DO 80 I=1,Nboucle
      IQ1(I,2)=IQ2(I,1)   
   80 IQ2(I,1)=10  
      END IF
!     Jusqu'à Nboucle ->il ne reste que IQ1(I,1) (1° colonne) et IQ1(I,2) (2° colonne)
!                     ->on fait un 1° passage avec la 1° colonne, puis un 2° avec la 2° colonne

!     PHASEi(I) sert comme incrément relatif lorsque l'on phase un paquet de pulses (FAM:0 ou PI)
      DO 90 I=1,9999
      IF(Ifasing.NE.0) PHASEi(I)=P(54+10*I)
   90 IF(Ifasing.EQ.0) Iref(I)=0
      RETURN                                                                   
      END 
!
!***********************************************************************
      SUBROUTINE PRINT(TL,P,Nphase,Level,IPRINT,Ifasing,IDNfase)
!***********************************************************************       
      USE common_module
      DIMENSION TL(61,30),P(*),Nphase(20),Level(0:20)
      IF(IPRINT.EQ.1) THEN
                        WRITE(6,1000) ((TL(I,J),J=1,30),P(I),I=1,2)
      IF(ISPEED.EQ.0)   WRITE(6,1000)  (TL(3,J),J=1,30),P(3)
                        WRITE(6,1000) ((TL(I,J),J=1,30),P(I),I=4,7)
      IF(ISPEED.EQ.0)   WRITE(6,1000) ((TL(I,J),J=1,30),P(I),I=8,9)
                        WRITE(6,1000) ((TL(I,J),J=1,30),P(I),I=10,17)
      IF(P(18).NE.0.0)  WRITE(6,1000) ((TL(I,J),J=1,30),P(I),I=18,22)
      Nec=INT(P(23))                      
      IF(Nec.GT.0) WRITE(6,1000) ((TL(I,J),J=1,30),P(I),I=23,23+5*Nec)  
      IF(ISPEED.EQ.0)   WRITE(6,1000)  (TL(39,J),J=1,30),P(39)
                        WRITE(6,1000)  (TL(40,J),J=1,30),P(40)
      IF(SPINi.GT.0.1)  WRITE(6,1000) ((TL(I,J),J=1,30),P(I),I=41,59)  
                        WRITE(6,1000)  (TL(60,J),J=1,30),P(60) 
      WRITE(6,1001) IPRINT,Nphasing,Nboucle,Ifasing,NCYCL
      IF(Ifasing.EQ.1) WRITE(6,1002) (Nphase(I),I=1,IDNfase+1)
      IF(Ifasing.EQ.1) WRITE(6,1003)  (Level(I),I=0,IDNfase+1)
      END IF 
 1000 FORMAT(30A2,F13.3)
 1001 FORMAT(13X,I1,15X,I1,11X,I2,16X,I1,9X,I4)
 1002 FORMAT(13X,20(1X,I2)) 
 1003 FORMAT(10X,21(1X,I2))
      RETURN                                                                   
      END 
!
!***********************************************************************
      SUBROUTINE PREPARATION(T2s,T2i,P,PI,PI2,RFstep)
!***********************************************************************       
      USE common_module
      DIMENSION P(*)
!           DONNEES CYCLES DE PREPARATION : PULSES ET ATTENTES
      DO 100 NC=1,NCYCL
      PULSEs(NC,1)=P(52+10*NC)
      PULSEs(NC,2)=P(53+10*NC)
      PULSEs(NC,3)=P(54+10*NC)*PI2    
      PULSEi(NC,1)=P(56+10*NC)
      PULSEi(NC,2)=P(57+10*NC)
      PULSEi(NC,3)=P(58+10*NC)*PI2
      PULSEs(NC,4)=P(59+10*NC)*0.000002*PI  
      NTP(NC)= INT(P(59+10*NC)/RFstep)+1 
      NNTP=NTP(NC)
      DTP(NC)=P(59+10*NC)*0.000002*PI/FLOAT(NNTP)
      DTP1(NC)=0.
      NTPint(NC)=0
!         Lorsque le pulse est <1 tour, on découpe chaque pulse en NTP petits pulses de durée DTP
!         Le 'reste' (DTP1 et NTPint)est annulé
      ITOUR(NC)=INT(P(59+10*NC)*WROT/1000000.)-1
      IF(ITOUR(NC).GE.0) THEN
          NTP(NC)=INT(1000000./WROT/RFstep)+1
          NNTP=NTP(NC)
          DTP(NC)=2.*PI/WROT/(FLOAT(NNTP)-1.)
          NTPint(NC)=INT(2.*PI*P(59+10*NC)/DTP(NC)/1000000.)-(ITOUR(NC)+1)*(NTP(NC)-1)
          DTP1(NC)=2.*PI*P(59+10*NC)/1000000.-DTP(NC)*((ITOUR(NC)+1)*(NTP(NC)-1)+NTPint(NC))
!         Lorsque le pulse est >1 tour, on le découpe en ITOUR+1 nombre entier de tours plus un 'reste'.
!         Chaque tour entier est découpé exactement en NTP pulses de durée DTP, calculée pour.
!         Le reste du dernier tour est découpé en NTPint pulses de meme durée DTP, plus un dernier
!                                                               pulse de durée inférieure: DTP1<DTP
      END IF
      DELAY(NC)=P(60+10*NC)*0.000002*PI
      T2ss(NC)=EXP(-(PULSEs(NC,4)+DELAY(NC))/T2s)
      T2ii(NC)=EXP(-(PULSEs(NC,4)+DELAY(NC))/T2i)
  100 T2si(NC)=(T2ss(NC)+T2ii(NC))/2.
      RETURN                                                                   
      END 
!
!***********************************************************************
      SUBROUTINE AMcalc(CA,SA,C2A,S2A,ZZ,SB,ARts,AIts,ARQi,AIQi,ARCi,AICi,ARQs,AIQs,ARCs,AICs,ARDs,AIDs,PR,PI,QR,QI)
!***********************************************************************       
      USE common_module

      real, parameter :: sq6=2.4494897427831780981972840747059

      DIMENSION ARQs(10,-4:4),AIQs(10,-4:4),ARCs(10,-4:4),AICs(10,-4:4) &
     ,ARDs(10,-4:4),AIDs(10,-4:4),ARts(-4:4),AIts(-4:4)                &
     ,ARQi(10,-4:4),AIQi(10,-4:4),ARCi(10,-4:4),AICi(10,-4:4)          &
     ,PR(11,-20:20),PI(11,-20:20),QR(11,-20:20),QI(11,-20:20)

      DO 10 I=-4,4
      ARts(I)=0.
      AIts(I)=0.
      DO 10 J=1,10
      ARQi(J,I)=0.
      AIQi(J,I)=0.
      ARCi(J,I)=0.
      AICi(J,I)=0.
      ARQs(J,I)=0.
      AIQs(J,I)=0.
      ARCs(J,I)=0.
      AICs(J,I)=0.
      ARDs(J,I)=0.
   10 AIDs(J,I)=0.


      DO 20 K=1,Ns-1
      DO 20 NSBs=-NSB,NSB
      PR(K,NSBs)=0.
      PI(K,NSBs)=0.
      QR(K,NSBs)=0.
   20 QI(K,NSBs)=0.

      IF(ABS(COEFPs).GE.0.1) CALL AM(C2A,S2A,ZZ,SB,ETAQs,1,ARQs,AIQs)
      IF(ABS(DELTACs).GE.0.1) THEN
          IF(A123cs.LT.0.01) CALL AM(C2A,S2A,ZZ,SB,ETACs,1,ARCs,AICs)
          IF(A123cs.GE.0.01) CALL AMdif(CA,SA,C2A,S2A,ZZ,SB,1,FFCs,GGCs,ARCs,AICs)
      END IF

      DO 30 L=0,2
      ARts(L)=SQ6*DELTACs*ARCs(1,L)
      AIts(L)=SQ6*DELTACs*AICs(1,L)
      DO 30 INs=1,NATsd
          IF(ABS(WdipIs(INs)).GE.0.1) THEN
          IF(A23Ds(INs).LT.0.01) CALL AM(C2A,S2A,ZZ,SB,0.,INs,ARDs,AIDs)
          IF(A23Ds(INs).GE.0.01) CALL AMdif(CA,SA,C2A,S2A,ZZ,SB,INs,FFDs,GGDs,ARDs,AIDs)
          END IF
      ARts(L)=ARts(L)+SQ6*WdipIs(INs)*ARDs(INs,L)
   30 AIts(L)=AIts(L)+SQ6*WdipIs(INs)*AIDs(INs,L)

        IF(ABS(Wdip(4)).GE.0.1) THEN
          IF(A23Ds(4).LT.0.01) CALL AM(C2A,S2A,ZZ,SB,0.,4,ARDs,AIDs)     
          IF(A23Ds(4).GE.0.01) CALL AMdif(CA,SA,C2A,S2A,ZZ,SB,4,FFDs,GGDs,ARDs,AIDs)
        END IF
       
      IF(ABS(COEFPi).GE.0.1) THEN                                              
         IF(A123qi.LT.0.01) CALL AM(C2A,S2A,ZZ,SB,ETAQi,1,ARQi,AIQi)
         IF(A123qi.GE.0.01) CALL AMdif(CA,SA,C2A,S2A,ZZ,SB,1,FFQi,GGQi,ARQi,AIQi)
      END IF

      IF(ABS(DELTACi).GE.0.1) THEN
         IF(A123ci.LT.0.01) CALL AM(C2A,S2A,ZZ,SB,ETACi,1,ARCi,AICi)
         IF(A123ci.GE.0.01) CALL AMdif(CA,SA,C2A,S2A,ZZ,SB,1,FFCi,GGCi,ARCi,AICi)
      END IF

      RETURN                                                                   
      END
!***********************************************************************
      SUBROUTINE AM(C2A,S2A,CB,SB,ETA,N,AR,AI)
!***********************************************************************
      REAL*4 AR(10,-4:4),AI(10,-4:4)
      AR(N,-2)=SB*SB/2.-ETA*(1.+CB*CB)*C2A/6.
      AI(N,-2)=-ETA*CB*S2A/3.
      AR(N,-1)=-SB*CB*(ETA*C2A/3.+1.) 
      AI(N,-1)=-ETA*SB*S2A/3.
      AR(N,0)=(3.*CB*CB-1.-ETA*SB*SB*C2A)/SQRT(6.)
      AR(N,2)=AR(N,-2)
      AR(N,1)=-AR(N,-1)
      AI(N,2)=-AI(N,-2)
      AI(N,1)=AI(N,-1)
      AI(N,0)=0.
      RETURN
      END
! 
!***********************************************************************
      SUBROUTINE ORIENTATION(A1,A2,A3,ETA,N,FF,GG)
!***********************************************************************
      REAL*4 FF(10,0:2),GG(10,0:2),D(-4:4,-4:4)
      V0=0.8164965809
      V2=-ETA/3.
      CALL WIGNER(COS(A2),SIN(A2),D)
      DO 10 M=1,2
      FF(N,M)=V0*D(0,M)*COS(M*A3)+V2*(D(2,M)*COS(M*A3+2.*A1)+D(-2,M)*COS(M*A3-2.*A1))
   10 GG(N,M)=V0*D(0,M)*SIN(M*A3)+V2*(D(2,M)*SIN(M*A3+2.*A1)+D(-2,M)*SIN(M*A3-2.*A1))
      FF(N,0)=V0*D(0,0)+2.*V2*D(2,0)*COS(2.*A1)
      GG(N,0)=0.
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE AMdif(CA,SA,C2A,S2A,C,S,N,FF,GG,AR,AI)
!***********************************************************************
      REAL*4 AR(10,-4:4),AI(10,-4:4),FF(10,0:2),GG(10,0:2)
      SQ15=SQRT(1.5)
      FGM1=CA*FF(N,1)-SA*GG(N,1)
      FGM2=C2A*FF(N,2)-S2A*GG(N,2)
      FGP1=CA*GG(N,1)+SA*FF(N,1)
      FGP2=C2A*GG(N,2)+S2A*FF(N,2)
      AR(N,0)=(3.*C*C-1.)*FF(N,0)/2.+SQ15*S*(S*FGM2-2.*C*FGM1)
      AR(N,1)=S*C*(SQ15*FF(N,0)-FGM2)+(2.*C*C-1.)*FGM1
      AR(N,2)=(SQ15*S*S*FF(N,0)+(1.+C*C)*FGM2)/2.+S*C*FGM1
      AI(N,1)=S*FGP2-C*FGP1
      AI(N,2)=-C*FGP2-S*FGP1
      AI(N,0)=0.
      AR(N,-2)=AR(N,2)
      AR(N,-1)=-AR(N,1)
      AI(N,-2)=-AI(N,2)
      AI(N,-1)=AI(N,1)
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE EDGEA()
!***********************************************************************
      USE common_module
      DO 10 JJ=0,NT
      DO 10 II=0,NT-JJ
      R=SQRT(FLOAT((NT-II-JJ)**2+II*II+JJ*JJ))
      R3=R*R*R
      XI=FLOAT(NT-II-JJ)/R
      YI=FLOAT(II)/R
      ZI=FLOAT(JJ)/R
   10 CALL LL(XI,YI,ZI,II,JJ,R3)
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE EDGEB()
!***********************************************************************
      USE common_module
      DO 40 JJ=0,NT-1
      DO 20 II=0,NT-JJ
      R=SQRT(FLOAT((NT-II-JJ)**2+II*II+JJ*JJ))
      R3=R*R*R
      XI=FLOAT(NT-II-JJ)/R
      YI=FLOAT(II)/R
      ZI=FLOAT(JJ)/R
   20 CALL LL(XI,YI,ZI,II,JJ,R3)
      DO 30 II=NT-JJ+1,NT
      R=SQRT(FLOAT((NT-II-JJ)**2+(NT-JJ)**2+(NT-II)**2))
      R3=R*R*R
      XI=FLOAT(NT-II-JJ)/R
      YI=FLOAT(NT-JJ)/R
      ZI=FLOAT(NT-II)/R
   30 CALL LL(XI,YI,ZI,II,JJ,R3)
   40 CONTINUE
      DO 70 JJ=NT,2*NT-1                    
      DO 50 II=JJ-NT+1,NT-1
      R=SQRT(FLOAT((JJ-NT-II)**2+(NT-JJ)**2+(NT-II)**2))
      R3=R*R*R
      XI=FLOAT(JJ-NT-II)/R
      YI=FLOAT(NT-JJ)/R
      ZI=FLOAT(NT-II)/R
   50 CALL LL(XI,YI,ZI,II,JJ,R3)
      DO 60 II=1,JJ-NT           
      R=SQRT(FLOAT((JJ-NT-II)**2+II*II+(2*NT-JJ)**2))
      R3=R*R*R
      XI=FLOAT(JJ-NT-II)/R
      YI=-FLOAT(II)/R
      ZI=FLOAT(2*NT-JJ)/R
   60 CALL LL(XI,YI,ZI,II,JJ,R3)
   70 CONTINUE

      DO 110 K=1,Ns-1
      DO 80 JJ=0,NT-1
      WR(0,2*NT-JJ,K)=WR(0,JJ,K)
      DO 80 NSBs=-NSB,NSB                            
      AMPX(0,2*NT-JJ,K,NSBs)=AMPX(0,JJ,K,NSBs)
   80 AMPY(0,2*NT-JJ,K,NSBs)=AMPY(0,JJ,K,NSBs)
      DO 90 II=0,NT
      WR(NT,NT+II,K)=WR(II,0,K)
      DO 90 NSBs=-NSB,NSB                    
      AMPX(NT,NT+II,K,NSBs)=AMPX(II,0,K,NSBs)
   90 AMPY(NT,NT+II,K,NSBs)=AMPY(II,0,K,NSBs)
      DO 100 JJ=1,NT-1
      WR(NT-JJ,2*NT,K)=WR(NT,JJ,K)
      DO 100 NSBS=-NSB,NSB                      
      AMPX(NT-JJ,2*NT,K,NSBs)=AMPX(NT,JJ,K,NSBs)
  100 AMPY(NT-JJ,2*NT,K,NSBs)=AMPY(NT,JJ,K,NSBs)
  110 CONTINUE

      R3=FLOAT(NT*NT*NT)
      CALL LL(0.,0.,1.,0,NT,R3)      

      RETURN
      END
!
!***********************************************************************
      SUBROUTINE SPECTREA()
!***********************************************************************
      USE common_module
      DO 31 K=1,Ns-1
      IF(NALL.EQ.0.AND.(2*K).NE.INT(2*SPINs+1)) GOTO 31
      DO 30 NSBs=-NSB,NSB
      DECAL=NSBs*WROT
      DO 10 I=0,NT-2
      DO 10 J=0,NT-2-I
      AMPT=AMPX(I+1,J,K,NSBs)+AMPX(I,J+1,K,NSBs)+AMPX(I,J,K,NSBs)
      CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I,J,K),AMPT,FST,DECAL,NPTS+1,SPX)
      AMPT=AMPY(I+1,J,K,NSBs)+AMPY(I,J+1,K,NSBs)+AMPY(I,J,K,NSBs)
      CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I,J,K),AMPT,FST,DECAL,NPTS+1,SPY)
      AMPT=AMPX(I+1,J,K,NSBs)+AMPX(I,J+1,K,NSBs)+AMPX(I+1,J+1,K,NSBs)
      CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I+1,J+1,K),AMPT,FST,DECAL,NPTS+1,SPX)
      AMPT=AMPY(I+1,J,K,NSBs)+AMPY(I,J+1,K,NSBs)+AMPY(I+1,J+1,K,NSBs)
   10 CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I+1,J+1,K),AMPT,FST,DECAL,NPTS+1,SPY)
      DO 20 I=0,NT-1
      J=NT-1-I
      AMPT=AMPX(I+1,J,K,NSBs)+AMPX(I,J+1,K,NSBs)+AMPX(I,J,K,NSBs)
      CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I,J,K),AMPT,FST,DECAL,NPTS+1,SPX)
      AMPT=AMPY(I+1,J,K,NSBs)+AMPY(I,J+1,K,NSBs)+AMPY(I,J,K,NSBs)
   20 CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I,J,K),AMPT,FST,DECAL,NPTS+1,SPY)
   30 CONTINUE
   31 CONTINUE
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE SPECTREB()
!***********************************************************************
      USE common_module
      DO 31 K=1,Ns-1
      IF(NALL.EQ.0.AND.(2*K).NE.INT(2*SPINs+1)) GOTO 31
      DO 30 NSBs=-NSB,NSB
      DECAL=NSBs*WROT
      DO 30 I=0,NT-1                             
      DO 10 J=0,NT-1                                                         
      AMPT=AMPX(I+1,J,K,NSBs)+AMPX(I,J+1,K,NSBs)+AMPX(I,J,K,NSBs)
      CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I,J,K),AMPT,FST,DECAL,NPTS+1,SPX)
      AMPT=AMPY(I+1,J,K,NSBs)+AMPY(I,J+1,K,NSBs)+AMPY(I,J,K,NSBs)
      CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I,J,K),AMPT,FST,DECAL,NPTS+1,SPY)
      AMPT=AMPX(I+1,J,K,NSBs)+AMPX(I,J+1,K,NSBs)+AMPX(I+1,J+1,K,NSBs)
      CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I+1,J+1,K),AMPT,FST,DECAL,NPTS+1,SPX)
      AMPT=AMPY(I+1,J,K,NSBs)+AMPY(I,J+1,K,NSBs)+AMPY(I+1,J+1,K,NSBs)
   10 CALL TENT(WR(I+1,J,K),WR(I,J+1,K),WR(I+1,J+1,K),AMPT,FST,DECAL,NPTS+1,SPY)
      DO 20 J=NT,2*NT-1                                               
      AMPT=AMPX(I,J,K,NSBs)+AMPX(I+1,J+1,K,NSBs)+AMPX(I+1,J,K,NSBs)
      CALL TENT(WR(I,J,K),WR(I+1,J+1,K),WR(I+1,J,K),AMPT,FST,DECAL,NPTS+1,SPX)
      AMPT=AMPY(I,J,K,NSBs)+AMPY(I+1,J+1,K,NSBs)+AMPY(I+1,J,K,NSBs)
      CALL TENT(WR(I,J,K),WR(I+1,J+1,K),WR(I+1,J,K),AMPT,FST,DECAL,NPTS+1,SPY)
      AMPT=AMPX(I,J,K,NSBs)+AMPX(I+1,J+1,K,NSBs)+AMPX(I,J+1,K,NSBs)
      CALL TENT(WR(I,J,K),WR(I+1,J+1,K),WR(I,J+1,K),AMPT,FST,DECAL,NPTS+1,SPX)
      AMPT=AMPY(I,J,K,NSBs)+AMPY(I+1,J+1,K,NSBs)+AMPY(I,J+1,K,NSBs)
   20 CALL TENT(WR(I,J,K),WR(I+1,J+1,K),WR(I,J+1,K),AMPT,FST,DECAL,NPTS+1,SPY)
   30 CONTINUE
   31 CONTINUE
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE DIAG(N,HTSRT,HTSIT,VP,ZR,ZI)
!***********************************************************************
      DIMENSION E(100),TAU(2,100),VP(100),ZR(100,100),ZI(100,100)
      DO 20 I=1,N
      DO 10 J=1,N
      ZR(I,J)=0.
   10 ZI(I,J)=0.
      VP(I)=0. 
   20 ZR(I,I)=1.
      CALL HTRIDI(N,HTSRT,HTSIT,VP,E,TAU)
      CALL TQL2(N,VP,E,ZR)
      CALL HTRIBK(N,HTSRT,HTSIT,TAU,N,ZR,ZI)
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE HTRIDI(N,HTSR,HTSI,VP,E,TAU)
!***********************************************************************
      REAL HTSR(100,100),HTSI(100,100),VP(100),E(100),TAU(2,100)
      TAU(1,N)=1.
      TAU(2,N)=0.
      DO 100 I=1,N
  100 VP(I)=HTSR(I,I)
      DO 300 II=1,N
      I=N+1-II
      L=I-1
      H=0.
      SCALE=0.
      IF(L.LT.1) GOTO 130
      DO 120 K=1,L
  120 SCALE=SCALE+ABS(HTSR(I,K))+ABS(HTSI(I,K))
      IF(SCALE.NE.0.0) GOTO 140
      TAU(1,L)=1.
      TAU(2,L)=0.
  130 E(I)=0.
      GOTO 290
  140 DO 150 K=1,L
          HTSR(I,K)=HTSR(I,K)/SCALE
          HTSI(I,K)=HTSI(I,K)/SCALE
  150 H=H+HTSR(I,K)*HTSR(I,K)+HTSI(I,K)*HTSI(I,K)
      G=SQRT(H)
      E(I)=SCALE*G
      F=CABS(CMPLX(HTSR(I,L),HTSI(I,L)))
      IF(F.EQ.0.0) GOTO 160
      TAU(1,L)=(HTSI(I,L)*TAU(2,I)-HTSR(I,L)*TAU(1,I))/F
      SI=(HTSR(I,L)*TAU(2,I)+HTSI(I,L)*TAU(1,I))/F
      H=H+F*G
      G=1.+G/F
      HTSR(I,L)=G*HTSR(I,L)
      HTSI(I,L)=G*HTSI(I,L)
      IF(L.EQ.1) GOTO 270
      GOTO 170
  160 TAU(1,L)=-TAU(1,I)
      SI=TAU(2,I)
      HTSR(I,L)=G
  170 F=0.
      DO 240 J=1,L
      G=0.
      GI=0.
      DO 180 K=1,J                                            
      G=G+HTSR(J,K)*HTSR(I,K)+HTSI(J,K)*HTSI(I,K)
  180 GI=GI-HTSR(J,K)*HTSI(I,K)+HTSI(J,K)*HTSR(I,K)
      JP1=J+1
      IF(L.LT.JP1) GOTO 220
      DO 200 K=JP1,L
      G=G+HTSR(K,J)*HTSR(I,K)-HTSI(K,J)*HTSI(I,K)
  200 GI=GI-HTSR(K,J)*HTSI(I,K)-HTSI(K,J)*HTSR(I,K)
  220 E(J)=G/H
      TAU(2,J)=GI/H
  240 F=F+E(J)*HTSR(I,J)-TAU(2,J)*HTSI(I,J)
      HH=F/(H+H)
      DO 260 J=1,L
      F=HTSR(I,J)
      G=E(J)-HH*F
      E(J)=G
      FI=-HTSI(I,J)
      GI=TAU(2,J)-HH*FI
      TAU(2,J)=-GI
      DO 260 K=1,J
      HTSR(J,K)=HTSR(J,K)-F*E(K)-G*HTSR(I,K)+FI*TAU(2,K)+GI*HTSI(I,K)
  260 HTSI(J,K)=HTSI(J,K)-F*TAU(2,K)-G*HTSI(I,K)-FI*E(K)-GI*HTSR(I,K)
  270 DO 280 K=1,L
      HTSR(I,K)=SCALE*HTSR(I,K)
  280 HTSI(I,K)=SCALE*HTSI(I,K)
      TAU(2,L)=-SI
  290 HH=VP(I)
      VP(I)=HTSR(I,I)
      HTSR(I,I)=HH
  300 HTSI(I,I)=SCALE*SQRT(H)
       RETURN
       END
!
!***********************************************************************
       SUBROUTINE TQL2(N,VP,E,Z)
!***********************************************************************       
       REAL VP(100),E(100),Z(100,100)
       ZACHEP=2.**(-26)
       DO 100 I=2,N
  100  E(I-1)=E(I)
       F=0.
       B=0.
       E(N)=0.
       DO 240 L=1,N
          J=0
          H=ZACHEP*(ABS(VP(L))+ABS(E(L)))
          IF(B.LT.H) B=H
          DO 110 M=L,N
  110     IF(ABS(E(M)).LE.B) GOTO 120
  120     IF(M.EQ.L) GOTO 220
  130     IF(J.EQ.30) WRITE(6,*) 'CA MERDE'
          J=J+1
          L1=L+1
          G=VP(L)
          P=(VP(L1)-G)/(2.*E(L))
          R=SQRT(P*P+1.)
          VP(L)=E(L)/(P+SIGN(R,P))
          H=G-VP(L)
          DO 140 I=L1,N
  140     VP(I)=VP(I)-H
          F=F+H
          P=VP(M)
          C=1.
          S=0.
          MML=M-L
          DO 200 II=1,MML
             I=M-II
             G=C*E(I)
             H=C*P
             IF(ABS(P).LT.ABS(E(I))) GOTO 150
             C=E(I)/P
             R=SQRT(C*C+1.)
             E(I+1)=S*P*R
             S=C/R
             C=1./R
             GOTO 160
  150        C=P/E(I)
             R=SQRT(C*C+1.)
             E(I+1)=S*E(I)*R
             S=1./R
             C=C*S
  160        P=C*VP(I)-S*G
             VP(I+1)=H+S*(C*G+S*VP(I))
             DO 180 K=1,N
                H=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*H
  180           Z(K,I)=C*Z(K,I)-S*H
  200     CONTINUE
          E(L)=S*P
          VP(L)=C*P
          IF(ABS(E(L)).GT.B) GOTO 130
  220     VP(L)=VP(L)+F
  240  CONTINUE
       DO 300 II=2,N
          I=II-1
          K=I
          P=VP(I)
          DO 260 J=II,N
             IF(VP(J).GE.P) GOTO 260
             K=J
             P=VP(J)
  260     CONTINUE 
          IF(K.EQ.I) GOTO 300
          VP(K)=VP(I)
          VP(I)=P
          DO 280 J=1,N
             P=Z(J,I)
             Z(J,I)=Z(J,K)
             Z(J,K)=P
  280     CONTINUE
  300  CONTINUE
       RETURN
       END
!
!***********************************************************************
       SUBROUTINE HTRIBK(N,HTSR,HTSI,TAU,M,ZR,ZI)
!***********************************************************************       
       REAL HTSR(100,100),HTSI(100,100),TAU(2,100)
       REAL ZR(100,100),ZI(100,100)
       DO 50 K=1,N
          DO 50 J=1,M
          ZI(K,J)=-ZR(K,J)*TAU(2,K)
  50      ZR(K,J)=ZR(K,J)*TAU(1,K)
       DO 140 I=2,N
          L=I-1
          H=HTSI(I,I)
          IF(H.EQ.0.0) GOTO 140
          DO 130 J=1,M
             S=0.
             SI=0.
             DO 110 K=1,L
                S=S+HTSR(I,K)*ZR(K,J)-HTSI(I,K)*ZI(K,J)
  110           SI=SI+HTSR(I,K)*ZI(K,J)+HTSI(I,K)*ZR(K,J)
             S=S/H/H
             SI=SI/H/H
             DO 120 K=1,L
                ZR(K,J)=ZR(K,J)-S*HTSR(I,K)-SI*HTSI(I,K)
  120           ZI(K,J)=ZI(K,J)-SI*HTSR(I,K)+S*HTSI(I,K)
  130     CONTINUE
  140  CONTINUE
       RETURN
       END
!
!***********************************************************************
      SUBROUTINE WIGNER(C,S,D)
!***********************************************************************
      REAL*4 D(-4:4,-4:4)
      D(-2,-2)=((1.+C)**2)/4.
      D(2,2)=D(-2,-2)
      D(2,1)=-(1.+C)*S/2.
      D(-2,-1)=-D(2,1)
      D(1,2)=-D(2,1)
      D(-1,-2)=D(2,1)
      D(2,0)=SQRT(3./8.)*S*S
      D(0,2)=D(2,0)
      D(-2,0)=D(2,0)
      D(0,-2)=D(2,0)
      D(2,-1)=-(1.-C)*S/2.
      D(-1,2)=-D(2,-1)
      D(-2,1)=-D(2,-1)
      D(1,-2)=D(2,-1)
      D(2,-2)=((1.-C)**2)/4.
      D(-2,2)=D(2,-2)
      D(1,1)=C*C-(1.-C)/2.
      D(-1,-1)=D(1,1)
      D(1,-1)=(1.+C)/2.-C*C
      D(-1,1)=D(1,-1)
      D(1,0)=-SQRT(3./8.)*2.*C*S
      D(0,1)=-D(1,0)
      D(-1,0)=-D(1,0)
      D(0,-1)=D(1,0)
      D(0,0)=(3.*C*C-1.)/2.
      RETURN
      END
!                                                      
!***********************************************************************
      SUBROUTINE FFT(NP,TR,SPCR,SPCI)
!***********************************************************************
      DIMENSION TR(*),SPCR(*),SPCI(*)
      NP2=NP/2
      CALL FAST(NP,TR)
      DO 10 I=1,NP2
      I2=2*I
      SPCR(I)=TR(I2-1)                   
      SPCR(NP+2-I)=TR(I2-1)
      SPCI(I)=TR(I2)
   10 SPCI(NP+2-I)=-TR(I2)
      SPCR(NP2+1)=TR(NP+1)
      SPCI(NP2+1)=TR(NP+2)
      RETURN
      END
!                                                      
!***********************************************************************
      SUBROUTINE FFT1(NP,TR,TI,SPCR,SPCI)
!***********************************************************************
      DIMENSION TR(*),TI(*),SPCR(*),SPCI(*)
      NP2=NP/2
      CALL FAST(NP,TR)
      CALL FAST(NP,TI)
      DO 10 I=1,NP2
      I2=2*I
      SPCR(I)=TR(I2-1)-TI(I2)                   
      SPCR(NP+2-I)=TR(I2-1)+TI(I2)
      SPCI(I)=TR(I2)+TI(I2-1)
   10 SPCI(NP+2-I)=-TR(I2)+TI(I2-1)
      SPCR(NP2+1)=TR(NP+1)-TI(NP+2)
      SPCI(NP2+1)=TR(NP+2)+TI(NP+1)
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE FAST(NP,TR)
!***********************************************************************
      REAL*4 TR(*)
      M=NINT(LOG(FLOAT(NP))/LOG(2.))       
      N4POW=M/2
      NN=1               
      IF((M-N4POW*2).GT.0) THEN    
      NN=2
      NNN=NP/NN
      DO 10 K=1,NNN
      T=TR(K)+TR(NNN+K)
      TR(NNN+K)=TR(K)-TR(NNN+K)
   10 TR(K)=T
      END IF 
      DO 20 I=1,N4POW  
      NN=NN*4
      NNN=NP/NN
   20 CALL FR4TR(NNN,NN,TR)                                     
      CALL FORD1(NP,TR)
      CALL FORD2(M,TR)
      T=TR(2)
      TR(2)=0.
      TR(NP+1)=T
      TR(NP+2)=0.
      DO 30 I=4,NP,2 
   30 TR(I)=-TR(I)    
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE FR4TR(NNN,NN,TR)
!***********************************************************************
      REAL*4 TR(1)
      INTEGER*4 L(15)
      N2=2*NNN
      N3=3*NNN
      DO 10 I=2,15
   10 L(I)=2
      L(1)=NN/4
      IF(L(1).LE.2) THEN
      L(1)=2
      ELSE
      K=1
   20 K=K+1
      L(K)=L(K-1)/2
      IF(L(K).NE.2) GOTO 20
      END IF
      L15=L(1)
      L14=L(2)
      L13=L(3)
      L12=L(4)
      L11=L(5)
      L10=L(6)
      L9=L(7)
      L8=L(8)
      L7=L(9)
      L6=L(10)
      L5=L(11)
      L4=L(12)
      L3=L(13)
      L2=L(14)
      L1=L(15)
      JI=3
      JL=2
      JR=2     
      P7=0.707106782
      PIOVN=3.141592654/FLOAT(NN)
      DO 70 J1=2,L1,2
      DO 70 J2=J1,L2,L1
      DO 70 J3=J2,L3,L2
      DO 70 J4=J3,L4,L3
      DO 70 J5=J4,L5,L4
      DO 70 J6=J5,L6,L5
      DO 70 J7=J6,L7,L6
      DO 70 J8=J7,L8,L7
      DO 70 J9=J8,L9,L8
      DO 70 J10=J9,L10,L9
      DO 70 J11=J10,L11,L10
      DO 70 J12=J11,L12,L11
      DO 70 J13=J12,L13,L12
      DO 70 J14=J13,L14,L13
      DO 70 JTHET=J14,L15,L14
      IF(JTHET.LE.2) THEN
      DO 30 K=1,NNN
      T0=TR(K)+TR(N2+K)
      T1=TR(NNN+K)+TR(N3+K)
      TR(N2+K)=TR(K)-TR(N2+K)
      TR(N3+K)=TR(NNN+K)-TR(N3+K)
      TR(K)=T0+T1
   30 TR(NNN+K)=T0-T1
      IF(NN.LE.4) GOTO 70
      K0=4*NNN+1
      K1=K0+NNN-1
      DO 40 K=K0,K1
      XR=P7*(TR(NNN+K)-TR(N3+K))
      XI=P7*(TR(NNN+K)+TR(N3+K))
      TR(N3+K)=TR(N2+K)+XI
      TR(NNN+K)=XI-TR(N2+K)
      TR(N2+K)=TR(K)-XR
   40 TR(K)=XR+TR(K)
      ELSE
      ARG=FLOAT(JTHET-2)*PIOVN
      C1=COS(ARG)
      S1=SIN(ARG)
      C2=C1*C1-S1*S1
      S2=C1*S1*2.
      C3=C1*C2-S1*S2
      S3=C2*S1+S2*C1
      J0=4*NNN*JR+1
      K0=4*NNN*JI+1
      JLAST=J0+NNN-1
      DO 50 J=J0,JLAST 
      K=K0+J-J0
      R1=TR(NNN+J)*C1-TR(NNN+K)*S1
      R5=TR(NNN+J)*S1+TR(NNN+K)*C1
      T2=TR(N2+J)*C2-TR(N2+K)*S2
      T6=TR(N2+J)*S2+TR(N2+K)*C2
      T3=TR(N3+J)*C3-TR(N3+K)*S3
      T7=TR(N3+J)*S3+TR(N3+K)*C3
      T0=TR(J)+T2
      T4=TR(K)+T6
      T2=TR(J)-T2
      T6=TR(K)-T6
      T1=R1+T3
      T5=R5+T7
      T3=R1-T3
      T7=R5-T7
      TR(J)=T0+T1
      TR(N3+K)=T4+T5
      TR(N2+K)=T0-T1
      TR(NNN+J)=T5-T4
      TR(N2+J)=T2-T7
      TR(NNN+K)=T6+T3
      TR(K)=T2+T7
   50 TR(N3+J)=T3-T6
      JR=JR+2
      JI=JI-2
      IF(JI.GT.JL) GOTO 60 
      JI=2*JR-1
      JL=JR
   60 END IF 
   70 CONTINUE
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE FORD1(NP,TR)
!***********************************************************************
      REAL*4 TR(1)
      K=4
      K1=2
      DO 10 J=4,NP,2
      IF(K.GT.J) THEN
      T=TR(J)
      TR(J)=TR(K)
      TR(K)=T
      END IF
      K=K-2
      IF(K.LE.K1) THEN
      K=2*J
      K1=J
      END IF
   10 CONTINUE
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE FORD2(M,TR)
!***********************************************************************
      REAL*4 TR(1)
      INTEGER*4 L(15)
      NP=2**M
      DO 10 I=2,15
   10 L(I)=2
      L(1)=NP
      DO 20 K=2,M
   20 L(K)=L(K-1)/2
      L15=L(1)
      L14=L(2)
      L13=L(3)
      L12=L(4)
      L11=L(5)
      L10=L(6)
      L9=L(7)
      L8=L(8)
      L7=L(9)
      L6=L(10)
      L5=L(11)
      L4=L(12)
      L3=L(13)
      L2=L(14)
      L1=L(15)
      IJ=2
      DO 30 J1=2,L1,2
      DO 30 J2=J1,L2,L1
      DO 30 J3=J2,L3,L2
      DO 30 J4=J3,L4,L3
      DO 30 J5=J4,L5,L4
      DO 30 J6=J5,L6,L5
      DO 30 J7=J6,L7,L6
      DO 30 J8=J7,L8,L7
      DO 30 J9=J8,L9,L8
      DO 30 J10=J9,L10,L9 
      DO 30 J11=J10,L11,L10
      DO 30 J12=J11,L12,L11
      DO 30 J13=J12,L13,L12
      DO 30 J14=J13,L14,L13
      DO 30 JI=J14,L15,L14
      IF(IJ.LT.JI) THEN
      T=TR(IJ-1)
      TR(IJ-1)=TR(JI-1)
      TR(JI-1)=T
      T=TR(IJ)
      TR(IJ)=TR(JI)
      TR(JI)=T
      END IF
      IJ=IJ+2
   30 CONTINUE
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE WIDTH(NPTS,DLB,DGB,FST,SPC,FID,SPCI)
!***********************************************************************
      REAL*4 SPC(1),SPCI(1),FID(1)
      DO 10 I=1,NPTS          
   10 SPC(2*NPTS+2-I)=SPC(I)
      CALL FFT(2*NPTS,SPC,FID,SPCI)               
      DNPTS=FLOAT(NPTS)
      DO 20 I=1,NPTS+1  
      WIDTHH=DLB-(I-1)*DGB*DGB*1.133/FST/4.
      FID(I)=FID(I)*EXP(3.1415926*(I-1)*WIDTHH/FST/4.)/DNPTS
   20 FID(2*NPTS+3-I)=0.       
      FID(1)=FID(1)/2.
      FID(NPTS+1)=FID(NPTS+1)/2. 
      RETURN
      END
!
!***********************************************************************
      SUBROUTINE TENT(W1,W2,W3,AMPT,FST,DEC,NPTS,SPC)
!***********************************************************************       
      REAL*4 SPC(1)
      FINC=-2.*FST/FLOAT(NPTS)                                                 
      FF1=AMAX1(W1,W2)                                                         
      FF2=AMAX1(W2,W3)                                                         
      FF3=AMAX1(W3,W1)                                                         
      FMIN=AMIN1(W1,W2,W3)+DEC                                                 
      FMID=AMIN1(FF1,FF2,FF3)+DEC                                          
      FMAX=AMAX1(W1,W2,W3)+DEC                                                
      FDN=FMID-FMIN+0.0000001                            
      FXD=FMAX-FMID+0.0000001                                           
      FMN=FMAX-FMIN+0.0000001                                               
      TOP=AMPT/FMN                                                             
      NP=INT((FMIN-FST)/FINC)+1                                                
      NPMID=INT((FMID-FST)/FINC)+1                                             
      NPMAX=INT((FMAX-FST)/FINC)+1                                             
      IF(NPMAX.GT.NPTS.OR.NP.LT.1) GOTO 80                                   
      IF(NP.NE.NPMID) GOTO 10                                                
      SPC(NP)=SPC(NP)+FDN*TOP                                                  
      GOTO 40                                                                  
   10 F2=FINC*FLOAT(NP)+FST                                                    
      SPC(NP)=SPC(NP)+(F2-FMIN)**2*TOP/FDN                                     
   20 NP=NP+1                                                                  
      F1=F2                                                                    
      IF(NP.EQ.NPMID) GOTO 30                                                 
      F2=FINC*FLOAT(NP)+FST                                                    
      SPC(NP)=SPC(NP)+FINC*(F2+F1-2.*FMIN)*TOP/FDN                             
      GOTO 20                                                                  
   30 SPC(NP)=SPC(NP)+(FMID-F1)*(FDN+F1-FMIN)*TOP/FDN                          
   40 IF(NP.NE.NPMAX) GOTO 50                                                 
      SPC(NP)=SPC(NP)+FXD*TOP                                                  
      GOTO 80                                                                  
   50 F2=FINC*FLOAT(NPMID)+FST                                                 
      SPC(NP)=SPC(NP)+(F2-FMID)*(FMAX-F2+FXD)*TOP/FXD                          
   60 NP=NP+1                                                                  
      F1=F2                                                                    
      IF(NP.EQ.NPMAX) GOTO 70                                                
      F2=FINC*FLOAT(NP)+FST                                                    
      SPC(NP)=SPC(NP)+FINC*(2.*FMAX-F1-F2)*TOP/FXD                             
      GOTO 60                                                                 
   70 SPC(NP)=SPC(NP)+(FMAX-F1)**2*TOP/FXD                                     
   80 CONTINUE                                                                 
      RETURN                                                                   
      END                                                                      
!
!***********************************************************************
      SUBROUTINE OPERATOR()
!***********************************************************************       
      USE common_module
      DIMENSION PRIPi(11,11),PRIPs(11,11),PRIZii(11,11),PRQZii(11,11)
      DIMENSION PRQS1ii(11,11),PRQS2ii(11,11)       
      DO 100 I=1,100
      DO 100 J=1,100
      PRD(I,J)=0.
      PRIZi(I,J)=0. 
      PRQZi(I,J)=0.
      PRIXi(I,J)=0.
      PIIYi(I,J)=0.
      PRQS1i(I,J)=0.
      PRQS2i(I,J)=0.
      PRIZs(I,J)=0.
      PRQZs(I,J)=0.
      PRIXs(I,J)=0.
      PIIYs(I,J)=0.
      PRQS1s(I,J)=0.
  100 PRQS2s(I,J)=0.   
      DO 110 I=1,11
      DO 110 J=1,11
      PRIPi(I,J)=0.
      PRIPs(I,J)=0.
      PRIZii(I,J)=0.
      PRQZii(I,J)=0.
      PRQS1ii(I,J)=0.
  110 PRQS2ii(I,J)=0.  

      DO 220 L=1,4
  220 WM(L)=FLOAT(L)*WROT

      SS=SPINs
      DO 230 I=1,Ns
      PRQS1s(I,I)=SS*(4.*SPINs*(SPINs+1.)-8.*SS*SS-1.)/2.
      PRQS2s(I,I)=SS*(2.*SPINs*(SPINs+1.)-2.*SS*SS-1.)/2.
      PRIZs(I,I)=SS
      PRQZs(I,I)=3.*SS*SS-SPINs*(SPINs+1.)
      IF(I.LT.Ns) PRIPs(I,I+1)=SQRT(SPINs*(SPINs+1.)-SS*(SS-1.)) 
  230 SS=SS-1.
      DO 240 I=1,Ns
      DO 240 KK=0,Ni-1
      PRQS1s(I+KK*Ns,I+KK*Ns)=PRQS1s(I,I)
      PRQS2s(I+KK*Ns,I+KK*Ns)=PRQS2s(I,I)
      PRIZs(I+KK*Ns,I+KK*Ns)=PRIZs(I,I)
      PRQZs(I+KK*Ns,I+KK*Ns)=PRQZs(I,I)
      PRIXs(I+KK*Ns,I+1+KK*Ns)=+PRIPs(I,I+1)/2.
      PIIYs(I+KK*Ns,I+1+KK*Ns)=-PRIPs(I,I+1)/2.
      PRIXs(I+1+KK*Ns,I+KK*Ns)=+PRIPs(I,I+1)/2.  
  240 PIIYs(I+1+KK*Ns,I+KK*Ns)=+PRIPs(I,I+1)/2. 

      SS=SPINi
      DO 250 I=1,Ni
      PRQS1ii(I,I)=SS*(4.*SPINi*(SPINi+1.)-8.*SS*SS-1.)/2.
      PRQS2ii(I,I)=SS*(2.*SPINi*(SPINi+1.)-2.*SS*SS-1.)/2.
      PRIZii(I,I)=SS
      PRQZii(I,I)=3.*SS*SS-SPINi*(SPINi+1.)
      IF(I.LT.Ni) PRIPi(I,I+1)=SQRT(SPINi*(SPINi+1.)-SS*(SS-1.))
  250 SS=SS-1. 
      DO 260 I=1,Ni
      DO 260 KK=1,Ns   
      PRD(KK+(I-1)*Ns,KK+(I-1)*Ns)=PRIZii(I,I)*PRIZs(KK+(I-1)*Ns,KK+(I-1)*Ns)
      PRQS1i(KK+(I-1)*Ns,KK+(I-1)*Ns)=PRQS1ii(I,I)
      PRQS2i(KK+(I-1)*Ns,KK+(I-1)*Ns)=PRQS2ii(I,I)
      PRIZi(KK+(I-1)*Ns,KK+(I-1)*Ns)=PRIZii(I,I)
      PRQZi(KK+(I-1)*Ns,KK+(I-1)*Ns)=PRQZii(I,I)
      PRIXi(KK+(I-1)*Ns,KK+I*Ns)=+PRIPi(I,I+1)/2.
      PIIYi(KK+(I-1)*Ns,KK+I*Ns)=-PRIPi(I,I+1)/2.
      PRIXi(KK+I*Ns,KK+(I-1)*Ns)=+PRIPi(I,I+1)/2.
  260 PIIYi(KK+I*Ns,KK+(I-1)*Ns)=+PRIPi(I,I+1)/2.
      RETURN 
      END
!
!***********************************************************************
      SUBROUTINE ZERO()
!***********************************************************************
      USE common_module
      DO 100 IT=1,NPTS
      SPX(IT)=0.
  100 SPY(IT)=0.
        
      DO 200 I=1,100
      DO 200 J=1,100
      RorTt(I,J)=0.
  200 RoiTt(I,J)=0. 
         
      DO 300 I=1,10
      DO 300 J=0,2
      FFCi(I,J)=0.
      GGCi(I,J)=0.
      FFQi(I,J)=0.
      GGQi(I,J)=0.
      FFDs(I,J)=0.
      GGDs(I,J)=0.  
      FFCs(I,J)=0.
 300  GGCs(I,J)=0.  
      RETURN
      END
