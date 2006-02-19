!     Last change:  CF   19 Feb 2006   10:55 am
! File: $Id: ll.f95, v 0.1
! ----------------------------------------------------------------------
! PULSAR Project
!        Copyright (C) 2006 Jean-Paul Amoureux, Christian Fernandez
!        JPA - Unité de Catalyse et Chimie du Solide, Lille.
!        CF  - Laboratoire Catalyse et Spectrochimie, Caen.
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
!  ' ll' subroutines
! ----------------------------------------------------------------------
! LAST MODIFICATIONS:
!   20060215 CF - Creating this file 'll.for'
!         19 CF - Modification parametres dans LLcp etc... (bien que ces functions ne soit pas encore utilisées)
!-----------------------------------------------------------------------

!***********************************************************************
      SUBROUTINE LLcp(ZZ)
!***********************************************************************
      USE common_module
      DIMENSION B(-2:2),BB(-2:2)
!                                          SUBROUTINE DE CALCUL DU CP
      eps=-2.*FST
      Wes=SQRT(PULSEs(1,1)**2+PULSEs(1,2)**2)
      Wei=SQRT(PULSEi(1,1)**2+PULSEi(1,2)**2)
      DWe=Wes+eps*Wei
      SB=SQRT(ABS(1.-ZZ*ZZ))
      B(0)=Wdip(4)*(1.-3.*ZZ*ZZ)/2.
      B(1)=Wdip(4)*SB*ZZ*0.7071067
      B(2)=Wdip(4)*SB*SB/4.

      DO 5 I=0,2
      BB(I)=(B(I)*PULSEi(1,1)*PULSEs(1,1)/Wei/Wes)**2
    5 BB(-I)=BB(I)

!ccc      DO 10 I=-2,2    à vérifier: celà a changé
!ccc   10    RorTt(1,I+4)=RorTt(1,I+4)+(1.-COS(B(I)*PULSEs(1,4)))/R3
!ccc     REN=SQRT((DWe-I*WROT)**2+BB(I)+0.00001)
!ccc   10 RorTt(1,I+4)=RorTt(1,I+4)-eps*(1.-COS(REN*PULSEs(1,4)))
!ccc     &*BB(I)/REN/REN/R3

!ccc      RorTt(1,1)=RorTs(1,1)+1./R3 
      RETURN
      END
!       
!***********************************************************************
      SUBROUTINE LLredor(YY,ZZ,R3)
!***********************************************************************
      USE common_module
!                            INTRODUIRE LLredor dans EDGEA pour la moyenne de poudre
      TPI=6.2831853072
      CST=2.828427125
      SB=SQRT(ABS(1.-ZZ*ZZ))
      S2B=2.*ZZ*SB
      SA=0.
      IF(SB.GT.0.0000001) SA=YY/SB
!     DNcs est le nombre total de périodes pendant lesquelles agit le déphasage dipolaire 
      DNcs=Pfit
      PHIss=Wdip(4)*CST*DNcs*S2B*SA/WROT
      PHIs=AMOD(PHIss,TPI)
      RorTt(1,1)=RorTt(1,1)+COS(PHIs)/R3      
      RorTt(2,2)=RorTt(2,2)+(COS(PHIs))**2/R3
      RETURN
      END
!              
!***********************************************************************
      SUBROUTINE LL3spinsRedor(XX,YY,ZZ,R3)
!***********************************************************************
      USE common_module
!                             TRES IMPORTANT: UTILISER EDGEB POUR LA MOYENNE DE POUDRE
!     SI ON INTRODUIT UNE SEULE PAIRE ISOLEE, LA METTRE SOIT EN Dip1, SOIT EN Dip2 avec Angle=0
!     L'ordre de Dip1 et Dip2 n'a pas d'importance
!     NE PAS OUBLIER DE CHANGER A123T>0  pour passer par EDGEB

      PI=3.141592654
      SB=SQRT(ABS(1.-ZZ*ZZ))
      IF(SB.GT.0.00000001) THEN
              CA=XX/SB
              SA=YY/SB
      ELSE
              CA=1.
              SA=0.
      END IF

!     DNcs est le nombre total de périodes pendant lesquelles agit le déphasage dipolaire 
      DNcs=Pfit
      CSTTMP=DNcs*5.656854249/WROT
      Ntgamma=100            !nombre de pas d'intégration en gamma sur 360°
      FNtgamma=FLOAT(Ntgamma) 
      Angle=81.06             !angle entre les 2 paires en degrés
      Dip1=1808.              !valeur en Hz du dipolaire de la paire 1 
      Dip2=1620.              !valeur en Hz du dipolaire de la paire 2
      Cangle=COS(PI*Angle/180.)
      Sangle=SIN(PI*Angle/180.)
      COS1=COS(CSTTMP*Dip1*ZZ*SB*SA)
      Sumgamma=0. 
      DO 10 Ng=1,Ntgamma
         gamma=2.*PI*Ng/FNtgamma
         Sgamma=SIN(gamma)
         Cgamma=COS(gamma)
         COS2=COS(CSTTMP*Dip2*(Sgamma*SB*Sangle+ZZ*Cangle)*(Sangle*(Cgamma*CA-Sgamma*ZZ*SA)+Cangle*SA*SB))
   10 Sumgamma=Sumgamma+COS2
      RorTt(1,1)=RorTt(1,1)+Sumgamma*COS1/R3/FNtgamma  
      RETURN
      END
!       
!***********************************************************************
      SUBROUTINE LLtedor(YY,ZZ,R3)
!***********************************************************************
      USE common_module
      TPI=6.2831853072
      CST=2.828427125
      SB=SQRT(ABS(1.-ZZ*ZZ))
      S2B=2.*ZZ*SB
      SA=0.
      IF(SB.GT.0.0000001) SA=YY/SB
      DNci=2.
      DNcs=4.
      PHIss=Wdip(4)*CST*DNcs*S2B*SA/WROT
      PHIii=Wdip(4)*CST*DNci*S2B*SA/WROT
      PHIs=AMOD(PHIss,TPI)
      PHIi=AMOD(PHIii,TPI)
      RorTt(1,2)=RorTt(1,2)+SIN(PHIs)*SIN(PHIi)/R3
      RorTt(2,2)=RorTt(2,2)+COS(PHIs)*COS(PHIi)/R3
      RETURN
      END
!


!***********************************************************************     
      SUBROUTINE LL(XX,YY,ZZ,II,JJ,R3)
!***********************************************************************     
      USE common_module

      real, parameter :: sq6=2.4494897427831780981972840747059

      DIMENSION Ror(100,100),Roi(100,100),YR(100,100),YI(100,100)    &
     ,TR(100,100),TI(100,100),RoTR(100,100),RoTI(100,100),ZL(0:4)    &
     ,ZR(100,100),ZI(100,100),VP(100),VVi(2),VVs(2),VVIi(2),VVIs(2)  &
     ,ARQs(10,-4:4),AIQs(10,-4:4),ARCs(10,-4:4),AICs(10,-4:4)        &
     ,ARDs(10,-4:4),AIDs(10,-4:4),ARts(-4:4),AIts(-4:4)              &
     ,ARQi(10,-4:4),AIQi(10,-4:4),ARCi(10,-4:4),AICi(10,-4:4)        &
     ,PR(11,-20:20),PI(11,-20:20),QR(11,-20:20),QI(11,-20:20)        &
     ,HTSR(100,100),HTSRT(100,100),HTSIT(100,100),AB(11,4),BB(11,4)  &
     ,ARQPi(0:4),AIQPi(0:4),ARQPs(0:4),AIQPs(0:4)                    &
     ,ARCPi(0:4),AICPi(0:4),ARDIPis(0:4),AIDIPis(0:4)
      DIMENSION ARDCPs(0:4),AIDCPs(0:4)                             &
     ,ARQSi(2,0:4),AIQSi(2,0:4),ARQSs(2,0:4),AIQSs(2,0:4)           &
     ,TRi(100,100),TIi(100,100),TRf(100,100),TIf(100,100)           &
     ,TRsi(100,100),TIsi(100,100),TRsf(100,100),TIsf(100,100)       &
     ,TRint(100,100),TIint(100,100),Rorsto(100,100),Roisto(100,100) &
     ,STRi(11,11),STIi(11,11)
      DATA ZL/1.,2.,2.,2.,2./
      TPI=6.2831853072 

      SB=SQRT(ABS(1.-ZZ*ZZ))
      IF(SB.NE.0.) THEN
              CA=XX/SB
              SA=YY/SB
              C2A=(XX*XX-YY*YY)/SB/SB
              S2A=2.*XX*YY/SB/SB
      ELSE
              CA=1.
              SA=0.
              C2A=1.
              S2A=0.
      END IF

      CALL AMcalc(CA,SA,C2A,S2A,ZZ,SB,ARts,AIts,ARQi,AIQi,ARCi,AICi,ARQs,AIQs,ARCs,AICs,ARDs,AIDs,PR,PI,QR,QI)

      DO 31 L=0,4
      IF((ISPEED.EQ.1).AND.(L.NE.0))  GOTO 31                        
      RRRR=ZL(L)*DE(L,0)
      ARQPi(L)=RRRR*ARQi(1,L)
      AIQPi(L)=RRRR*AIQi(1,L)
      ARCPi(L)=RRRR*ARCi(1,L)*SQ6*DELTACi  
      AICPi(L)=RRRR*AICi(1,L)*SQ6*DELTACi
      ARDIPis(L)=SQ6*RRRR*ARDs(4,L)*Wdip(4) 
      AIDIPis(L)=SQ6*RRRR*AIDs(4,L)*Wdip(4) 
      ARQPs(L)=RRRR*ARQs(1,L)
      AIQPs(L)=RRRR*AIQs(1,L)
      ARDCPs(L)=RRRR*ARts(L)
      AIDCPs(L)=RRRR*AIts(L)
      DO 30 KS=1,2
      ARQSi(KS,L)=0.
      AIQSi(KS,L)=0.
      ARQSs(KS,L)=0.
      AIQSs(KS,L)=0.
      DO 30 I=L-2,2
      ARQSi(KS,L)=ARQSi(KS,L)+ZL(L)*DK(KS,L,I)*(ARQi(1,I)*ARQi(1,L-I)-AIQi(1,I)*AIQi(1,L-I))
      AIQSi(KS,L)=AIQSi(KS,L)+ZL(L)*DK(KS,L,I)*(ARQi(1,I)*AIQi(1,L-I)+AIQi(1,I)*ARQi(1,L-I))
      ARQSs(KS,L)=ARQSs(KS,L)+ZL(L)*DK(KS,L,I)*(ARQs(1,I)*ARQs(1,L-I)-AIQs(1,I)*AIQs(1,L-I))
   30 AIQSs(KS,L)=AIQSs(KS,L)+ZL(L)*DK(KS,L,I)*(ARQs(1,I)*AIQs(1,L-I)+AIQs(1,I)*ARQs(1,L-I))
   31 CONTINUE

      DO 41 K=1,Ns-1
      DDM=2.*FLOAT(K)-2.*SPINs-1.
      WR(II,JJ,K)=U1s(K)*ARQSs(1,0)+U2s(K)*ARQSs(2,0)+ARDCPs(0)-3.*DDM*COEFPs*ARQPs(0)+DISOst
      IF(ISPEED.EQ.1)    GOTO 41                    
      DO 40 L=1,4                          
      AB(K,L)=(U1s(K)*ARQSs(1,L)+U2s(K)*ARQSs(2,L)+ARDCPs(L)-3.*DDM*COEFPs*ARQPs(L))/WM(L)
   40 BB(K,L)=(U1s(K)*AIQSs(1,L)+U2s(K)*AIQSs(2,L)+AIDCPs(L)-3.*DDM*COEFPs*AIQPs(L))/WM(L)
   41 CONTINUE

!     WR(II,JJ,K),PR(K,NSBs)et PI(K,NSBs) sont indépendants de gamma et du cyclage de phase.
!     on les calcule donc une fois pour toute pour chaque cristallite.

      U=-DU
      DO 70 IU=1,NU
      U=U+DU
      DO 70 K=1,Ns-1
         ARG=0.     
         DO 60 L=1,4                       
   60    ARG=ARG+AB(K,L)*SIN(L*U)-BB(K,L)*COS(L*U)
      DO 70 NSBs=-NSB,NSB
      PR(K,NSBs)=PR(K,NSBs)+COS(ARG-NSBs*U)
   70 PI(K,NSBs)=PI(K,NSBs)+SIN(ARG-NSBs*U)  
!*****************************************BOUCLE SUR GAMMA**************************************************
! 'T' matrices (T,Tf,Tint,Tsi,Tsf,Ti) are in fact U evolution operator matrices
      NPHI=1
      G=-DU

      DO 700 IG=1,NU
      G=G+DU

!                        matrice densité réduite initiale 
  100 DO 110 I=1,Nsi
      DO 110 J=1,Nsi
      Ror(I,J)=PRIZs(I,J)+PRIZi(I,J)*WLi/WLs
  110 Roi(I,J)=0.

      Tat=0.
!**************** BOUCLE SUR LES PERIODES DE PREPARATION **********
      DO 590 NC=1,NCYCL
!     ------period with rf-fields---------------------------------
      TP=-DTP(NC)/2.
       
!     lorsque le pulse est court (<1 tour), TRi et TIi ne servent à rien : matrices identitées
!     lorsque le pulse est long  (>1 tour), Tri et TIi ont été ré-initialisées (280) et servent
!     de point de départ pour le dernier petit pas de pulse.
      DO 150 I=1,Nsi
      DO 150 J=1,Nsi
      TRi(I,J)=0.
      TIi(I,J)=0. 
  150 IF(I.EQ.J) TRi(I,J)=1.
!     --------pulse-decomposition in small rf-steps                   
      DO 300 ITP=1,NTP(NC)
      DTPP=DTP(NC)
      TP=TP+DTPP
      IF(ITOUR(NC).GE.0.AND.ITP.EQ.NTP(NC)) THEN 
!                  dans le cas où le pulse est plus long qu'un tour, on traite à la  fin le dernier petit pulse
                   DTPP=DTP1(NC)
                   IF(DTPP.LT.0.0000001) GOTO 300
                   TP=PULSEs(NC,4)-DTPP/2.
      END IF
                 
      V0Qi=ARQPi(0)
      V0Qs=ARQPs(0)
      V0Ci=ARCPi(0)
      V0DCs=ARDCPs(0)
      V0DIP=ARDIPis(0)
      VVi(1)=ARQSi(1,0)
      VVs(1)=ARQSs(1,0)
      VVi(2)=ARQSi(2,0)
      VVs(2)=ARQSs(2,0)
      IF(ISPEED.EQ.1)    GOTO 170

      DO 160 L=1,4
      TETAF=WM(L)*(Tat+TP)+L*G
      CTETAF=COS(TETAF)
      STETAF=SIN(TETAF)
      V0Qi=V0Qi+ARQPi(L)*CTETAF+AIQPi(L)*STETAF      
      V0Qs=V0Qs+ARQPs(L)*CTETAF+AIQPs(L)*STETAF
      V0Ci=V0Ci+ARCPi(L)*CTETAF+AICPi(L)*STETAF
      V0DCs=V0DCs+ARDCPs(L)*CTETAF+AIDCPs(L)*STETAF
      V0DIP=V0DIP+ARDIPis(L)*CTETAF+AIDIPis(L)*STETAF
      VVi(1)=VVi(1)+ARQSi(1,L)*CTETAF+AIQSi(1,L)*STETAF
      VVs(1)=VVs(1)+ARQSs(1,L)*CTETAF+AIQSs(1,L)*STETAF
      VVi(2)=VVi(2)+ARQSi(2,L)*CTETAF+AIQSi(2,L)*STETAF  
  160 VVs(2)=VVs(2)+ARQSs(2,L)*CTETAF+AIQSs(2,L)*STETAF 
              
  170 DO 180 I=1,Nsi
      DO 180 J=1,Nsi
      HTSRT(I,J)=(DISOi-PULSEi(NC,2)+V0Ci)*PRIZi(I,J)              &
               +(DISOst-PULSEs(NC,2)+V0DCs)*PRIZs(I,J)            &
               +(V0DIP+DJ(4))*PRD(I,J)                            &
               +COEFPi*V0Qi*PRQZi(I,J)+COEFPs*V0Qs*PRQZs(I,J)     &
               +COEFSi*(VVi(1)*PRQS1i(I,J)+VVi(2)*PRQS2i(I,J))    &
               +COEFSs*(VVs(1)*PRQS1s(I,J)+VVs(2)*PRQS2s(I,J))    &
               +PULSEi(NC,1)*PRIXi(I,J)*COS(PULSEi(NC,3))         &
               +PULSEs(NC,1)*PRIXs(I,J)*COS(PULSEs(NC,3))
  180 HTSIT(I,J)=PULSEi(NC,1)*PIIYi(I,J)*SIN(PULSEi(NC,3))         &
               +PULSEs(NC,1)*PIIYs(I,J)*SIN(PULSEs(NC,3))
   
      CALL DIAG(Nsi,HTSRT,HTSIT,VP,ZR,ZI)
     
      DO 190 J=1,Nsi
      CVP=COS(VP(J)*DTPP)
      SVP=SIN(VP(J)*DTPP)
      DO 190 I=1,Nsi
      YR(I,J)=ZR(I,J)*CVP-ZI(I,J)*SVP
  190 YI(I,J)=ZR(I,J)*SVP+ZI(I,J)*CVP 

      DO 200 I=1,Nsi
      DO 200 J=1,Nsi
         TR(I,J)=0.
         TI(I,J)=0.
         DO 200 K=1,Nsi
            TR(I,J)=TR(I,J)+YR(I,K)*ZR(J,K)+YI(I,K)*ZI(J,K)
  200       TI(I,J)=TI(I,J)-YR(I,K)*ZI(J,K)+YI(I,K)*ZR(J,K)       
            
      DO 210 I=1,Nsi
      DO 210 J=1,Nsi
         TRf(I,J)=0.
         TIf(I,J)=0.
         DO 210 K=1,Nsi
            TRf(I,J)=TRf(I,J)+TRi(I,K)*TR(K,J)-TIi(I,K)*TI(K,J)
  210       TIf(I,J)=TIf(I,J)+TRi(I,K)*TI(K,J)+TIi(I,K)*TR(K,J)
                                             
      IF(ITP.EQ.NTPint(NC)) THEN
!            storage de l'opérateur d'évolution intermédiaire dans le cas où le pulse>1 tour
      DO 220 I=1,Nsi
      DO 220 J=1,Nsi
      TRint(I,J)=TRf(I,J)
  220 TIint(I,J)=TIf(I,J) 
      END IF
!     ++++++++++++++++++++++++lorsque pulse>1 tour+++++++++++++++++++++       
      IF(ITOUR(NC).GE.0.AND.ITP.EQ.NTP(NC)-1) THEN
      DO 230 I=1,Nsi
      DO 230 J=1,Nsi
      TRsi(I,J)=TRf(I,J)
      TIsi(I,J)=TIf(I,J)  
      TRsf(I,J)=TRf(I,J)
  230 TIsf(I,J)=TIf(I,J) 
!     à la fin du 1° tour,on calcule (de 230 à 260) la matrice densité correspondant à ITOUR+1 tours 
!     entiers: c'est la puissance ITOUR+1 de celle correspondant à 1 tour.
!     on doit donc remultiplier ITOUR fois celle correspondant à 1 tour.
      DO 260 IB=1,ITOUR(NC)
      DO 240 I=1,Nsi 
      DO 240 J=1,Nsi
         TRf(I,J)=0.
         TIf(I,J)=0.
         DO 240 K=1,Nsi
            TRf(I,J)=TRf(I,J)+TRsf(I,K)*TRsi(K,J)-TIsf(I,K)*TIsi(K,J) 
  240       TIf(I,J)=TIf(I,J)+TRsf(I,K)*TIsi(K,J)+TIsf(I,K)*TRsi(K,J)        
      DO 250 I=1,Nsi
      DO 250 J=1,Nsi
         TRsf(I,J)=TRf(I,J)
  250    TIsf(I,J)=TIf(I,J) 
  260 CONTINUE 
                                           
      IF(NTPint(NC).GT.0) THEN
!     aprés ITOUR+1 il y a un reste, dont la 1° partie correspond à NTPint pas de durée DTP.
!     l'opérateur d'évolution correspondant à ces NTPint pas, a été stoqué au passage en 220.
      DO 270 I=1,Nsi
      DO 270 J=1,Nsi
         TRf(I,J)=0.
         TIf(I,J)=0.
         DO 270 K=1,Nsi
            TRf(I,J)=TRf(I,J)+TRsf(I,K)*TRint(K,J)-TIsf(I,K)*TIint(K,J)
  270       TIf(I,J)=TIf(I,J)+TRsf(I,K)*TIint(K,J)+TIsf(I,K)*TRint(K,J)        
      END IF
      END IF
!     +++++++++++++++on a traité tout sauf le dernier petit pulse++++++++++++  
      DO 280 I=1,Nsi
      DO 280 J=1,Nsi
      TRi(I,J)=TRf(I,J)
  280 TIi(I,J)=TIf(I,J)  
!     on redémarre de là pour traiter ce dernier petit pulse: on repart à 150

  300 CONTINUE 
             
      DO 310 I=1,Nsi 
      DO 310 J=1,Nsi
         RoTR(I,J)=0.
         RoTI(I,J)=0.
         DO 310 K=1,Nsi
            RoTR(I,J)=RoTR(I,J)+Ror(I,K)*TRf(K,J)-Roi(I,K)*TIf(K,J)
  310       RoTI(I,J)=RoTI(I,J)+Ror(I,K)*TIf(K,J)+Roi(I,K)*TRf(K,J) 
           
      DO 320 I=1,Nsi
      DO 320 J=1,Nsi
         Ror(I,J)=0.
         Roi(I,J)=0.
         DO 320 K=1,Nsi
            Ror(I,J)=Ror(I,J)+TRf(K,I)*RoTR(K,J)+TIf(K,I)*RoTI(K,J)
  320       Roi(I,J)=Roi(I,J)+TRf(K,I)*RoTI(K,J)-TIf(K,I)*RoTR(K,J) 
!     ---------period without rf-field---------------------------- 
      IF(DELAY(NC).EQ.0.) GOTO 490

      V0QIs=    ARQPs(0)*DELAY(NC)                                        
      V0QIi=    ARQPi(0)*DELAY(NC)
      V0DCIs=  ARDCPs(0)*DELAY(NC)                                        
      V0CIi=    ARCPi(0)*DELAY(NC)
      V0DIP=  ARDIPis(0)*DELAY(NC)      
      VVIs(1)=ARQSs(1,0)*DELAY(NC)
      VVIi(1)=ARQSi(1,0)*DELAY(NC)
      VVIs(2)=ARQSs(2,0)*DELAY(NC)
      VVIi(2)=ARQSi(2,0)*DELAY(NC)
      IF(ISPEED.EQ.1)                      GOTO 420               

      DO 410 L=1,4                  
      TETAI=WM(L)*(Tat+PULSEs(NC,4))+L*G
      TETAF=TETAI+WM(L)*DELAY(NC)
      DTETAI=AMOD(TETAI,TPI)
      DTETAF=AMOD(TETAF,TPI)
      CTETA=(COS(DTETAF)-COS(DTETAI))/WM(L)                     
      STETA=(SIN(DTETAF)-SIN(DTETAI))/WM(L)                    
      V0QIs=V0QIs+ARQPs(L)*STETA-AIQPs(L)*CTETA
      V0QIi=V0QIi+ARQPi(L)*STETA-AIQPi(L)*CTETA
      V0DCIs=V0DCIs+ARDCPs(L)*STETA-AIDCPs(L)*CTETA
      V0CIi=V0CIi+ARCPi(L)*STETA-AICPi(L)*CTETA
      V0DIP=V0DIP+ARDIPis(L)*STETA-AIDIPis(L)*CTETA      
      VVIs(1)=VVIs(1)+ARQSs(1,L)*STETA-AIQSs(1,L)*CTETA
      VVIi(1)=VVIi(1)+ARQSi(1,L)*STETA-AIQSi(1,L)*CTETA
      VVIs(2)=VVIs(2)+ARQSs(2,L)*STETA-AIQSs(2,L)*CTETA
  410 VVIi(2)=VVIi(2)+ARQSi(2,L)*STETA-AIQSi(2,L)*CTETA 

  420 DO 430 I=1,Nsi
  430 HTSR(I,I)=                                             &
       ((DISOst-PULSEs(NC,2))*DELAY(NC)+V0DCIs)*PRIZs(I,I)  &
       +((DISOi-PULSEi(NC,2))*DELAY(NC)+V0CIi)*PRIZi(I,I)   &
       +Decouple(NC)*(V0DIP+DJ(4)*DELAY(NC))*PRD(I,I)       &
       +COEFPs*V0QIs*PRQZs(I,I)+COEFPi*V0QIi*PRQZi(I,I)     &
       +COEFSs*(VVIs(1)*PRQS1s(I,I)+VVIs(2)*PRQS2s(I,I))    &
       +COEFSi*(VVIi(1)*PRQS1i(I,I)+VVIi(2)*PRQS2i(I,I))

      DO 440 I=1,Nsi
      DO 440 J=1,Nsi
      VPIJ=HTSR(I,I)-HTSR(J,J)
      DVPIJ=AMOD(VPIJ,TPI)
      CVP=COS(DVPIJ)
      SVP=SIN(DVPIJ)
      RoRR=    Ror(I,J)*CVP+Roi(I,J)*SVP
      Roi(I,J)=Roi(I,J)*CVP-Ror(I,J)*SVP
  440 Ror(I,J)=RoRR    
                                            
  490 IF((T2ss(NC)+T2ii(NC)).GT.1.999) GOTO 510 
!     -------------- introduction des T2 pour chaque cohérence --------
      DO 500 I=1,Nsi
      DO 500 J=1,Nsi  
      INDx=INT((I-1)/Ns) 
      INDy=INT((J-1)/Ns)      
      IF(INDx.EQ.INDy.AND.I.NE.J) THEN 
                   Ror(I,J)=Ror(I,J)*T2ss(NC)
                   Roi(I,J)=Roi(I,J)*T2ss(NC)
      END IF       
      IF(INDx.NE.INDy) THEN 
               IF((I-J).EQ.(Ns*((I-J)/Ns))) THEN    
                   Ror(I,J)=Ror(I,J)*T2ii(NC)
                   Roi(I,J)=Roi(I,J)*T2ii(NC)
               ELSE
                   Ror(I,J)=Ror(I,J)*T2si(NC)
                   Roi(I,J)=Roi(I,J)*T2si(NC)    
               END IF
      END IF        
  500 CONTINUE        
  510 IF(IQ1(NC,1).GE.10) GOTO 580  
!     ---------------- sélection des niveaux de quanta ----------------
!     On sélectionne les niveaux de quanta désirés sur S, tout en réintroduisant les termes de la matrice
!     densité sur I avant sélection. Divise ne sert que lorsque Nboucle>0 et lorsque l'on a fini le 1° passage.
!     On calcule la matrice densité de I (projection de Ror et Roi): chacun de ses termes doit etre
!     conservé à tout moment lorsque l'on sélectionne des niveaux de quanta sur S.
!     Lorsqu'un seul niveau de quanta est sélectionné, il doit l'etre dans la 1° colonne (IQ1).

      DO 520 I=1,Ni     
      DO 520 J=1,Ni        
      STRi(I,J)=0.         
      STIi(I,J)=0.         
      DO 520 KK=1,Ns
      STRi(I,J)=STRi(I,J)+Ror(KK+(I-1)*Ns,KK+(J-1)*Ns)/FLOAT(Ns)
  520 STIi(I,J)=STIi(I,J)+Roi(KK+(I-1)*Ns,KK+(J-1)*Ns)/FLOAT(Ns) 

      Divise=1.          
      IF(NC.EQ.Nboucle) Divise=2.                                                    
      DO 530 I=1,Nsi                 
      DO 530 J=1,Nsi              
      INDx=INT((I-1)/Ns)          
      INDy=INT((J-1)/Ns)      
      IF(INDx.EQ.INDy.AND.(J-I).EQ.IQ1(NC,NPHI).OR.(J-I).EQ.IQ2(NC,NPHI)) GOTO 530
      Ror(I,J)=0.             
      Roi(I,J)=0.              
      IF((I-J).EQ.(Ns*((I-J)/Ns))) THEN   
                   Ror(I,J)=STRi(INDx+1,INDy+1)/Divise
                   Roi(I,J)=STIi(INDx+1,INDy+1)/Divise
      END IF                                 
  530 CONTINUE
                                
      IF(NC.EQ.Nboucle) THEN
            IF(NPHI.EQ.1) THEN
!                      1° passage: on stoque la moitié (Divise=2) de la matrice densité totale
!                      puis on repart au début: goto 100
                       DO 540 I=1,Nsi
                       DO 540 J=1,Nsi
                       Rorsto(I,J)=Ror(I,J)
  540                  Roisto(I,J)=Roi(I,J)               
                   NPHI=2
                   GOTO 100 
            ELSE         
!                      2° passage: on additionne les matrices densités totales des 2 passages (/2)
                       DO 550 I=1,Nsi
                       DO 550 J=1,Nsi
                       Ror(I,J)=Ror(I,J)+Rorsto(I,J)
  550                  Roi(I,J)=Roi(I,J)+Roisto(I,J) 
                   NPHI=1
            END IF
      END IF

  580 CONTINUE

  590 Tat=Tat+PULSEs(NC,4)+DELAY(NC)
      
!     on calcule la matrice densité totale pondérée (R3) par l'angle polaire.
      DO 595 I=1,Nsi
      DO 595 J=1,Nsi
      RorTt(I,J)=RorTt(I,J)+Ror(I,J)/R3
  595 RoiTt(I,J)=RoiTt(I,J)+Roi(I,J)/R3   
!*****INTEGRALES LIEES A LA FIN DE LA PERIODE DE PREPARATION******       
             Gt=WROT*Tat+G
             DO 610 K=1,Ns-1
                      ARG2=0.
                      DO 600 L=1,4         
  600                 ARG2=ARG2+BB(K,L)*COS(L*Gt)-AB(K,L)*SIN(L*Gt)
             DO 610 NSBs=-NSB,NSB
             CARG=COS(ARG2+NSBs*Gt)
             SARG=SIN(ARG2+NSBs*Gt)
                       SIGNALR=0.
                       SIGNALI=0.
                       DO 605 KK=0,Ni-1 
                       SIGNALR=SIGNALR+Ror(K+1+KK*Ns,K+KK*Ns)
  605                  SIGNALI=SIGNALI-Roi(K+1+KK*Ns,K+KK*Ns) 
             QR(K,NSBs)=QR(K,NSBs)+SIGNALR*CARG+SIGNALI*SARG
  610        QI(K,NSBs)=QI(K,NSBs)+SIGNALR*SARG-SIGNALI*CARG     
!============= FIN DE BOUCLE SUR GAMMA ========================================
  700 CONTINUE
!============= PERIODE D'OBSERVATION : Sx et S-y ==============================
             DO 810 K=1,Ns-1
             DO 810 NSBs=-NSB,NSB
        DEM=-DNORM*PRIXs(K,K+1)/(1.+QQ*(WR(II,JJ,K)+NSBs*WROT)**2)/R3
        QRT(II,JJ,K,NSBs)=QRT(II,JJ,K,NSBs)+QR(K,NSBs)  
        QIT(II,JJ,K,NSBs)=QIT(II,JJ,K,NSBs)+QI(K,NSBs)
        AMPX(II,JJ,K,NSBs)=(PI(K,NSBs)*QIT(II,JJ,K,NSBs)       &
                         -PR(K,NSBs)*QRT(II,JJ,K,NSBs))*DEM
 810    AMPY(II,JJ,K,NSBs)=(PI(K,NSBs)*QRT(II,JJ,K,NSBs)       &
                          +PR(K,NSBs)*QIT(II,JJ,K,NSBs))*DEM
 
      RETURN                                        
      END                                                                      
!

