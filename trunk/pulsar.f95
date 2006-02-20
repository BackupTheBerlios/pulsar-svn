!     Last change:  CF   20 Feb 2006    7:44 pm
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
!
! ----------------------------------------------------------------------
! List of subroutines in the various files
!
!PULSAR.f95
!----------
!   PULSAR  Main program
!   Main    Main subroutine
!
! LL.f95
!-------
!   LLcp
!   LLredor
!   LL3spinsRedor
!   LLtedor
!   LL
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

    REAL:: STARTTIME,ENDTIME

    Call CPU_TIME(STARTTIME)

    CALL GetArg(1, cl)            ! Get command line (input filename)

    OPEN (10,file=cl,status="old") ! Open this file for input
    OPEN (20,file="test.out",status="replace") ! Open file for writing

    CALL Main(cl)                 ! Call the main program

    CLOSE(10)                     ! Close the input file
    CLOSE(20)                     ! Close the output file

    Call CPU_TIME(ENDTIME)

    PRINT *, "Execution time: ", ENDTIME-STARTTIME

END PROGRAM


SUBROUTINE Main()
!================

  USE common_module
  USE constantes_module

  CHARACTER(LEN=256) :: fmt0,fmt1,fmt2,fmt3,fmt4,fmtRho,fmtSpec
  CHARACTER*32 :: STR(65)

  REAL :: NISPs(10)

  REAL :: FIDX(32800) ! 32768+2? suffisant à vérifier
  REAL :: FIDY(32800)
  REAL :: SPIX(32800)
  REAL :: SPIY(32800)
  REAL :: P(9999)
  INTEGER :: Idec(9999)
  REAL :: PHASEi(9999)
  INTEGER :: Iref(9999)
  REAL :: RorTs(11,11)
  REAL :: RoiTs(11,11)
  REAL :: RorTi(11,11)
  REAL :: RoiTi(11,11)
  INTEGER :: JJ(20)
  INTEGER :: Nphase(20)
  REAL :: FASE(20)
  INTEGER :: Level(0:20)

  INTEGER :: Inorm1(16)
  INTEGER :: Inorm2(16)
  DATA Inorm1/2154752,933774,589820,430020,338069,278417,236611,205703, &
              181927,163068,147749,135059,124378,115254,107382,100517/
  DATA Inorm2/6122925,3121593,2088467,1568265,1255355,1046463,897112, &
              785115,697918,628170,571061,523524,483263,448748,418821,392639/

  !---------------------!
  ! Reading pulsar file !
  !---------------------!

  READ(10,*) STR(1), P(1)
  Verbose=.false.
  if (P(1).EQ.1) then
    Verbose=.true.
  end if

  DO I=2,59
    READ(10,*) STR(I),P(I)
  END DO

  READ(10,*) STR(60),Nphasing
  READ(10,*) STR(61),Nboucle
  READ(10,*) STR(62),Ifasing

  READ(10,*) (Nphase(I),I=1,20)
  READ(10,*) (Level(J),J=0,20)

  READ(10,*) STR(63),NCYCL

  DO NC=1,NCYCL
    READ(10,*) (P(51+J+10*NC),J=1,3),Iref(NC),IQ1(NC,1),IQ2(NC,1),(P(51+J+10*NC),J=5,9),Idec(NC)
  END DO

  IF (Verbose) THEN
    PRINT *,"Reading Parameters..."
    fmt0="(I3,2X,A10,F13.3)"

    WRITE (*,fmt0) 1,STR(1), P(1)
    DO I=2,59
      WRITE (*,fmt0) I,STR(I), P(I)
    END DO
    fmt1="(A10,I3)"
    WRITE (*,fmt1) STR(60),Nphasing
    WRITE (*,fmt1) STR(61),Nboucle
    WRITE (*,fmt1) STR(62),Ifasing
    fmt2="(A6,3X,20I3)"
    WRITE (*,fmt2) "NPhase",(Nphase(I),I=1,20)
    fmt3="(A6,21I3)"
    WRITE (*,fmt3) "Level",(Level(J),J=0,20)
    WRITE (*,fmt1) STR(63),NCYCL
    fmt4="(A6,I3,F9.0,F7.1,F8.3,3I3,F9.0,F7.1,2F8.3,F10.3,I2)"
    DO NC=1,NCYCL
      WRITE (*,fmt4) &
             "PULSE",NC,(P(51+J+10*NC),J=1,3),Iref(NC),IQ1(NC,1),IQ2(NC,1),(P(51+J+10*NC),J=5,9),Idec(NC)
    END DO
    PRINT *,"End of reading"
  END IF



      ISPEED=0


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

  boucle_t1: DO IUY=10,10,1

    Pfit=0.5*FLOAT(IUY)

    if (Verbose) CALL INFO(STR,P,Nphase,Level,Ifasing,IDNfase)

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
    CALL DIPOLAR(3,P,NISPs,COEFdip)
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

           ! Write some information about the sequence
                                                      IF(Verbose) THEN
           	                                            PRINT *,"SEQUENCE"
                                                        PRINT *,"========"
                                                        DO NC=1,NCYCL
                                                          WRITE (*,fmt4) &
  "PULSE",NC,(P(51+J+10*NC),J=1,2),PULSEs(NC,3)/PI2,Iref(NC),IQ1(NC,1),IQ2(NC,1),(P(51+J+10*NC),J=5,9),Idec(NC)
                                                        END DO
                                                      END IF

    IF(A123T.LT.0.1) THEN
           DNORM=FLOAT(NT)*COEFdip/FLOAT(NU)**2/FLOAT(Ns)/DNscan
           CALL EDGEA()
    ELSE
           DNORM=FLOAT(NT)*COEFdip/FLOAT(NU)**2/FLOAT(Ns)/DNscan/4.
           CALL EDGEB()
    END IF

500 CONTINUE

    FNEXP=FLOAT(NU)*DNscan  ! warning message from gfortran??? : 'fnexp' may be used uninitialized in this function
                            ! Possibly a problem with loops
    IF(A123T.LT.0.1) THEN
           FNEXP=FNEXP*powder1/COEFdip
           CALL SPECTREA()
    ELSE
           FNEXP=FNEXP*powder2/COEFdip
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

    fmtSpec= "(F17.5,2X,F17.5,2X,F17.5)"

    WRITE (20,*) "<spectre>"
    DO I=1,NPTS
      WRITE (20,fmtSpec) FST+(I-1)*(-2*FST/NPTS),1000000*SPX(I),1000000*SPY(I)
    END DO
    WRITE (20,*) "</spectre>"

    ! writing density matrices

    fmtRho= "(100(2X,F7.1,1X,F7.1))"

    WRITE (20,*) "<SPINs>",SPINs,"</SPINs>"
    WRITE (20,*) "<MDs>"
    DO I=1,Ns
      WRITE(20,fmtRho) (RorTs(I,J),RoiTs(I,J),J=1,Ns)
    END DO
    WRITE (20,*) "</MDs>"

    IF(SPINi.GT.0.1) THEN

      WRITE (20,*) "<SPINi>",SPINi,"</SPINi>"
      WRITE (20,*) "<MDi>"
      DO I=1,Ni
      	WRITE (20,fmtRho) (RorTi(I,J),RoiTi(I,J),J=1,Ni)
      END DO
      WRITE (20,*) "</MDi>"

      WRITE (20,*) "<MDis>"
      DO I=1,Nsi
    	WRITE(20,fmtRho) (RorTt(I,J)/FNEXP,RoiTt(I,J)/FNEXP,J=1,Nsi)
      END DO
      WRITE (20,*) "</MDis>"

    END IF

  ENDDO boucle_t1

END SUBROUTINE

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
      SUBROUTINE DIPOLAR(NATsmax,P,NISPs,COEFdip)
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
SUBROUTINE INFO(STR,P,Nphase,Level,Ifasing,IDNfase)
!***********************************************************************       
USE common_module
REAL :: P(*)
INTEGER :: Nphase(20)
INTEGER :: Level(0:20)
CHARACTER*32 :: STR(*)
CHARACTER*256 :: fmt0,fmt1

  fmt0="(A20,F13.3)"
  WRITE(*,fmt0) (STR(I),P(I),I=1,2)
  IF(ISPEED.EQ.0) WRITE(*,fmt0)  STR(3),P(3)
  WRITE(*,fmt0) (STR(I),P(I),I=4,7)
  IF(ISPEED.EQ.0) WRITE(*,fmt0) (STR(I),P(I),I=8,9)
  WRITE(*,fmt0) (STR(I),P(I),I=10,17)
  IF(P(18).NE.0.0) WRITE(*,fmt0) (STR(I),P(I),I=18,22)
  Nec=INT(P(23))
  IF(Nec.GT.0) WRITE(*,fmt0) (STR(I),P(I),I=23,23+5*Nec)
  IF(ISPEED.EQ.0) WRITE(*,fmt0)  STR(39),P(39)
  WRITE(*,fmt0)  STR(40),P(40)
  IF(SPINi.GT.0.1) WRITE(*,fmt0)  (STR(I),P(I),I=41,59)
  WRITE(*,*) "Phase setting"
  WRITE(*,*) "============="
  fmt1="(13X,I1,15X,I1,11X,I2,16X,I1,9X,I4)"
  WRITE(*,fmt1) Nphasing,Nboucle,Ifasing,NCYCL
  IF(Ifasing.EQ.1) WRITE(6,"(A7,3X,20(1X,I2))") "PHASE",(Nphase(I),I=1,IDNfase+1)
  IF(Ifasing.EQ.1) WRITE(6,"(A7,21(1X,I2))") "LEVEL",(Level(I), I=0,IDNfase+1)

RETURN
END!
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





