!     Last change:  CF   20 Feb 2006    5:32 am
! File: $Id: Module_common.f95,v 1.0
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
! Original Author of file: JP Amoureux and C Fernandez
! Purpose of file:
!     Common declaration block for pulsar
! ----------------------------------------------------------------------

module common_module
! ******************

  implicit none

      LOGICAL*4 :: Verbose

      REAL :: SPX(32800)
      REAL :: SPY(32800)
      REAL :: WR(0:40,0:80,11)
      REAL :: RorTt(100,100)
      REAL :: RoiTt(100,100)
      REAL :: AMPX(0:40,0:80,11,-20:20)
      REAL :: AMPY(0:40,0:80,11,-20:20)
      REAL :: QRT(0:40,0:80,11,-20:20)
      REAL :: QIT(0:40,0:80,11,-20:20)
      REAL :: FFCs(10,0:2)
      REAL :: GGCs(10,0:2)
      REAL :: FFDs(10,0:2)
      REAL :: GGDs(10,0:2)
      REAL :: FFCi(10,0:2)
      REAL :: GGCi(10,0:2)
      REAL :: FFQi(10,0:2)
      REAL :: GGQi(10,0:2)
      REAL :: A23Ds(10)
      REAL :: WdipIs(10)
      REAL :: Wdip(10)
      REAL :: DJ(10)
      REAL :: WM(4)
      REAL :: DK(2,0:4,-2:2)
      REAL :: DE(-4:4,-4:4)
      REAL :: U1s(11)
      REAL :: U2s(11)
      REAL :: PULSEs(9999,4)
      REAL :: PULSEi(9999,4)
      REAL :: T2ss(9999)
      REAL :: T2ii(9999)
      REAL :: DELAY(9999)
      REAL :: DTP(9999)
      REAL :: DTP1(9999)
      REAL :: T2si(9999)
      REAL :: Decouple(9999)
      INTEGER :: NTP(9999)
      INTEGER :: NTPint(9999)
      INTEGER :: ITOUR(9999)
      INTEGER :: IQ1(9999,2)
      INTEGER :: IQ2(9999,2)
      REAL :: PRIXs(100,100)
      REAL :: PIIYs(100,100)
      REAL :: PRIZs(100,100)
      REAL :: PRQZs(100,100)
      REAL :: PRQS1s(100,100)
      REAL :: PRQS2s(100,100)
      REAL :: PRIXi(100,100)
      REAL :: PIIYi(100,100)
      REAL :: PRIZi(100,100)
      REAL :: PRD(100,100)
      REAL :: PRQZi(100,100)
      REAL :: PRQS1i(100,100)
      REAL :: PRQS2i(100,100)

      INTEGER :: NCYCL
      INTEGER :: Nphasing
      INTEGER :: Nboucle
      INTEGER :: Ns
      INTEGER :: Ni
      INTEGER :: Nsi
      INTEGER :: NT
      INTEGER :: NPTS
      INTEGER :: NALL
      INTEGER :: NU
      INTEGER :: ISPEED
      INTEGER :: NSB
      INTEGER :: NATsd
      REAL :: COEFPs
      REAL :: COEFSs
      REAL :: ETAQs
      REAL :: DISOst
      REAL :: DELTACs
      REAL :: ETACs
      REAL :: A123cs
      REAL :: COEFPi
      REAL :: COEFSi
      REAL :: ETAQi
      REAL :: DISOi
      REAL :: DELTACi
      REAL :: ETACi
      REAL :: A123ci
      REAL :: A123qi
      REAL :: SPINs
      REAL :: SPINi
      REAL :: FST
      REAL :: QQ
      REAL :: DNORM
      REAL :: WROT
      REAL :: DU
      REAL :: WLs
      REAL :: WLi
      REAL :: Pfit


end module common_module













































