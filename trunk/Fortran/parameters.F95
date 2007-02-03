!     Last change:  CF    3 Feb 2007    5:36 pm
!
!$parameters.f95
!
!author = "Jean-Paul Amoureux <jean-paul.amoureux@univ-lille1.fr>" 
!contributors = "C. Fernandez <christian.fernandez@ensicaen.fr>"
!revision = "1" 
!
!PULSAR Project:
!
! Copyright (C) 2006-2007 Jean-Paul Amoureux, Christian Fernandez
! JPA - Unite de Catalyse et Chimie du Solide, Lille, France.
! CF  - Laboratoire Catalyse et Spectrochimie, Caen, France.
!
!LICENSE:
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License (GPL)
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! To read the license please visit http://www.gnu.org/copyleft/gpl.html
!
!PURPOSE OF THIS FILE:
!
!  "Parameters" Module for Pulsar
!_______________________________________________________________________

module parameters
!----------------
  implicit none

  logical*8 :: verbose=.false.

                                                 ! info debbugage
                                                 !---------------
  logical*8 :: debug=.false.

                                                 !---------------------------------------------------------
                                                 ! VERY IMPORTANT: To be used in conjunction with python,
                                                 ! all static variables must be initialized...
                                                 ! for an unknown reason,
                                                 ! if the variable are not initialized or initialised to zero,
                                                 ! the passing of values between python and fortran doesn't work
                                                 ! (at least on my window system: to check on linux)
                                                 !---------------------------------------------------------

                                                 ! useful constants
                                                 ! ----------------
  real*8 :: pi
  real*8 :: sq6
  real*8 :: eps
  real*8 :: twopi
  real*8 :: twopim6
  CHARACTER(LEN=2)  :: CRLF=CHAR(13)//CHAR(10)
  CHARACTER(LEN=2)  :: TAB=CHAR(09)
  CHARACTER(LEN=32) :: UNDERLINE

  !PARAMETER (CRLF = CHAR(13)//CHAR(10))     ! Not reconised by F2PY
  PARAMETER (UNDERLINE="--------------------------------")

  parameter (pi=3.1415926535897932384626433832795d0)
  parameter (sq6=2.4494897427831780981972840747059d0)
  parameter (eps=1.0d-30)
  parameter (twopi=6.283185307179586476925286766559d0)
  parameter (twopim6=0.000006283185307179586476925286766559d0)


  CHARACTER (LEN=4096) :: message
  CHARACTER (LEN=255)  :: message_fmt


                                                 !limits for array allocations
                                                 !----------------------------
  integer*8 :: nptsmax=16384
  integer*8 :: accuracymax=16

                                                 !spectrometer
                                                 !------------
  real*8 :: spectrometerfield=eps                !spectrometer field
  real*8 :: protonfrequency=eps                  !proton frequency

                                                 !probehead
                                                 !---------
  real*8 :: spinningangle=eps                    !spinner  angle / bo
  real*8 :: spinningspeed=eps                    !spinner frequency
  real*8 :: qfactor=eps                          !probe quality q-factor (f.w.h.m=wo/q_factor)

                                                 !spin_system
                                                 !-----------
  real*8, allocatable :: nucleus(:,:)            !nucleus 'x_n'
  real*8, allocatable :: chemicalshift(:,:)      !iso, csa, eta, euler angles
  real*8, allocatable :: indirect(:,:)           !i, j, jcoupling  (note: jcoupling is a reserved keyword in f2py cannot be used here)
  real*8, allocatable :: dipole(:,:)             !i, j, dipole, euler angles
  real*8, allocatable :: quadrupole(:,:)         !quadrupolar constant, eta and euler angle
  real*8, allocatable :: t2(:,:)                 !t2 homogeneous relaxation time

                                                 !simulation_parameters
                                                 !---------------------
  integer*8 :: accuracy=1                        !accuracy on euler angles (1-16)
  real*8    :: rfstep=eps                        !rf integration step
  integer*8 :: nsb=10                            !number of sidebands
  integer*8 :: npts=512                          !spectrum points number
  real*8    :: sw=eps                            !total spectral width
  integer*8 :: ncycles=-1                        !number of period (pulse+delay) in the sequence
  integer*8 :: nall=999                          !undefined=999, all transitions :0, central: 1, satellite: 2
                                                 !pulses
                                                 !------
  REAL*8, ALLOCATABLE :: pulseangle(:)           !pulse angle
  real*8, allocatable :: pulse(:,:)              !pulse definition
  real*8, allocatable :: delay(:,:)              !delay definition
  real*8, allocatable :: coher(:,:)              !phase selection and coherence jump
  real*8, allocatable :: ctp(:,:)                !pathways selection list

                                                 !final simulation arrays
                                                 !-----------------------
  complex*8, allocatable :: spec(:)              !spectrum
  complex*8, allocatable :: ros(:,:)             !density matrix for the spin s
  complex*8, allocatable :: roi(:,:)             !density matrix for the spin i

                                                 !phases  (TO BE SUPRESSED LATER)
                                                 !------
  integer*8 :: nboucle=-1                        !index of the pulse where a CTP merging occurs
  integer*8 :: ifasing=1                         !flag to say if a phase cycling is used or not (by default phase cycling is used)

  integer*8, allocatable :: nphase(:)            !number of phase for each pulse  (2pi/nphase phase cycling)
  integer*8, allocatable :: level(:)             !selected coherence level after each pulse (during a phase cycling)
  integer*8, allocatable :: qu1(:,:)             !selected coherence in a direct selection in the density matrix
  integer*8, allocatable :: qu2(:,:)             !a second coherence with may be also retained
  integer*8, allocatable :: iref(:)

  LOGICAL*8 :: idealpulse=.true.

  contains

    subroutine controls(pydebug_msg)
    !--------------------------------------------------------------------------------!
    ! this routine is used to display some data from the parameter module and check  !
    ! their validity during the debugging process                                    !
    !--------------------------------------------------------------------------------!
      implicit none

      integer k,i,j
      character*32 ::fmts="(10(g10.2,2x))"
      character*32 ::fmts2="(10(I10,2x))"
      character*32 ::fmt1="(A32,g15.3)"
      character*32 ::fmt2="(A32,I15)"
      EXTERNAL pydebug_msg

      if (.not.debug) return

      message =  "Parameter control:"
      call pydebug_msg()
      message =  UNDERLINE // UNDERLINE
      call pydebug_msg()

      write (message,fmt1) "spectrometerfield", spectrometerfield
      call pydebug_msg()
      write (message,fmt1) "  protonfrequency", protonfrequency
      call pydebug_msg()
      write (message,fmt1) "    spinningangle", spinningangle
      call pydebug_msg()
      write (message,fmt1) "    spinningspeed", spinningspeed
      call pydebug_msg()
      write (message,fmt1) "          qfactor", qfactor
      call pydebug_msg()
      write (message,fmt2) "         accuracy", accuracy
      call pydebug_msg()
      write (message,fmt1) "           rfstep", rfstep
      call pydebug_msg()
      write (message,fmt2) "              nsb", nsb
      call pydebug_msg()
      write (message,fmt2) "             npts", npts
      call pydebug_msg()
      write (message,fmt1) "               sw", sw
      call pydebug_msg()
      write (message,fmt2) "          ncycles", ncycles
      call pydebug_msg()
      write (message,fmt2) "             nall", nall
      call pydebug_msg()
      message=""
      call pydebug_msg()

      if (allocated(nucleus)) then
        message = "Array 'nucleus' : "
        call pydebug_msg()
        do k = 1,size(nucleus,1)
          write(message,fmts)  nucleus(k,1:size(nucleus,2))
          call pydebug_msg()
        enddo
      else
        message = "Array 'nucleus' is not allocated"
        call pydebug_msg()
      endif

      if (allocated(chemicalshift)) then
       message = "Array 'chemicalshift' : "
       call pydebug_msg()
       do k = 1,size(nucleus,1)
          write(message,fmts)  chemicalshift(k,1:size(chemicalshift,2))
          call pydebug_msg()
       enddo
      else
        message = "Array 'chemicalshift' is not allocated"
        call pydebug_msg()
      endif

      if (allocated(indirect)) then
        message = "Array 'indirect' : "
        call pydebug_msg()
        do k = 1,size(indirect,1)
          write(message,fmts)  indirect(k,1:size(indirect,2))
          call pydebug_msg()
        enddo
      else
        message = "Array 'indirect' is not allocated"
        call pydebug_msg()
      endif

      if (allocated(dipole)) then
       message = "Array 'dipole' : "
       call pydebug_msg()
       do k = 1,size(dipole,1)
          write(message,fmts)  dipole(k,1:size(dipole,2))
          call pydebug_msg()
       enddo
      else
        message =  "Array 'dipole' is not allocated"
        call pydebug_msg()
      endif

      if (allocated(quadrupole)) then
        message = "Array 'quadrupole' :"
        call pydebug_msg()
        do k = 1,size(quadrupole,1)
          write(message,fmts)  quadrupole(k,1:size(quadrupole,2))
          call pydebug_msg()
        enddo
      else
        message =  "Array 'quadrupole' is not allocated"
        call pydebug_msg()
      endif

      if (allocated(t2)) then
        message = "Array 't2' : "
        call pydebug_msg()
        do k = 1,size(t2,1)
          write(message,fmts)  t2(k,1:size(t2,2))
          call pydebug_msg()
        enddo
      else
        message =  "Array 't2' is not allocated"
        call pydebug_msg()
      endif

      if (allocated(pulse)) then
        message =  "Array 'pulse' : "
        call pydebug_msg()
        do k = 1,size(pulse,1)
           write(message,fmts)  pulse(k,1:size(pulse,2))
           call pydebug_msg()
        enddo
      else
        message =  "Array 'pulse' is not allocated"
        call pydebug_msg()
      endif

      if (allocated(delay)) then
       message = "Array 'delay' : "
       call pydebug_msg()
       do k = 1,size(delay,1)
          write(message,fmts)  delay(k,1:size(delay,2))
          call pydebug_msg()
       enddo
      else
       message = "Array 'delay' is not allocated"
       call pydebug_msg()
      endif

      if (allocated(coher)) then
        message = "Array 'coherences' : "
        call pydebug_msg()
        do k=1,SIZE(coher,1)
          WRITE(message,fmts2) int(coher(k,1:SIZE(coher,2)))
          call pydebug_msg()
        end do
      else
        message =  "Array 'coher' is not allocated"
        call pydebug_msg()
      end if

      if (allocated(ctp)) then
        message = "Array 'Pathway' (ctp) :"
        call pydebug_msg()
        do k=1,SIZE(ctp,1)
          WRITE(message,fmts2) int(ctp(k,1:SIZE(ctp,2)))
          call pydebug_msg()
        end do
      else
        message = "Array 'ctp' is not allocated"
        call pydebug_msg()
      end if

    end subroutine controls

end module parameters



