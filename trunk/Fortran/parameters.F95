!     Last change:  CF   14 Feb 2007    6:30 pm
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
!  Module for Pulsar
!
! modules:
!   parameters
!   operators
!   share
!
! subroutines:
!   show_operators
!   reset_operators
!   control_parameters
!   reset_parameters
!_______________________________________________________________________

module parameters
!----------------
  implicit none

  logical*8 :: verbose=.false.

                                                 ! info debbugage
                                                 !---------------
  logical*8 :: debug=.false.

                                                 !---------------------------------------------------------
                                                 ! very important: to be used in conjunction with python,
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
  character(len=2)  :: crlf=char(13)//char(10)
  character(len=2)  :: tab=char(09)
  character(len=32) :: underline

  !parameter (crlf = char(13)//char(10))     ! not reconised by f2py
  parameter (underline="----------------------------------------------")

  parameter (pi=3.1415926535897932384626433832795d0)
  parameter (sq6=2.4494897427831780981972840747059d0)
  parameter (eps=1.0d-30)
  parameter (twopi=6.283185307179586476925286766559d0)
  parameter (twopim6=0.000006283185307179586476925286766559d0)


  character (len=32768) :: message
  character (len=255)  :: message_fmt


                                                 !limits for array allocations
                                                 !----------------------------
  integer*8 :: nptsmax=16384
  integer*8 :: accuracymax=16

                                                 !spectrometer
                                                 !------------
  real*8 :: spectrometerfield=2.3d0              !spectrometer field
  real*8 :: protonfrequency=100.d0               !proton frequency

                                                 !probehead
                                                 !---------
  real*8 :: spinningangle=eps                    !spinner  angle / bo
  real*8 :: spinningspeed=eps                    !spinner frequency
  real*8 :: qfactor=0.001                        !probe quality q-factor (f.w.h.m=wo/q_factor)

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
  real*8    :: rfstep=5.d0                        !rf integration step
  integer*8 :: nsb=10                            !number of sidebands
  integer*8 :: npts=512                          !spectrum points number
  real*8    :: sw=10000.d0                       !total spectral width
  integer*8 :: ncycles=1                         !number of period (pulse+delay) in the sequence
  integer*8 :: nall=999                          !undefined=999, all transitions :0, central: 1, satellite: 2
                                                 !pulses
                                                 !------
  real*8, allocatable :: pulseangle(:)           !pulse angle
  real*8, allocatable :: pulse(:,:)              !pulse definition
  real*8, allocatable :: delay(:,:)              !delay definition
  real*8, allocatable :: coher(:,:)              !phase selection and coherence jump
  real*8, allocatable :: ctp(:,:)                !pathways selection list

                                                 !final simulation arrays
                                                 !-----------------------
  complex*8, allocatable :: spec(:)              !spectrum
  complex*8, allocatable :: ros(:,:)             !density matrix for the spin s
  complex*8, allocatable :: roi(:,:)             !density matrix for the spin i

                                                 !phases  (to be supressed later)
                                                 !------
  integer*8 :: nboucle=-1                        !index of the pulse where a ctp merging occurs
  integer*8 :: ifasing=1                         !flag to say if a phase cycling is used or not (by default phase cycling is used)

  integer*8, allocatable :: nphase(:)            !number of phase for each pulse  (2pi/nphase phase cycling)
  integer*8, allocatable :: level(:)             !selected coherence level after each pulse (during a phase cycling)
  integer*8, allocatable :: qu1(:,:)             !selected coherence in a direct selection in the density matrix
  integer*8, allocatable :: qu2(:,:)             !a second coherence with may be also retained
  integer*8, allocatable :: iref(:)

  logical*8 :: idealpulse=.true.
  LOGICAL*8 :: zerospectrum=.TRUE.                         ! the spectrum memo
  REAL*8 		:: rcph = 0.0d0                                ! receiver phase

end module parameters

!------------------------
MODULE share
!------------------------
! other shared parameters
!------------------------
  USE parameters


  integer*8 :: nt=-1
  integer*8 :: ng=-1
  logical*2 :: staticsample=.false.              ! (replace ispeed in the original version)
  real*8    :: dg=eps

  real*8    :: de(-4:4,-4:4)=eps                 ! wigner matrix
  real*8    :: dk(2,0:4,-2:2)=eps                !

  real*8    :: spins=eps                         ! spin i
  real*8    :: spini=eps                         ! spin s

  real*8    :: vls=eps                           ! larmor of spin s
  real*8    :: cqs=eps
  REAL*8    :: abundance=1.0d0
  real*8    :: etaqs=eps
  real*8    :: isos=eps
  real*8    :: csas=eps
  real*8    :: etacs=eps
  real*8    :: t2s=eps
  real*8    :: qq=eps                            ! effect of quality factor
  real*8    :: wm(4)=eps
  real*8    :: fst=eps
  integer*8 :: icoupled=-1                       ! if i exist
  real*8    :: vli=eps                           ! larmor of spin i
  real*8    :: cqi=eps
  real*8    :: etaqi=eps
  real*8    :: isoi=eps
  real*8    :: csai=eps
  real*8    :: etaci=eps
  real*8    :: t2i=eps
  real*8    :: jsi=eps                           ! j coupling between s and i
  real*8    :: dsi=eps                           ! dipolar coupling i-s
  integer*8 :: kcoupled=-1                       ! number of coupled k nucleus to s
  real*8    :: coefps=eps
  real*8    :: coefpi=eps
  real*8    :: coefss=eps
  real*8    :: coefsi=eps
  real*8    :: isost=eps

  complex*8 :: r2mdsi(0:2)                       ! r2m component for is dipolar tensor
  complex*8 :: r2mqi(0:2)
  complex*8 :: r2mcs(0:2)
  complex*8 :: r2mci(0:2)

  integer*8 :: ns=-1
  integer*8 :: ni=-1
  integer*8 :: nsi=-1
  real*8    :: fnexp=eps

  real*8 :: vvi(2)
  real*8 :: vvs(2)
  real*8 :: vvii(2)
  real*8 :: vvis(2)
  complex*8 :: ats(-4:4)                   ! arts and aits
  complex*8 :: aqs(-4:4)                  ! arqs(10,-4:4) and aiqs(10,-4:4)
  complex*8 :: acs(-4:4)                  ! arcs and aics
  complex*8 :: ads(-4:4)                     ! ads(4,:) correspondant  dipolaire i-s
  complex*8 :: aqi(-4:4)                     ! arqi and aiqi
  complex*8 :: aci(-4:4)                     ! arci and aici
  complex*8 :: aqpi(0:4)
  complex*8 :: aqps(0:4)
  complex*8 :: adcps(0:4)
  complex*8 :: aqsi(2,0:4)
  complex*8 :: aqss(2,0:4)
  complex*8 :: acspi(0:4)
  complex*8 :: adipolis(0:4)

  real*8 :: zl(0:4)=(/1.0d0,2.0d0,2.0d0,2.0d0,2.0d0/)

  real*8    :: qsangle(3)
  real*8    :: csangle(3)
  real*8    :: qiangle(3)
  real*8    :: ciangle(3)
  real*8    :: dsiangle(3)                   ! dipolar orientation

  real*8, allocatable    :: spink(:)             ! non -excited and non-observed nuclei
  real*8, allocatable    :: vlk(:)               ! larmor of spins k
  real*8, allocatable    :: jsk(:)               ! j coupling between s-k
  real*8, allocatable    :: dsk(:)               ! dipolar coupling s-k
  real*8, allocatable    :: dskangle(:,:)        ! dipolar orientation of s-k

  real*8, allocatable    :: rfpowers(:)
  real*8, allocatable    :: rfoffsets(:)
  real*8, allocatable    :: rfphases(:)
  real*8, allocatable    :: rfpoweri(:)
  real*8, allocatable    :: rfoffseti(:)
  real*8, allocatable    :: rfphasei(:)
  real*8, allocatable    :: rflength(:)
  real*8, allocatable    :: rfdelay(:)
  integer*8, allocatable :: decouple(:)

  integer*8, allocatable :: itour(:)
  integer*8, allocatable :: ntp(:)
  real*8, allocatable    :: deltatp(:)
  integer*8, allocatable :: ntpadd(:)
  real*8, allocatable    :: deltatpadd(:)
  real*8, allocatable    :: t2ss(:)
  real*8, allocatable    :: t2ii(:)
  real*8, allocatable    :: t2si(:)

  complex*8, allocatable :: rott(:,:)            ! rortt and roitt
  complex*8, allocatable :: qt(:,:,:,:)          ! qrt and qit(0:40,0:80,11,-20:20)

  real*8, allocatable    :: u1s(:)
  real*8, allocatable    :: u2s(:)

  complex*8, allocatable :: r2mdsk(:,:)          ! tensor components ks
  integer*8, allocatable :: nk(:)

  real*8, allocatable    :: dipk(:)
  integer*8, allocatable :: mk(:)

  real*8, allocatable    :: wr(:,:,:)
  complex*8, allocatable :: amp(:,:,:,:)
  real*8                 :: dnorm=eps


  complex*8, allocatable :: adks(:,:)            ! ards and aids   correspondant au dipolaire k-s
  complex*8, allocatable :: a(:,:)               ! ab(11,4),bb(11,4)
  complex*8, allocatable :: p(:,:)               ! pr(11,-20:20) and pi(11,-20:20)
  complex*8, allocatable :: u(:,:)               ! tr and ti
  complex*8, allocatable :: ui(:,:)              ! tri(100,100),tii(100,100)
  complex*8, allocatable :: uf(:,:)              ! trf(100,100),tif(100,100)
  complex*8, allocatable :: usi(:,:)             ! trsi(100,100),tisi(100,100)
  complex*8, allocatable :: usf(:,:)             ! trsf(100,100),tisf(100,100)
  complex*8, allocatable :: uint(:,:)            ! trint(100,100),tiint(100,100)
  complex*8, allocatable :: ro(:,:)              ! ror and roi
  complex*8, allocatable :: rot(:,:)             ! rotr and roti
  complex*8, allocatable :: y(:,:)               ! yr and yi
  complex*8, allocatable :: htst(:,:)            ! htsrt(100,100) and htsit(100,100)
  complex*8, allocatable :: z(:,:)               ! zr and zi
  real*8, allocatable    :: vp(:)
  complex*8, allocatable :: expvp(:,:)           ! exp(ivp)
  real*8, allocatable    :: htsr(:,:)            ! htsr(100,100)
  complex*8, allocatable :: sti(:,:)             ! stri(11,11),stii(11,11)
  complex*8, allocatable :: rosto(:,:)           !
  complex*8, allocatable :: q(:,:)               ! qr(11,-20:20) and qi(11,-20:20)
  real*8, allocatable    :: phasei(:)

  COMPLEX*8, ALLOCATABLE :: roctp(:,:,:)
  INTEGER*8 :: nctp=-1

END MODULE share

MODULE operators
!---------------
                                                 ! todo: make this module able to handle more than two nuclei
                                                 !
  implicit none

  real*8, allocatable    :: unity(:,:)             ! identity matrix

  real*8, allocatable    :: s2(:,:)              ! s^2
  real*8, allocatable    :: sx(:,:)              ! sx
  real*8, allocatable    :: sy(:,:)              ! sy
  real*8, allocatable    :: sz(:,:)              ! sz
  real*8, allocatable    :: sq(:,:)              ! 3sz^2-s^2
  real*8, allocatable    :: sq1(:,:)             ! sz(4s^2-8sz^2-1)/2
  real*8, allocatable    :: sq2(:,:)             ! sz(2s^2-2sz^2-1)/2
  real*8, allocatable    :: sp(:,:)              ! s+
  real*8, allocatable    :: sm(:,:)              ! s-

  real*8, allocatable    :: i2(:,:)              ! i^2
  real*8, allocatable    :: ix(:,:)              ! ix
  real*8, allocatable    :: iy(:,:)              ! iy
  real*8, allocatable    :: iz(:,:)              ! iz
  real*8, allocatable    :: iq(:,:)              ! 3iz^2-i^2
  real*8, allocatable    :: iq1(:,:)             ! iz(4i^2-8iz^2-1)/2
  real*8, allocatable    :: iq2(:,:)             ! iz(2i^2-2iz^2-1)/2
  real*8, allocatable    :: ip(:,:)              ! i+
  real*8, allocatable    :: im(:,:)              ! i-

  real*8, allocatable    :: sziz(:,:)            ! product iz*sz

END MODULE

subroutine reset_operators()
!--------------------------

use share
USE operators

implicit none

                                             !loop counters
  integer*8 :: i
  integer*8 :: j

  ! first deallocate the operators array if they are still allocated
  if (allocated(sz)) deallocate (sz)
  if (allocated(s2)) deallocate (s2)
  if (allocated(sx)) deallocate (sx)
  if (allocated(sy)) deallocate (sy)
  if (allocated(sq)) deallocate (sq)
  if (allocated(sq1)) deallocate (sq1)
  if (allocated(sq2)) deallocate (sq2)

  if (allocated(iz)) deallocate (iz)
  if (allocated(i2)) deallocate (i2)
  if (allocated(ix)) deallocate (ix)
  if (allocated(iy)) deallocate (iy)
  if (allocated(iq)) deallocate (iq)
  if (allocated(iq1)) deallocate (iq1)
  if (allocated(iq2)) deallocate (iq2)

  if (allocated(sziz)) deallocate (sziz)
  if (allocated(unity)) deallocate (unity)

  if (allocated(sp)) deallocate (sp)
  if (allocated(sp)) deallocate (sm)
  if (allocated(im)) deallocate (im)
  if (allocated(ip)) deallocate (ip)
  ns=nint(2.0d0*spins+1.0d0)
  ni=nint(2.0d0*spini+1.0d0)
  nsi=ns*ni

  allocate(unity(nsi,nsi))
  unity=0.0d0
  do j=1,ni
   do i=1,ns
     unity(i+(j-1)*ns,i+(j-1)*ns)=1.0d0
   end do
  end do

  allocate(sz(nsi,nsi))
  allocate(s2(nsi,nsi))
  allocate(iz(nsi,nsi))
  allocate(i2(nsi,nsi))
  sz=0.0d0
  s2=0.0d0
  iz=0.0d0
  i2=0.0d0
  do j=1,ni
    do i=1,ns
      s2(i+(j-1)*ns,i+(j-1)*ns)=spins*(spins+1.0d0)
      sz(i+(j-1)*ns,i+(j-1)*ns)=spins-dble(i-1)
      i2(i+(j-1)*ns,i+(j-1)*ns)=spini*(spini+1.0d0)
      iz(i+(j-1)*ns,i+(j-1)*ns)=spini-dble(j-1)
    end do
  end do
  allocate(sp(nsi,nsi))
  allocate(sm(nsi,nsi))
  allocate(sx(nsi,nsi))
  allocate(sy(nsi,nsi))
  allocate(sq(nsi,nsi))
  allocate(sq1(nsi,nsi))
  allocate(sq2(nsi,nsi))
  sp=0.0D0
  sm=0.0D0
  sp=cshift(sqrt(s2-matmul(sz,sz-unity)),-1,2) ! circular shift of unity row or one column to create s+ or s-
  sm=transpose(sp)
  sx=(sp+sm)/2.0d0
  sy=-(sp-sm)/2.0d0
  sq=3.0d0*matmul(sz,sz)-s2
  sq1=matmul(sz,4.0d0*s2-8.0d0*matmul(sz,sz)-unity)/2.0d0
  sq2=matmul(sz,2.0d0*s2-2.0d0*matmul(sz,sz)-unity)/2.0d0
  deallocate(sp)
  deallocate(sm)

  allocate(ip(nsi,nsi))
  allocate(im(nsi,nsi))
  allocate(ix(nsi,nsi))
  allocate(iy(nsi,nsi))
  allocate(iq(nsi,nsi))
  allocate(iq1(nsi,nsi))
  allocate(iq2(nsi,nsi))
  ip=0.0D0
  im=0.0D0
  ip=cshift(sqrt(i2-matmul(iz,iz-unity)),-ns,2)! circular shift of ns row or ns column to create i+ or i-
  im=transpose(ip)
  ix=(ip+im)/2.0d0
  iy=-(ip-im)/2.0d0
  iq=3.0d0*matmul(iz,iz)-i2
  iq1=matmul(iz,4.0d0*i2-8.0d0*matmul(iz,iz)-unity)/2.0d0
  iq2=matmul(iz,2.0d0*i2-2.0d0*matmul(iz,iz)-unity)/2.0d0
  deallocate(ip)
  deallocate(im)

  allocate (sziz(nsi,nsi))
  sziz=0.0D0
  sziz=matmul(sz,iz)

end subroutine reset_operators


!--------------------------------------------
subroutine show_operators(pydebug_msg)
!--------------------------------------------
! shows the initials operators in python log
!--------------------------------------------

use share, ONLY: message,nsi
USE operators

implicit none

 EXTERNAL pydebug_msg
                                             ! local
                                             !------
  character (len=3) :: str
  character (len=32) :: fmt1

                                             ! loop counters
                                             !--------------
  integer*8 :: i,j

  write (str,"(I3)") nsi
  if (allocated(sz)) then

    fmt1="(" // TRIM(adjustl(str)) // "(f10.4,1x))"
    message="size of the operator matrices: " // TRIM(adjustr(str)) // "x" // TRIM(adjustl(str))
    call pydebug_msg()
    message="sz"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (sz(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "sx"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (sx(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message ="isy"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (sy(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "sq"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (sq(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "sq1"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (sq1(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "sq2"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (sq2(j,i),i=1,nsi)
      call pydebug_msg()
    end do

    message= "----"
    call pydebug_msg()

  end if

  if (allocated(iz)) then

    message= "iz"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (iz(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "ix"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (ix(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "iiy"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (iy(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "iq"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (iq(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "iq1"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (iq1(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "iq2"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (iq2(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "------"
    call pydebug_msg()
    message= "sziz"
    call pydebug_msg()
    do j=1,nsi
      write( message, fmt1)  (sziz(j,i),i=1,nsi)
      call pydebug_msg()
    end do
    message= "------"
    call pydebug_msg()

  end if

end subroutine show_operators


!-------------------------------
subroutine control_parameters(pydebug_msg)
!--------------------------------------------------------------------------------!
! this routine is used to display some data from the parameter module and check  !
! their validity during the debugging process                                    !
!--------------------------------------------------------------------------------!
  USE parameters

  implicit none

  integer i,j,k
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
    i=size(nucleus,1)
    j=size(nucleus,2)
    do k = 1,i
      write(message,fmts)  nucleus(k,1:j)
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

end subroutine

!-------------------------------------------------
subroutine reset_parameters()
!-------------------------------------------------
! Initialisation of the main parameters and arrays
!-------------------------------------------------
  USE parameters

  verbose=.false.
  debug=.false.
  message=''
  message_fmt=''
  nptsmax=16384
  accuracymax=16
  spectrometerfield=2.3d0
  protonfrequency=100.0d0
  spinningangle=0.0d0
  spinningspeed=0.0d0
  qfactor=0.001
  if (ALLOCATED(nucleus)) DEALLOCATE(nucleus)
  if (ALLOCATED(chemicalshift)) deallocate (chemicalshift)
  if (ALLOCATED(indirect)) DEALLOCATE (indirect)
  if (allocated(dipole)) deallocate (dipole)
  if (ALLOCATED(quadrupole)) DEALLOCATE (quadrupole)
  if (ALLOCATED(t2)) deallocate (t2)
  accuracy=1
  rfstep=5.0
  nsb=0.0D0
  npts=512
  sw=10000.0
  ncycles=1
  nall=0
  if (allocated (pulseangle)) DEALLOCATE(pulseangle)
  if (allocated (pulse)) DEALLOCATE(pulse)
  if (ALLOCATED (delay)) DEALLOCATE(delay)
  if (allocated (coher)) DEALLOCATE(coher)
  if (allocated (ctp)) DEALLOCATE(ctp)
  if (allocated (spec)) DEALLOCATE (spec)
  if (allocated (ros)) DEALLOCATE (ros)
  if (allocated (roi)) DEALLOCATE (roi)
  nboucle=-1
  ifasing=1
  if (allocated (nphase)) DEALLOCATE (nphase)
  if (allocated (level)) DEALLOCATE (level)
  if (allocated (qu1)) deallocate (qu1)
  if (allocated (qu2)) deallocate (qu2)
  if (allocated (iref)) deallocate (iref)
  idealpulse=.true.

end subroutine reset_parameters

!-----------------------
subroutine reset_share()
!-----------------------

USE share

      nt=-1
      ng=-1
      staticsample=.false.              ! (replace ispeed in the original version)
      dg=eps
      de=eps                            ! wigner matrix
      dk=eps                            !
      spins=eps                         ! spin i
      spini=eps                         ! spin s

      vls=eps                           ! larmor of spin s
      cqs=eps
      abundance=1.0d0
      etaqs=eps
      isos=eps
      csas=eps
      etacs=eps
      t2s=eps
      qq=eps                            ! effect of quality factor
      wm=eps
      fst=eps
      icoupled=-1                       ! if i exist
      vli=eps                           ! larmor of spin i
      cqi=eps
      etaqi=eps
      isoi=eps
      csai=eps
      etaci=eps
      t2i=eps
      jsi=eps                           ! j coupling between s and i
      dsi=eps                           ! dipolar coupling i-s
      kcoupled=-1                       ! number of coupled k nucleus to s
      coefps=eps
      coefpi=eps
      coefss=eps
      coefsi=eps
      isost=eps

      r2mdsi=CMPLX(0.0d0,0.0d0)                       ! r2m component for is dipolar tensor
      r2mqi=CMPLX(0.0d0,0.0d0)
      r2mcs=CMPLX(0.0d0,0.0d0)
      r2mci=CMPLX(0.0d0,0.0d0)

      ns=-1
      ni=-1
      nsi=-1
      fnexp=eps

      vvi=0.0d0
      vvs=0.0d0
      vvii=0.0d0
      vvis=0.0D0
      ats=CMPLX(0.0d0,0.0d0)                   ! arts and aits
      aqs=CMPLX(0.0d0,0.0d0)                  ! arqs(10,-4:4) and aiqs(10,-4:4)
      acs=CMPLX(0.0d0,0.0d0)                  ! arcs and aics
      ads=CMPLX(0.0d0,0.0d0)                     ! ads(4,:) correspondant  dipolaire i-s
      aqi=CMPLX(0.0d0,0.0d0)                     ! arqi and aiqi
      aci=CMPLX(0.0d0,0.0d0)                     ! arci and aici
      aqpi=CMPLX(0.0d0,0.0d0)
      aqps=CMPLX(0.0d0,0.0d0)
      adcps=CMPLX(0.0d0,0.0d0)
      aqsi=CMPLX(0.0d0,0.0d0)
      aqss=CMPLX(0.0d0,0.0d0)
      acspi=CMPLX(0.0d0,0.0d0)
      adipolis=CMPLX(0.0d0,0.0d0)

      zl(0:4)=(/1.0d0,2.0d0,2.0d0,2.0d0,2.0d0/)

      qsangle=0.0d0
      csangle=0.0d0
      qiangle=0.0d0
      ciangle=0.0d0
      dsiangle=0.0d0

      if (ALLOCATED(spink)) DEALLOCATE (spink)
      if (ALLOCATED(vlk)) DEALLOCATE (vlk)
      if (ALLOCATED(jsk)) DEALLOCATE (jsk)
      if (ALLOCATED(dsk)) deallocate (dsk)
      if (ALLOCATED(dskangle)) deallocate (dskangle)

      if (ALLOCATED(rfpowers)) deallocate (rfpowers)
      if (ALLOCATED(rfoffsets)) deallocate (rfoffsets)
      if (ALLOCATED(rfphases)) deallocate (rfphases)
      if (ALLOCATED(rfpoweri)) deallocate (rfpoweri)
      if (ALLOCATED(rfoffseti)) deallocate (rfoffseti)
      if (ALLOCATED(rfphasei)) deallocate (rfphasei)
      if (ALLOCATED(rflength)) deallocate (rflength)
      if (ALLOCATED(rfdelay)) deallocate (rfdelay)
      if (ALLOCATED(decouple)) deallocate (decouple)

      if (ALLOCATED(itour)) deallocate (itour)
      if (ALLOCATED(ntp)) deallocate (ntp)
      if (ALLOCATED(deltatp)) deallocate (deltatp)
      if (ALLOCATED(ntpadd)) deallocate (ntpadd)
      if (ALLOCATED(deltatpadd)) deallocate (deltatpadd)
      if (ALLOCATED(t2ss)) deallocate (t2ss)
      if (ALLOCATED(t2ii)) deallocate (t2ii)
      if (ALLOCATED(t2si)) deallocate (t2si)

      if (ALLOCATED(rott)) DEALLOCATE(rott)
      if (ALLOCATED(qt)) DEALLOCATE (qt)

      if (ALLOCATED(u1s)) deallocate (u1s)
      if (ALLOCATED(u2s)) deallocate (u2s)

      if (ALLOCATED(r2mdsk)) DEALLOCATE (r2mdsk)          ! tensor components ks
      if (ALLOCATED(nk)) deallocate (nk)

      if (ALLOCATED(dipk)) DEALLOCATE(dipk)
      if (ALLOCATED(mk)) DEALLOCATE(mk)

      if (ALLOCATED(wr)) deallocate (wr)
      if (ALLOCATED(amp)) deallocate (amp)
      dnorm=0.0d0


      if (ALLOCATED(adks)) deallocate (adks)
      if (ALLOCATED(a)) deallocate (a)
      if (ALLOCATED(p)) deallocate (p)
      if (ALLOCATED(u)) deallocate (u)
      if (ALLOCATED(ui)) deallocate (ui)
      if (ALLOCATED(uf)) deallocate (uf)
      if (ALLOCATED(usi)) deallocate (usi)
      if (ALLOCATED(usf)) deallocate (usf)
      if (ALLOCATED(uint)) deallocate (uint)
      if (ALLOCATED(ro)) deallocate (ro)
      if (ALLOCATED(rot)) deallocate (rot)
      if (ALLOCATED(y)) deallocate (y)
      if (ALLOCATED(htst)) deallocate (htst)
      if (ALLOCATED(z)) deallocate (z)
      if (ALLOCATED(vp)) deallocate (vp)
      if (ALLOCATED( expvp)) deallocate (expvp)
      if (ALLOCATED(htsr)) deallocate (htsr)
      if (ALLOCATED(sti)) deallocate (sti)
      if (ALLOCATED(rosto)) deallocate (rosto)
      if (ALLOCATED(q)) deallocate (q)
      if (ALLOCATED(phasei)) deallocate (phasei)
      if (ALLOCATED(roctp)) DEALLOCATE(roctp)
      nctp=-1

end subroutine

subroutine reset()
!--------------------------------------
! Reset and deallocation of all tables
!--------------------------------------
      USE diagonalize
      call dealloc()
      call reset_share()
      call reset_operators()
      call reset_parameters()

end subroutine reset







