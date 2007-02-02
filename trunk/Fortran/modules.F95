!     Last change:  CF    2 Feb 2007    9:37 am
!
!$modules.f95
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
! Shared parameters and operators
!

MODULE share
!-----------
                                                 !---------------------------------------------------------
                                                 ! VERY IMPORTANT: To be used in conjunction with python,
                                                 ! all static variables must be initialized...
                                                 ! or incorrect results happen!
                                                 !---------------------------------------------------------

  use parameters

  implicit none

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
  REAL*8    :: abundance=1
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
  real*8                 :: dnorm


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


  contains

    subroutine init_operators()
    !--------------------------

      use share

      implicit none

                                                 !loop counters
      integer*8 :: i
      integer*8 :: j

      ns=nint(2.0d0*spins+1.0d0)
      ni=nint(2.0d0*spini+1.0d0)
      nsi=ns*ni
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

    end subroutine init_operators

    subroutine show_operators(pydebug_msg)
    !-------------------------------------
      use share

      implicit none

      EXTERNAL pydebug_msg
                                                 ! local
                                                 !------
      character (len=2) :: str
      character (len=8) :: fmt1

                                                 ! loop counters
                                                 !--------------
      integer*8 :: i,j

      write (str,"(I5)") nsi
      if (allocated(sz)) then

        fmt1="(" // adjustr(str) // "f5.2)"
        message="size of the operator matrices: " // adjustr(str) // "x" // adjustl(str)
        call pydebug_msg()
        print *, "sz"
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

      end if

    end subroutine show_operators

END MODULE operators


