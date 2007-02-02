!     Last change:  CF    1 Feb 2007   10:16 pm
!
!$fourier.f95
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
!  fft routines
!

module fourier
!-------------

  contains

    subroutine fft(np,tr,spcr,spci)
    !------------------------------
      implicit none

      integer*8 :: np
      real*8    :: tr(:),spcr(:),spci(:)

      integer*8 :: np2
      integer*8 :: i, i2

      np2=np/2

      call fast(np,tr)
      spcr=0.0d0
      spci=0.0d0
      do i=1,np2
        i2=2*i
        spcr(i)=tr(i2-1)
        spcr(np+1-i)=tr(i2-1)
 !        spcr(np+2-i)=tr(i2-1)
        spci(i)=tr(i2)
        spci(np+1-i)=-tr(i2)
 !        spci(np+2-i)=-tr(i2)
      end do
      spcr(np2+1)=tr(np+1)
      spci(np2+1)=tr(np+2)
      spcr(1)=spcr(1)/2.0d0
      spci(1)=spci(1)/2.0d0
      spcr(np2+1)=spcr(np2+1)/2.0d0
      spci(np2+1)=spci(np2+1)/2.0d0

    end subroutine fft

    subroutine fast(np,tr)
    !---------------------
      implicit none

      integer*8 :: np
      real*8 :: tr(:)

      integer*8 :: nnn, nn, m
      integer*8 :: n4pow
      integer*8 :: i,k
      real*8 :: t

      m=nint(log(dble(np))/log(2.0d0))
      n4pow=m/2
      nn=1
      if((m-n4pow*2).gt.0) then
      nn=2
      nnn=np/nn
      do k=1,nnn
        t=tr(k)+tr(nnn+k)
        tr(nnn+k)=tr(k)-tr(nnn+k)
        tr(k)=t
      end do
      end if
      do i=1,n4pow
        nn=nn*4
        nnn=np/nn
        call fr4tr(nnn,nn,tr)
      end do
      call ford1(np,tr)
      call ford2(m,tr)
      t=tr(2)
      tr(2)=0.0d0
      tr(np+1)=t
      tr(np+2)=0.0d0
      do i=4,np,2
        tr(i)=-tr(i)
      end do

    end subroutine fast

    subroutine fr4tr(nnn,nn,tr)
    !--------------------------
      implicit none

      integer*8 :: nnn, nn
      real*8    :: tr(:)
      integer*8 :: l(15)
      real*8 :: arg
      real*8 :: p7, pi
      parameter (p7=0.70710678118654752440084436210485d0)
      parameter (pi=3.1415926535897932384626433832795d0)

      integer*8 :: n2,n3, i, j, k
      integer*8 :: l15, l14, l13, l12, l11, l10, l9, l8, l7, l6, l5, l4, l3, l2, l1
      integer*8 :: ji, jl, jr
      real*8 :: piovn
      integer*8 :: j0, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14, jthet
      integer*8 :: k0, k1
      real*8 :: c1,s1,c2,s2,c3,s3
      real*8 :: t0, t1, t2, t3, t4, t5, t6, t7
      real*8 :: xr, xi
      integer*8 :: jlast
      real*8 :: r1, r5

      n2=2*nnn
      n3=3*nnn
      do i=2,15
        l(i)=2
      enddo
      l(1)=nn/4
      if(l(1).le.2) then
        l(1)=2
      else
        k=1
    20  k=k+1
        l(k)=l(k-1)/2
        if(l(k).ne.2) goto 20
      end if
      l15=l(1)
      l14=l(2)
      l13=l(3)
      l12=l(4)
      l11=l(5)
      l10=l(6)
      l9=l(7)
      l8=l(8)
      l7=l(9)
      l6=l(10)
      l5=l(11)
      l4=l(12)
      l3=l(13)
      l2=l(14)
      l1=l(15)
      ji=3
      jl=2
      jr=2
      piovn=pi/dble(nn)
      do j1=2,l1,2
        do j2=j1,l2,l1
          do j3=j2,l3,l2
            do j4=j3,l4,l3
              do j5=j4,l5,l4
                do j6=j5,l6,l5
                  do j7=j6,l7,l6
                    do j8=j7,l8,l7
                      do j9=j8,l9,l8
                        do j10=j9,l10,l9
                          do j11=j10,l11,l10
                            do j12=j11,l12,l11
                              do j13=j12,l13,l12
                                do j14=j13,l14,l13
                                  do jthet=j14,l15,l14
                                    if(jthet.le.2) then
                                      do k=1,nnn
                                        t0=tr(k)+tr(n2+k)
                                        t1=tr(nnn+k)+tr(n3+k)
                                        tr(n2+k)=tr(k)-tr(n2+k)
                                        tr(n3+k)=tr(nnn+k)-tr(n3+k)
                                        tr(k)=t0+t1
                                        tr(nnn+k)=t0-t1
                                      enddo
                                      if(nn.le.4) cycle
                                      k0=4*nnn+1
                                      k1=k0+nnn-1
                                      do k=k0,k1
                                        xr=p7*(tr(nnn+k)-tr(n3+k))
                                        xi=p7*(tr(nnn+k)+tr(n3+k))
                                        tr(n3+k)=tr(n2+k)+xi
                                        tr(nnn+k)=xi-tr(n2+k)
                                        tr(n2+k)=tr(k)-xr
                                        tr(k)=xr+tr(k)
                                      end do
                                    else
                                      arg=dble(jthet-2)*piovn
                                      c1=dcos(arg)
                                      s1=dsin(arg)
                                      c2=c1*c1-s1*s1
                                      s2=c1*s1*2.
                                      c3=c1*c2-s1*s2
                                      s3=c2*s1+s2*c1
                                      j0=4*nnn*jr+1
                                      k0=4*nnn*ji+1
                                      jlast=j0+nnn-1
                                      do j=j0,jlast
                                      k=k0+j-j0
                                        r1=tr(nnn+j)*c1-tr(nnn+k)*s1
                                        r5=tr(nnn+j)*s1+tr(nnn+k)*c1
                                        t2=tr(n2+j)*c2-tr(n2+k)*s2
                                        t6=tr(n2+j)*s2+tr(n2+k)*c2
                                        t3=tr(n3+j)*c3-tr(n3+k)*s3
                                        t7=tr(n3+j)*s3+tr(n3+k)*c3
                                        t0=tr(j)+t2
                                        t4=tr(k)+t6
                                        t2=tr(j)-t2
                                        t6=tr(k)-t6
                                        t1=r1+t3
                                        t5=r5+t7
                                        t3=r1-t3
                                        t7=r5-t7
                                        tr(j)=t0+t1
                                        tr(n3+k)=t4+t5
                                        tr(n2+k)=t0-t1
                                        tr(nnn+j)=t5-t4
                                        tr(n2+j)=t2-t7
                                        tr(nnn+k)=t6+t3
                                        tr(k)=t2+t7
                                        tr(n3+j)=t3-t6
                                      enddo
                                      jr=jr+2
                                      ji=ji-2
                                      if(ji.le.jl) then
                                        ji=2*jr-1
                                        jl=jr
                                      end if
                                    end if
                                  enddo
                                enddo
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

    end subroutine fr4tr

    subroutine ford1(np,tr)
    !----------------------
      implicit none

      integer*8 :: np
      real*8    :: tr(:)

      integer*8 :: k, k1, j
      real*8 :: t

      k=4
      k1=2
      do j=4,np,2
        if(k.gt.j) then
          t=tr(j)
          tr(j)=tr(k)
          tr(k)=t
        end if
        k=k-2
        if(k.le.k1) then
          k=2*j
          k1=j
        end if
      enddo

    end subroutine ford1

    subroutine ford2(m,tr)
    !---------------------
      implicit none

      integer*8 :: m
      real*8    :: tr(:)
      integer*8 :: l(15)

      integer*8 :: i,k, ij, ji
      integer*8 :: l15, l14, l13, l12, l11, l10, l9, l8, l7, l6, l5, l4, l3, l2, l1
      integer*8 :: j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14
      integer*8 :: np
      real*8 :: t

      np=2**m
      do i=2,15
        l(i)=2
      end do
      l(1)=np
      do k=2,m
        l(k)=l(k-1)/2
      end do
      l15=l(1)
      l14=l(2)
      l13=l(3)
      l12=l(4)
      l11=l(5)
      l10=l(6)
      l9=l(7)
      l8=l(8)
      l7=l(9)
      l6=l(10)
      l5=l(11)
      l4=l(12)
      l3=l(13)
      l2=l(14)
      l1=l(15)
      ij=2
      do j1=2,l1,2
        do j2=j1,l2,l1
          do j3=j2,l3,l2
            do j4=j3,l4,l3
              do j5=j4,l5,l4
                do j6=j5,l6,l5
                  do j7=j6,l7,l6
                    do j8=j7,l8,l7
                      do j9=j8,l9,l8
                        do j10=j9,l10,l9
                          do j11=j10,l11,l10
                            do j12=j11,l12,l11
                              do j13=j12,l13,l12
                                do j14=j13,l14,l13
                                  do ji=j14,l15,l14
                                    if(ij.lt.ji) then
                                      t=tr(ij-1)
                                      tr(ij-1)=tr(ji-1)
                                      tr(ji-1)=t
                                      t=tr(ij)
                                      tr(ij)=tr(ji)
                                      tr(ji)=t
                                    end if
                                    ij=ij+2
                                  enddo
                                enddo
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

    end subroutine ford2

end module

