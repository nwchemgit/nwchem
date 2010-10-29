c
c $Id$
c

      SUBROUTINE smd_util_rebox(n,latt,rlatt,aaa)

      implicit none

      integer n
      double precision rlatt(3,3),latt(3,3)
      double precision  aaa(n,3)
c
      integer i
      double precision  ssx,ssy,ssz,xss,yss,zss
      logical oprint


      if(n.eq.1) then
       oprint =.true.
      else
       oprint = .false.
      end if
      oprint = .false.
      do i=1,n

       if(oprint) 
     >          write(*,*) "rebox",aaa(i,1),aaa(i,2),aaa(i,3)
       ssx=(rlatt(1,1)*aaa(i,1)+rlatt(1,2)*aaa(i,2)+rlatt(1,3)*aaa(i,3))
       ssy=(rlatt(2,1)*aaa(i,1)+rlatt(2,2)*aaa(i,2)+rlatt(2,3)*aaa(i,3))
       ssz=(rlatt(3,1)*aaa(i,1)+rlatt(3,2)*aaa(i,2)+rlatt(3,3)*aaa(i,3))

       xss=ssx-nint(ssx)
       yss=ssy-nint(ssy)
       zss=ssz-nint(ssz)

       aaa(i,1)=(latt(1,1)*xss+latt(1,2)*yss+latt(1,3)*zss)
       aaa(i,2)=(latt(2,1)*xss+latt(2,2)*yss+latt(2,3)*zss)
       aaa(i,3)=(latt(3,1)*xss+latt(3,2)*yss+latt(3,3)*zss)

      enddo

      return

      END

      subroutine smd_util_print_force_array(un,na,
     >                           fff)

      implicit none
      integer un
      integer na
      double precision fff(na,3)
c
      integer i
      do i=1,na
         write(un,*) i, fff(i,1),fff(i,2),fff(i,3)
      end do
      return

      end

      subroutine smd_util_print_charge_array(un,na,lo,hi,
     >                           fff)

      implicit none
      integer un
      integer na
      double precision fff(na)
c
      integer i,lo,hi
      do i=lo,hi
         write(un,*) i, fff(i)
      end do
      return

      end

      SUBROUTINE tool_randm(iseed,x)

      implicit none

      integer i,iseed,itoz,itozz,mz,mult

      real*8 x,add,dimax,ddimax
c     real*8 rand

      logical newjob

      dimension mz(250)

      save newjob,itoz,dimax,ddimax

      data newjob/.true./

*     **** external functions ****
      real*8   util_random
      external util_random

      if(newjob)then
       if(mod(iseed,2).eq.0)iseed=iseed+1
       mult=65539
       add=2147483648.0d00
       dimax=1.0d00/add
       ddimax=0.50d00*dimax
       do i=1,250
        x=util_random(iseed)
        mz(i)=x*iseed
       enddo
       itoz=1
       newjob=.false.
      else
       itoz=itoz+1
       if(itoz.gt.250)itoz=itoz-250
       itozz=itoz+103
       if(itozz.gt.250)itozz=itozz-250
       mz(itoz)=ieor(mz(itoz),mz(itozz))
       x=mz(itoz)*dimax+ddimax
       x=2.0d00*x
      endif

      return

      END

       FUNCTION erfxc(x)

      implicit none

      double precision erfxc
      real*8 a1,a2,a3,a4,a5,p

      parameter ( a1 = 0.254829592, a2 = -0.284496736 )
      parameter ( a3 = 1.421413741, a4 = -1.453152027 )
      parameter ( a5 = 1.061405429, p  =  0.327591100 )

      real*8 t,x,xsq,tp

      t=1.0/(1.0+p*x)
      xsq=x*x

      tp=t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))

      erfxc=tp*exp(-xsq)

      return

      END

      subroutine smd_ma_char_print(nt,ns,t,un)
      implicit none
#include "smd_const_data.fh"
c
      integer nt,un,ns
      character(1) t(nt*ns)

      integer i,j

      do i=1,nt
          write(un,'(72A1)')
     >      (t(ns*(i-1)+j),j=1,ns)

      end do

      end
      subroutine smd_set_char_to_ma(nt,myname,t)
      implicit none
c
      integer nt
      character*(32) t(nt)
      character*(*) myname(nt)

      integer i

      do i=1,nt
        t(i) = myname(i)
      end do

      end
      subroutine smd_set_ma_to_char(nt,myname,t)
      implicit none
c
      integer nt
      character*(*) t(nt)
      character*(32) myname(nt)

      integer i

      do i=1,nt
        t(i) = " "
        t(i) = myname(i)
      end do

      end
      subroutine smd_set_ma_char(nt,myname,t)
      implicit none
c
      integer nt
      character(1) t(nt)
      character(1) myname(nt)

      integer i

      do i=1,nt
        t(i) = myname(i)
      end do

      end

      subroutine smd_ma_char_count(s,t,nt)
      implicit none
c
#include "inp.fh"
      integer nt
      character*(*) s
      character*(*) t

      integer istart,iend

      istart = 0
      nt = 0
  10  if (inp_strtok(t,s, istart, iend)) then
      nt = nt +1
      goto 10
      endif

      end

