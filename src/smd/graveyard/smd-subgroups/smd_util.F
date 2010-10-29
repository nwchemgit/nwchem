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
