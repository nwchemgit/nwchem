c
c $Id$
c

      subroutine smd_util_global_sync()
      implicit none
c
#include "global.fh"      
      call ga_sync(ga_pgroup_get_world())

      return
 
      end

      subroutine smd_util_sort(n,a)
      implicit none
      integer n
      integer a(n)
c
c     local variables:
      integer i
      integer pass  
      integer sorted 
      integer temp
      character*32 pname
      
      pass = 1
      sorted = 0
      do while(sorted .eq. 0) 
        sorted = 1
        do 2 i = 1,n-pass
          if(a(i) .gt. a(i+1)) then
            temp = a(i)
            a(i) = a(i+1)
            a(i+1) = temp
            sorted = 0
          endif
 2      continue
        pass = pass +1
      end do

      return
 
      end

      subroutine smd_print_pdb(nt,myname,c,q,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      double precision c(3,nt)
      double precision q(nt)
      character*(smd_string_size) t(nt)
      character*(*) myname
      character*(nw_max_path_len) filename

      integer i
      integer n
      integer nf
      integer ns
      integer un

      call util_file_prefix(" ",filename)
      nf=inp_strlen(filename)-1
      ns=inp_strlen(myname)
      filename = filename(1:nf)//myname(1:ns)
c
      if(.not.util_get_io_unit(un)) 
     >   call errquit("cannot get file number",0,0)
c

      open(unit=un,status="unknown",form="formatted",file=filename)
      write(*,*) "pdb is in "//filename
      do i=1,nt
       write(un,FMT=9000)
     >           i,t(i),c(1,i),
     >            c(2,i),
     >            c(3,i), q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      call util_flush(un)
      close(un)


      end

      subroutine smd_file_open(un,myname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
#include "util.fh"
#include "inp.fh"
      integer un
      character*(*) myname
c
      if(.not.util_get_io_unit(un)) 
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="unknown",form="formatted",file=myname)

      end

      SUBROUTINE smd_util_rebox_alt(n,latt,rlatt,aaa)

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

      SUBROUTINE smd_util_rebox(n,latt,rlatt,aaa)

      implicit none

      integer n
      double precision rlatt(3,3),latt(3,3)
      double precision  aaa(3,n)
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
     >          write(*,*) "rebox",aaa(1,i),aaa(2,i),aaa(3,i)
       ssx=(rlatt(1,1)*aaa(1,i)+rlatt(1,2)*aaa(2,i)+rlatt(1,3)*aaa(3,i))
       ssy=(rlatt(2,1)*aaa(1,i)+rlatt(2,2)*aaa(2,i)+rlatt(2,3)*aaa(3,i))
       ssz=(rlatt(3,1)*aaa(1,i)+rlatt(3,2)*aaa(2,i)+rlatt(3,3)*aaa(3,i))

       xss=ssx-nint(ssx)
       yss=ssy-nint(ssy)
       zss=ssz-nint(ssz)

       aaa(1,i)=(latt(1,1)*xss+latt(1,2)*yss+latt(1,3)*zss)
       aaa(2,i)=(latt(2,1)*xss+latt(2,2)*yss+latt(2,3)*zss)
       aaa(3,i)=(latt(3,1)*xss+latt(3,2)*yss+latt(3,3)*zss)

      enddo

      return

      END
C>    Folds coordinates back into the box
C>    based on center of coordinates
C>    so that "molecule" remains intact
C>    @param[in] n number of atoms
C>    @param[in] latt(3,3) lattice vectors
C>    @param[in] rlatt(3,3) inverse lattice vectors
C>    @param[in,out] aaa[3,n] coordinates
      SUBROUTINE smd_util_rebox_by_cog(n,latt,rlatt,aaa)

      implicit none

      integer n
      double precision rlatt(3,3),latt(3,3)
      double precision  aaa(3,n)
c
      integer i,j
      integer issx,issy,issz
      double precision  ssx,ssy,ssz,xss,yss,zss
      double precision cog(3),shift(3)
      logical oprint


      cog(1) = 0.0d0
      cog(2) = 0.0d0
      cog(3) = 0.0d0
      do i=1,n
      do j=1,3
         cog(j)=cog(j)+aaa(j,i)
      end do
      end do
c      
      cog(1) = cog(1)/n
      cog(2) = cog(2)/n
      cog(3) = cog(3)/n
c       
      ssx=(rlatt(1,1)*cog(1)+rlatt(1,2)*cog(2)+rlatt(1,3)*cog(3))
      ssy=(rlatt(2,1)*cog(1)+rlatt(2,2)*cog(2)+rlatt(2,3)*cog(3))
      ssz=(rlatt(3,1)*cog(1)+rlatt(3,2)*cog(2)+rlatt(3,3)*cog(3))
c       
      issx=nint(ssx)
      issy=nint(ssy)
      issz=nint(ssz)

      shift(1)=(latt(1,1)*issx+latt(1,2)*issy+latt(1,3)*issz)
      shift(2)=(latt(2,1)*issx+latt(2,2)*issy+latt(2,3)*issz)
      shift(3)=(latt(3,1)*issx+latt(3,2)*issy+latt(3,3)*issz)
c
      do i=1,n

       aaa(1,i)= aaa(1,i)-shift(1)
       aaa(2,i)= aaa(2,i)-shift(2)
       aaa(3,i)= aaa(3,i)-shift(3)

      enddo

      return

      END

C>    Folds coordinates back into the box
C>    while preserving molecules/segments
C>    @param[in] nt number of atoms
C>    @param[in] n number of segments
C>    @param[in] in[n] starting index of given segment in coordinate array
C>    @param[in] nn[n] number of atoms in given segment
C>    @param[in] latt(3,3) lattice vectors
C>    @param[in] rlatt(3,3) inverse lattice vectors
C>    @param[in,out] aaa[3,nt] coordinates

      SUBROUTINE smd_util_rebox_segments(nt,n,in,nn,latt,rlatt,aaa)

      implicit none
      integer nt
      integer n
      integer in(n)
      integer nn(n)
      double precision rlatt(3,3),latt(3,3)
      double precision  aaa(3,nt)
c
      integer i,j
      integer issx,issy,issz
      double precision  ssx,ssy,ssz,xss,yss,zss
      double precision cog(3),shift(3)
      logical oprint
      integer ilo,ihi,ir


      do ir=1,n
        ilo=in(ir)
        ihi=in(ir)+nn(ir)-1
        cog(1) = 0.0d0
        cog(2) = 0.0d0
        cog(3) = 0.0d0
        do i=ilo,ihi
        do j=1,3
           cog(j)=cog(j)+aaa(j,i)
        end do
        end do
c        
        cog(1) = cog(1)/nn(ir)
        cog(2) = cog(2)/nn(ir)
        cog(3) = cog(3)/nn(ir)
c         
        ssx=(rlatt(1,1)*cog(1)+rlatt(1,2)*cog(2)+rlatt(1,3)*cog(3))
        ssy=(rlatt(2,1)*cog(1)+rlatt(2,2)*cog(2)+rlatt(2,3)*cog(3))
        ssz=(rlatt(3,1)*cog(1)+rlatt(3,2)*cog(2)+rlatt(3,3)*cog(3))
c         
        issx=nint(ssx)
        issy=nint(ssy)
        issz=nint(ssz)

        shift(1)=(latt(1,1)*issx+latt(1,2)*issy+latt(1,3)*issz)
        shift(2)=(latt(2,1)*issx+latt(2,2)*issy+latt(2,3)*issz)
        shift(3)=(latt(3,1)*issx+latt(3,2)*issy+latt(3,3)*issz)
c
        do i=ilo,ihi

         aaa(1,i)= aaa(1,i)-shift(1)
         aaa(2,i)= aaa(2,i)-shift(2)
         aaa(3,i)= aaa(3,i)-shift(3)

        enddo
      end do

      return

      END

      subroutine smd_util_print_force_array1(un,na,
     >                           fff)

      implicit none
      integer un
      integer na
      double precision fff(3,na)
c
      integer i
      do i=1,na
         write(un,'(I4,2X,3(E16.8,2X))') i, fff(1,i),fff(2,i),fff(3,i)
      end do
      return

      end

      subroutine smd_util_test_array(na,
     >                           fff1,fff2)

      implicit none
      integer na
      double precision fff1(na,3)
      double precision fff2(3,na)
c
      integer i
      do i=1,na
         if(fff1(i,3).ne.fff2(3,i)) then
           write(*,*) "failed array test"
           stop
         end if
      end do
      write(*,*) "arrays are the same"
      return

      end

      subroutine smd_util_equal_array(na,
     >                           fff1,fff2)

      implicit none
      integer na
      double precision fff1(na,3)
      double precision fff2(3,na)
c
      integer i,k
      do i=1,na
      do k=1,3
         fff2(k,i)=fff1(i,k)  
      end do
      end do
      return

      end

      subroutine smd_util_swap_array(n1,n2,f)
      implicit none
#include "mafdecls.fh"
      integer n1,n2
      double precision f(n1,n2)
c
      integer i_c,h_c

      if(.not.ma_push_get(mt_dbl,n1*n2,'util tmp',h_c,i_c))
     + call errquit('Failed to allocate memory for util tmp',
     + 0, 0)
      call smd_util_swap_array1(n1,n2,f,dbl_mb(i_c))
      if(.not.ma_pop_stack(h_c))
     & call errquit('Failed to deallocate stack h_c',0,
     &       0)
      return

      end
      subroutine smd_util_swap_array1(n1,n2,f1,f2)

      implicit none
      integer n1,n2
      double precision f1(n1,n2)
      double precision f2(n2,n1)
c
      integer i1,i2

      do i1=1,n1
      do i2=1,n2
         f2(n2,n1)=f1(n1,n2)
      end do
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
      subroutine smd_set_ma_to_char1(myname,ilo,ihi,nt,t)
      implicit none
c
      integer nt,ilo,ihi
      character*(32) t(nt)
      character*(*) myname(ihi-ilo+1)

      integer i,j

      do i=1,ihi-ilo+1
        j=ilo+i-1
        myname(i)=t(j)
      end do

      end
      subroutine smd_set_char_to_ma1(myname,ilo,ihi,nt,t)
      implicit none
c
      integer nt,ilo,ihi
      character*(32) t(nt)
      character*(*) myname(ihi-ilo+1)

      integer i,j

      do i=1,ihi-ilo+1
        j=ilo+i-1
        t(j) = myname(i)
      end do

      end
      subroutine smd_set_char_to_ma2(myname,nt,t)
      implicit none
c
      integer nt
      character*(32) t(nt)
      character*(*) myname

      integer i,j

      do i=1,nt
        t(i) = myname
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
      subroutine smd_wipe_ma_char(nt,t)
      implicit none
c
      integer nt
      character*(32) t(nt)

      integer i

      do i=1,nt
        t(i) = " "
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

