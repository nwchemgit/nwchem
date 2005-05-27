c
c $Id: qmmm_utils.F,v 1.8 2005-05-27 17:59:22 marat Exp $
c

      function qmmm_geom_distance(r1,r2)

      implicit none
      double precision r1(3),r2(3)
      double precision qmmm_geom_distance
c 
      integer k
      double precision r

      r=0.0d0
      do k=1,3
        r=r+(r1(k)-r2(k))*(r1(k)-r2(k))
      end do
      r=sqrt(r)

      qmmm_geom_distance=r

      end

      subroutine qmmm_sort(n,a)
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
      do i=1,n-1
       if(a(i).eq.a(i+1)) a(i)=-1 
      end do

      return
 
      end

      function qmmm_map(nr,ar,
     >                     n,a,map)
      implicit none
      integer nr
      integer ar(nr)
      integer n
      integer a(n)
      integer map(n)
      logical qmmm_map
c
c     local variables:
      integer i
      integer ir
      integer im

      qmmm_map = .true.
      im = 0
      do ir=1,nr
        do i=1,n
          if(a(i).eq.ar(ir)) then
            map(i)=ir
            im = im + 1 
          end if
        end do 
        if(im.eq.n) return
      end do

      qmmm_map = .false.
      return
      end

      function qmmm_map1(nr,ar,
     >                     n,a,nm,map)
      implicit none
      integer nr
      integer ar(nr)
      integer n
      integer a(n)
      integer nm
      integer map(n)
      integer qmmm_map1
c
c     local variables:
      integer i
      integer ir

      nm = 0
      do ir=1,nr
        do i=1,n
          if(a(i).eq.ar(ir)) then
            map(i)=ir
            nm = nm + 1 
          end if
        end do 
        if(nm.eq.n) goto 1
      end do

1     qmmm_map1 = nm

      return
      end

      subroutine qmmm_print_bq(nt,c,q)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      double precision c(3,nt)
      double precision q(nt)
      character*16 t

      integer i


      t = "Bq"
      do i=1,nt
       write(6,FMT=9000)
     >           i,t,c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  format(i5,2x,a16,1x,3f15.8,3x,"charge",3x,f15.8)


      end

      subroutine qmmm_print_pdb(nt,myname,c,q,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      double precision c(3,nt)
      double precision q(nt)
      character*16 t(nt)
      character*(*) myname
      character*(nw_max_path_len) filename

      integer i
      integer n
      integer nf
      integer ns

      call util_file_prefix("",filename)
      nf=inp_strlen(filename)-1
      ns=inp_strlen(myname)
      filename = filename(1:nf)//myname(1:ns)

      open(unit=46,status="unknown",form="formatted",file=filename)

      do i=1,nt
       write(46,FMT=9000)
     >           i,t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      call util_flush(46)
      close(46)


      end

      subroutine qmmm_set_ma_char(nt,myname,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      character*16 t(nt)
      character*(*) myname

      integer i

      do i=1,nt
       t(i) = myname
      end do

      end

