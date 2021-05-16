c
c $Id$
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

      subroutine qmmm_sort_unique(nr,ar,
     >                     n,a)
      implicit none
      integer nr
      integer ar(nr)
      integer n
      integer a(nr)
c
c     local variables:
c
      integer i
      integer ir
      integer im

      n = 0
      do ir=1,nr
        do i=1,im
          if(ar(ir).eq.a(i)) goto 1
        end do 
        n = n + 1
        a(n) = ar(ir)
1       continue
      end do

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
      integer un

      call util_file_prefix(" ",filename)
      nf=inp_strlen(filename)-1
      ns=inp_strlen(myname)
      filename = filename(1:nf)//myname(1:ns)
c
      if(.not.qmmm_get_io_unit(un)) 
     >   call errquit("cannot get file number",0,0)
c

      open(unit=un,status="unknown",form="formatted",file=filename)

      do i=1,nt
       write(un,FMT=9000)
     >           i,t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      call util_flush(un)
      close(un)


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

      subroutine qmmm_print_pdb_bq(nt,myname,c,q)
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
      character*(*) myname
      character*(nw_max_path_len) filename

      integer i
      integer n
      integer nf
      integer ns

      call util_file_prefix(" ",filename)
      nf=inp_strlen(filename)-1
      ns=inp_strlen(myname)
      filename = filename(1:nf)//myname(1:ns)

      open(unit=46,status="unknown",form="formatted",file=filename)

      t = "Bq"
      do i=1,nt
       write(46,FMT=9000)
     >           i,t,c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      call util_flush(46)
      close(46)


      end

      subroutine qmmm_print_pdb_forces(nt,myname,c,f)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      double precision c(3,nt)
      double precision f(3,nt)
      character*16 t
      character*(*) myname
      character*(nw_max_path_len) filename

      integer i
      integer n
      integer nf
      integer ns

      call util_file_prefix(" ",filename)
      nf=inp_strlen(filename)-1
      ns=inp_strlen(myname)
      filename = filename(1:nf)//myname(1:ns)

      open(unit=46,status="unknown",form="formatted",file=filename)

      t = "Bq"
      do i=1,nt
       write(46,FMT=9000)
     >           i,t,c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, 
     >            f(1,i),f(2,i),f(3,i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,3F12.6)

      call util_flush(46)
      close(46)


      end

      subroutine qmmm_print_pdbi(un,nt,ai,c,q,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer un
      integer nt
      integer ai(nt)
      double precision c(3,nt)
      double precision q(nt)
      character*16 t(nt)

      integer i
      integer n
      integer nf
      integer ns

      do i=1,nt
       write(un,FMT=9000)
     >           ai(i),t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      end

      subroutine qmmm_print_xyz(un,nt,c,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer un
      integer nt
      double precision c(3,nt)
      character*16 t(nt)

      integer i
      integer n
      integer nf
      integer ns

      write(un,FMT=1) nt

      do i=1,nt
       write(un,FMT=2)
     >            t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00

      end do
1     FORMAT(1X,I9,/, "   ")
2     FORMAT(1X,A5,6X,3(F12.6,4X))

      end

      subroutine qmmm_print_xyzi(un,nt,ai,c,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer un
      integer nt
      integer ai(nt)
      double precision c(3,nt)
      character*16 t(nt)
      character*4 t1

      integer i
      integer n
      integer nf
      integer ns

      write(un,FMT=1) nt

      do i=1,nt
        t1 = t(i)
c       write(un,FMT=2)
       write(un,2)
     >            t1,c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00,ai(i)

      end do
1     FORMAT(1X,I9,/, "   ")
2     FORMAT(1X,A5,6X,3(F12.6,4X),I6)

      end

      function qmmm_get_io_unit(fn)

      implicit none
      integer fn
      logical qmmm_get_io_unit
c 
      integer k
      logical ostatus
c
      do k=80,90
        INQUIRE(UNIT=k,OPENED=ostatus)
        ostatus = .not.ostatus
        if(ostatus) 
     >    INQUIRE(UNIT=k,EXIST=ostatus)
        if(ostatus) then
          fn = k
          qmmm_get_io_unit = .true.
          return
        end if 
      end do
      qmmm_get_io_unit = .false.
      return
      end

      function qmmm_file_exist(filename,perm)

      implicit none
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
      character*(*) filename
      logical perm
      logical qmmm_file_exist
c 
      character*255 filename0
      logical ofile
c
      filename0 = filename
      if(perm)
     >    call util_file_name_resolve(filename0,.false.)
      if(ga_nodeid().eq.0) then 
        inquire(file=filename0,exist=ofile)
      end if
      call ga_brdcst(msg_qmmm_misc, ofile,
     $     MA_sizeof(MT_INT,1,MT_BYTE),0)
      call ga_sync()
c
      qmmm_file_exist = ofile
      return
      end

      subroutine qmmm_interp_xyzi_file(aname1,aname2,aname3,lambda)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) aname1    
      character*(*) aname2    
      character*(*) aname3    
      double precision lambda
c     local variables
      integer ns
      integer i
      integer k
      logical title


      integer i_itmp,h_itmp
      integer i_ctmp,h_ctmp
      character*32 pname
      character*30 buf
      character*72 message
      character*255 filename1
      character*255 filename2
      character*4 tag(1000)
      double precision r1(3),r2(3)
      integer ai1,ai2
 
      integer fn1,fn2
      integer ns1,ns2

      pname = "qmmm_interp_xyzi_file"

      filename1 = aname1
      call util_file_name_resolve(filename1,.false.)
      if(.not.qmmm_get_io_unit(fn1)) 
     >       call errquit("cannot get file number",0,0)
c
      message = "opening filename1 for reading"
      open(fn1,file=filename1,form='formatted',status='old',
     $          err=133)
c
      filename2 = aname2
      call util_file_name_resolve(filename2,.false.)
      if(.not.qmmm_get_io_unit(fn2)) 
     >       call errquit("cannot get file number",0,0)
c
      message = "opening filename2 for reading"
      open(fn2,file=filename2,form='formatted',status='old',
     $          err=133)
       
c     get number of atoms
      read(fn1,*) ns1
      read(fn2,*) ns2
      if(ns1.ne.ns2) 
     >  call errquit( pname//'different number of atoms',0,0)
      ns = ns1
      if(ns.gt.1000)
     > call errquit(pname//"increase memory for atom tag",ns,0)

      if(.not.ma_push_get(mt_int,ns,'itmp',h_itmp,i_itmp))
     + call errquit( pname//'Failed to allocate memory for itmp1',
     + ns, MA_ERR)
      call ifill(ns,0,int_mb(i_itmp),1)

      if(.not.ma_push_get(mt_dbl,3*ns,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp1',
     + 3*ns, MA_ERR)
      call dfill(3*ns,0.0d0,dbl_mb(i_ctmp),1)

c
c     forward to cooordinates
c     -----------------------
      rewind(fn1)
      title =.false.
      do i=1,ns+2
        read(fn1,*,end=1) buf
      end do
      title = .true.
1     continue
      rewind(fn1) 
      read(fn1,*) buf
      if(title) 
     +   read(fn1,*) buf

      rewind(fn2)
      title =.false.
      do i=1,ns+2
        read(fn2,*,end=2) buf
      end do
      title = .true.
2     continue
      rewind(fn2) 
      read(fn2,*) buf
      if(title) 
     +   read(fn2,*) buf

      do i=1,ns
       read(fn1,*) tag(i), (r1(k),k=1,3),
     +             ai1
       read(fn2,*) buf, (r2(k),k=1,3),
     +             ai2

       if(ai2.ne.ai1)
     >   call errquit( pname//'different global index',0,0) 
       int_mb(i_itmp+i-1) = ai1
       do k=1,3
        dbl_mb(i_ctmp+3*(i-1)+k-1) = r1(k)+lambda*(r2(k)-r1(k))
       end do
      end do
c
c
      filename1 = aname3
      call util_file_name_resolve(filename1,.false.)
c
      message = "opening filename1 for writing"
      open(fn1,file=filename1,form='formatted',status='unknown',
     $          err=133)

      write(fn1,3) ns
      do i=1,ns
       write(fn1,4) tag(i), (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3),
     +             int_mb(i_itmp+i-1)

      end do
3     FORMAT(1X,I9,/, "   ")
4     FORMAT(1X,A5,6X,3(F12.6,4X),I6)
      close(fn1)
      close(fn2)
 
c
      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_itmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack i_itmp',ns,
     &       MA_ERR)


      return

 133  call errquit(pname//'error '//message,0, 0)

      end

      subroutine qmmm_interp_esp_file(aname1,aname2,aname3,lambda)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) aname1    
      character*(*) aname2    
      character*(*) aname3    
      double precision lambda
c     local variables
      integer ns
      integer i
      integer k
      logical title


      integer i_itmp,h_itmp
      integer i_ctmp,h_ctmp
      character*32 pname
      character*30 buf
      character*72 message
      character*255 filename1
      character*255 filename2
      double precision q1,q2,q3(1000)
      integer ai1,ai2,ai3(1000)
 
      integer fn1,fn2
      integer ns1,ns2

      pname = "qmmm_interp_esp_file"

      filename1 = aname1
      call util_file_name_resolve(filename1,.false.)
      if(.not.qmmm_get_io_unit(fn1)) 
     >       call errquit("cannot get file number",0,0)
c
      message = "opening filename1 for reading"
      open(fn1,file=filename1,form='formatted',status='old',
     $          err=133)
c
      filename2 = aname2
      call util_file_name_resolve(filename2,.false.)
      if(.not.qmmm_get_io_unit(fn2)) 
     >       call errquit("cannot get file number",0,0)
c
      message = "opening filename2 for reading"
      open(fn2,file=filename2,form='formatted',status='old',
     $          err=133)
       

      ns = 0
      do i=1,1000
       read(fn1,*,end=1) ai1,q1
       read(fn2,*,end=1) ai2,q2

       ns = ns + 1
       if(ai2.ne.ai1)
     >   call errquit( pname//'different global index',0,0) 
       if(ns.gt.1000)
     > call errquit(pname//"increase memory for atom tag",ns,0)
       ai3(i) = ai1
       q3(i) = q1+lambda*(q2-q1)
      end do
c
 1    continue
c
c
      filename1 = aname3
      call util_file_name_resolve(filename1,.false.)
c
      message = "opening filename1 for writing"
      open(fn1,file=filename1,form='formatted',status='unknown',
     $          err=133)

      do i=1,ns
       write(fn1,*) ai3(i),q3(i)
      end do

      close(fn1)
      close(fn2)

 
      return

 133  call errquit(pname//'error '//message,0, 0)

      end
