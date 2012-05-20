      subroutine smd_paramfile_natoms(na)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer na
      logical result
c
      character*(180) buffer
      character*(nw_max_path_len) filename
      integer un
      character*30 pname
c
      pname = "smd_paramfile_natoms"
c
      call smd_paramfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no parameter file found',0, 0)

c
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)


      na = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        na = na + 1
      end if
      goto 100
200   continue
      close(un)

      end

      subroutine smd_paramfile_nshake(na)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer na
c
      character*(180) buffer
      integer un
c
      character*30 pname
      logical result
      character*(nw_max_path_len) filename
c
      pname = "smd_paramfile_nshake"
c      write(*,*) "in "//pname
c
      call smd_paramfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no parameter file found',0, 0)
c
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)


      na = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:5).eq."SHAKE") then
        na = na + 1
      end if
      goto 100
200   continue
      close(un)

      end

      subroutine smd_paramfile_nbond(na)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer na
c
      character*(180) buffer
      integer un
      logical result
      character*(nw_max_path_len) filename
c
      character*30 pname

      pname = "smd_paramfile_nbond"
c      write(*,*) "in "//pname
c
      call smd_paramfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no parameter file found',0, 0)
c
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)


      na = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."BOND") then
        na = na + 1
      end if
      goto 100
200   continue
      close(un)

      end

      subroutine smd_paramfile_nres(na)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer na
c
      character*(180) buffer
      integer un
      character*30 pname
      logical result
      character*(nw_max_path_len) filename
c
      pname = "smd_paramfile_nres"

c
      call smd_paramfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no parameter file found',0, 0)
c

      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)


      na = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:3).eq."MOL") then
        na = na + 1
      end if
      goto 100
200   continue
      close(un)

      end

      subroutine smd_paramfile_read_atom(na,ta,tt,tr,mass,q)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer na
      character*(smd_string_size) ta(na)
      character*(smd_string_size) tt(na)
      character*(smd_string_size) tr(na)
      double precision mass(na)
      double precision q(na)
c
      character*180 buffer
      character*4 tag
      character*30 pname
      character*(smd_string_size) rtag
      integer i
      integer un
      logical result
      character*(nw_max_path_len) filename

      pname = "smd_paramfile_read_params"
c
      call smd_paramfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no parameter file found',0, 0)
c

      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)


      i=0
100   continue
      read(un,'(A180)',end=150) buffer
      if(buffer(1:3).eq."MOL") then
        read(buffer,*) tag,rtag
        goto 100
      end if
      if(buffer(1:4).eq."ATOM") then
         i = i +1
         if(i.gt.na)
     >  call errquit(pname//"insufficient array size",0,0)
         read(buffer,*) tag,tag,ta(i),tt(i),mass(i),q(i)
         tr(i)=rtag
         goto 100
      end if
      goto 100
150   continue
      close(un)
      return
      end

      subroutine smd_paramfile_read_pot(np,pt,
     +                                  p1,p2,p3,
     +                                  nt,tt,tr)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer np
      integer pt(np)
      double precision p1(np)
      double precision p2(np)
      double precision p3(np)
      integer nt
      character*(smd_string_size) tt(nt)
      character*(smd_string_size) tr(nt)
c
      character*180 buffer
      character*4 tag
      character*(smd_string_size) tag1,tag2
      character*30 pname
      integer i,i1,i2,npr
      integer un

      logical result
      character*(nw_max_path_len) filename

      pname = "smd_paramfile_read_pot"

c      write(*,*) "in "//pname
c
      call smd_paramfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no parameter file found',0, 0)
c

      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)

      i = 0
      npr = 0
      do i1=1,nt
       do i2=i1,nt
        i = i+1
         rewind(un)
100      continue
         read(un,'(A180)',end=150) buffer
         if(buffer(1:3).eq."POT") then
          read(buffer,*) tag,tag1,tag2
          if(tt(i1).eq.tag1) then
           if(tt(i2).eq.tag2) goto 110
          end if
          if(tt(i1).eq.tag2) then
           if(tt(i2).eq.tag1) goto 110
          end if
         end if
         goto 100
110      continue
         read(buffer,*) tag,tag1,tag2,pt(i),
     +                  p1(i),p2(i),p3(i)
         npr = npr + 1
150      continue
       end do
      end do
      
      np = npr

c      i = 0
c      do i1=1,nt
c       do i2=i1,nt
c        i = i+1
c        write(*,*) "index",i,
c     +             ((i1-1)*(2*nt+2-i1))/2+i2-i1+1,
cc     +             (i1*(i1+1))/2+i2-i1,
c     +             i1,i2,nt
c        write(*,*) pt(i),p1(i),p2(i),p3(i)
c       end do
c      end do
c 
      close(un)
      return
      end

      subroutine smd_paramfile_read_shake(ns,nt,s1,s2,d,tt)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer nt
      integer ns
      integer s1(ns)
      integer s2(ns)
      double precision d(ns)
      character*(smd_string_size) tt(nt)
c
      character*180 buffer
      character*4 tag
      character*(smd_string_size) tag1,tag2
      character*30 pname
      integer i,j,i1,i2
      double precision dtmp
      integer un

      logical result
      character*(nw_max_path_len) filename

      pname = "smd_paramfile_read_shake"

c      write(*,*) "in "//pname
c
      call smd_paramfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no parameter file found',0, 0)
c

      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)

      i = 0
      do i1=1,nt
       do i2=i1+1,nt
         rewind(un)
100      continue
         read(un,'(A180)',end=150) buffer
         if(buffer(1:5).eq."SHAKE") then
          read(buffer,*) tag,tag1,tag2,dtmp
          if(tt(i1).eq.tag1) then
           if(tt(i2).eq.tag2) goto 110
          end if
          if(tt(i1).eq.tag2) then
           if(tt(i2).eq.tag1) goto 110
          end if
         end if
         goto 100
110      continue
         i = i + 1
         s1(i) = i1
         s2(i) = i2
         d(i)  =dtmp
150      continue
       end do
      end do
       
      ns = i

      close(un)
      return
      end

      subroutine smd_paramfile_read_bond(nb,nt,s1,s2,
     >                                   it,d,k,ta,tr)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer nt
      integer nb
      integer s1(nb)
      integer s2(nb)
      integer it(nb)
      double precision d(nb)
      double precision k(nb)
      character*(smd_string_size) ta(nt)
      character*(smd_string_size) tr(nt)
c
      character*180 buffer
      character*4 tag
      character*(smd_string_size) atag1,atag2,rtag
      character*30 pname
      integer i,j,i1,i2
      integer s1t,s2t,itt
      double precision dt,kt
      logical lr,la1,la2
      integer un

      logical result
      character*(nw_max_path_len) filename

      pname = "smd_paramfile_read_bond"

c      write(*,*) "in "//pname
c
      call smd_paramfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no parameter file found',0, 0)
c

      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)


      i=0
100   continue
      read(un,'(A180)',end=150) buffer
      if(buffer(1:3).eq."MOL") then
        read(buffer,*) tag,rtag
        goto 100
      end if
      if(buffer(1:4).eq."BOND") then
         if(i+1.gt.nb)
     >  call errquit(pname//"insufficient array size",0,0)
         read(buffer,*) tag,atag1,atag2,itt,dt,kt
         s1t = 0
         s2t = 0
         do j=1,nt
           lr = rtag.eq.tr(j)
           la1 = atag1.eq.ta(j)
           la2 = atag2.eq.ta(j)
           if(lr.and.la1)  s1t=j
           if(lr.and.la2)  s2t=j
           if(s1t*s2t.ne.0) then
             i = i +1
             s1(i) = min(s1t,s2t)
             s2(i) = max(s1t,s2t)
             it(i) = itt
             d(i) = dt
             k(i) = kt
             goto 100
           end if
         end do 
      end if
      goto 100
150   continue

      close(un)
      return
      end

      subroutine smd_paramfile(filename,result)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "global.fh"

c
      character*(*) filename
      logical result
c
      character*30 pname

      pname = "smd_paramfile"

      result = .true.
      call smd_rtdb_get_string("smd:paramfile",1,
     >                           filename,result)

      end
c $Id$
