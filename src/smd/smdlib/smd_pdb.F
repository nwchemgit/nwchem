      subroutine smd_pdb_natoms(filename,nt)
      implicit none
#include "errquit.fh"
#include "util.fh"
      character*(*) filename
      integer nt
c
      character*(4) buffer
      integer un
c
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)

      nt = 0
100   continue
      read(un,'(A4)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        nt = nt + 1
      end if
      goto 100
200   continue
      close(un)

      end

      subroutine smd_pdb_read_coords(filename,nt,c)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
      character*(*) filename
      integer nt
      double precision c(nt,3)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname

      pname = "smd_pdb_read_coords"
c
      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)
      i = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        i = i +1
        if(i.gt.nt)
     >  call errquit(pname//"insufficient array size",0,0)
        read(buffer,*) tag,tag,tag,tag,tag,
     >                 c(i,1),c(i,2),c(i,3)
       end if
      goto 100
200   continue
      close(un)

      if(i.ne.nt)
     >  call errquit(pname//"incompatible array sizes",0,0)
      return

      end

      subroutine smd_pdb_read_coords1(filename,nt,c)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
      character*(*) filename
      integer nt
      double precision c(3,nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname

      pname = "smd_pdb_read_coords"
c
      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)
      i = 0
100   continue
      read(un,'(A180)',end=200) buffer
      call util_flush(34)
      if(buffer(1:4).eq."ATOM") then
        i = i +1
        if(i.gt.nt)
     >  call errquit(pname//"insufficient array size",0,0)
        read(buffer,*) tag,tag,tag,tag,tag,
     >                 c(1,i),c(2,i),c(3,i)
       end if
      goto 100
200   continue
      close(un)

      if(i.ne.nt)
     >  call errquit(pname//"incompatible array sizes",0,0)
      return

      end
      subroutine smd_pdb_read_atomres(filename,nt,ta,tr,ir)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      character*(*) filename
      integer nt
      character*(smd_string_size) ta(nt)
      character*(smd_string_size) tr(nt)
      integer ir(nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname

      pname = "smd_pdb_read_atomres"
c
      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)
      i = 0
100   continue
      read(un,'(A180)',end=200) buffer
      call util_flush(34)
      if(buffer(1:4).eq."ATOM") then
        i = i +1
        if(i.gt.nt)
     >  call errquit(pname//"insufficient array size",0,0)
        read(buffer,*) tag,tag,ta(i),tr(i),ir(i),
     >                 tag,tag,tag
        call util_flush(6)
       end if
      goto 100
200   continue
      close(un)

      if(i.ne.nt)
     >  call errquit(pname//"incompatible array sizes",0,0)
      return

      end

      subroutine smd_pdb_read(filename,nt,ta,tr,ir,c)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      character*(*) filename
      integer nt
      character*(smd_string_size) ta(nt)
      character*(smd_string_size) tr(nt)
      integer ir(nt)
      double precision c(3,nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname

      pname = "smd_pdb_read"
c
      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)
      i = 0
100   continue
      read(un,'(A180)',end=200) buffer
      call util_flush(34)
      if(buffer(1:4).eq."ATOM") then
        i = i +1
        if(i.gt.nt)
     >  call errquit(pname//"insufficient array size",0,0)
        read(buffer,*) tag,tag,ta(i),tr(i),ir(i),
     >                 c(1,i),c(2,i),c(3,i)
        call util_flush(6)
       end if
      goto 100
200   continue
      close(un)

      if(i.ne.nt)
     >  call errquit(pname//"incompatible array sizes",0,0)
      return

      end

      subroutine smd_pdb_nres(filename,nr)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
      character*(*) filename
      integer nr
c
      character*(180) buffer
      character*(4)  tag
      integer ir0,ir
      integer un
c
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)

c      reset residue arrays to be the size of number of residues only
      nr = 0
      ir0 = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        read(buffer,*) tag,tag,tag,tag,ir
        if(ir0.ne.ir) then
          nr = nr + 1
          ir0=ir
        end if
      end if
      goto 100
200   continue
      close(un)

      end

      subroutine smd_pdb_read_res(filename,nt,nr,tr,ir,nm)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      character*(*) filename
      integer nt,nr,nc
      character*(smd_string_size) tr(nr)
      integer ir(nt)
      integer nm(nr)
c
      character*(30) pname
      character*(180) buffer
      character*(4)  tag
      character*(smd_string_size)  rtag,rtag0
      integer ir0,nr0
      integer ncenter
      integer un
c
      pname = "sg_read_res"
c
      if(.not.util_get_io_unit(un))
     >   call errquit(pname//"cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)

      ncenter = 0
      nr0 = 0
      rtag0 = " "
      ir0 = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        ncenter = ncenter + 1
        read(buffer,*) tag,tag,tag,rtag,ir(ncenter)
        write(*,*) "ir",ncenter,ir(ncenter)
        if(ir0.ne.ir(ncenter)) then
          ir0=ir(ncenter)
          nr0 = nr0 + 1
          if(ncenter.gt.nt)
     >  call errquit(pname//"insufficient array size",0,0)
          if(nr0.gt.nr)
     >  call errquit(pname//"insufficient array size",0,0)
          tr(nr0) = rtag
          rtag0=rtag 
        end if
        ir(ncenter) = nr0
        nm(nr0) = nm(nr0) + 1
      end if
      goto 100
200   continue

      close(un)

      write(*,*) "ncenter,nt",ncenter,nt
      if(ncenter.ne.nt)
     >  call errquit(pname//"incompatible atom array sizes",
     >             ncenter,nt)

      if(nr.ne.nr0)
     >  call errquit(pname//"incompatible residue array sizes",nr0,nr)

      return
      end
c $Id$
