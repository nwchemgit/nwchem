C $Id$
      logical function grid_grinfoinit(fdg,grinfo,enuc_in,what)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "eaf.fh"
#include "stdio.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
      integer fdg
      character*(*) grinfo, what
      double precision enuc_in
c
      integer ierr
      character*80 errmsg
      double precision offset
      integer uvect(4)
      logical otrace
c
      otrace = util_print('grinfodebug',print_debug)
      if(otrace) then
         write(luout,*) '1fdg in ',fdg
         call util_flush(luout)
      endif
      if(fdg.eq.-1) then
cold      ierr=eaf_open(grinfo, eaf_rw, fdg)
c now we open a MA of 1k since eaftype=10000+1024
         ierr=eaf_open(grinfo, 1001024, fdg)
         if (ierr .ne. 0) then
            call eaf_errmsg(ierr,errmsg)
            write(LuOut,*) ga_nodeid(),errmsg
            call errquit('grideaf: eaf_open grinfo failed', ierr,
     &           DISK_ERR)
         endif
      endif
cgrinfo     enuc in first record
      offset=0
      if(what.eq.'write') then
         ierr = eaf_write(fdg, offset, enuc_in,
     .        8)
cdump garbage on next record to invalidate fields
         uvect(1)=-123
         uvect(2)=-123
         uvect(3)=-123
         uvect(4)=-123
         ierr = eaf_write(fdg, 8d0, uvect,
     .        4*ma_sizeof(mt_int,1,mt_byte))
         else
         ierr = eaf_read(fdg, offset, enuc_in,
     .        8)
      endif
      grid_grinfoinit=ierr.eq.0
      if (otrace) then
         write(LuOut,*) '1INFO ' // what(1:4)// ' from ',grinfo
         write(LuOut,*) '      enuc',enuc_in
         call util_flush(LuOut)
      endif
      return
      end
      logical function grid_grinfo2(fdg,grinfo,what,
     I     max_file_rec, max_cache_rec,n_rec_in_file,
     L     file_opened)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "eaf.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
      integer fdg
      character*(*) grinfo, what
      integer max_file_rec, max_cache_rec,n_rec_in_file
      logical file_opened
c
      integer foplog,intlgt,ierr
      integer uvect(4)
      double precision offset
      logical otrace
      character*80 errmsg
c
      intlgt=ma_sizeof(mt_int,1,mt_byte)
      otrace = util_print('grinfodebug',print_debug)
      if(otrace) then
         write(luout,*) '2fdg in ',fdg
         call util_flush(luout)
      endif
      offset=8d0
cgrinfo     max_file_rec in 2nd record
cgrinfo     max_cache_rec in 3nd record
cgrinfo     n_rec_in_file in 4th record
cgrinfo     foplog in 5th record
      if(what.eq.'write') then
         foplog=1
         if(.not.file_opened) foplog=0
         uvect(1)=max_file_rec
         uvect(2)=max_cache_rec 
         uvect(3)=n_rec_in_file 
         uvect(4)=foplog
         ierr = eaf_write(fdg, offset, uvect,
     .        4*intlgt)
      else
         ierr = eaf_read(fdg, offset, uvect,
     .        4*intlgt)
         max_file_rec =uvect(1)
         max_cache_rec=uvect(2) 
         n_rec_in_file=uvect(3) 
         foplog       =uvect(4)
         file_opened=foplog.eq.1
      endif
      if (ierr .ne. 0) then
         call eaf_errmsg(ierr,errmsg)
         write(LuOut,*) ga_nodeid(),errmsg
         call util_flush(LuOut)         
         call errquit('gridfileopen: eaf_' // what(1:4) //' failed', 
     I        ierr, DISK_ERR)
      endif
      if (otrace) then
         write(LuOut,*) '2INFO ' // what(1:4) // ' from ',grinfo
         write(LuOut,*) ' ierr ',ierr
         write(LuOut,*) '      mxcrec',max_cache_rec
         write(LuOut,*) '      nrecfil',n_rec_in_file
         write(LuOut,*) '      flopend',file_opened
         call util_flush(LuOut)
      endif
      grid_grinfo2=ierr.eq.0
      return
      end
