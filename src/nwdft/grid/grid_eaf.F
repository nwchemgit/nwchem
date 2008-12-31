C $Id: grid_eaf.F,v 1.2 2008-12-31 17:14:21 edo Exp $
      logical function grid_grinfoinit(fdg,grinfo,enuc_in,what)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "eaf.fh"
#include "stdio.fh"
#include "util.fh"
#include "inp.fh"
      double precision fdg
      character*(*) grinfo, what
      double precision enuc_in
c
      integer ierr
      character*80 errmsg
      double precision offset
      logical otrace
c
      otrace = util_print('griddebug',print_debug)
      ierr=eaf_open(grinfo, eaf_rw, fdg)
      if (ierr .ne. 0) then
         call eaf_errmsg(ierr,errmsg)
         write(LuOut,*) ga_nodeid(),errmsg
         call errquit('gridfileopen: eaf_open grinfo failed', ierr,
     &        DISK_ERR)
      endif
cgrinfo     enuc in first record
      offset=0
      if(what.eq.'write') then
         ierr = eaf_write(fdg, offset, enuc_in,
     .        8)
         else
         ierr = eaf_read(fdg, offset, enuc_in,
     .        8)
      endif
      grid_grinfoinit=ierr.eq.0
      if (otrace) then
         write(LuOut,*) ' INFO ' // what(1:4)// ' from ',grinfo
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
      double precision fdg
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
      otrace = util_print('griddebug',print_debug)
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
         write(LuOut,*) ' INFO ' // what(1:4) // ' from ',grinfo
         write(LuOut,*) ' ierr ',ierr
         write(LuOut,*) '      mxcrec',max_cache_rec
         write(LuOut,*) '      nrecfil',n_rec_in_file
         write(LuOut,*) '      flopend',file_opened
         call util_flush(LuOut)
      endif
      grid_grinfo2=ierr.eq.0
      return
      end
