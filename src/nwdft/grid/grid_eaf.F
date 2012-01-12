C $Id$
      logical function grid_grinfoinit(fdg,grinfo,enuc_in,what)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "eaf.fh"
#include "stdio.fh"
#include "util.fh"
      integer fdg
      character*(*) grinfo, what
      double precision enuc_in
c
      double  precision uenuc   ! grinfo enuc
      integer uvect(4) ! grinfo contents
      common /grgr/ uenuc,uvect
      integer ierr
      logical otrace
c
      otrace = util_print('grinfodebug',print_debug)
      if(otrace) then
         write(luout,*) '1fdg in ',fdg
         call util_flush(luout)
      endif
cgrinfo     enuc in first record
      if(what.eq.'write') then
         uenuc=enuc_in
cdump garbage on next record to invalidate fields
         uvect(1)=-123
         uvect(2)=-123
         uvect(3)=-123
         uvect(4)=-123
         else
         enuc_in=uenuc
      endif
      ierr=0
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
#include "stdio.fh"
#include "util.fh"
      integer fdg
      character*(*) grinfo, what
      integer max_file_rec, max_cache_rec,n_rec_in_file
      logical file_opened
c
      double  precision uenuc   ! grinfo enuc
      integer uvect(4) ! grinfo contents
      common /grgr/ uenuc,uvect
      integer foplog,ierr
      logical otrace
c
      otrace = util_print('grinfodebug',print_debug)
      if(otrace) then
         write(luout,*) '2fdg in ',fdg
         call util_flush(luout)
      endif
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
      else
         max_file_rec =uvect(1)
         max_cache_rec=uvect(2) 
         n_rec_in_file=uvect(3) 
         foplog       =uvect(4)
         file_opened=foplog.eq.1
      endif
      ierr=0
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
