C $Id$
      block data gridptsdata
      implicit none
#include "grid_eaf.fh"
      data grpsize /-1/
      end
      logical function grid_grinfoinit(fdg,grinfo,enuc_in,what)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "grid_eaf.fh"
#include "inp.fh"
      integer fdg
      character*(*) grinfo, what
      double precision enuc_in
c
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
#if 1
cdump garbage on next record to invalidate fields
         uvect(1)=-123
         uvect(2)=-123
         uvect(3)=-123
         uvect(4)=-123
#endif
         else
         enuc_in=uenuc
      endif
      ierr=0
      grid_grinfoinit=ierr.eq.0
      if (otrace) then
         write(LuOut,125) '1INFO ' // what(1:4)// ' from ',
     A        grinfo(1:inp_strlen(grinfo)), ' enuc ',enuc_in
 125     format(A,A,A,e15.8)
      endif
      return
      end
      logical function grid_grinfo2(fdg,grinfo,what,
     I     max_file_rec, max_cache_rec,n_rec_in_file,
     L     file_opened)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "inp.fh"
#include "grid_eaf.fh"
      integer fdg
      character*(*) grinfo, what
      integer max_file_rec, max_cache_rec,n_rec_in_file
      logical file_opened
c
      integer foplog,ierr
      logical otrace
      ierr=0
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
          fdg=1
      else
         max_file_rec =uvect(1)
         max_cache_rec=uvect(2) 
         n_rec_in_file=uvect(3) 
         foplog       =uvect(4)
         file_opened=foplog.eq.1
c         if(fdg.eq.-1) then
c            ierr=-1
c            file_opened=.false.
c            n_rec_in_file=0
c         endif
      endif
      if (otrace) then
         write(LuOut,*) '2INFO ' // what(1:4) // ' from ',
     R        grinfo(1:inp_strlen(grinfo)),
     A        '  mxcrec',max_cache_rec,
     A        '  nrecfil',n_rec_in_file,
     A        '  flopend ',file_opened
 123     format(3A,I3,A,I5,A,I5,A,I1)
         call util_flush(LuOut)
      endif
      grid_grinfo2=ierr.eq.0
      return
      end
      logical function grid_ptsalloc(fd,eaf_size,
     S     eafname)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "grid_eaf.fh"
#include "errquit.fh"
      integer fd
      character*(*) eafname
      integer eaf_size
c
      integer ierr
      logical otrace,grid_ptsdel
      external grid_ptsdel
c
      otrace = util_print('grinfodebug',print_debug)
      if(otrace) then
         write(luout,*) ' in: alloc fd ',fd,
     S        ' size in dbl ',eaf_size
         call util_flush(luout)
      endif
      ierr=0
      if(fd.ne.-1.and.eaf_size.ne.grpsize) then
chere we want to change the alloc size
         fd=-1
      if(.not.grid_ptsdel(fd)) call errquit(
     E        'gridptsall: ptsdel failed',0,MA_ERR)
      endif
      if(fd.eq.-1) then
         if(.not.MA_Alloc_Get(MT_dbl,eaf_size,
     A        eafname,l_gridp,k_gridp)) then
            call
     C        errquit(' ptsalloc: increase memory in input line',
     C           0,MA_ERR)
         endif
      fd=k_gridp
      grpsize=eaf_size
      endif
      grid_ptsalloc=.true.
      if(otrace) then
         write(luout,*) ' out: alloc fd ',fd,
     S        ' size in dbl ',eaf_size
         call util_flush(luout)
      endif
      return
      end
      logical function grid_ptsdel(fd)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"
#include "grid_eaf.fh"
      integer fd
c
      integer ierr
      logical otrace
c
      otrace = util_print('grinfodebug',print_debug)
      if(otrace) then
         write(luout,*) 'ptswrite fd in ',fd,
     A        ' kgridp ',k_gridp
         call util_flush(luout)
      endif
c     return because gridpts not allocated
      if(fd.eq.-1) then
      if(otrace)  write(luout,*) ' quick return for fd= ',fd
         grid_ptsdel=.true.
         return
      endif

csanity check
      if(fd.ne.k_gridp) call errquit(
     A     'gridptsdel: fd ne kgridp',0,MA_ERR)
      ierr=0
      if (.not.MA_Free_Heap(l_gridp))
     &   ierr=-1
      fd=-1
      grid_ptsdel=ierr.eq.0
      return
      end
      integer function grid_ptswrite(fd,offset,buf_in,
     S     lngth)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"
      integer fd
      double precision offset,buf_in(*)
      integer lngth
c
      integer l_gridp,k_gridp,grpsize
      common /grpgrp/ l_gridp,k_gridp,grpsize
      integer offset_int,lngth_dbl
      logical otrace
c
      otrace = util_print('grinfodebug',print_debug)
      if(otrace) then
         write(luout,*) '3fd in ',fd
         call util_flush(luout)
      endif
csanity check
      if(fd.ne.k_gridp) call errquit(
     A     'gridptwrite: fd ne kgridp',0,MA_ERR)
      offset_int=offset
      lngth_dbl=lngth/MA_sizeof(MT_DBL,1,MT_BYTE)
      if(offset_int+lngth_dbl.gt.grpsize) then
         write(luout,*) ' add input line set grid:eaf_size_in_dbl',
     A        (offset_int+lngth_dbl)
         call
     E     errquit('grpwrite: insuff eafsize ',grpsize,MA_ERR)
      endif
      call dcopy(lngth_dbl,
     B     buf_in,1,dbl_mb(k_gridp+offset_int),1)
      grid_ptswrite=0
      return
      end
      integer function grid_ptsread(fd,offset,buf_out,
     S     lngth)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"
      integer fd
      double precision offset,buf_out(*)
      integer lngth
c
      integer l_gridp,k_gridp,grpsize
      common /grpgrp/ l_gridp,k_gridp,grpsize
      integer offset_int,lngth_dbl
      logical otrace
c
      otrace = util_print('grinfodebug',print_debug)
      if(otrace) then
         write(luout,*) '4fd in ',fd
         call util_flush(luout)
      endif
csanity check
      if(fd.ne.k_gridp) call errquit(
     A     'gridptread: fd ne kgridp',0,MA_ERR)
      offset_int=offset
      lngth_dbl=lngth/MA_sizeof(MT_DBL,1,MT_BYTE)
      if(offset_int+lngth_dbl.gt.grpsize) call
     E     errquit('grpwrite: insuff eafsize ',grpsize,MA_ERR)
      call dcopy(lngth_dbl,
     B     dbl_mb(k_gridp+offset_int),1,buf_out,1)
      grid_ptsread=0
      return
      end
