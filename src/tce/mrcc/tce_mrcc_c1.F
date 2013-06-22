c      SUBROUTINE tce_mrcc_c1(size_t1m,k_t1_offsetm,d_t1m,nref,iref,
c     1 d_r1m,k_r1_offsetm)
      SUBROUTINE tce_mrcc_c1(size_t1m,k_t1_offsetm,d_t1m,iref,
     1 d_r1m,k_r1_offsetm)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "util.fh"
#include "msgids.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_mrcc.fh"
#include "stdio.fh"

      integer inoab,invab
      integer size_t1m(maxref)
      integer k_t1_offsetm(maxref)
      integer d_t1m(maxref)
c      integer iref,nref 
      integer iref
      integer size
      integer dima,dimb
      integer p1b,p2b,h3b,h4b
      integer l_t2,k_t2
      integer l_t1a1,k_t1a1
      integer l_t1b1,k_t1b1
      integer l_t1a2,k_t1a2
      integer l_t1b2,k_t1b2
      integer l_tmp,k_tmp
      integer d_r1m(maxref)
      integer k_t2_2,l_t2_2
      integer l_r2,k_r2
      integer l_c,k_c
      integer irefnew
      integer p1,p2,h3,h4
      integer i
      character*255 filename
      integer d_t1c
      integer k_r1_offsetm(maxref)
      double precision coef
      integer d_t1col

          call tce_filename('t1c',filename)
       if(lusesub) then
          call pgcreatefile(filename,d_t1c,size_t1m(iref),mypgid)
       else
          call createfile(filename,d_t1c,size_t1m(iref))
       endif
          call tce_filename('t1col',filename)
       if(lusesub) then
          call pgcreatefile(filename,d_t1col,size_t1m(iref),mypgid)
       else
          call createfile(filename,d_t1col,size_t1m(iref))
       endif
c          call ga_zero(d_t1col)

c      write(6,*)'IREF ',iref
c      call ga_print(d_t1m(iref))
      call util_flush(LuOut)

      do irefnew=1,nref
      if(irefnew.ne.iref) then

       coef = dbl_mb(k_heff+iref-1+(irefnew-1)*nref)*
     1 dbl_mb(k_sqc+(irefnew-1)*nref+mkroot-1)

       if(dbl_mb(k_sqc+(iref-1)*nref+mkroot-1).lt.0.0d0)
     1 coef=-1.0d0*coef

c     2 dbl_mb(k_sqc+(mkroot-1)*nref+iref-1)
c      write(6,"('USING',2F16.12)")coef
c      write(6,"('USING',2F16.12,2F16.12)")
c     1dbl_mb(k_heff+iref-1+(irefnew-1)*nref),
c     2dbl_mb(k_sqc+(irefnew-1)*nref+mkroot-1)


      call ga_zero(d_t1c)
       if(lusesub) then
         call ga_pgroup_sync(mypgid)
       else
      call reconcilefile(d_t1c,size_t1m(iref))
       endif

      call tce_mrcc_c1_offset(k_t1_offsetm,d_t1m,iref,irefnew,d_t1c)
       if(lusesub) then
         call ga_pgroup_sync(mypgid)
       else
         call reconcilefile(d_t1c,size_t1m(iref))
       endif

c      write(6,*)'T1T1T1T1T1T1'
c      call tce_mrcc_print_t1type(d_t1m(iref),k_t1_offsetm(iref),iref,
c     1irefnew)

c      write(6,*)'SELECTED T1'
c      call tce_mrcc_print_t1type(d_t1c,k_t1_offsetm(iref),iref,irefnew)

c      write(6,*)'R1R1R1R1R1R1'
c      call tce_mrcc_print_t1type(d_r1m(iref),k_t1_offsetm(iref),
c     2 iref,irefnew)

c      call ga_add(1.0d0,d_t1c,-1.0d0,d_t1m(iref),d_t1c)
      call tce_mrcc_mul_1(1.0d0,d_t1c,-1.0d0,d_t1m(iref),
     1 k_t1_offsetm(iref),d_t1c,iref)

       if(lusesub) then
         call ga_pgroup_sync(mypgid)
       else
      call reconcilefile(d_t1c,size_t1m(iref))
       endif

c      write(6,"(I4,'COEF : ',2F16.12,2F16.12,2F16.12)")
c     1ga_nodeid(),dbl_mb(k_heff+iref-1+(irefnew-1)*nref),
c     1 dbl_mb(k_sqc+(irefnew-1)*nref+mkroot-1),coef
c      call ga_add(coef,d_t1c,1.0d0,d_t1col,d_t1col)
      call tce_mrcc_mul_1(coef,d_t1c,1.0d0,d_t1col,
     1 k_t1_offsetm(iref),d_t1col,iref)
       if(lusesub) then
         call ga_pgroup_sync(mypgid)
       else
      call reconcilefile(d_t1col,size_t1m(iref))
       endif

      endif 
      enddo

      call util_flush(LuOut)

c abs(dbl_mb(k_sqc+(mkroot-1)*nref+iref-1))
      call tce_mrcc_mul_1(abs(dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)),
     1 d_r1m(iref),1.0d0,d_t1col,k_t1_offsetm(iref),d_r1m(iref),iref)

c      call ga_add(abs(dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)),
c     1 d_r1m(iref),1.0d0,d_t1col,d_r1m(iref))
       if(lusesub) then
         call ga_pgroup_sync(mypgid)
       else
         call reconcilefile(d_r1m(iref),size_t1m(iref))
       endif

c       write(6,*)'FINAL R1'
c       call tce_mrcc_print_t1type(d_r1m(iref),k_t1_offsetm(iref),
c     2 iref,irefnew)

      call deletefile(d_t1c) 
      call deletefile (d_t1col)

      return
      end
c $Id$
