c $Id$
      subroutine et_det(op,handle,n,detval)
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
c
c returns the determinant of the square matrix with file handle 'handle'.
c n is the dimension of the array, detval is the returned value of the
c determinant. If op='N', a=a. If op='T', a=a'.
c
      character*1 op
      integer handle
      integer n,i,j,info
      integer l_ipvt,k_ipvt,l_a,k_a,l_work,k_work
      double precision detval,det(2),ten
c
      ten = 10.0d0
      if (.not.ma_push_get(mt_int,n,'ipvt',
     1  l_ipvt,k_ipvt)) call errquit
     2  ('et_det: failed to allocate ipvt',0,0)
      if (.not.ma_push_get(mt_dbl,n*n,'a',
     1  l_a,k_a)) call errquit
     2  ('et_det: failed to allocate a',0,0)
      if (.not.ma_push_get(mt_dbl,n,'work',
     1  l_work,k_work)) call errquit
     2  ('et_det: failed to allocate work',0,0)
c
      do i = 1,n
      do j = 1,n
       if(op.eq.'N') call ga_get(handle,i,i,j,j,
     1 dbl_mb(k_a+i-1+(j-1)*n),1)
       if(op.eq.'T') call ga_get(handle,j,j,i,i,
     1 dbl_mb(k_a+i-1+(j-1)*n),1)
      enddo
      enddo
c
      call dgefa(dbl_mb(k_a),n,n,int_mb(k_ipvt),info)
      call dgedi(dbl_mb(k_a),n,n,int_mb(k_ipvt),det,dbl_mb(k_work),10)
      detval =  det(1) * ten**det(2)
c
      if (.not.ma_pop_stack(l_work)) call errquit
     1  ('et_det: failed to deallocate work',0,0)
      if (.not.ma_pop_stack(l_a)) call errquit
     1  ('et_det: failed to deallocate a',0,0)
      if (.not.ma_pop_stack(l_ipvt)) call errquit
     1  ('et_det: failed to deallocate ipvt',0,0)
      call ga_sync()
      return
      end

