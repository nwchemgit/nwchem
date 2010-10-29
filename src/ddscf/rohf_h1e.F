      subroutine rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $     g_fcv, g_fpv, g_fcp,
     $     g_x, g_ax)
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c     
      integer basis, geom
      integer nmo, nclosed, nopen
      integer g_fcv, g_fpv, g_fcp
      integer g_x
      integer g_ax
c     
      integer nvir, voff, ooff, oend, xoff, xend
      integer g_tmp1, g_tmp2, g_tmp3
      double precision four, mfour, one, mone, two, mtwo, zero
      parameter (four=4.0d0, mfour=-4.0d0, one=1.0d0, mone=-1.0d0,
     $     two=2.0d0, mtwo=-2.0d0, zero=0.0d0)
c
      integer ivec, nvec, gtype, vlen
c
      call ga_inquire(g_x,gtype,vlen,nvec)
c     
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
      ooff = nclosed + 1
      oend = nclosed + nopen
c     
c     Reshape argument vector by copying patches into matrix
c     
c       cv      cv  cv       cv  cv
c     S   =  4.F  .x   -  4.x  .F
c     
      if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_hv1e: tmp1',
     $     32, 32, g_tmp1)) call errquit('rohf_hv1e: tmp1', 0, GA_ERR)
      if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_hv1e: tmp2',
     $     32, 32, g_tmp2)) call errquit('rohf_hv1e: tmp2', 0, GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_hv1e: tmp3',
     $        32, 32, g_tmp3)) call errquit('rohf_hv1e: tmp3', 0,
     &       GA_ERR)
      endif
c
      do ivec = 1, nvec
         call ga_zero(g_tmp1)
         call ga_vec_to_mat(g_tmp1, voff, nmo, 1, nclosed, 
     $                      g_x, 1, ivec)
*         call ga_copy_patch('n', g_x, 1, (nvir*nclosed), ivec, ivec,
*     $        g_tmp1, voff, nmo, 1, nclosed )
         call ga_zero(g_tmp2)
         call ga_matmul_patch('n','n',four,zero,
     $        g_fcv, voff, nmo, voff, nmo,
     $        g_tmp1, voff, nmo, 1, nclosed,
     $        g_tmp2, voff, nmo, 1, nclosed )
         call ga_matmul_patch('n','n',mfour,one,
     $        g_tmp1, voff, nmo, 1, nclosed, 
     $        g_fcv, 1, nclosed, 1, nclosed,
     $        g_tmp2, voff, nmo, 1, nclosed )
c
         if (nopen .gt. 0) then
c     
            xoff = nclosed*nvir + 1
            xend = (nclosed + nopen)*nvir
            call ga_copy_patch('n', g_x, xoff, xend, ivec, ivec,
     $           g_tmp1, voff, nmo, ooff, oend )
            xoff = nvir*(nclosed+nopen) + 1
            xend = nvir*(nclosed+nopen) + nclosed*nopen
            call ga_copy_patch('n', g_x, xoff, xend, ivec, ivec,
     $           g_tmp1, ooff, oend, 1, nclosed )
c     
c       pv       pv  pv       pv  pv
c     S   =  2.F  .x   -  2.x  .F
c     
            call ga_matmul_patch('n','n',two,one,
     $           g_fpv, voff, nmo, voff, nmo,
     $           g_tmp1, voff, nmo, ooff, oend,
     $           g_tmp2, voff, nmo, ooff, oend )
            call ga_matmul_patch('n','n', mtwo, one,
     $           g_tmp1, voff, nmo, ooff, oend,
     $           g_fpv, ooff, oend, ooff, oend,
     $           g_tmp2, voff, nmo, ooff, oend )
c     
c       cp       cp  cp       cp  cp
c     S   =  2.F  .x   -  2.x  .F
c     
            call ga_matmul_patch('n','n',two,one,
     $           g_fcp, ooff, oend, ooff, oend,
     $           g_tmp1, ooff, oend, 1, nclosed,
     $           g_tmp2, ooff, oend, 1, nclosed )
            call ga_matmul_patch('n','n',mtwo,one,
     $           g_tmp1, ooff, oend, 1, nclosed,
     $           g_fcp, 1, nclosed, 1, nclosed,
     $           g_tmp2, ooff, oend, 1, nclosed )
c     
c       cv         cp     pv   cp
c     S    += (2.F   +  F  ).x
c     
c       cp         cp     pv   cv
c     S    += (2.F   +  F  ).x
c     
            call ga_dadd(two,g_fcp,one,g_fpv,g_tmp3)
            call ga_matmul_patch('n', 'n', one, one,
     $           g_tmp3, voff, nmo, ooff, oend,
     $           g_tmp1, ooff, oend, 1, nclosed,
     $           g_tmp2, voff, nmo, 1, nclosed )
            call ga_matmul_patch('t', 'n', one, one,
     $           g_tmp3, ooff, oend, voff, nmo,
     $           g_tmp1, voff, nmo, 1, nclosed,
     $           g_tmp2, ooff, oend, 1, nclosed )
c     
c       cv         pv     cp   pv
c     S    -= (2.F   +  F  ).x
c     
c       pv         pv     cp   cv
c     S    -= (2.F   +  F  ).x
c     
            call ga_dadd(one,g_fcp,two,g_fpv,g_tmp3)
            call ga_matmul_patch('n','n',mone,one,
     $           g_tmp1, voff, nmo, ooff, oend,
     $           g_tmp3, ooff, oend, 1, nclosed,
     $           g_tmp2, voff, nmo, 1, nclosed )
            call ga_matmul_patch('n','t',mone,one,
     $           g_tmp1, voff, nmo, 1, nclosed,
     $           g_tmp3, 1, nclosed, ooff, oend,
     $           g_tmp2, voff, nmo, ooff, oend )
c     
c       pv       cp     pv   cp
c     S    += (F   -  F  ).x
c     
c       cp       cp     pv   pv
c     S    += (F   -  F  ).x
c     
            call ga_dadd(one,g_fcp,mone,g_fpv,g_tmp3)
            call ga_matmul_patch('n', 't', one, one,
     $           g_tmp3, voff, nmo, 1, nclosed,
     $           g_tmp1, 1, nclosed, ooff, oend,
     $           g_tmp2, voff, nmo, ooff, oend )
            call ga_matmul_patch('t', 'n', one, one,
     $           g_tmp1, ooff, oend, voff, nmo,
     $           g_tmp3, voff, nmo, 1, nclosed,
     $           g_tmp2, ooff, oend, 1, nclosed )
c
         endif
c     
         xend = nclosed*nvir
c
         call ga_mat_to_vec(g_tmp2, voff, nmo, 1, nclosed,
     $       g_ax, 1, ivec, 1.0d0, '+')
*         call ga_dadd_patch( one, g_tmp2, voff, nmo, 1, nclosed,
*     $        one, g_ax, 1, xend, ivec, ivec,
*     $        g_ax, 1, xend, ivec, ivec)
c
         if (nopen .gt. 0) then
            xoff = xend + 1
            xend = xend + nopen*nvir
*
            call ga_mat_to_vec(g_tmp2, voff, nmo, ooff, oend,
     $         g_ax, xoff, ivec, 1.0d0, '+')
*            call ga_dadd_patch( one, g_tmp2, voff, nmo, ooff, oend,
*     $           one, g_ax, xoff, xend, ivec, ivec,
*     $           g_ax, xoff, xend, ivec, ivec)
            xoff = xend + 1
            xend = xend + nopen*nclosed
            call ga_mat_to_vec(g_tmp2, ooff, oend, 1, nclosed,
     $         g_ax, xoff, ivec, 1.0d0, '+')
*            call ga_dadd_patch( one, g_tmp2, ooff, oend, 1, nclosed,
*     $           one, g_ax, xoff, xend, ivec, ivec,
*     $           g_ax, xoff, xend, ivec, ivec)
         endif
      enddo
c
      if (.not. ga_destroy(g_tmp1)) call errquit('rohf_hv1?',0, GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('rohf_hv1?',0, GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_tmp3)) call errquit('rohf_hv1?',0,
     &       GA_ERR)
      endif
c
      end
