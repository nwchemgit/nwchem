      subroutine rohf_hdiag_scale( nmo, nclosed, nopen, lshift,
     $     g_fcv, g_fpv, g_fcp, g_x )
C     $Id: rohf_hscale.F,v 1.3 1996-06-27 22:42:41 d3g681 Exp $
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
c     
      integer nmo, nclosed, nopen
      double precision lshift
      integer g_fcv, g_fpv, g_fcp
      integer g_x
c     
c     
c     Scale the vector x by the ROHF 1e-Hessian diagonal
c     
c               x                    x
c                pq                   ia
c     x   =  --------           = ---------    [ closed-virtual ]
c      pq      (1e)                 cv    cv
c             H                   4F  - 4F
c              pq,pq                aa    ii
c     
c               x
c                pa
c         = ---------    [ partial-virtual ]
c             pv    pv
c           2F  - 2F
c             aa   pp
c     
c               x
c                ip
c         = ---------    [ closed-partial ]
c             cp    cp
c           2F  - 2F
c             pp    ii
c     
      integer nvir
      integer i, j, ioff, cbase, obase, vbase
      integer l_fcv, k_fcv
      integer l_fcp, k_fcp
      integer l_fpv, k_fpv
      integer l_x, k_x
      double precision fcv, fcp, fpv, denominator
      fcv(i) = dbl_mb(k_fcv + i - 1)
      fpv(i) = dbl_mb(k_fpv + i - 1)
      fcp(i) = dbl_mb(k_fcp + i - 1)
c     
      nvir = nmo - nclosed - nopen
      cbase = 0                 ! Offsets into diagonals
      obase = nclosed
      vbase = nclosed + nopen
c     
c     Get diagonals of each matrix
c
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_fcv,k_fcv))
     $     call errquit('rohf_hdiag_scale: cannot allocate',0)
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_x,k_x))
     $     call errquit('rohf_hdiag_scale: cannot allocate',0)
      if (nopen .gt. 0) then
         if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_fcp,k_fcp))
     $        call errquit('rohf_hdiag_scale: cannot allocate',0)
         if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_fpv,k_fpv))
     $        call errquit('rohf_hdiag_scale: cannot allocate',0)
      endif
c     
      call ga_get_diagonal(g_fcv, dbl_mb(k_fcv))
      if (nopen .gt. 0) then
         call ga_get_diagonal(g_fcp, dbl_mb(k_fcp))
         call ga_get_diagonal(g_fpv, dbl_mb(k_fpv))
      endif
c
c     Scale closed-virtual piece
c
      if (nvir*nclosed .gt. 0) then
         do i = ga_nodeid()+1, nclosed, ga_nnodes()
            ioff = (i-1)*nvir + 1
            call ga_get(g_x, ioff, ioff+nvir-1, 1, 1, dbl_mb(k_x),1)
            do j=1,nvir
               denominator = 4.d0*(fcv(j+vbase) - fcv(i+cbase)) + lshift
               if (denominator .lt. 0.1d0) denominator = 0.1d0
               dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
            enddo
            call ga_put(g_x, ioff, ioff+nvir-1, 1, 1, dbl_mb(k_x),1)
         enddo
      endif
c
      if (nvir*nopen .gt. 0) then
c     
c     Scale open-virtual piece
c     
         do i = ga_nodeid()+1, nopen, ga_nnodes()
            ioff = (nclosed+i-1)*nvir + 1
            call ga_get(g_x, ioff, ioff+nvir-1, 1, 1, dbl_mb(k_x),1)
            do j=1,nvir
               denominator = 2.d0*(fpv(j+vbase) - fpv(i+obase)) + lshift
               if (denominator .lt. 0.1d0) denominator = 0.1d0
               dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
            enddo
            call ga_put(g_x, ioff, ioff+nvir-1, 1, 1, dbl_mb(k_x),1)
         enddo
      endif
c     
c     Scale closed-open piece
c     
      if (nclosed*nopen .gt. 0) then
         do i = ga_nodeid()+1, nclosed, ga_nnodes()
            ioff = (nclosed+nopen)*nvir + (i-1)*nopen + 1
            call ga_get(g_x, ioff, ioff+nopen-1, 1, 1, dbl_mb(k_x),1)
            do j=1,nopen
               denominator = 2.d0*(fcp(j+obase) - fcp(i+cbase)) + lshift
               if (denominator .lt. 0.1d0) denominator = 0.1d0
               dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
            enddo
            call ga_put(g_x, ioff, ioff+nopen-1, 1, 1, dbl_mb(k_x),1)
         enddo
      endif
c     
      if (nopen .gt. 0) then
         if (.not.ma_pop_stack(l_fpv))
     $        call errquit('rohf_hdiag_scale: cannot pop stack',0)
         if (.not.ma_pop_stack(l_fcp))
     $        call errquit('rohf_hdiag_scale: cannot pop stack',0)
      endif
      if (.not.ma_pop_stack(l_x))
     $     call errquit('rohf_hdiag_scale: cannot pop stack',0)
      if (.not.ma_pop_stack(l_fcv))
     $     call errquit('rohf_hdiag_scale: cannot pop stack',0)
c
      end

