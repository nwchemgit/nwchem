      subroutine rohf_hdiag_scale( nmo, nclosed, nopen, lshift,
     $     g_fcv, g_fpv, g_fcp, g_x )
C     $Id$
c
c ... jochen: modified version. below was the check: 
c
c     denominator = 4.d0*(fcv(j+vbase) - fcv(i+cbase)) +
c    $                 lshift
c     if (denominator .lt. 0.1d0) denominator = 0.1d0
c     dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
c
c     where I changed the number to 0.00001 instead, and to print
c     a warning! (ONLY C-V PIECE)
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c     
      integer nmo, nclosed, nopen
      double precision lshift, small, sign
      parameter (small = 1d-5)
      integer g_fcv, g_fpv, g_fcp
      integer g_x
      logical warning
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
c
      integer ivec, nvec, gtype, vlen
c
      fcv(i) = dbl_mb(k_fcv + i - 1) ! Statement function definitions
      fpv(i) = dbl_mb(k_fpv + i - 1)
      fcp(i) = dbl_mb(k_fcp + i - 1)
c
c     =================================================================
c
      warning = .false.
c
      call ga_inquire(g_x,gtype,vlen,nvec)

c     
      nvir = nmo - nclosed - nopen
      cbase = 0                 ! Offsets into diagonals
      obase = nclosed
      vbase = nclosed + nopen
c     
c     Get diagonals of each matrix
c
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_fcv,k_fcv))
     $     call errquit('rohf_hdiag_scale: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_x,k_x))
     $     call errquit('rohf_hdiag_scale: cannot allocate',0, MA_ERR)
      if (nopen .gt. 0) then
         if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_fcp,k_fcp))
     $        call errquit('rohf_hdiag_scale: cannot allocate',0,
     &       MA_ERR)
         if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_fpv,k_fpv))
     $        call errquit('rohf_hdiag_scale: cannot allocate',0,
     &       MA_ERR)
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
      do ivec = 1, nvec
         if (nvir*nclosed .gt. 0) then
            do i = ga_nodeid()+1, nclosed, ga_nnodes()
               ioff = (i-1)*nvir + 1
               call ga_get(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x),nvir)
               do j=1,nvir
                  denominator = 4.d0*(fcv(j+vbase) - fcv(i+cbase)) +
     $                 lshift
                  if (abs(denominator).lt.small) then
                    sign = 1d0
                    if (denominator.lt.0d0) sign = -1d0
                    denominator = small * sign
                    if (ga_nodeid().eq.0 .and. .not.warning)
     &                 write (LuOut,9000) 
                    warning = .true.
                  endif
                  dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
               enddo
               call ga_put(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x),nvir)
            enddo
         endif
c     
         if (nvir*nopen .gt. 0) then
c     
c     Scale open-virtual piece
c     
            do i = ga_nodeid()+1, nopen, ga_nnodes()
               ioff = (nclosed+i-1)*nvir + 1
               call ga_get(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x),nvir)
               do j=1,nvir
                  denominator = 2.d0*(fpv(j+vbase) - fpv(i+obase)) +
     $                 lshift
                  if (denominator .lt. 0.1d0) denominator = 0.1d0
                  dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
               enddo
               call ga_put(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x),nvir)
            enddo
         endif
c     
c     Scale closed-open piece
c     
         if (nclosed*nopen .gt. 0) then
            do i = ga_nodeid()+1, nclosed, ga_nnodes()
               ioff = (nclosed+nopen)*nvir + (i-1)*nopen + 1
               call ga_get(g_x, ioff, ioff+nopen-1, ivec, ivec,
     $              dbl_mb(k_x),nopen)
               do j=1,nopen
                  denominator = 2.d0*(fcp(j+obase) - fcp(i+cbase)) +
     $                 lshift
                  if (denominator .lt. 0.1d0) denominator = 0.1d0
                  dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
               enddo
               call ga_put(g_x, ioff, ioff+nopen-1, ivec, ivec,
     $              dbl_mb(k_x),nopen)
            enddo
         endif
      enddo
c
      call ga_sync
c     
      if (nopen .gt. 0) then
         if (.not.ma_pop_stack(l_fpv))
     $        call errquit('rohf_hdiag_scale: cannot pop stack',0,
     &       MA_ERR)
         if (.not.ma_pop_stack(l_fcp))
     $        call errquit('rohf_hdiag_scale: cannot pop stack',0,
     &       MA_ERR)
      endif
      if (.not.ma_pop_stack(l_x))
     $     call errquit('rohf_hdiag_scale: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_fcv))
     $     call errquit('rohf_hdiag_scale: cannot pop stack',0, MA_ERR)
c
      return
c
 9000 format (/' *** WARNING: SMALL OCC-VIRT ORBITAL ENERGY DIFFERENCE',
     &   ' DETECTED'/' (Including frequency factor for response)'/
     &   ' This is going to lead to diverging respose properties'/
     &   ' I have adjusted the denominator to a value of 1E-5')
c
      end

c     =======================================================================

      subroutine rohf_hdiag_scale_damp( 
     &                 nmo, 
     &                 nclosed, 
     &                 nopen, 
     &                 lshift,
     &                 gamwidth,
     $                 g_fcv, 
     &                 g_fpv, 
     &                 g_fcp, 
     &                 g_x, 
     &                 g_x_im)
C     $Id$
c
c ... jochen: modified version of rohf_hdiag_scale to allow for
c     a situation where we have damping and a coupled set of equations.
c 
c     further, the original routine contained the check: 
c
c     denominator = 4.d0*(fcv(j+vbase) - fcv(i+cbase)) +
c    $                 lshift
c     if (denominator .lt. 0.1d0) denominator = 0.1d0
c     dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
c
c     where I changed the number to 0.00001 instead, and to print
c     a warning! (ONLY C-V PIECE)
c
c     NOTE: THIS ROUTINE PRESENTLY WORKS ONLY FOR A CLOSED SHELL SYSTEM
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c     
      integer nmo, nclosed, nopen
      double precision lshift, gamwidth, small, sign
      parameter (small = 1d-5)
      integer g_fcv, g_fpv, g_fcp
      integer g_x, g_x_im
      logical warning
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
c ... jochen: we need some memory for the imaginary part, too
      integer l_x_im, k_x_im
      double precision fcv, fcp, fpv, denominator, factor_re,
     &   factor_im, tmpre, tmpim,den
c
      integer ivec, nvec, gtype, vlen
c
      fcv(i) = dbl_mb(k_fcv + i - 1) ! Statement function definitions
      fpv(i) = dbl_mb(k_fpv + i - 1)
      fcp(i) = dbl_mb(k_fcp + i - 1)
c
c     =================================================================
c
      warning = .false.
c
      call ga_inquire(g_x,gtype,vlen,nvec)

c     
      nvir = nmo - nclosed - nopen
      cbase = 0                 ! Offsets into diagonals
      obase = nclosed
      vbase = nclosed + nopen
c     
c     Get diagonals of each matrix
c
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: fcv',l_fcv,k_fcv))
     $     call errquit('hdiag_scale_3im: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: x',l_x,k_x))
     $     call errquit('hdiag_scale_3im: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: x_im',l_x_im,k_x_im))
     $     call errquit('hdiag_scale_3im: cannot allocate',0, MA_ERR)

      if (nopen .gt. 0) then
         if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_fcp,k_fcp))
     $        call errquit('hdiag_scale_3im: cannot allocate',0,
     &       MA_ERR)
         if (.not.ma_push_get(MT_DBL,nmo,'rohf: tmp',l_fpv,k_fpv))
     $        call errquit('hdiag_scale_3im: cannot allocate',0,
     &       MA_ERR)
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
      do ivec = 1, nvec
         if (nvir*nclosed .gt. 0) then
            do i = ga_nodeid()+1, nclosed, ga_nnodes()
               ioff = (i-1)*nvir + 1
c
c              put elements of g_x, g_x_im to local memory:
               call ga_get(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x),nvir)
               call ga_get(g_x_im, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x_im),nvir)
c
               do j=1,nvir
c                orbital energy differences:
c                denominator = (fcv(j+vbase) - fcv(i+cbase)) + lshift ! FA-03-10-14
                 denominator = 4.d0*(fcv(j+vbase) - fcv(i+cbase)) +
     &              lshift

c                check for small HOMO-LUMO gap:
                 if (abs(denominator).lt.small) then
                   sign = 1d0
                   if (denominator.lt.0d0) sign = -1d0
                   denominator = small * sign
                   if (ga_nodeid().eq.0 .and. .not.warning)
     &                write (LuOut,9000) 
                   warning = .true.
                 endif
c                
c                take care of the damping here: set up scale_factor
c                and multiply instead of dividing by denominator                
c                FA-03-12-14 ------- START
c                 den=(denominator**2 + gamwidth**2)
c                 factor_re =4.0d0*denominator/den
c                 factor_im =4.0d0*gamwidth/den
c                FA-03-12-14 ------- END

                 factor_re =
     &             denominator / (denominator**2 + gamwidth**2)
                 factor_im = gamwidth / (denominator**2 + gamwidth**2)

                 
                 tmpre = dbl_mb(k_x   +j-1)
                 tmpim = dbl_mb(k_x_im+j-1)

                 dbl_mb(k_x   +j-1) =
     &             tmpre * factor_re +
     &             tmpim * factor_im
                 dbl_mb(k_x_im+j-1) =
     &             tmpim * factor_re -
     &             tmpre * factor_im
                 
c                 dbl_mb(k_x   +j-1)=tmpre*factor_re+tmpim*factor_im
c FA-03-12-14 ------- START
c                 dbl_mb(k_x_im+j-1)=tmpim*factor_re-tmpre*factor_im
c FA-03-12-14 ------- END
               enddo
               
               call ga_put(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x),nvir)
               call ga_put(g_x_im, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x_im),nvir)
            enddo
         endif
c     
         if (nvir*nopen .gt. 0) then
c          
           call errquit ('rohf_hscale_damp open shell not implemented',
     &        0, UNKNOWN_ERR)
c          THIS PIECE OF CODE SHOULD NOT BE EXECUTED IN RESPONSE CALC
C          BECAUSE THIS HAS NOT BEEN ADAPTED TO FREQUENCY AND DAMPING
c          
c     Scale open-virtual piece
c     
            do i = ga_nodeid()+1, nopen, ga_nnodes()
               ioff = (nclosed+i-1)*nvir + 1
               call ga_get(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x),nvir)
               do j=1,nvir
                  denominator = 2.d0*(fpv(j+vbase) - fpv(i+obase)) +
     $                 lshift
                  if (denominator .lt. 0.1d0) denominator = 0.1d0
                  dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
               enddo
               call ga_put(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $              dbl_mb(k_x),nvir)
            enddo
         endif
c     
c     Scale closed-open piece
c     
         if (nclosed*nopen .gt. 0) then
           call errquit ('rohf_hscale_damp open shell not implemented',
     &        0, UNKNOWN_ERR)
c          THIS PIECE OF CODE SHOULD NOT BE EXECUTED IN RESPONSE CALC
C          BECAUSE THIS HAS NOT BEEN ADAPTED TO FREQUENCY AND DAMPING
           do i = ga_nodeid()+1, nclosed, ga_nnodes()
             ioff = (nclosed+nopen)*nvir + (i-1)*nopen + 1
             call ga_get(g_x, ioff, ioff+nopen-1, ivec, ivec,
     $          dbl_mb(k_x),nopen)
             do j=1,nopen
               denominator = 2.d0*(fcp(j+obase) - fcp(i+cbase)) +
     $            lshift
               if (denominator .lt. 0.1d0) denominator = 0.1d0
               dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
             enddo
             call ga_put(g_x, ioff, ioff+nopen-1, ivec, ivec,
     $          dbl_mb(k_x),nopen)
           enddo
         endif
       enddo
c
      call ga_sync
c     
      if (nopen .gt. 0) then
         if (.not.ma_pop_stack(l_fpv))
     $        call errquit('rohf_hdiag_scale_im3: cannot pop stack fpv',
     &     0, MA_ERR)
         if (.not.ma_pop_stack(l_fcp))
     $        call errquit('rohf_hdiag_scale_im3: cannot pop stack fcp',
     &      0, MA_ERR)
      endif
      if (.not.ma_pop_stack(l_x_im))
     $     call errquit('rohf_hdiag_scale_im3: cannot pop stack xim',
     &   0, MA_ERR)
      if (.not.ma_pop_stack(l_x))
     $     call errquit('rohf_hdiag_scale_im3: cannot pop stack lx',
     &   0, MA_ERR)
      if (.not.ma_pop_stack(l_fcv))
     $     call errquit('rohf_hdiag_scale_im3: cannot pop stack fcv',
     &   0, MA_ERR)
c
      return
c
 9000 format (/' *** WARNING: SMALL OCC-VIRT ORBITAL ENERGY DIFFERENCE',
     &   ' DETECTED'/' (Including frequency factor for response)'/
     &   ' This is going to lead to diverging respose properties'/
     &   ' I have adjusted the denominator to a value of 1E-5')
c
      end

