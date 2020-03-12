      subroutine rohf_hessv2( acc, g_x, g_ax )
      implicit none
#include "global.fh"
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "stdio.fh"
#include "util.fh"
c     
c     $Id$
c
      integer g_x, g_ax
      double precision acc
c     
      integer gtype,grow,gcol,growp,gcolp
      integer icol,ilo(2),ihi(2)
      logical oprint, olprint
      double precision dnrm
c     
c     Check
c     
      oprint = util_print('hessv',print_high)
      olprint = oprint.and.(ga_nodeid().eq.0)
      if (crohf_init_flag.ne.1)
     $     call errquit('rohf_hessv: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_ax,gtype,growp,gcolp)
      if (growp.ne.crohf_vlen)
     $     call errquit('rohf_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      if (gcol.ne.gcolp)
     $     call errquit('rohf_hessv: invalid no. of vectors',0,
     &       UNKNOWN_ERR)
c
      if (oprint) then
        do icol = 1, gcol
          ilo(1) = 1
          ilo(2) = icol
          ihi(1) = grow
          ihi(2) = icol
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": g_x = ",i4,f24.8)')
     +      'rohf_hessv2',icol,dnrm
          endif
        enddo
      endif
c     
c     Call internal routine
c     
      call rohf_hessv_xx2( basis, geom, nbf, nmo,
     $     nclosed, nopen,
     $     pflg, g_movecs, oskel, noskew,
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $     acc, lshift, g_x, g_ax )
c
c     Zap numbers much smaller than acc to ensure hard zeroes 
c     remain unpolluted ... cannot use a threshold larger than the
c     integral accuracy since can break symmetry in non-abelian groups
c     Also must ensure that the threshold tends to zero to permit
c     tight convergence.
c
*      call ga_screen(g_ax, acc*0.01d0)
      call ga_screen(g_ax, min(acc*acc,acc*0.01d0,1d-12))
c
      if (oprint) then
        do icol = 1, gcolp
          ilo(1) = 1
          ilo(2) = icol
          ihi(1) = growp
          ihi(2) = icol
          call nga_normf_patch(g_ax,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": g_ax = ",i4,f24.8)')
     +      'rohf_hessv2',icol,dnrm
          endif
        enddo
      endif
c
      end
      subroutine rohf_hessv_xx2( basis, geom, nbf, nmo, nclosed, nopen, 
     $     pflg,
     $     g_movecs, oskel, noskew, g_fcv, g_fpv, g_fcp,
     $     acc, lshift, g_x, g_ax )
C     $Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
c     
      integer basis, geom
      integer nbf, nmo, nclosed, nopen
      integer pflg
      integer g_movecs(*)
      logical oskel, noskew
      integer g_fcv, g_fpv, g_fcp
      double precision acc
      double precision lshift
      integer g_x, g_ax, vlen, nvec, g_tmp, gtype
      call ga_zero(g_ax)
      if (pflg.gt.2 .or. pflg.le.0) then
         call errquit('rohf_hessv_xx: pflg invalid ', pflg,
     &       UNKNOWN_ERR)
      endif
      if (pflg .gt. 0) then
         call ga_dadd( lshift, g_x, 0.d0, g_ax, g_ax )
         call rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $        g_fcv, g_fpv, g_fcp,
     $        g_x, g_ax )
      endif
 200  continue
      if (pflg .gt. 1) then
         call rohf_hessv_2e2( basis, geom, nbf, nmo, nclosed, nopen,
     $        g_movecs, oskel, noskew, g_x, acc, ! was min(1d-6,acc)
     $        g_ax )
      endif
c
c     
      end
      subroutine rohf_hessv_2e2( basis, geom, nbf, nmo, nclosed, nopen,
     $     g_movec, oskel, noskew, g_x, acc, g_ax )
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "cscfps.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "case.fh"
c     
c     Return the ROHF orbital 2e-Hessian vector product, g_ax = A * g_x
c     
ccccccccccccccc This code does NOT work for open shell!!!!!ccccccccccccccccc
c
c
      integer basis, geom       ! basis & geom handle
      integer nbf, nclosed, nopen ! Basis size and occupation
      integer nmo               ! No. of linearly dependent MOs
      integer g_movec(*)        ! MO coefficients
      logical oskel
      integer g_x               ! Argument
      double precision acc      ! Accuracy of "Fock" construction
      integer g_ax              ! Hessian product
c     
      integer nvir, voff, xoff, xend
      integer ooff, oend
c     
      integer nfock, ifock
      integer g_tmp1, g_tmp2, g_tmp3
c
      integer ivec, nvec, gtype, vlen
      integer k_jfac,l_jfac,k_kfac,l_kfac
c
      integer g_dcv,  g_dcp,  g_dpv, g_dens, g_fock
      double precision tol2e
      logical odebug
      logical oprint
      logical noskew
c
      integer g_work
c
      logical ldimqm
      integer g_dim
c
      integer dims(3), chunk(3), alo(3), ahi(3), blo(2), bhi(2)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision one, zero, mone, four, half, mhalf, two, mtwo
      double precision itol_floor, itol_ceil
      double precision xc_xfac1
      logical docoul
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0, mtwo=-2.0d0)
      oprint= util_print('rohf_hessv2',print_debug)
      xc_xfac1=bgj_kfac()
c
c     This get's cleaned up a lot when the Fock build accepts
c     a single multiple dimension GA for its input.
c
c     Lots more to do on memory management.
c
      if (nopen.ne.0) call errquit
     $     ('rohf_h2e2: does not work for open shells',nopen,
     &       UNKNOWN_ERR)
c
      call ga_inquire(g_x,gtype,vlen,nvec)
      odebug = util_print('rohf_hessv', print_debug)
c     
      tol2e = min(max(acc,itol_floor),itol_ceil)
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
      ooff = nclosed + 1
      oend = nclosed + nopen
c
c     allocate Mas for jfacs
c
      if (.not.MA_Push_Get(mt_dbl,nvec*2,'hessv jfacs',
     &     l_jfac,k_jfac))
     &     call errquit('hessv2: cannot allocate jfac',nvec*2,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,nvec*2,'hessv kfacs',
     &     l_kfac,k_kfac))
     &     call errquit('hessv2: cannot allocate kfac',nvec*2,
     &       MA_ERR)

c     Allocate global for work space
c     
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp2')
      if (nopen .gt. 0) then
         g_tmp3 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp3')
      endif
c
c     Allocate arrays and generate the AO densities
c
      dims(1) = nvec
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
c
      if (.not. nga_create (MT_DBL, 3, dims, 'CPKS dens',chunk, g_dens))
     & call errquit('rohf_h2e: could not allocate g_dens',555,GA_ERR)
      call ga_zero(g_dens)
c
      if (.not. nga_create (MT_DBL, 3, dims, 'CPKS dens',chunk, g_work))
     & call errquit('rohf_h2e: could not allocate g_work',555,GA_ERR)
      call ga_zero(g_work)
c
      g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dcv')
      if (nopen .gt. 0) then
         g_dcp = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dcp')
         g_dpv = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dpv')
      endif
c
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      nfock = 0
      do ivec = 1, nvec
         nfock          = nfock + 1
         if (noskew) then
            dbl_mb(k_jfac+nfock-1)    =  1.0d0   ! CV
         else
            dbl_mb(k_jfac+nfock-1)    =  0.0d0   ! CV
         endif
         dbl_mb(k_kfac+nfock-1)    = -0.5d0*xc_xfac1   ! CV
c        jfac(nfock)    =  1.0d0   ! CV
c        kfac(nfock)    = 1.0d0   ! CV
c        jfac(nfock)    =  4.0d0   ! CV
c        kfac(nfock)    = -1.0d0   ! CV
         if (nopen .gt. 0) then
            nfock          = nfock + 1
            if (noskew) then
               dbl_mb(k_jfac+nfock-1)    =  1.0d0   ! CP
            else
               dbl_mb(k_jfac+nfock-1)    =  0.0d0   ! CP
            endif
            dbl_mb(k_kfac+nfock-1)    = -0.5d0*xc_xfac1 ! CP
            nfock          = nfock + 1
            dbl_mb(k_jfac+nfock-1) =  0.0d0 ! PV
            dbl_mb(k_kfac+nfock-1) =  0.5d0*xc_xfac1 ! PV
         endif
c     
c     Compute CV, PV & CP "densities" from argument vector
c     
         call ga_zero(g_tmp1)
         call ga_zero(g_tmp2)
         call ga_zero(g_dcv)
         call ga_vec_to_mat(g_tmp1, 1, nvir, 1, nclosed,
     $                      g_x, 1, ivec)
         call ga_matmul_patch('n','t',one,zero,
     $        g_tmp1, 1, nvir, 1, nclosed,
     $        g_movec, 1, nclosed, 1, nbf,
     $        g_tmp2, 1, nvir, 1, nbf )
         call ga_matmul_patch('n','n',one,zero,
     $        g_movec, 1, nbf, voff, nmo,
     $        g_tmp2, 1, nvir, 1, nbf,
     $        g_dcv, 1, nbf, 1, nbf)
c
         if (nopen .gt. 0) then
c     
c     PV density
c     
            xoff = nclosed*nvir + 1
            xend = nclosed*nvir + nopen*nvir
            call ga_zero(g_tmp1)
            call ga_zero(g_tmp2)
            call ga_zero(g_dpv)
            call ga_copy_patch('n',
     $           g_x, xoff, xend, ivec, ivec,
     $           g_tmp1, 1, nvir, 1, nopen )
            call ga_matmul_patch('n', 't', one, zero,
     $           g_tmp1, 1, nvir, 1, nopen,
     $           g_movec, ooff, oend, 1, nbf,
     $           g_tmp2, 1, nvir, 1, nbf )
            call ga_matmul_patch('n', 'n', one, zero,
     $           g_movec, 1, nbf, voff, nmo,
     $           g_tmp2, 1, nvir, 1, nbf,
     $           g_dpv, 1, nbf, 1, nbf )
c
c     CP density
c     
            call ga_zero(g_tmp1)
            call ga_zero(g_tmp2)
            call ga_zero(g_dcp)
            xoff = (nclosed+nopen)*nvir + 1
            xend = (nclosed+nopen)*nvir + nclosed*nopen
            call ga_copy_patch('n', g_x, xoff, xend, ivec, ivec,
     $           g_tmp1, 1, nopen, 1, nclosed )
            call ga_matmul_patch('n', 't', one, zero,
     $           g_tmp1, 1, nopen, 1, nclosed,
     $           g_movec, 1, nclosed, 1, nbf,
     $           g_tmp2, 1, nopen, 1, nbf )
            call ga_matmul_patch('n', 'n', one, zero,
     $           g_movec, 1, nbf, ooff, oend,
     $           g_tmp2, 1, nopen, 1, nbf,
     $           g_dcp, 1, nbf, 1, nbf )
c
         endif
c     
c     We actually use dcp+dpv and dcp-dpv to build fock matrix
c     
         if (nopen .gt. 0) then
            call ga_dadd(one, g_dcp, one, g_dpv, g_tmp1)
            call ga_dadd(one, g_dcp,mone, g_dpv, g_tmp2)
            call ga_copy(g_tmp1, g_dcp)
            call ga_copy(g_tmp2, g_dpv)
         endif
c
         
         if (noskew) then
            call ga_symmetrize(g_dcv)
         else
            call ga_antisymmetrize(g_dcv)
         endif
         if (nopen .gt. 0) then
            if (noskew) then
               call ga_symmetrize(g_dpv)
               call ga_symmetrize(g_dcp)
            else
               call ga_antisymmetrize(g_dpv)
               call ga_antisymmetrize(g_dcp)
            endif
         endif
c
         alo(1) = ivec
         ahi(1) = ivec
         call nga_copy_patch('N',g_dcv,blo,bhi,g_dens,alo,ahi)
         if (nopen .gt. 0) then
c Note that this is NOT correct!
            call nga_copy_patch('N',g_dpv,blo,bhi,g_dens,alo,ahi)
            call nga_copy_patch('N',g_dcp,blo,bhi,g_dens,alo,ahi)
         endif
      enddo
c
         if(oprint) then
            call ga_print(g_dens)
         endif
      if (.not. ga_destroy(g_tmp1)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_dcv)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_tmp3)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_dcp)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_dpv)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      endif
c
c     Fock build
c     
      if (.not. nga_create (MT_DBL, 3, dims, 'Fockv',chunk,
     &     g_fock))
     &   call errquit('rohf_h2e: could not allocate g_fock',555,
     &       GA_ERR)
      call ga_zero(g_fock)
c
c     noskew is true, except for nmr perturbations (skew symmetric)
c
      if (.not.cam_exch) then
        call shell_fock_build(geom, basis, 0, nfock, 
     $   dbl_mb(k_jfac), dbl_mb(k_kfac),
     $   tol2e, g_dens, g_fock, noskew)
      else   ! for attenuated calculations
        call shell_fock_build_cam(geom, basis, 0, nfock, 
     $   dbl_mb(k_jfac), dbl_mb(k_kfac),
     $   tol2e, g_dens, g_fock, noskew)
      end if
c     
      if (oskel) 
     $   call errquit('rohf_h2e2: no way',0, UNKNOWN_ERR)
c
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp2')
      g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dcv')
      if (nopen .gt. 0) then
         g_tmp3 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp3')
         g_dcp = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dcp')
         g_dpv = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dpv')
      endif
c
c     DIM/QM JEM
      if (.not.rtdb_get(bgj_get_rtdb_handle(), 'dimqm:lrsp', mt_log,
     $                  1, ldimqm)) ldimqm = .false.
c     Calculate DIM dipoles and DIM potential
      if (ldimqm) then
c  jbecca START: added false variable to subroutine call
        call dimqm_rohf_calcDIM(bgj_get_rtdb_handle(), geom, basis,
     $                          nbf, nfock, g_dens, g_dim, .false.)
c jbecca END
      end if

      ifock = 0
      do ivec = 1, nvec
         ifock = ifock + 1      ! History ... copy into tmp arrays
         alo(1) = ivec
         ahi(1) = ivec
         call nga_copy_patch('N',g_fock,alo,ahi,g_tmp1,blo,bhi)
         call nga_copy_patch('N',g_dens,alo,ahi,g_dcv,blo,bhi)
         if (nopen .gt. 0) then 
            ifock = ifock + 1
            alo(1) = ifock
            ahi(1) = ifock
            call nga_copy_patch('N',g_fock,alo,ahi,g_tmp2,blo,bhi)
            call nga_copy_patch('N',g_dens,alo,ahi,g_dcp,blo,bhi)
            ifock = ifock + 1
            alo(1) = ifock
            ahi(1) = ifock
            call nga_copy_patch('N',g_fock,alo,ahi,g_tmp3,blo,bhi)
            call nga_copy_patch('N',g_dens,alo,ahi,g_dpv,blo,bhi)
         endif
c     
c     P      =  4(ij|kl) - (ik|jl) - (il|kj)
c     ij,kl
c     
c     K      =  (ik|jl) + (il|kj)
c     ij,kl
c     
c     cv         cv          pv   cp
c     Z   =  2P.[D  ]  +  P.[D  + D  ]
c     
c     pv          cv           cp   pv
c     Z   =  0.5d0*Z   + 0.5*K.[D  - D  ]
c     
c     cp          cv           cp   pv
c     Z   =  0.5d0*Z   - 0.5*K.[D  - D  ]
c     
c     Add the Fock matrices together overwriting the density
c     matrices to form the results above
c     
         call ga_copy(g_tmp1, g_dcv)
         call ga_dscal(g_dcv, two)
         if (nopen .gt. 0) then
            call ga_dadd(one,  g_dcv,  one, g_tmp2, g_dcv)
            call ga_dadd(half, g_dcv, half, g_tmp3, g_dpv)
            call ga_dadd(half, g_dcv,mhalf, g_tmp3, g_dcp)
         endif
c     
c     Closed-Virtual bit
c     
         call ga_matmul_patch('n','n',two,zero,
     $        g_dcv, 1, nbf, 1, nbf,
     $        g_movec, 1, nbf, 1, nclosed,
     $        g_tmp1, 1, nbf, 1, nclosed )
         call ga_matmul_patch('t','n',one,zero,
     $        g_movec, voff, nmo, 1, nbf, 
     $        g_tmp1, 1, nbf, 1, nclosed,
     $        g_tmp2, 1, nvir, 1, nclosed )
         xoff = 1
         xend = nvir*nclosed
         call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed,
     $        g_ax, xoff, ivec, four, '+')
*         call ga_dadd_patch( four, g_tmp2, 1, nvir, 1, nclosed,
*     $        one, g_ax, xoff, xend, ivec, ivec,
*     $        g_ax, xoff, xend, ivec, ivec )
c
c     DIM/QM JEM
c     Add DIM potential to the Fock Matrix
c        if (ldimqm) then
c          call dimqm_rohf_addDIM(rtdb, g_ax, ivec, xoff, voff,
c     $                           nclosed, nvir, nmo, nbf, 
c     $                           g_dcv, g_movec, g_tmp1, g_tmp2,
c     $                           alo, ahi, blo, bhi)
c        end if
         if(ldimqm) then
           call nga_copy_patch('N',g_dim,alo,ahi,g_dcv,blo,bhi)
           call ga_scale(g_dcv, TWO)
           call ga_matmul_patch('n', 'n', TWO, zero,
     $                           g_dcv,   1, nbf, 1, nbf,
     $                           g_movec, 1, nbf, 1, nclosed,
     $                           g_tmp1,  1, nbf, 1, nclosed)
           call ga_matmul_patch('t', 'n', one, zero,
     $                           g_movec, voff, nmo, 1, nbf,
     $                           g_tmp1, 1, nbf,  1, nclosed,
     $                           g_tmp2, 1, nvir, 1, nclosed)
           call ga_sync()
           call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_ax,
     $                         xoff, ivec, four, '+')
         end if
c
c     
         if (nopen .gt. 0) then
c     
c     Open-Virtual bit
c     
            call ga_matmul_patch('n','n',two,zero,
     $           g_dpv, 1, nbf, 1, nbf,
     $           g_movec, 1, nbf, ooff, oend,
     $           g_tmp1, 1, nbf, 1, nopen )
            call ga_matmul_patch('t','n',one,zero,
     $           g_movec, voff, nmo, 1, nbf,
     $           g_tmp1, 1, nbf, 1, nopen,
     $           g_tmp2, 1, nvir, 1, nopen )
            xoff = xend + 1
            xend = xend + nvir*nopen
            call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nopen,
     $          g_ax, xoff, ivec, four, '+')
*            call ga_dadd_patch( four, g_tmp2, 1, nvir, 1, nopen,
*     $           one, g_ax, xoff, xend, ivec, ivec,
*     $           g_ax, xoff, xend, ivec, ivec)
c     
c     Closed-Open bit
c     
            call ga_matmul_patch('n','n',two,zero,
     $           g_dcp, 1, nbf, 1, nbf,
     $           g_movec, 1, nbf, 1, nclosed,
     $           g_tmp1, 1, nbf, 1, nclosed )
            call ga_zero(g_tmp2)
            call ga_matmul_patch('t','n',one,zero,
     $           g_movec, ooff, oend, 1, nbf,
     $           g_tmp1, 1, nbf, 1, nclosed,
     $           g_tmp2, 1, nopen, 1, nclosed )
            xoff = xend + 1
            xend = xend + nopen*nclosed
            call ga_mat_to_vec(g_tmp2, 1, nopen, 1, nclosed,
     $          g_ax, xoff, ivec, four, '+')
*            call ga_dadd_patch( four, g_tmp2, 1, nopen, 1, nclosed,
*     $           one, g_ax, xoff, xend, ivec, ivec,
*     $           g_ax, xoff, xend, ivec, ivec)
         endif
      enddo
c    DIM/QM JEM
      if (ldimqm) then
        if (.not.ga_destroy(g_dim))
     &    call errquit('hessv_2e2: problem destroying g_dim array',1,
     &                 GA_ERR)
      end if
      if (.not.ma_chop_stack(l_jfac))
     &   call errquit('hessv2: cannot chop stack',0, MA_ERR)
c
      if (.not. ga_destroy(g_dens)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_work)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp1)) call errquit('rohf:destroy?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('rohf:destroy?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_dcv)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_tmp3)) call errquit('rohf:destroy?',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_dcp)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_dpv)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      endif
c     
      end
