C
C $Id$
C
      subroutine twodd_cont(rtdb)
c
c  Add the two electron derivative (first and seconds) contributions to 
c  the hessian matrix and the RHS.  Part of this is stolen from grad_force
c  to set up the explicit hessian contributions.
c
      implicit none
c
#include "errquit.fh"
#include "bas.fh"
#include "hess_info.fh"
#include "bgj.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "case.fh"
c
      integer rtdb, i
c
      integer nshblocks
c
      integer blen              ! bf-blocking size for shell ordering
      integer maxblen           ! maximum value for blen = max no. bf in group
      integer maxsh             ! maximum no. of shells in a group
      integer maxq              ! max quartets in a request
      parameter (maxblen=36, maxsh=10)
      parameter (maxq=maxsh**4)
c
      integer
     $     lh_ij,  ld_ij, lh_kl, ld_kl, lh_ik, ld_ik,
     $     lh_jl,  ld_jl, lh_il, ld_il,
     $     lh_jk,  ld_jk,
     $     lh_ij2, ld_ij2, lh_kl2, ld_kl2, lh_ik2, ld_ik2,
     $     lh_jl2, ld_jl2, lh_il2, ld_il2, lh_jk2, ld_jk2,
     $     ld_ij3, ld_kl3, ld_ik3,
     $     ld_jl3, ld_il3, ld_jk3,
     $     ld_ij4, ld_kl4, ld_ik4,
     $     ld_jl4, ld_il4, ld_jk4,
     $     ld_ik5, ld_jl5, ld_il5, ld_jk5,
     $     ld_ik6, ld_jl6, ld_il6, ld_jk6,
     $     ld_ik7, ld_jl7, ld_il7, ld_jk7,
     $     ld_ik8, ld_jl8, ld_il8, ld_jk8
      integer max_at_bf, max_sh_bf, lsqa, nsh,
     $ max1e, max2e, mscratch_1e, mscratch_2e, lbuf, lscratch
      integer l_shmap, l_shglo, l_shghi, l_bfmap, l_rbfmap,
     $     l_bfglo, l_bfghi, l_labels, l_list, l_q4, l_bftoat
      integer k_shmap, k_shglo, k_shghi, k_bfmap, k_rbfmap,
     $     k_bfglo, k_bfghi, k_labels, k_list, k_q4, k_bftoat
      integer l_shbflo, l_shbfhi
      integer k_shbflo, k_shbfhi
      integer g_dens2
      integer nfock
      integer dims(3), chunk(3), alo(3), ahi(3), blo(2), bhi(2)
c
c dummy variable needed for grad2 call, but not used here
c
      integer g_force
c
      integer l_buf, k_buf, l_scr, k_scr, iat, nproc
c     pdm2d is now used by all methods.  pdm2/a/b/c only by MCSCF
      integer l_pdm2d ! MCSCF 2-pdm
      integer k_pdm2, k_pdm2a, k_pdm2b, k_pdm2c, k_pdm2d ! MCSCF 2-pdm
      integer k_coeff  ! MCSCF local copy of Active space MOs
c
      logical omp2, oprint, otwoprint,oprintsome
c
      double precision tol2e,t0
      double precision one, mone
      parameter ( mone = -1.D0 )
      parameter ( one  =  1.D0 )
c
c Set up printing
c
      oprint = util_print('hess_follow',print_high)
      otwoprint = util_print('twodd_cont',print_debug)
      oprintsome = util_print('hess-follow',print_low)
      oprintsome = oprintsome .and. (ga_nodeid().eq.0)
c
c Need to do a bit of cleanup so DFT and SCF work OK
c
      call schwarz_tidy()    ! need to get rid of temp buffers
      call intdd_terminate()
c
c  First let's do the DFT part.
c  For now we send in seperate memory, k_hess_xc and g_rhs_xc, to make sure
c  that the xc contributions are coming out O.K.  Eventually, the xc routines
c  will need to modified to add their contribution to k_exy and g_rhs and
c  not zero out the matrices when they enter the routines.
c  I will add the results into the actual "permanent" storage for now so
c  the rest of the code will work properly.
c 
      if (oprint.and.(ga_nodeid().eq.0)) 
     *     write (LuOut,*) 'HESS: DFT contribution'
      if (theory.eq.'dft') then
         call int_init(rtdb, nbasis, basis_list)
         call xc_hessian(geom,rtdb,nbf,basis,
     *        g_dens, nactive, log_mb(k_act), 
     *        dbl_mb(k_hess_xc),scftype)
         call xc_cpks_rhs(geom,rtdb,nbf,basis,
     *        g_dens,  nactive, log_mb(k_act),g_rhs_xc,scftype)
         call int_terminate()
         if (ocdfit) then
           call int_init(rtdb, nbasis, basis_list) !should this be intd?
           call schwarz_init(geom, basis)
           call j_cpks_rhs(g_dens, log_mb(k_act), nactive,
     *          g_rhs_xc)
           call schwarz_tidy()   ! need to get rid of temp buffers
           call int_terminate()  !should this be intd?
         endif
         if (otwoprint) then
           if(ga_nodeid().eq.0) then
             write(LuOut,*)
             write(LuOut,*) 'The Hessian dft contribution'
             call hess_hssout(dbl_mb(k_hess_xc),n3xyz,n3xyz,n3xyz)
             if(skelsym) then
               call sym_hess_symmetrize(geom,dbl_mb(k_hess_xc))  
               write(LuOut,*) 'The Hessian dft contribution symmetrized'
               call hess_hssout(dbl_mb(k_hess_xc),n3xyz,n3xyz,n3xyz)
             endif
             write(LuOut,*)
             write(LuOut,*) 'The g_rhs dft contribution'
             call util_flush(LuOut)
           endif
           call ga_print(g_rhs_xc)
           if(ga_nodeid().eq.0)write(LuOut,*)
           if(skelsym) then
             if(ga_nodeid().eq.0) then
               write(LuOut,*) 'symmetrized'
               call util_flush(LuOut)
             endif
             call  sym_symmetrize(geom,basis,.false.,g_rhs_xc)
             call ga_print(g_rhs_xc)
           endif
         endif
c
c        The next statements are only temporarily in until we are
c        confident the dft code is giving the correct contributions.
c
         call ga_add(1.0d+00,g_rhs_xc,1.0d+00,g_rhs,g_rhs)
         call daxpy(n3xyz2,1.0d+00,dbl_mb(k_hess_xc),1,
     *        dbl_mb(k_exy),1)
      endif
c
c     Release the global fock derivative matrices
c
      if (.not.ga_destroy(g_rhs_xc)) 
     *  call errquit('twodd_cont: problem destroying g_rhs_xc',1,
     &       GA_ERR)
c
c    Print out the hessian before adding anything else
c
      if (otwoprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'The full Hessian after the dft contribution'
          call hess_hssout(dbl_mb(k_exy),n3xyz,n3xyz,n3xyz)
          write(LuOut,*)
          write(LuOut,*) 'The full g_rhs after the dft contribution'
          call util_flush(LuOut)
        endif
        call ga_print(g_rhs)
        if(ga_nodeid().eq.0)write(LuOut,*)
      endif
c
c     Coulomb and exact exchange second derivative contributions to
c     hessian. Need to do a lot of setup for this.
c
c     Get some basis set and shell information.
c
      if (oprint.and.(ga_nodeid().eq.0)) 
     .     write(LuOut,*) 'HESS: setup and block order'
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('twodd_cont: could not get nsh',0, BASIS_ERR)
      if (.not. bas_nbf_ce_max(basis,max_at_bf))
     $     call errquit('twodd_cont: could not get max_at_bf',0,
     &       BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('twodd_cont: could not get max_sh_bf',0,
     &       BASIS_ERR)
c
C     define threshold for Schwarz screening(same as in SCF)
c
      if (.not. rtdb_get(rtdb, 'scf:tol2e', mt_dbl, 1, tol2e)) then
         if (rtdb_get(rtdb,'scf:thresh',mt_dbl,1,tol2e)) then
            tol2e = min(1d-7,tol2e * 1d-2)
         else
            tol2e = 1.0d-7
         end if
      end if
c
c     initialize for integrals
c
      call int_init(rtdb, nbasis, basis_list)
      call schwarz_init(geom, basis)
      call int_terminate()
      call intd_init(rtdb, nbasis, basis_list)
      t0=util_wallsec()
      
c
c     For now I am not dealing with the active atom stuff, just setting
c     every atom to active.
c
c     if (.not. ma_push_get(MT_LOG,ncent,'active atoms',l_act,k_act))
c    $     call errquit('twodd_cont: could not allocate l_act',ncent) 
      do i = 1, ncent
        log_mb(k_act+i-1) = .true.
      end do
c
c     Take care of the first derivative contributions to the rhs
c
      if(otwoprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*) 'Original Density'
          call util_flush(LuOut)
        endif
        do i = 1, ndens
          call ga_print(g_dens(i))
        enddo
      endif
C
c     This is temporary until I get the rest of the code working with
c     multi-dimensional GA for the density
c
      dims(1) = ndens*n3xyz
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'CPKS density',chunk,
     &     g_dens2))
     &   call errquit('twodd_cont: could not allocate g_dens2',555,
     &       GA_ERR)
      call ga_zero(g_dens2)
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      do i = 1, ndens
        alo(1) = i
        ahi(1) = i
        call nga_copy_patch('N',g_dens(i),blo,bhi,g_dens2,alo,ahi)
      enddo
c
      nfock = ndens*n3xyz
      if (scftype .eq. 'UHF') call uhf_fock_setup(g_dens2,g_rhs,
     *   jfac,kfac,ndens,nfock,nbf,otwoprint)
      call shell_fock_build (geom,basis, 1, ndens, jfac, kfac,
     *   tol2e, g_dens2, g_rhs,.true.)
      if (scftype .eq. 'UHF') call uhf_fock_finish(g_dens2,g_rhs,
     *   jfac,kfac,ndens,nfock,nbf,otwoprint)
      call intd_terminate()
      if(oprintsome) write(LuOut,12) util_wallsec()-t0
 12     format(/' HESSIAN: 2-el 1st deriv.',
     $    ' term done in ',14x,f9.1,'s'/)
      t0=util_wallsec()
      call intdd_init(rtdb, nbasis, basis_list)
c will eventually bring back - log_mb(k_act))
c
c     Get rid of temporary memory
c
      if (.not. ga_destroy(g_dens2)) call errquit
     *   ('twodd_cont: failed to free g_dens2', ndens, GA_ERR)
      if (otwoprint) then
         if(ga_nodeid().eq.0) then
           write(LuOut,*) 'after going into shell_fock_build:'
           write(LuOut,*) 'geom, basis, n3xyz, ndens, tol2e',geom,basis,
     *                     n3xyz, ndens, tol2e
           call util_flush(LuOut)
         endif
         call ga_print(g_rhs)
         do i = 1, ndens
           call ga_print(g_dens(i))
         enddo
      endif
c
c     Set up some variables that would be set using mp2, etc... for the
c     the call to twodd_coul_ex.
c
      omp2 = .false.
      k_pdm2  = 1            ! To avoid SEGV
      k_pdm2a = 1
      k_pdm2b = 1
      k_pdm2c = 1
      k_coeff = 1
      ld_ij2  = 1
      ld_kl2  = 1
      ld_ik2  = 1
      ld_jl2  = 1
      ld_il2  = 1
      ld_jk2  = 1
      ld_ij3  = 1
      ld_kl3  = 1
      ld_ik3  = 1
      ld_jl3  = 1
      ld_il3  = 1
      ld_jk3  = 1
      ld_ij4  = 1
      ld_kl4  = 1
      ld_ik4  = 1
      ld_jl4  = 1
      ld_il4  = 1
      ld_jk4  = 1
      ld_ik5  = 1    ! Should be 0?
      ld_jl5  = 1
      ld_il5  = 1
      ld_jk5  = 1
      ld_ik6  = 1
      ld_jl6  = 1
      ld_il6  = 1
      ld_jk6  = 1
      ld_ik7  = 1
      ld_jl7  = 1
      ld_il7  = 1
      ld_jk7  = 1
      ld_ik8  = 1
      ld_jl8  = 1
      ld_il8  = 1
      ld_jk8  = 1
c
      blen = min(nbf,maxblen,6*max_sh_bf) ! d(6)*6=36,6**4=1296 quartets
c
c     Need to put the density in the D- and D+ formalism required for
c     2pdm formalism in twodd_coul_ex.  g_wdens is used as temporary
c     space.
c
      if (scftype .eq. 'UHF') then
        call ga_dadd (one, g_dens(1), mone, g_dens(2), g_wdens)
        call ga_dadd (one, g_dens(1), one,  g_dens(2), g_dens(1))
        call ga_copy (g_wdens, g_dens(2))
        call ga_symmetrize(g_dens(1))  ! Need this?
        call ga_symmetrize(g_dens(2))
      endif
      if (.not. ga_destroy(g_wdens))
     *  call errquit('twodd_cont: failed to free g_wdens',555, GA_ERR)
c
c     Set up memory needs
c
      lsqa = blen*blen
 
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ij', lh_ij, ld_ij))
     $     call errquit('twodd_cont:could not allocate ld_ij',lsqa,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_kl', lh_kl, ld_kl))
     $     call errquit('twodd_cont:could not allocate ld_kl',lsqa,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik', lh_ik, ld_ik))
     $     call errquit('twodd_cont:could not allocate ld_ik',lsqa,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl', lh_jl, ld_jl))
     $     call errquit('twodd_cont:could not allocate ld_jl',lsqa,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il', lh_il, ld_il))
     $     call errquit('twodd_cont:could not allocate ld_il',lsqa,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk', lh_jk, ld_jk))
     $     call errquit('twodd_cont:could not allocate ld_jk',lsqa,
     &       MA_ERR)
      if (scftype .eq. 'UHF' .or.
     $     scftype .eq. 'ROHF' .or. omp2) then ! UHF or ROHF or MP2
         if (scftype .eq. 'ROHF' .or. omp2) then
            if (.not. ma_push_get(mt_dbl,lsqa,'ld_ij2',lh_ij2,ld_ij2))
     $           call errquit(
     $            'twodd_cont:could not allocate ld_ij2',lsqa,
     &       MA_ERR)
            if (.not. ma_push_get(mt_dbl,lsqa,'ld_kl2',lh_kl2,ld_kl2))
     $           call errquit('twodd_cont:could not allocate ld_kl2',
     $            lsqa, MA_ERR)
         end if
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik2', lh_ik2, ld_ik2))
     $        call errquit('twodd_cont:could not allocate ld_ik2',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl2', lh_jl2, ld_jl2))
     $        call errquit('twodd_cont:could not allocate ld_jl2',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il2', lh_il2, ld_il2))
     $        call errquit('twodd_cont:could not allocate ld_il2',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk2', lh_jk2, ld_jk2))
     $        call errquit('twodd_cont:could not allocate ld_jk2',lsqa,
     &       MA_ERR)
      end if
c
      if (.not. ma_push_get(mt_dbl, blen**4, 'pdm2d',
     $     l_pdm2d, k_pdm2d)) call errquit
     $     ('twodd_cont:failed allocating pdm2d',blen**4, MA_ERR)
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call intb_mem_2e4c(max2e, mscratch_2e) ! blocking algorithm
      max2e = max2e/12
      max1e = max1e/12
c     max2e = max(max2e,1296*100)          ! 100 D quartets
      lbuf = max(max1e, max2e)
      lscratch = max(mscratch_1e, mscratch_2e)
c
c     probably need the 12 to be 78  CHECK
c     if (.not. ma_push_get(mt_dbl,12*lbuf,'deriv buffer',l_buf,
      if (.not. ma_push_get(mt_dbl,78*lbuf,'deriv buffer',l_buf,
     $    k_buf)) call errquit('twodd_cont:could not allocate buffer'
     $    ,78*lbuf, MA_ERR)
*      call ma_summarize_allocated_blocks()
      if (.not. ma_push_get(mt_dbl,lscratch,'deriv scratch',
     $     l_scr,k_scr))call errquit
     $     ('twodd_cont: scratch alloc failed', lscratch, MA_ERR)
      if (.not. ma_push_get(mt_int,4*lbuf,'labels',l_labels,k_labels))
     $     call errquit('twodd_cont: could not allocate labels',
     $     4*lbuf, MA_ERR)
c
c     Block the shells for Texas
c
      if (.not. ma_push_get(mt_int, nsh, 'shmap', l_shmap, k_shmap))
     $     call errquit('twodd_cont:could not allocate shmap',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     $     call errquit('twodd_cont:could not allocate blo',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     $     call errquit('twodd_cont:could not allocate bhi',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     $     call errquit('twodd_cont:could not allocate bflo',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     $     call errquit('twodd_cont:could not allocate bfhi',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfglo', l_bfglo, k_bfglo))
     $     call errquit('twodd_cont:could not allocate blo',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfghi', l_bfghi, k_bfghi))
     $     call errquit('twodd_cont:could not allocate bhi',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfmap', l_bfmap, k_bfmap))
     $     call errquit('twodd_cont:could not allocate bfmap',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'rbfmap', l_rbfmap, k_rbfmap))
     $     call errquit('twodd_cont:could not allocate rbfmap',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bftoat', l_bftoat, k_bftoat))
     $     call errquit('twodd_cont:could not allocate bftoat',nbf,
     &       MA_ERR)
 
      call grad_shorder(basis, nsh, nbf, maxsh, blen,
     $     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     &     int_mb(k_shmap),
     $     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     $     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))
 
      do i = 1, nbf
         if (.not. bas_bf2ce(basis, i, iat)) call errquit('bf2ce',i,
     &       BASIS_ERR)
         int_mb(k_bftoat+i-1) = iat
      end do
c
c     ==================================================================
c
c     BEWARE:  below here all AO indices have been put in texas order!
c
c     ==================================================================
c
c     Now reorder all of the arrays with AO indices according to the map
c
      do i = 1, ndens
         call ga_reorder(g_dens(i), .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))
         if (otwoprint) then
           if(ga_nodeid().eq.0) then
             write(LuOut,*) 'Density after reordering'
             call util_flush(LuOut)
           endif
           call ga_print(g_dens(i))
         endif
      end do
c
      if (.not. ma_push_get(mt_int, 4*maxq, 'list', l_list, k_list))
     $     call errquit('twodd_cont:could not allocate list',4*maxq,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('twodd_cont:could not allocate q4',maxq,
     &       MA_ERR)
      nproc = ga_nnodes()
c
c  Actually do the deed
c
      if (oprint.and.(ga_nodeid().eq.0)) 
     $     write(LuOut,*) 'HESS: before twodd_coul_ex call'
c
      if (jfac(1).ne.0.0d0 .or. kfac(1).ne.0.0d0) call twodd_coul_ex(
     $     dbl_mb(ld_ij), dbl_mb(ld_kl), dbl_mb(ld_ik),
     $     dbl_mb(ld_jl), dbl_mb(ld_il), dbl_mb(ld_jk),
     $     dbl_mb(ld_ij2),dbl_mb(ld_kl2),dbl_mb(ld_ik2),
     $     dbl_mb(ld_jl2),dbl_mb(ld_il2),dbl_mb(ld_jk2),
     $     dbl_mb(ld_ij3),dbl_mb(ld_kl3),dbl_mb(ld_ik3),
     $     dbl_mb(ld_jl3),dbl_mb(ld_il3),dbl_mb(ld_jk3),
     $     dbl_mb(ld_ij4),dbl_mb(ld_kl4),dbl_mb(ld_ik4),
     $     dbl_mb(ld_jl4),dbl_mb(ld_il4),dbl_mb(ld_jk4),
     $     dbl_mb(ld_ik5),dbl_mb(ld_jl5),dbl_mb(ld_il5),
     $     dbl_mb(ld_jk5),
     $     dbl_mb(ld_ik6),dbl_mb(ld_jl6),dbl_mb(ld_il6),
     $     dbl_mb(ld_jk6),
     $     dbl_mb(ld_ik7),dbl_mb(ld_jl7),dbl_mb(ld_il7),
     $     dbl_mb(ld_jk7),
     $     dbl_mb(ld_ik8),dbl_mb(ld_jl8),dbl_mb(ld_il8),
     $     dbl_mb(ld_jk8),
     $     dbl_mb(k_exy), g_dens, g_force, blen,
     $     geom, basis, nproc, ncent,
     $     lscratch, dbl_mb(k_scr), lbuf, dbl_mb(k_buf),
     $     int_mb(k_labels), maxq, int_mb(k_list), dbl_mb(k_q4),
     $     tol2e, nsh,
     $     log_mb(k_act), skelsym, scftype, omp2, nopen, nbf,
     $     dbl_mb(k_pdm2),dbl_mb(k_pdm2a),dbl_mb(k_pdm2b), ! MCSCF
     $     dbl_mb(k_pdm2c),dbl_mb(k_pdm2d), dbl_mb(k_coeff), ! MCSCF
     $     nshblocks,
     $     int_mb(k_shmap), int_mb(k_shglo), int_mb(k_shghi),
     $     int_mb(k_bfglo), int_mb(k_bfghi),
     $     int_mb(k_bfmap), int_mb(k_rbfmap),
     $     int_mb(k_bftoat), int_mb(k_shbflo), int_mb(k_shbfhi),
     $     jfac, kfac) 
c    $     jfac, kfac,.false.) 
c
c     Print out the hessian before adding anything else
c
      if (otwoprint.and.(ga_nodeid().eq.0)) then
        write(LuOut,*)
        write(LuOut,*) 'The Hessian after the twodd contribution'
        call hess_hssout(dbl_mb(k_exy),n3xyz,n3xyz,n3xyz)
        write(LuOut,*)
      endif
c
c     Get rid of the g_dens arrays since they are no longer needed.
c
      do i = 1, ndens
        if (.not. ga_destroy(g_dens(i))) call errquit
     *     ('twodd_cont: failed to free g_dens()', ndens, GA_ERR)
      end do
c
c     Need to get rid of extra memory
c
      if (.not. ma_chop_stack(lh_ij))
     $  call errquit('twodd_cont: failed chopping MA stack',555, MA_ERR)
c
c     call terminate since cphf will call it's own integral
c     initialization
c
      call schwarz_tidy()    ! need to get rid of temp buffers
      call intdd_terminate()
      if(oprintsome) write(LuOut,13) util_wallsec()-t0
 13   format(/' HESSIAN: 2-el 2nd deriv.',
     $    ' term done in ',14x,f9.1,'s'/)
c
      return
      end
