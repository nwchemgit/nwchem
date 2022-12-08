c
c     Compute the right-hand-sides R of the linear system of equations.
c     The current implementation exploits subroutine tddft_nga_cont
c     to avoid reimplementing the integral processing. This comes at a
c     price though. The algorithm below needs 3 integral evaluations
c     to complete R, whereas the optimal solution would need only
c     1 integral evaluation.
c     HvD 11/2007, NG 11/2012
c
      subroutine tddft_grad_compute_r(rtdb,ihdl_geom,ihdl_bfao,tol2e,
     +                             tda,ipol,nao,nfc,naoc,nocc,nav,nfv,
     +                             nroot,g_mo,g_xpy,g_xmy,g_tp,g_r,
     +                             kfac,lhashf,otriplet,lhascd)
c
c     Eq. 19-22 (Furche & Ahlrichs)
c     The algorithm is
c     1. Do the T contributions:
c       1.1. Transform T MO->AO
c       1.2. Compute (A+B)T in AO basis
c       1.3. Transform (A+B)T AO->MO(ia) and store in R
c     2. Do the (X+Y) and (X-Y) contributions:
c       2.1. Transform (X+Y) and (X-Y) in MO basis -> X and Y in AO basis
c       2.2. Compute (A+B)X and (A-B)X
c       2.3. Compute (A+B)Y and (A-B)Y
c       2.4. Transform (A+B)X, (A-B)X, (A+B)Y and (A-B)Y from AO basis to (A+B)(X+Y) and (A-B)(X-Y) in MO basis
c       2.5. Contract with (X+Y) and (X-Y) to generate contributions to R.
c     3. Do the Gxc contributions:
c       3.1. Compute the electron density matrices in AO basis
c       3.2. Transform (X+Y) to AO basis
c       3.3. Call fock_xc to evaluate gxc in AO basis
c       3.4. Transform gxc from AO basis to MO(ia) and add onto R.
c
c     Note:  H+[X+Y]ij, H+[X+Y]ab are symmetric,  H-[X-Y]ij, H-[X-Y]ab are anti-symmetric.
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tddft_grad_util.fh"
#include "stdio.fh"
c
c     Input:
c
      integer rtdb      ! the runtime database handle
      integer ihdl_geom ! the geometry handle
      integer ihdl_bfao ! the AO basis set handle
      logical tda       ! True if Tamm-Dancoff approximation
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer nao       ! the number of AOs
      integer nfc(2)    ! the number of frozen cores
      integer naoc(2)   ! the number of active occupied orbitals
      integer nocc(2)   ! the number of occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nfv(2)    ! the number of frozen virtuals
      integer nroot     ! the number of roots to be treated
      integer g_mo(2)   ! global arrays to the MOs
      integer g_xpy(2)  ! global arrays to hold X+Y
      integer g_xmy(2)  ! global arrays to hold X-Y
      integer g_tp(2)   ! global arrays to hold T, and later P
c
      logical lhashf    ! =.true.  hybrid functionals
                        ! =.false. otherwise
      logical otriplet  ! =.true.  triplet excited states
                        ! =.false. singlet excited states
                        ! value does not matter for TDUDFT
      logical lhascd
c
      double precision kfac  ! the weight of the Hartree-Fock exchange
                             ! contributions
      double precision tol2e ! 2-electron integral tolerance
c
c     Output:
c
      integer g_r(2)    ! global arrays holding the right-hand-sides
c
c     Local for now...
c
      logical oskel 
      parameter (oskel=.false.)
c
c     Functions:
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  xc_gotxc
      external xc_gotxc
c
c     Local
c
      integer g_tuv         ! work space array
      integer g_apbt        ! work space array for (A+B)T
      integer g_ambt        ! work space array for (A-B)T
      integer g_x           ! work space array for X (AO basis)
      integer g_y           ! work space array for Y (AO basis)
      integer g_apbx        ! work space array for (A+B)X
      integer g_ambx        ! work space array for (A-B)X
      integer g_apby        ! work space array for (A+B)Y
      integer g_amby        ! work space array for (A-B)Y
      integer g_hij(2)      ! work space arrays
      integer g_hab(2)      ! work space arrays
      integer idim(3)       ! dimensions
      integer ichnk(3)      ! chunking
      integer alo(3),ahi(3) ! lower and upper limits for GA patches
      integer blo(3),bhi(3) ! lower and upper limits for GA patches
      integer clo(3),chi(3) ! lower and upper limits for GA patches
      integer ir,ip         ! counters for roots and spins
      integer calc_type     ! the calculation type for fock_xc
      integer l_dens,k_dens ! memory for array of density handles
      integer l_den2,k_den2 ! memory for array of density handles
      integer ndens         ! the number of density matrices
      integer l_gxc ,k_gxc  ! memory for array of Gxc matrix handles
      integer ngxc          ! the number of Gxc matrices
      integer i             ! counter for abitrary purposes
      logical oroot
      integer iwhich
c
      character*32 pname
c
c Daniel (2-8-13): This is a local TDA variable, so that we don't need
c to change the global tda variable.
      logical tdaloc, tddft_grad_tdaloc
      external tddft_grad_tdaloc
      logical oprint
      integer g_r_x(2),g_r_org(2)
c
      double precision Exc(2)    ! Exchange-correlation energy
c
      pname="tddft_grad_compute_r: "
      iwhich = 0  ! call to tddft_nga_cont()
      write(6,*) ' 0 tda ',tda
c
c     1. Do the T contributions
c
c     1.1. Create global array to hold Tuv
c
      oprint=util_print('tddft_grad_r',print_debug)
      idim(1) = nroot*ipol
      idim(2) = nao
      idim(3) = nao
      ichnk(1) = nroot*ipol
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not.nga_create(mt_dbl,3,idim,'vectors Tuv',ichnk,g_tuv))
     +  call errquit(pname//'failed to create g_tuv',0,GA_ERR)
c
c     1.2. Tuv = sum_ij Cui*Tij*Cvj
c
      call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ij",g_mo,g_tp,"pq",g_tuv)
c
c     1.3. Tuv = Tuv + sum_ab Cua*Tab*Cvb
c
      call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,1.0d0,"ab",g_mo,g_tp,"pq",g_tuv)
c
c     1.4. Create global array for (A+B)T in AO basis
c
      if (.not.nga_create(mt_dbl,3,idim,'vec (A+B)T',ichnk,g_apbt))
     +  call errquit(pname//'failed to create g_apbt',0,GA_ERR)
      call ga_zero(g_apbt)
c
c     1.5. Compute (A+B)T in AO basis (currently we have to compute
c          (A-B)T as well although we do not need it...)
c Daniel (12/1/12): Technically we are calculating H^+[T] (and H^-[T]),
c here (e.g. what some people call the coupling matrices contracted with
c a density matrix).  (A+B) actually contains the diagonal part of the
c orbital rotation Hessian (MO energy part), which is absent here.
c Daniel (12-17-12): Note that H^+[T] =/= 0 while H^-[T] = 0, since
c T is a symmetric matrix.  In tddft_contract routines, the density
c matrix is first symmetrized then anti-symmetrized (for HF exchange),
c however the latter is zero for a symmetric matrix.
c
      if (.not.nga_create(mt_dbl,3,idim,'vec (A-B)T',ichnk,g_ambt))
     +    call errquit(pname//'failed to create g_ambt',0,GA_ERR)
      call ga_zero(g_ambt)
c
c Daniel (2-26-13): It was not obvious that we need to unset
c fock_xc:triplet here for restricted triplet calculations to work.
      if (otriplet) then
        if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1    call errquit(pname//'failed to set triplet',0,RTDB_ERR)
      endif
c
      call tddft_nga_cont(rtdb,ihdl_geom,ihdl_bfao,g_tuv,g_apbt,g_ambt,
     +nao,ipol,tol2e,tda,oskel,kfac,lhashf,.false.,nroot,iwhich,lhascd)
      write(6,*) ' 2 tda ',tda
c Daniel (2-26-13): Reset fock_xc:triplet here for restricted triplet 
c calculations to work.
      if (otriplet) then
        if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.true.))
     1    call errquit(pname//'failed to set triplet',0,RTDB_ERR)
      endif
c
      if (.not.ga_destroy(g_ambt)) 
     +    call errquit(pname//'failed to destroy g_ambt',0,GA_ERR)
c
c     1.6. Destroy global array for Tuv
c
      if (.not.ga_destroy(g_tuv))
     +    call errquit(pname//'failed to destroy g_tuv',0,GA_ERR)
c
c     1.7. Transform (A+B)T from AO basis to Ria in MO basis
c
c Daniel (1-5-13): Here, the R vector is half of what it should be for
c a CIS calculation.  When we do RPA with B = 0, the value here is
c double what comes from the CIS code, which is a consequence of the
c doubling of g_apbt (and g_ambt) for RPA to account for the B
c contribution in tddft_nga_cont.
      if (tda) then
        call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +       nroot,2.0d0,0.0d0,"ib",g_mo,g_apbt,g_r,"ib")
      else
        call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +       nroot,1.0d0,0.0d0,"ib",g_mo,g_apbt,g_r,"ib")
      endif
c
c     1.8. Destroy global array for (A+B)T in AO basis
c
      if (.not.ga_destroy(g_apbt)) 
     +    call errquit(pname//'failed to destroy g_apbt',0,GA_ERR)
c
c     1. Done T contributions
cDEBUG
      if (oprint) then
       oroot = ga_nodeid().eq.0
       if (oroot) write(LuOut,*)'DEBUG A: tddft_grad_comp_r: R_T'
       call tddft_grad_print_array(ipol,nroot,g_r,dble(ipol))
      endif
cDEBUG
c
c     2. Do the (X+Y) and (X-Y) contributions
c
c     2.1 Create global arrays for X and Y in AO basis
c
c Daniel (12-4-12): nga_create constructs an initialized global array,
c so we don't need to zero it ourselves.
      if (.not.nga_create(mt_dbl,3,idim,'vectors Xuv',ichnk,g_x))
     +    call errquit(pname//'failed to create g_x',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors Yuv',ichnk,g_y))
     +    call errquit(pname//'failed to create g_y',0,GA_ERR)
c
c     2.2 Transform (X+Y) to AO basis in g_x and (X-Y) to AO basis
c         in g_y
c
      call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ib",g_mo,g_xpy,"ib",g_x)
      if (.not.tda) then
        call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ib",g_mo,g_xmy,"ib",g_y)
c Daniel (1-5-13): In order to ensure Y = 0 for later parts of this
c code, we put g_xpy into g_y such that the lines in 2.3 work the same
c way they do for RPA.  This is needed so that CIS and RPA (B = 0) are
c identical.
      else
        call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ib",g_mo,g_xpy,"ib",g_y)
      endif
c
c     2.3 Compute X and Y from (X+Y) and (X-Y) in place
c
      alo(1) = 1
      alo(2) = 1
      alo(3) = 1
      ahi(1) = nroot*ipol
      ahi(2) = nao
      ahi(3) = nao
c Daniel (1-5-13): With the modification above, the lines here behave
c like you'd expect (i.e. g_x = X and g_y = Y = 0 for CIS, rather than
c what was here before which made g_x = 0.50*X and g_y = 0.50*X).
      call nga_add_patch(0.5d0,g_x,alo,ahi,0.5d0,g_y,alo,ahi,g_x,
     +                   alo,ahi)
      call nga_add_patch(1.0d0,g_x,alo,ahi,-1.0d0,g_y,alo,ahi,g_y,
     +                   alo,ahi)
c
c     2.4 Allocate various workspace arrays
c
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A+B)X',ichnk,g_apbx))
     +    call errquit(pname//'failed to create g_apbx',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A-B)X',ichnk,g_ambx))
     +    call errquit(pname//'failed to create g_ambx',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A+B)Y',ichnk,g_apby))
     +    call errquit(pname//'failed to create g_apby',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A-B)Y',ichnk,g_amby))
     +    call errquit(pname//'failed to create g_amby',0,GA_ERR)
      tdaloc=tddft_grad_tdaloc(rtdb,tda,lhashf)
c     
c     2.5 Compute (A+B)X, (A-B)X, (A+B)Y and (A-B)Y
c
      call tddft_nga_cont(rtdb,ihdl_geom,ihdl_bfao,g_x,g_apbx,g_ambx,
     +nao,ipol,tol2e,tdaloc,oskel,kfac,lhashf,otriplet,nroot,iwhich,
     +lhascd)
      call tddft_nga_cont(rtdb,ihdl_geom,ihdl_bfao,g_y,g_apby,g_amby,
     +    nao,ipol,tol2e,tdaloc,oskel,kfac,lhashf,otriplet,nroot,iwhich,
     +    lhascd)
c
c     2.6 Dispose of X and Y
c
      if (.not.ga_destroy(g_x)) 
     +  call errquit(pname//'failed to destroy g_x',0,GA_ERR)
      if (.not.ga_destroy(g_y)) 
     +  call errquit(pname//'failed to destroy g_y',0,GA_ERR)
c
c     2.7 Compute (A+B)(X+Y) and (A-B)(X-Y)
c
      call nga_add_patch(1.0d0,g_apbx,alo,ahi,1.0d0,g_apby,alo,ahi,
     +                   g_apbx,alo,ahi)
      call nga_add_patch(1.0d0,g_ambx,alo,ahi,-1.0d0,g_amby,alo,ahi,
     +                   g_ambx,alo,ahi)
c
c     2.8 Dispose of (A+B)Y and (A-B)Y
c
      if (.not.ga_destroy(g_apby)) 
     +  call errquit(pname//'failed to destroy g_apby',0,GA_ERR)
      if (.not.ga_destroy(g_amby)) 
     +  call errquit(pname//'failed to destroy g_amby',0,GA_ERR)
c
c     2.9 Allocate (A+-B)(X+-Y)ij
c
      do ip = 1, ipol
        idim(1) = nroot
        idim(2) = naoc(ip)
        idim(3) = naoc(ip)
        ichnk(1) = nroot
        ichnk(2) = -1
        ichnk(3) = -1
        if (.not.nga_create(mt_dbl,3,idim,'vec hij',ichnk,g_hij(ip)))
     +    call errquit(pname//'failed to create g_hij',0,GA_ERR)
      enddo
c
c     2.10 Transform (A+B)(X+Y) to MO basis occupied-occupied block only
c
      call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ij",g_mo,g_apbx,g_hij,"ij")
      if(oprint)then
c     extra ga to keep track of X
         if (.not. ga_duplicate(g_r(1), g_r_x(1), 'gtmp2'))
     c        call errquit('tdgcr: gadupl failed',0,0)
         call ga_zero(g_r_x(1))
         g_r_org(1)=g_r(1)
         g_r(1)=g_r_x(1)
         if(ipol.eq.2) then
         if (.not. ga_duplicate(g_r(2), g_r_x(2), 'gtmp2'))
     c           call errquit('tdgcr: gadupl failed',0,0)
         call ga_zero(g_r_x(2))
         g_r_org(2)=g_r(2)
         g_r(2)=g_r_x(2)
         endif
      endif
c
c     2.11 Add -sum_j (X+Y)ja [(A+B)(X+Y)ji] to Ria
c
      do ip = 1, ipol
        do ir=1,nroot
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = naoc(ip)
          alo(3) = 1
          ahi(3) = naoc(ip)
          blo(1) = ir
          bhi(1) = ir
          blo(2) = 1
          bhi(2) = naoc(ip)
          blo(3) = 1
          bhi(3) = nav(ip)
          clo(1) = ir
          chi(1) = ir
          clo(2) = 1
          chi(2) = naoc(ip)
          clo(3) = 1
          chi(3) = nav(ip)
          call nga_matmul_patch('n','n',-1.0d0,1.0d0,g_hij(ip),alo,ahi,
     +         g_xpy(ip),blo,bhi,g_r(ip),clo,chi)
        enddo
      enddo
c
c Daniel (1-7-13): We can definitely avoid this part of the routine if
c we don't use exact exchange, since the linear transformation H^-[V]
c is zero in that case. 
      if (lhashf) then
c
c     2.12 Transform (A-B)(X-Y) to MO basis occupied-occupied block only
c
        call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +       nroot,1.0d0,0.0d0,"ij",g_mo,g_ambx,g_hij,"ij")
c
c     2.13 Add -sum_j (X-Y)ja [(A-B)(X-Y)ji] to Ria
c
        do ip = 1, ipol
          do ir=1,nroot
            alo(1) = ir
            ahi(1) = ir
            alo(2) = 1
            ahi(2) = naoc(ip)
            alo(3) = 1
            ahi(3) = naoc(ip)
            blo(1) = ir
            bhi(1) = ir
            blo(2) = 1
            bhi(2) = naoc(ip)
            blo(3) = 1
            bhi(3) = nav(ip)
            clo(1) = ir
            chi(1) = ir
            clo(2) = 1
            chi(2) = naoc(ip)
            clo(3) = 1
            chi(3) = nav(ip)
c Daniel (1-5-13): Manipulate the code here for CIS to use g_xpy here
c since (X+Y) = (X-Y) = X.  This is a consequence of not allocating 
c a g_xmy array for CIS.  The linear transformation H^-[X] still exists
c in CIS.
            if (.not.tda) then
              call nga_matmul_patch('n','n',+1.0d0,1.0d0,
     +             g_hij(ip),alo,ahi,
     +             g_xmy(ip),blo,bhi,g_r(ip),clo,chi)
            else
              call nga_matmul_patch('n','n',+1.0d0,1.0d0,
     +             g_hij(ip),alo,ahi,
     +             g_xpy(ip),blo,bhi,g_r(ip),clo,chi)
            endif
          enddo
        enddo
c
      endif ! lhashf 
c
c     2.14 Dispose of g_hij
c
      do ip = 1, ipol
        if (.not.ga_destroy(g_hij(ip))) 
     +   call errquit(pname//'failed to destroy g_hij',0,GA_ERR)
      enddo
c
c     2.15 Allocate (A+-B)(X+-Y)ab
c
      do ip = 1, ipol
        idim(1) = nroot
        idim(2) = nav(ip)
        idim(3) = nav(ip)
        ichnk(1) = nroot
        ichnk(2) = -1
        ichnk(3) = -1
        if (.not.nga_create(mt_dbl,3,idim,'vec hab',ichnk,g_hab(ip)))
     +   call errquit(pname//'failed to create g_hab',0,GA_ERR)
      enddo
c
c     2.16 Transform (A+B)(X+Y) to MO basis virtual-virtual block only
c
      call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ab",g_mo,g_apbx,g_hab,"ab")
c
c     2.17 Add +sum_b (X+Y)ib [(A+B)(X+Y)ab] to Ria
c
      do ip = 1, ipol
        do ir=1,nroot
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = naoc(ip)
          alo(3) = 1
          ahi(3) = nav(ip)
          blo(1) = ir
          bhi(1) = ir
          blo(2) = 1
          bhi(2) = nav(ip)
          blo(3) = 1
          bhi(3) = nav(ip)
          clo(1) = ir
          chi(1) = ir
          clo(2) = 1
          chi(2) = naoc(ip)
          clo(3) = 1
          chi(3) = nav(ip)
          call nga_matmul_patch('n','n', 1.0d0,1.0d0,g_xpy(ip),alo,ahi,
     +         g_hab(ip),blo,bhi,g_r(ip),clo,chi)
        enddo
      enddo
c
c Daniel (1-7-13): We can definitely avoid this part of the routine if
c we don't use exact exchange, since the linear transformation H^-[V]
c is zero in that case. 
      if (lhashf) then
c
c     2.18 Transform (A-B)(X-Y) to MO basis virtual-virtual block only
c
        call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +       nroot,1.0d0,0.0d0,"ab",g_mo,g_ambx,g_hab,"ab")
c
c     2.19 Add +sum_b (X-Y)ib [(A-B)(X-Y)ab] to Ria
c
        do ip = 1, ipol
          do ir=1,nroot
            alo(1) = ir
            ahi(1) = ir
            alo(2) = 1
            ahi(2) = naoc(ip)
            alo(3) = 1
            ahi(3) = nav(ip)
            blo(1) = ir
            bhi(1) = ir
            blo(2) = 1
            bhi(2) = nav(ip)
            blo(3) = 1
            bhi(3) = nav(ip)
            clo(1) = ir
            chi(1) = ir
            clo(2) = 1
            chi(2) = naoc(ip)
            clo(3) = 1
            chi(3) = nav(ip)
c Daniel (1-5-13): Manipulate the code here for CIS to use g_xpy here
c since (X+Y) = (X-Y) = X.  This is a consequence of not allocating 
c a g_xmy array for CIS.  The linear transformation H^-[X] still exists
c in CIS.
            if (.not.tda) then
              call nga_matmul_patch('n','n',-1.0d0,1.0d0,
     +             g_xmy(ip),alo,ahi,
     +             g_hab(ip),blo,bhi,g_r(ip),clo,chi)
            else
              call nga_matmul_patch('n','n',-1.0d0,1.0d0,
     +             g_xpy(ip),alo,ahi,
     +             g_hab(ip),blo,bhi,g_r(ip),clo,chi)
            endif
          enddo
        enddo
c
      endif ! lhashf
c
c     2.20 Dispose of g_hab
c
      do ip = 1, ipol
        if (.not.ga_destroy(g_hab(ip))) 
     +    call errquit(pname//'failed to destroy g_hab',0, GA_ERR)
      enddo
c
c     2.21 Dispose of (A+B)(X+Y) and (A+B)(X-Y)
c
      if (.not.ga_destroy(g_apbx)) 
     +   call errquit(pname//'failed to destroy g_apbx',0, GA_ERR)
      if (.not.ga_destroy(g_ambx)) 
     +   call errquit(pname//'failed to destroy g_ambx',0, GA_ERR)
c
c     2. Done (X+Y) and (X-Y) contributions
c
cDEBUG
      if (oprint) then
       oroot = ga_nodeid().eq.0
       if (oroot) write(LuOut,*)'DEBUG B: tddft_grad_comp_r: R_x'
       call tddft_grad_print_array(ipol,nroot,g_r_x,dble(ipol))
       call ga_add(1d0,g_r_x(1),1d0,g_r_org(1),g_r_org(1))
       g_r(1)=g_r_org(1)
       if (.not.ga_destroy(g_r_x(1))) 
     +     call errquit(pname//'failed to destroy grx',1,GA_ERR)
       if(ipol.eq.2) then
       call ga_add(1d0,g_r_x(2),1d0,g_r_org(2),g_r_org(2))
       g_r(2)=g_r_org(2)
       if (.not.ga_destroy(g_r_x(2))) 
     +     call errquit(pname//'failed to destroy grx',2,GA_ERR)
       endif
      endif
cDEBUG
c
c     3. Do the Gxc contributions
c
c     The DFT code in NWChem expects an array of global array handles
c     for the densities and the Gxc matrices involved.
c     Also we need to set fock_xc:calc_type to trigger the correct
c     matrices to the be computed. The values for calc_type are defined
c     in grid_quadv0.F .
c
c     The densities as expected by NWChem come in the following order
c
c     1. Nroot alpha (X+Y)_ao matrices
c     2. Nroot beta  (X+Y)_ao matrices
c     3. Alpha ground state density matrix
c     4. Beta  ground state density matrix
c
c     By the time we get to tddft_grad_gxc these densities have been
c     reodered (see xc_rhogen). We then have them in the order
c
c     A. For each root
c       1. Alpha (X+Y) quantities
c       2. Beta  (X+Y) quantities
c
      if (xc_gotxc()) then ! check if DFT xc is specified
c Daniel (1-10-13): There are nroot (X+Y) density matrices and 1 GS
c density matrix per spin.
        ndens = ipol*(nroot+1)
        ngxc  = ipol*nroot
        if (.not.ma_push_get(mt_int,ngxc,'gxc-s',l_gxc,k_gxc))
     +    call errquit(pname//'failed to allocate l_gxc',0,MA_ERR)
        if (.not.ma_push_get(mt_int,ndens,'densities',l_dens,k_dens))
     +    call errquit(pname//'failed to allocate l_dens',0,MA_ERR)
        if (.not.ma_push_get(mt_int,ndens,'dens-tmp',l_den2,k_den2))
     +    call errquit(pname//'failed to allocate l_den2',0,MA_ERR)
c
c       3.1. Create and calculate the AO basis density matrices
c
c Daniel (12/1/12): This is the ground-state density matrix.
        do ip = 0, ipol-1
          int_mb(k_dens+nroot*ipol+ip) = 
     +       ga_create_atom_blocked(ihdl_geom,ihdl_bfao,"d_ao")
        enddo
        call tddft_grad_compute_dao(ipol,nao,nocc,g_mo,
     +                           int_mb(k_dens+nroot*ipol))
        if (ipol.eq.1) call ga_scale(int_mb(k_dens+nroot*ipol),2.0d0)
c
c       3.2. Create and calculate the AO basis representation of (X+Y)
c
        do ip = 0, ipol-1
          do ir = 0, nroot-1
            int_mb(k_dens+ip*nroot+ir) = ga_create_atom_blocked(
     +                                   ihdl_geom,ihdl_bfao,"xpy_ao")
            call ga_zero(int_mb(k_dens+ip*nroot+ir))
          enddo
        enddo
        do ip = 1, ipol
          do ir = 1, nroot
            call tddft_grad_trans_mo2ao(1,nao,nfc(ip),naoc(ip),nocc(ip),
     +           nav(ip),nfv(ip),ir,1.0d0,0.0d0,"ib",g_mo(ip),g_xpy(ip),
     +           "ib",int_mb(k_dens+(ip-1)*nroot+ir-1))
            call ga_symmetrize(int_mb(k_dens+(ip-1)*nroot+ir-1))
c Daniel (2-16-13): This line is needed to get matching results from the
c unrestricted code compared to the restricted one.
           if (ipol.eq.1) then
             call ga_scale(int_mb(k_dens+(ip-1)*nroot+ir-1),2.0d0)
           endif
          enddo
        enddo
c
c       3.3. Create and calculate Gxc in AO basis using fock_xc
c
        do i = 0, ngxc-1
          int_mb(k_gxc+i) = ga_create_atom_blocked(ihdl_geom,ihdl_bfao,
     +                                             "gxc_ao")
          call ga_zero(int_mb(k_gxc+i))
        enddo
        if (.not.rtdb_get(rtdb,'fock_xc:calc_type',mt_int,1,calc_type))
     +    calc_type=0
c
c       Set calc_type = 5 for third derivatives
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,5))
     +    call errquit(pname//'failed to set calc_type 5',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type_tddft_r',mt_int,1,
     +      calc_type)) 
     +   call errquit(pname//'failed to set calc_type_tddft_r',0,
     +      RTDB_ERR)
c Daniel (2-16-13): There is a line in fock_xc, involving lcgmin, that
c will set l3d to true, even though we want it to be false for the 
c gradients.  We need to set the RTDB such that dft:cgmin is true so
c that the l3d is set correctly in fock_xc.
c **********************************************************************
c WE NEED TO FIX THIS BECAUSE TDDFT GRADIENTS DON'T WORK with CGMIN SET.
c **********************************************************************
c Daniel (2-18-13): This fix will not work for optimizations
c        if (.not.rtdb_put(rtdb,'dft:cgmin',mt_log,1,.true.))
c     1    call errquit(pname//'failed to set cgmin',0,RTDB_ERR)
c
c       Piss about with the order of the density matrices. 
c       We currently have the density grouped per state and then
c       for each state we have alpha and beta components. Subroutine
c       xc_rhogen with what.eq.1 however expects all alpha components
c       first and then all beta components. The densities expressed
c       on the grid however will be grouped per state once more (sigh).
c       This shit only matters for the unrestricted case.
c
        do i = 0, nroot
          do ip = 0, ipol-1
            int_mb(k_den2+i+ip*(nroot+1)) = int_mb(k_dens+ip+ipol*i)
          enddo
        enddo
        call ga_sync()
cDEBUG
        if (oprint) then
           oroot = ga_nodeid().eq.0
           if (oroot) write(LuOut,*)'DEBUG C: tddft_grad_comp_r: den'
       call tddft_grad_print_array(ipol,nroot,int_mb(k_den2),
     +          dble(ipol))
       call tddft_grad_print_array(ipol,nroot,int_mb(k_den2+1),
     +      dble(ipol))
      endif
cDEBUG
c Daniel (1-10-13): For a single restricted root ngxc = 1 and there are 2
c density matrices fed to fock_xc, the GS density matrix D and the
c (X+Y) density matrix.  The former is need to construct the 
c functional derivatives in the adiabatic approximation, while the
c latter is needed for the contraction with the third-order functional
c derivative. 
c Daniel (2-12-13): Need to feed a triplet flag through the RTDB for
c fock_xc.
        if (otriplet) then
          if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,otriplet))
     1      call errquit(pname//'failed to set triplet',0,RTDB_ERR)
        endif
        call fock_xc(ihdl_geom, nao, ihdl_bfao, ngxc, int_mb(k_den2),
     +               int_mb(k_gxc), Exc, ipol, .false.)
        call ga_sync()
c DEBUG
      if (oprint) then
        oroot = ga_nodeid().eq.0
        if (oroot) write(LuOut,*)'DEBUG: '//pname//'gxc'
        call tddft_grad_print_array(ipol,nroot,int_mb(k_gxc),
     +                              dble(ipol))
      endif
c      call errquit(pname//'quit after fock_xc for Gxc', 0, UNKNOWN_ERR)
c DEBUG
        if (.not.rtdb_get(rtdb,'fock_xc:calc_type_tddft_r',mt_int,1,
     +    calc_type)) 
     +   call errquit(pname//'failed to get calc_type_tddft_r',0,
     +    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,calc_type))
     +   call errquit(pname//'failed to reset calc_type',0,RTDB_ERR)
        if (.not.rtdb_delete(rtdb,'fock_xc:calc_type_tddft_r'))
     +   call errquit(pname//'failed to delete calc_type_tddft_r',
     +    0, RTDB_ERR)
c
c       3.4. Transform the Gxc matrices to MO basis and add them to Ria
c
        do i = 0, ndens-1
          if (.not.ga_destroy(int_mb(k_dens+i))) 
     +     call errquit(pname//'failed to destroy densities',0,GA_ERR)
        enddo
        if (.not.ma_pop_stack(l_den2)) 
     +    call errquit(pname//'failed to deallocate l_den2', 0,MA_ERR)
        if (.not.ma_pop_stack(l_dens)) 
     +    call errquit(pname//'failed to deallocate l_dens', 0, MA_ERR)
        do ip = 1, ipol
          do ir = 1, nroot
            call tddft_grad_trans_ao2mo(1,nao,nfc(ip),naoc(ip),nocc(ip),
     +           nav(ip),nfv(ip),ir,2.0d0,1.0d0,"ib",g_mo(ip),
     +           int_mb(k_gxc+(ip-1)*nroot+ir-1),g_r(ip),"ib")
          enddo
        enddo
        do i = 0, ngxc-1
          if (.not.ga_destroy(int_mb(k_gxc+i))) 
     +      call errquit(pname//'failed to destroy gxc_ao', 0, GA_ERR)
        enddo
        if (.not.ma_pop_stack(l_gxc)) 
     +      call errquit(pname//'failed to deallocate l_gxc', 0, MA_ERR)
c
      endif  ! xc_gotxc()
c
c     Debug print
      if (oprint) then
        oroot = ga_nodeid().eq.0
        if (oroot) write(LuOut,*)'DEBUG: '//pname//'R_all'
        call tddft_grad_print_array(ipol,nroot,g_r,dble(ipol))
      endif
c
      end
c $Id$
