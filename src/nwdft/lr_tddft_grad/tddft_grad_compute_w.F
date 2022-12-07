c
c     Evaluate the quantity W. Eq. 24, Furche & Ahlrichs
c     There are three parts that need evaluating separately. 
c     Occupied-occupied, Occupied-Virtual, Virtual-Virtual parts.
c     HvD 11/2007, NG 11/2012
c
      subroutine tddft_grad_compute_w(rtdb,ihdl_geom,ihdl_bfao,tol2e,
     +           tda,ipol,nroots,nfc,naoc,nocc,nav,
     +           nfv,nao,g_mo,g_p,g_z,g_xpy,g_xmy,eps,omega,g_w,
     +           kfac,lhashf,otriplet,lhascd)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tddft_grad_util.fh"
#include "stdio.fh"
c
c     Input:
c
      integer rtdb      ! runtime database handle
      integer ihdl_geom ! geometry handle
      integer ihdl_bfao ! basis set handle
      logical tda       ! .true. if Tamm-Dancoff approximation is used
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer nfc(2)    ! the number of frozen core orbitals
      integer naoc(2)   ! the number of active occupied orbitals
      integer nocc(2)   ! the total number of occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nfv(2)    ! the number of frozen virtual orbitals
      integer nao       ! the number of basis functions
      integer nroots    ! the number of states to consider
      integer g_mo(2)   ! global array handle for MO coefficients
      integer g_p(2)    ! global array handle for Ppq
      integer g_z(2)    ! global array handle for Zia
      integer g_xpy(2)  ! global array handle for (X+Y)
      integer g_xmy(2)  ! global array handle for (X-Y)
c
      logical lhashf    ! =.true.  hybrid functionals
                        ! =.false. otherwise
      logical otriplet  ! =.true.  triplet excited states
                        ! =.false. singlet excited states
                        ! value does not matter for TDUDFT
      logical lhascd
c
      double precision kfac          ! the weight of the Hartree-Fock
                                     ! exchange contributions
      double precision eps(nao,2)    ! orbital eigenvalues
      double precision omega(nroots) ! the excitation energies
      double precision tol2e         ! integral tolerance
c
c     Output:
c
      integer g_w(2)    ! global array handle for W
c
c     Functions:
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  xc_gotxc
      external xc_gotxc
c
c     Local:
c
      integer g_puv     ! global array handle for Puv
      integer g_apbp ! global array handle for (A+B)Puv
      integer g_ambp ! global array handle for (A-B)Puv
      integer g_x    ! global array handle for X
      integer g_y    ! global array handle for Y
      integer g_apbx ! global array handle for (A+B)X
      integer g_ambx ! global array handle for (A-B)X
      integer g_apby ! global array handle for (A+B)Y
      integer g_amby ! global array handle for (A-B)Y
      integer g_hij(2)  ! global array handle for Hij+-
      integer g_t       ! global array handle for transposer
c
      integer alo(3) ! lower chunck limits on A
      integer ahi(3) ! upper chunck limits on A
      integer blo(3) ! lower chunck limits on B
      integer bhi(3) ! upper chunck limits on B
      integer clo(3) ! lower chunck limits on C
      integer chi(3) ! upper chunck limits on C
      integer ld(3)  ! leading dimensions
c
      integer idim(3)  ! dimensions for global arrays
      integer ichnk(3) ! chunk sizes for distribution
c
      integer nproc  ! the number of processors
      integer iproc  ! the rank of the current processor
      integer jproc  ! the rank of some other processor
c
      integer ip     ! counter over spin components
      integer ir     ! counter over roots
      integer ic     ! counter for arbitrary purposes
      integer icount ! counter for arbitrary purposes
c
      integer mini   ! minimum value of i-label
      integer minj   ! minimum value of j-label
      integer mina   ! minimum value of a-label
      integer minb   ! minimum value of b-label
      integer maxi   ! maximum value of i-label
      integer maxj   ! maximum value of j-label
      integer maxa   ! maximum value of a-label
      integer maxb   ! maximum value of b-label
      integer numi   ! the number of i-labels
      integer numj   ! the number of j-labels
      integer numa   ! the number of a-labels
      integer numb   ! the number of b-labels
c
      integer i      ! counter over i-labels
      integer j      ! counter over j-labels
      integer a      ! counter over a-labels
      integer b      ! counter over b-labels
c
      integer imo    ! molecular orbital label for accessing eps
c
      integer ihdl_a ! the memory handle for block A
      integer ihdl_b ! the memory handle for block B
      integer ihdl_r ! the memory handle for block R
      integer ihdl_v ! the memory handle for block V
c
      integer iptr_a ! the memory index for block A
      integer iptr_b ! the memory index for block B
      integer iptr_r ! the memory index for block R
      integer iptr_v ! the memory index for block V
c
      integer calc_type      ! calculation type for fock_xc
      integer l_dens, k_dens ! memory for array of density handles
      integer l_den2, k_den2 ! memory for array of density handles
      integer ndens          ! the number of density matrices
      integer l_gxc , k_gxc  ! memory for array of Gxc matrix handles
      integer ngxc           ! the number of Gxc matrices
      logical oroot
      integer iwhich
c
      character*32 pname
c
      logical oskel
      parameter (oskel=.false.)
      double precision Exc(2)    ! Exchange-correlation energy
c
      logical tdaloc, tddft_grad_tdaloc
      external tddft_grad_tdaloc
      logical doitw,doitz
      logical doitxpy1,doitxpy2
      logical doitxmy1,doitxmy2
c
      pname = "tddft_grad_compute_w: "
      iwhich = 0   ! call to tddft_nga_cont()
c
c     General stuff
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      do ip = 1, ipol
        call ga_zero(g_w(ip))
      enddo
c Daniel (12/1/12): Here, we have that the occupied-occupied block is
c the sum of the ground state energy-weighted matrix with the
c energy-weighted difference density matrix.  The former gives the
c first term in the expression.  The ground state density matrix makes
c no other contributions to W.
c
c     Compute the occupied-occupied block
c
c     Wijs = delta_ij eps_i
c          + Omega Sum_a [(X+Y)ias(X-Y)jas+(X-Y)ias(X+Y)jas]/2
c          - Sum_a eps_a [(X+Y)ias(X+Y)jas+(X-Y)ias(X-Y)jas]/2
c          + Aijs'[P]
c          + Sum_kat Sum_ldr Gxc_ijskatldr(X+Y)kat(X+Y)ldr
c
c     Do the first term, the parallelisation is driven off the 
c     distribution of g_w.
c
      do ip = 1, ipol
        call nga_distribution(g_w(ip),iproc,alo,ahi)
        doitw = .not.(alo(1).ne.1.or.ahi(1).ne.nroots)
        if (doitw) then
         mini = max(alo(2),alo(3))
         maxi = min(ahi(2),ahi(3),naoc(ip))
         do i = mini, maxi
          imo = nfc(ip)+i
          do ir = 1, nroots
            blo(1) = ir
            bhi(1) = ir
            blo(2) = i
            bhi(2) = i
            blo(3) = i
            bhi(3) = i
            ld(1)  = 1
            ld(2)  = 1
            call nga_acc(g_w(ip),blo,bhi,(3-ipol)*eps(imo,ip),ld,1.0d0)
          enddo
         enddo
        endif ! doitw
      enddo ! ipol
c
c     Do the second term
c
      do ip = 1, ipol
          do ir = 1, nroots
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
            bhi(3) = naoc(ip)
            clo(1) = ir
            chi(1) = ir
            clo(2) = 1
            chi(2) = naoc(ip)
            clo(3) = 1
            chi(3) = naoc(ip)
            if (.not.tda) then
               call tddft_patch3mxm('n','t',0.5d0*omega(ir),1.0d0,
     +                                g_xpy(ip),alo,ahi,
     +                                g_xmy(ip),blo,bhi,
     +                                g_w(ip),clo,chi)
               call tddft_patch3mxm('n','t',0.5d0*omega(ir),1.0d0,
     +                                g_xmy(ip),alo,ahi,
     +                                g_xpy(ip),blo,bhi,
     +                                g_w(ip),clo,chi)
            else
c              g_xpy: X with CIS, Y = 0
c              g_xmy: not created with CIS
               call tddft_patch3mxm('n','t',omega(ir),1.0d0,
     +                                g_xpy(ip),alo,ahi,
     +                                g_xpy(ip),blo,bhi,
     +                                g_w(ip),clo,chi)
            endif
          enddo  ! ir
c
c         Do the third term, drive the parallelization off the 
c         distribution of (X+-Y)ib. The reason for this choice is that
c         the occ-occ block will be localized on a (small) sub-set of
c         the processors. So driving the parallelization off the 
c         distribution of W will result in dreadful load imbalance
c         problems. This choice will lead to a (slightly) less bad
c         load imbalance.
c
c         Do the (X+Y) terms first
c
          call nga_distribution(g_xpy(ip),iproc,blo,bhi)
          doitxpy1 = .not.(blo(1).ne.1.or.bhi(1).ne.nroots)
c
          if (doitxpy1) then
c
           minj = blo(2)
           maxj = bhi(2)
           numj = maxj-minj+1
           mina = blo(3)
           maxa = bhi(3)
           numa = maxa-mina+1
           if (.not.ma_push_get(mt_dbl,nroots*numa*numj,'block b',
     +                         ihdl_b,iptr_b))
     +      call errquit(pname//'failed to allocate blockb',0,MA_ERR)
           ld(1)  = nroots
           ld(2)  = numj
           call nga_get(g_xpy(ip),blo,bhi,dbl_mb(iptr_b),ld)
           ic = 0
           do a = 1, numa
             do j = 1, numj
              imo = nocc(ip)+mina-1+a
              do ir = 1, nroots
                dbl_mb(iptr_b+ic) = 0.5d0*eps(imo,ip)*dbl_mb(iptr_b+ic)
                ic = ic + 1
              enddo
            enddo
           enddo
c
           do icount = 0, nproc-1
           jproc = mod(iproc+icount,nproc)
c
           call nga_distribution(g_xpy(ip),jproc,alo,ahi)
           doitxpy2 = .not.(alo(1).ne.1.or.ahi(1).ne.nroots)
c
           if (doitxpy2.and.doitxpy1) then
            if (alo(3).eq.mina.and.ahi(3).eq.maxa) then
              mini = alo(2)
              maxi = ahi(2)
              numi = maxi-mini+1
c             alo(3) = mina
c             ahi(3) = maxa
              if (.not.ma_push_get(mt_dbl,nroots*numa*numi,'block a',
     +                             ihdl_a,iptr_a)) 
     +         call errquit(pname//'failed to allocate blockb',0,UERR)
              if (.not.ma_push_get(mt_dbl,nroots*numi*numj,'block r',
     +                             ihdl_r,iptr_r)) 
     +         call errquit(pname//'failed to allocate blockr',0,UERR)
              ld(1)  = nroots
              ld(2)  = numi
              call nga_get(g_xpy(ip),alo,ahi,dbl_mb(iptr_a),ld)
              call dfill(nroots*numi*numj,0.0d0,dbl_mb(iptr_r),1)
              do a = 0, numa-1
                do j = 0, numj-1
                  do i = 0, numi-1
                    do ir = 0, nroots-1
                      dbl_mb(iptr_r+ir+nroots*i+nroots*numi*j) 
     +                = dbl_mb(iptr_r+ir+nroots*i+nroots*numi*j) 
     +                - dbl_mb(iptr_a+ir+nroots*i+nroots*numi*a)
     +                * dbl_mb(iptr_b+ir+nroots*j+nroots*numj*a)
                    enddo
                  enddo
                enddo
              enddo
c
c Daniel (12-3-12): We need to account for (X-Y) = X in TDA here.
c We must either double the contribution from (X+Y) here or redo the
c procedure in this block of code.  I prefer scaling here since it's
c less messy.
              if (tda) then
                call dscal(nroots*numi*numj,2.0d0,dbl_mb(iptr_r),1)
              endif
              clo(1) = 1
              chi(1) = nroots
              clo(2) = mini
              chi(2) = maxi
              clo(3) = minj
              chi(3) = maxj
              ld(1) = nroots
              ld(2) = numi
              call nga_acc(g_w(ip),clo,chi,dbl_mb(iptr_r),ld,1.0d0)
              if (.not.ma_pop_stack(ihdl_r)) 
     +         call errquit(pname//'failed to deallocate blockr',0,UERR)
              if (.not.ma_pop_stack(ihdl_a)) 
     +         call errquit(pname//'failed to deallocate blocka',0,UERR)
            endif
           endif !doitxpy1 and doitxpy2
           enddo  ! icount
c
             if (.not.ma_pop_stack(ihdl_b)) 
     +         call errquit(pname//'failed to deallocate blockb',0,UERR)
           endif ! doitxpy1 

          call ga_sync()
c
          if (.not.tda) then
c
c         Do the (X-Y) terms next
c
            call nga_distribution(g_xmy(ip),iproc,blo,bhi)
            doitxmy1 = .not.(blo(1).ne.1.or.bhi(1).ne.nroots)
c
            if (doitxmy1) then
c
             minj = blo(2)
             maxj = bhi(2)
             numj = maxj-minj+1
             mina = blo(3)
             maxa = bhi(3)
             numa = maxa-mina+1
             if (.not.ma_push_get(mt_dbl,nroots*numa*numj,'block b',
     +                           ihdl_b,iptr_b))
     +        call errquit(pname//'failed to allocate blockb',0,UERR)
             ld(1)  = nroots
             ld(2)  = numj
             call nga_get(g_xmy(ip),blo,bhi,dbl_mb(iptr_b),ld)
             ic = 0
             do a = 1, numa
              do i = 1, numj
                imo = nocc(ip)+mina-1+a
                do ir = 1, nroots
                  dbl_mb(iptr_b+ic) = 0.5d0*eps(imo,ip)*
     +                                      dbl_mb(iptr_b+ic)
                  ic = ic + 1
                enddo
              enddo
             enddo
c
            do icount = 0, nproc-1
              jproc = mod(iproc+icount,nproc)
              call nga_distribution(g_xmy(ip),jproc,alo,ahi)
              doitxmy2 = .not.(alo(1).ne.1.or.ahi(1).ne.nroots)
c
               if (doitxmy2.and.doitxmy1) then
               if (alo(3).eq.mina.and.ahi(3).eq.maxa) then
                mini = alo(2)
                maxi = ahi(2)
                numi = maxi-mini+1
c               alo(3) = mina
c               ahi(3) = maxa
                if (.not.ma_push_get(mt_dbl,nroots*numa*numi,'block a',
     +                               ihdl_a,iptr_a)) 
     +           call errquit(pname//'failed to allocate blockb',0,UERR)
                if (.not.ma_push_get(mt_dbl,nroots*numi*numj,'block r',
     +                               ihdl_r,iptr_r)) 
     +           call errquit(pname//'failed to allocate blockr',0,UERR)
                ld(1)  = nroots
                ld(2)  = numi
                call nga_get(g_xmy(ip),alo,ahi,dbl_mb(iptr_a),ld)
                call dfill(nroots*numi*numj,0.0d0,dbl_mb(iptr_r),1)
                do a = 0, numa-1
                  do i = 0, numi-1
                    do j = 0, numj-1
                      do ir = 0, nroots-1
                        dbl_mb(iptr_r+ir+nroots*i+nroots*numi*j) 
     +                  = dbl_mb(iptr_r+ir+nroots*i+nroots*numi*j) 
     +                  - dbl_mb(iptr_a+ir+nroots*i+nroots*numi*a)
     +                  * dbl_mb(iptr_b+ir+nroots*j+nroots*numj*a)
                      enddo
                    enddo
                  enddo
                enddo
                clo(1) = 1
                chi(1) = nroots
                clo(2) = mini
                chi(2) = maxi
                clo(3) = minj
                chi(3) = maxj
                ld(1) = nroots
                ld(2) = numi
                call nga_acc(g_w(ip),clo,chi,dbl_mb(iptr_r),ld,1.0d0)
                if (.not.ma_pop_stack(ihdl_r)) 
     +            call errquit(pname//'failed to deallocate blockr',0,
     +              UERR)
                if (.not.ma_pop_stack(ihdl_a)) 
     +            call errquit(pname//'failed to deallocate blocka',0,
     +              UERR)
              endif
              endif ! doitxmy1 and doitxmy2
            enddo ! icount
c
             if (.not.ma_pop_stack(ihdl_b)) 
     +        call errquit(pname//'failed to deallocate blockb',0,UERR)
            endif ! doitxmy1
c
          call ga_sync()
c
          endif ! tda
c
      enddo ! ipol
c
c     Now do the fourth term
c
c     - Create Puv
c
      idim(1) = nroots*ipol
      idim(2) = nao
      idim(3) = nao
      ichnk(1) = nroots*ipol
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not.nga_create(mt_dbl,3,idim,'vec Puv',ichnk,g_puv))
     +     call errquit(pname//'failed to create g_puv',0,GA_ERR)
c
c     - Puv = sum_pq Cup*Ppq*Cvq
c
      call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroots,1.0d0,0.0d0,"pq",g_mo,g_p,"pq",g_puv)
c
c     - Create global array for (A+B)P in AO basis
c     - Compute (A+B)P in AO basis (currently we have to compute
c       (A-B)P as well although we do not need it...)
c
      if (.not.nga_create(mt_dbl,3,idim,'vec (A+B)P',ichnk,g_apbp))
     +    call errquit(pname//'failed to create g_apbp',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vec (A-B)P',ichnk,g_ambp))
     +    call errquit(pname//'failed to create g_ambp',0,GA_ERR)
      call ga_zero(g_ambp)
      call ga_zero(g_apbp)
c Daniel (2-26-13): It was not obvious that we need to unset
c fock_xc:triplet here for restricted triplet calculations to work.
      if (otriplet) then
        if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1    call errquit(pname//'failed to set triplet',0,RTDB_ERR)
      endif
      call tddft_nga_cont(rtdb,ihdl_geom,ihdl_bfao,g_puv,g_apbp,g_ambp,
     +     nao,ipol,tol2e,tda,oskel,kfac,lhashf,.false.,nroots,iwhich,
     +     lhascd)
c Daniel (2-26-13): Reset fock_xc:triplet here for restricted triplet 
c calculations to work.
      if (otriplet) then
        if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.true.))
     1    call errquit(pname//'failed to set triplet',0,RTDB_ERR)
      endif
c Daniel (1-5-13): This line accidentally halves the contribution of
c this term for CIS, but functions correctly for RPA. 
      if (.not.tda) then
        call ga_add(+0.5d0,g_apbp,+0.5d0,g_ambp,g_apbp)
      endif
c
c     - Destroy global array for Puv
c
      if (.not.ga_destroy(g_puv))
     +    call errquit(pname//'failed to destroy g_puv',0,GA_ERR)
c
c     - Transform (A+B)P from AO basis to Wij in MO basis
c
      call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroots,1.0d0,1.0d0,"ij",g_mo,g_apbp,g_w,"pq")
c
c     - Destroy (A+B)P and (A-B)P
c
      if (.not.ga_destroy(g_apbp))
     +    call errquit(pname//'failed to destroy g_apbp',0,GA_ERR)
      if (.not.ga_destroy(g_ambp))
     +    call errquit(pname//'failed to destroy g_apbp',0,GA_ERR)
c
c     Now do the fifth term (this follows the same approach as in
c     tddft_grad_comp_r)
c
      if (xc_gotxc()) then
        ndens = ipol*(nroots+1)
        ngxc  = ipol*nroots
        if (.not.ma_push_get(mt_int,ngxc,'gxc-s',l_gxc,k_gxc))
     +    call errquit(pname//'failed to allocate l_gxc',0,MA_ERR)
        if (.not.ma_push_get(mt_int,ndens,'densities',l_dens,k_dens))
     +    call errquit(pname//'failed to allocate l_dens',0,MA_ERR)
        if (.not.ma_push_get(mt_int,ndens,'dens-tmps',l_den2,k_den2))
     +    call errquit(pname//'failed to allocate l_den2',0,MA_ERR)
c
c       - Create and calculate the AO basis density matrices
c
        do ip = 0, ipol-1
          int_mb(k_dens+nroots*ipol+ip) = 
     +       ga_create_atom_blocked(ihdl_geom,ihdl_bfao,"d_ao")
        enddo
        call tddft_grad_compute_dao(ipol,nao,nocc,g_mo,
     +                           int_mb(k_dens+nroots*ipol))
        if (ipol.eq.1) call ga_scale(int_mb(k_dens+nroots*ipol),2.0d0)
c
c       - Create and calculate the AO basis representation of (X+Y)
c
        do ip = 0, ipol-1
          do ir = 0, nroots-1
            int_mb(k_dens+ip*nroots+ir) = ga_create_atom_blocked(
     +                                   ihdl_geom,ihdl_bfao,"xpy_ao")
c Daniel (1-14-13): Added this call for consistency with
c tddft_grad_compute_r
            call ga_zero(int_mb(k_dens+ip*nroots+ir))
          enddo
        enddo
        do ip = 1, ipol
          do ir = 1, nroots
            call tddft_grad_trans_mo2ao(1,nao,nfc(ip),naoc(ip),nocc(ip),
     +           nav(ip),nfv(ip),ir,1.0d0,0.0d0,"ib",g_mo(ip),g_xpy(ip),
     +           "ib",int_mb(k_dens+(ip-1)*nroots+ir-1))
            call ga_symmetrize(int_mb(k_dens+(ip-1)*nroots+ir-1))
c Daniel (2-16-13): This line is needed to get matching results from the
c unrestricted code compared to the restricted one.
            if (ipol.eq.1) then
              call ga_scale(int_mb(k_dens+(ip-1)*nroots+ir-1),2.0d0)
            endif
          enddo
        enddo
c
c       - Create and calculate Gxc in AO basis using fock_xc
c
        do i = 0, ngxc-1
          int_mb(k_gxc+i) = ga_create_atom_blocked(ihdl_geom,ihdl_bfao,
     +                                             "gxc_ao")
          call ga_zero(int_mb(k_gxc+i))
        enddo
        if (.not.rtdb_get(rtdb,'fock_xc:calc_type',mt_int,1,calc_type))
     +    calc_type=0
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,5))
     +    call errquit(pname//'failed to set calc_type 5',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type_tddft_w',mt_int,1,
     +      calc_type)) 
     +    call errquit(pname//'failed to set calc_type_tddft_w',
     +       0,RTDB_ERR)
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c Daniel (2-16-13): There is a line in fock_xc, involving lcgmin, that
c will set l3d to true, even though we want it to be false for the 
c gradients.  We need to set the RTDB such that dft:cgmin is true so
c that the l3d is set correctly in fock_xc.
c **********************************************************************
c WE NEED TO FIX THIS BECAUSE TDDFT GRADIENTS DON'T WORK with CGMIN SET.
c **********************************************************************
c Daniel (2-18-13): This fix will not work for optimizations.
c        if (.not.rtdb_put(rtdb,'dft:cgmin',mt_log,1,.true.))
c     1    call errquit(pname//'failed to set cgmin',0,RTDB_ERR)
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
c       Reorder the densities to match the expectation in fock_xc
c       See tddft_grad_comp_r for details.
c
        do i = 0, nroots
          do ip = 0, ipol-1
            int_mb(k_den2+i+ip*(nroots+1)) = int_mb(k_dens+ip+ipol*i)
          enddo
        enddo
        call ga_sync()
c Daniel (2-12-13): We need to set triplet here for fock_xc.
        if (otriplet) then
          if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,otriplet))
     1      call errquit(pname//'failed to set triplet',0,RTDB_ERR)
        endif
        call fock_xc(ihdl_geom, nao, ihdl_bfao, ngxc, int_mb(k_den2),
     +               int_mb(k_gxc), Exc, ipol, .false.)
        call ga_sync()
c DEBUG
      if (util_print('tddft grad w',print_debug)) then
        oroot = ga_nodeid().eq.0
        if (oroot) write(LuOut,*)'DEBUG: '//pname//'gxc'
        call tddft_grad_print_array(ipol,nroots,int_mb(k_gxc),
     +                              dble(ipol))
      endif
c DEBUG
c
        if (.not.rtdb_get(rtdb,'fock_xc:calc_type_tddft_w',mt_int,1,
     +    calc_type)) 
     +    call errquit(pname//'failed to get calc_type_tddft_w',0,
     +        RTDB_ERR)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,calc_type))
     +    call errquit(pname//'failed to reset calc_type',0,RTDB_ERR)
        if (.not.rtdb_delete(rtdb,'fock_xc:calc_type_tddft_w'))
     +    call errquit(pname//'failed to delete calc_type_tddft_w',0,
     +        RTDB_ERR)
c
c       - Transform the Gxc matrices to MO basis and add them to Wij
c
        do i = 0, ndens-1
          if (.not.ga_destroy(int_mb(k_dens+i))) 
     +     call errquit(pname//'failed to destroy densities',0,GA_ERR)
        enddo
        if (.not.ma_pop_stack(l_den2)) 
     +    call errquit(pname//'failed to deallocate l_den2', 0, MA_ERR)
        if (.not.ma_pop_stack(l_dens)) 
     +    call errquit(pname//'failed to deallocate l_dens', 0, MA_ERR)
        do ip = 1, ipol
          do ir = 1, nroots
            call tddft_grad_trans_ao2mo(1,nao,nfc(ip),naoc(ip),nocc(ip),
     +           nav(ip),nfv(ip),ir,1.0d0,1.0d0,"ij",g_mo(ip),
     +           int_mb(k_gxc+(ip-1)*nroots+ir-1),g_w(ip),"ij")
          enddo
        enddo
        do i = 0, ngxc-1
          if (.not.ga_destroy(int_mb(k_gxc+i))) 
     +     call errquit(pname//'failed to destroy gxc_ao', 0, GA_ERR)
        enddo
        if (.not.ma_pop_stack(l_gxc)) 
     +     call errquit(pname//'failed to deallocate l_gxc', 0, MA_ERR)
c
      endif
c
c     Compute the occupied-virtual block
c
c     Wias = eps_is Zias
c          + Sum_j {(X+Y)jasHji+[X+Y]+(X-Y)jasHjis-[X-Y]}
c
c     Do the first term, drive the parallelization off the distribution
c     of Zias
c
      do ip = 1, ipol
c
        call nga_distribution(g_z(ip),iproc,alo,ahi)
        doitz = .not.(alo(1).ne.1.or.ahi(1).ne.nroots)
c
        if (doitz) then
         mini = alo(2)
         maxi = ahi(2)
         numi = maxi-mini+1
         mina = alo(3)
         maxa = ahi(3)
         numa = maxa-mina+1
         if (.not.ma_push_get(mt_dbl,nroots*numi*numa,'block v',
     +                       ihdl_v,iptr_v)) 
     +   call errquit(pname//'failed to allocate block v',0, UERR)
         ld(1) = nroots
         ld(2) = numi
         call nga_get(g_z(ip),alo,ahi,dbl_mb(iptr_v),ld)
         do a = 0, numa-1
          do i = 0, numi-1
            imo = nfc(ip)+mini+i
            do ir = 0, nroots-1
              dbl_mb(iptr_v+ir+nroots*i+nroots*numi*a)
     +        = dbl_mb(iptr_v+ir+nroots*i+nroots*numi*a)
     +        * eps(imo,ip)
            enddo
          enddo
         enddo
        alo(3) = naoc(ip)+mina
        ahi(3) = naoc(ip)+maxa
        call nga_put(g_w(ip),alo,ahi,dbl_mb(iptr_v),ld)
        if (.not.ma_pop_stack(ihdl_v)) 
     +    call errquit(pname//'failed to deallocate block v',0, UERR)
       endif ! doitz 
      enddo  ! ip
c
c     Do the second term
c
c     - Do the (X+Y) and (X-Y) contributions
c
c     - Create global arrays for X and Y in AO basis
c
      idim(1) = nroots*ipol
      idim(2) = nao
      idim(3) = nao
      ichnk(1) = nroots*ipol
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not.nga_create(mt_dbl,3,idim,'vectors Xuv',ichnk,g_x))
     +    call errquit(pname//'failed to create g_x',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors Yuv',ichnk,g_y))
     +    call errquit(pname//'failed to create g_y',0,GA_ERR)
c
c     - Transform (X+Y) to AO basis in g_x and (X-Y) to AO basis
c       in g_y
c
      call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroots,1.0d0,0.0d0,"ib",g_mo,g_xpy,"ib",g_x)
      if (.not.tda) then
        call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +       nroots,1.0d0,0.0d0,"ib",g_mo,g_xmy,"ib",g_y)
      else
        call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +       nroots,1.0d0,0.0d0,"ib",g_mo,g_xpy,"ib",g_y)
      endif
c
c     - Compute X and Y from (X+Y) and (X-Y) in place
c
      alo(1) = 1
      alo(2) = 1
      alo(3) = 1
      ahi(1) = nroots*ipol
      ahi(2) = nao
      ahi(3) = nao
c Daniel (1-7-13): With the modification above, the lines here behave
c like you'd expect (i.e. g_x = X and g_y = Y = 0 for CIS, rather than
c what was here before which made g_x = 0.50*X and g_y = 0.50*X).
      call nga_add_patch(0.5d0,g_x,alo,ahi,0.5d0,g_y,alo,ahi,g_x,
     +                   alo,ahi)
      call nga_add_patch(1.0d0,g_x,alo,ahi,-1.0d0,g_y,alo,ahi,g_y,
     +                   alo,ahi)
c
c     - Allocate various workspace arrays
c
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A+B)X',ichnk,g_apbx))
     +    call errquit(pname//'failed to create g_apbx',0, GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A-B)X',ichnk,g_ambx))
     +    call errquit(pname//'failed to create g_ambx',0, GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A+B)Y',ichnk,g_apby))
     +    call errquit(pname//'failed to create g_apby',0, GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A-B)Y',ichnk,g_amby))
     +    call errquit(pname//'failed to create g_amby',0, GA_ERR)
c
c     - Compute (A+B)X, (A-B)X, (A+B)Y and (A-B)Y
c
      tdaloc=tddft_grad_tdaloc(rtdb,tda,lhashf)
      call tddft_nga_cont(rtdb,ihdl_geom,ihdl_bfao,g_x,g_apbx,g_ambx,
     +   nao,ipol,tol2e,tdaloc,oskel,kfac,lhashf,otriplet,nroots,iwhich,
     +   lhascd)
      call tddft_nga_cont(rtdb,ihdl_geom,ihdl_bfao,g_y,g_apby,g_amby,
     +   nao,ipol,tol2e,tdaloc,oskel,kfac,lhashf,otriplet,nroots,iwhich,
     +   lhascd)
c
c     - Dispose of X and Y
c
      if (.not.ga_destroy(g_x)) 
     +    call errquit(pname//'failed to destroy g_x',0, GA_ERR)
      if (.not.ga_destroy(g_y)) 
     +    call errquit(pname//'failed to destroy g_y',0, GA_ERR)
c
c     - Compute (A+B)(X+Y) and (A-B)(X-Y)
cdbg            (A+B)(X-Y)     (A-B)(X+Y)
c
      call nga_add_patch(1.0d0,g_apbx,alo,ahi,+1.0d0,g_apby,alo,ahi,
     +                   g_apbx,alo,ahi)
      call nga_add_patch(1.0d0,g_ambx,alo,ahi,-1.0d0,g_amby,alo,ahi,
     +                   g_ambx,alo,ahi)
c
c     - Dispose of (A+B)Y and (A-B)Y
c
      if (.not.ga_destroy(g_apby)) 
     +   call errquit(pname//'failed to destroy g_apby',0, GA_ERR)
      if (.not.ga_destroy(g_amby)) 
     +   call errquit(pname//'failed to destroy g_amby',0, GA_ERR)
c
c     - Allocate (A+-B)(X+-Y)ij
c
      do ip = 1, ipol
        idim(1) = nroots
        idim(2) = naoc(ip)
        idim(3) = naoc(ip)
        ichnk(1) = nroots
        ichnk(2) = -1
        ichnk(3) = -1
        if (.not.nga_create(mt_dbl,3,idim,'vec hij',ichnk,g_hij(ip)))
     +    call errquit(pname//'failed to create g_hij',0, GA_ERR)
      enddo
c
c     - Transform (A+B)(X+Y) to MO basis occupied-occupied block only
cdbg              (A+B)(X-Y)
c
      call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroots,1.0d0,0.0d0,"ij",g_mo,g_apbx,g_hij,"ij")
c
c     - Add sum_j (X+Y)ja [(A+B)(X+Y)ji] to Wia
cdbg              (X-Y)ja [(A+B)(X-Y)ji]
c       use symmetry of [(A+B)(X+Y)ji]
c
      do ip = 1, ipol
        do ir=1,nroots
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
          clo(3) = naoc(ip)+1
          chi(3) = naoc(ip)+nav(ip)
          call tddft_patch3mxm('n','n',1.0d0,1.0d0,g_hij(ip),alo,ahi,
     +         g_xpy(ip),blo,bhi,g_w(ip),clo,chi)
        enddo
      enddo
c
c     - Transform (A-B)(X-Y) to MO basis occupied-occupied block only
cdbg              (A-B)(X+Y)
c
      call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroots,1.0d0,0.0d0,"ij",g_mo,g_ambx,g_hij,"ij")
c
c     - Add sum_j (X-Y)ja [(A-B)(X-Y)ji] to Wia
cdbg              (X+Y)js [(A-B)(X+Y)ji]
c       use anti-symmetry of [(A-B)(X-Y)ji]
c
c Daniel (1-7-13): We can definitely avoid this part of the routine if
c we don't use exact exchange, since the linear transformation H^-[V]
c is zero in that case.  
      if (lhashf) then
        do ip = 1, ipol
          do ir = 1, nroots
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
            clo(3) = naoc(ip)+1
            chi(3) = naoc(ip)+nav(ip)
c Daniel (1-7-13): Manipulate the code here for CIS to use g_xpy here
c since (X+Y) = (X-Y) = X.  This is a consequence of not allocating 
c a g_xmy array for CIS.  The linear transformation H^-[X] still exists
c in CIS.
            if (.not.tda) then
              call tddft_patch3mxm('n','n',-1.0d0,1.0d0,
     +             g_hij(ip),alo,ahi,
     +             g_xmy(ip),blo,bhi,g_w(ip),clo,chi)
            else
              call tddft_patch3mxm('n','n',-1.0d0,1.0d0,
     +             g_hij(ip),alo,ahi,
     +             g_xpy(ip),blo,bhi,g_w(ip),clo,chi)
            endif
          enddo
        enddo
      endif ! lhashf
c
c     - Dispose of g_hij
c
      do ip = 1, ipol
        if (.not.ga_destroy(g_hij(ip))) 
     +   call errquit(pname//'failed to destroy g_hij',0, GA_ERR)
      enddo
c
c     - Dispose of (A+B)(X+Y) and (A+B)(X-Y)
c
      if (.not.ga_destroy(g_apbx)) 
     +   call errquit(pname//'failed to destroy g_apbx',0, GA_ERR)
      if (.not.ga_destroy(g_ambx)) 
     +   call errquit(pname//'failed to destroy g_ambx',0, GA_ERR)
c
c
c     Compute the virtual-virtual block
c
c     Wabs = Omega Sum_i  [(X+Y)ias(X-Y)ibs+(X-Y)ias(X+Y)ibs]/2
c          + Sum_i eps_is [(X+Y)ias(X+Y)ibs+(X-Y)ias(X-Y)ibs]/2
c
      do ip = 1, ipol
          do ir = 1, nroots
            alo(1) = ir
            ahi(1) = ir
            alo(2) = 1
            ahi(2) = nav(ip)
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
            clo(2) = naoc(ip)+1
            chi(2) = naoc(ip)+nav(ip)
            clo(3) = naoc(ip)+1
            chi(3) = naoc(ip)+nav(ip)
c
            if (.not.tda) then
            call tddft_patch3mxm('t','n',0.5d0*omega(ir),0.0d0,
     +                              g_xpy(ip),alo,ahi,
     +                              g_xmy(ip),blo,bhi,
     +                              g_w(ip),clo,chi)
            call tddft_patch3mxm('t','n',0.5d0*omega(ir),1.0d0,
     +                              g_xmy(ip),alo,ahi,
     +                              g_xpy(ip),blo,bhi,
     +                              g_w(ip),clo,chi)
            else
c            g_xpy: X with CIS, Y = 0
c            g_xmy: not created with CIS
             call tddft_patch3mxm('t','n',omega(ir),0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_w(ip),clo,chi)
            endif ! tda
          enddo ! ir
c
c         Do the second term, drive the parallelization off the 
c         distribution of (X+-Y)ib
c
c         Do the (X+Y) terms first
c
          call nga_distribution(g_xpy(ip),iproc,blo,bhi)
          doitxpy1 = .not.(blo(1).ne.1.or.bhi(1).ne.nroots)
c
          if (doitxpy1) then
c
           mini = blo(2)
           maxi = bhi(2)
           numi = maxi-mini+1
           minb = blo(3)
           maxb = bhi(3)
           numb = maxb-minb+1
           if (.not.ma_push_get(mt_dbl,nroots*numb*numi,'block b',
     +                         ihdl_b,iptr_b))
     +      call errquit(pname//'failed to allocate blockb',0,UERR)
           ld(1)  = nroots
           ld(2)  = numi
           call nga_get(g_xpy(ip),blo,bhi,dbl_mb(iptr_b),ld)
           ic = 0
           do b = 1, numb
            do i = 1, numi
              imo = nfc(ip)+mini-1+i
              do ir = 1, nroots
                dbl_mb(iptr_b+ic) = 0.5d0*eps(imo,ip)*dbl_mb(iptr_b+ic)
                ic = ic + 1
              enddo
            enddo
           enddo
c
           do icount = 0, nproc-1
           jproc = mod(iproc+icount,nproc)
c
           call nga_distribution(g_xpy(ip),jproc,alo,ahi)
           doitxpy2 = .not.(alo(1).ne.1.or.ahi(1).ne.nroots)
c
           if (doitxpy2.and.doitxpy1) then
           if (alo(2).eq.mini.and.ahi(2).eq.maxi) then
           mina = alo(3)
           maxa = ahi(3)
           numa = maxa-mina+1
           alo(2) = mini
           ahi(2) = maxi
           if (.not.ma_push_get(mt_dbl,nroots*numa*numi,'block a',
     +                           ihdl_a,iptr_a)) 
     +      call errquit(pname//'failed to allocate blockb',0, MA_ERR)
            if (.not.ma_push_get(mt_dbl,nroots*numa*numb,'block r',
     +                           ihdl_r,iptr_r)) 
     +      call errquit(pname//'failed to allocate blockr',0, MA_ERR)
            ld(1)  = nroots
            ld(2)  = numi
            call nga_get(g_xpy(ip),alo,ahi,dbl_mb(iptr_a),ld)
            call dfill(nroots*numa*numb,0.0d0,dbl_mb(iptr_r),1)
            do b = 0, numb-1
              do a = 0, numa-1
                do i = 0, numi-1
                  do ir = 0, nroots-1
                    dbl_mb(iptr_r+ir+nroots*a+nroots*numa*b) 
     +              = dbl_mb(iptr_r+ir+nroots*a+nroots*numa*b) 
     +              + dbl_mb(iptr_a+ir+nroots*i+nroots*numi*a)
     +              * dbl_mb(iptr_b+ir+nroots*i+nroots*numi*b)
                  enddo
                enddo
              enddo
            enddo
c Daniel (12-3-12): We need to account for (X-Y) = X in TDA here.
c We must either double the contribution from (X+Y) here or redo the
c procedure in this block of code.  I prefer scaling here since it's
c less messy.
            if (tda) then
              call dscal(nroots*numa*numb,2.0d0,dbl_mb(iptr_r),1)
            endif
            clo(1) = 1
            chi(1) = nroots
            clo(2) = naoc(ip)+mina
            chi(2) = naoc(ip)+maxa
            clo(3) = naoc(ip)+minb
            chi(3) = naoc(ip)+maxb
            ld(1) = nroots
            ld(2) = numa
            call nga_acc(g_w(ip),clo,chi,dbl_mb(iptr_r),ld,1.0d0)
            if (.not.ma_pop_stack(ihdl_r)) 
     +       call errquit(pname//'failed to deallocate blockr',0,MA_ERR)
            if (.not.ma_pop_stack(ihdl_a)) 
     +       call errquit(pname//'failed to deallocate blocka',0,MA_ERR)
            endif
           endif ! doitxpy1 and doitxpy2
          enddo ! icount
c
           if (.not.ma_pop_stack(ihdl_b)) 
     +      call errquit(pname//'failed to deallocate blockb',0,MA_ERR)
          endif  ! doitxpy1
c
          call ga_sync()
c
          if (.not.tda) then
c
c         Do the (X-Y) terms next
c
          call nga_distribution(g_xmy(ip),iproc,blo,bhi)
          doitxmy1 = .not.(blo(1).ne.1.or.bhi(1).ne.nroots)
c
          if (doitxmy1) then
c
           mini = blo(2)
           maxi = bhi(2)
           numi = maxi-mini+1
           minb = blo(3)
           maxb = bhi(3)
           numb = maxb-minb+1
           if (.not.ma_push_get(mt_dbl,nroots*numb*numi,'block b',
     +                         ihdl_b,iptr_b))
     +      call errquit(pname//'failed to allocate blockb',0,UERR)
           ld(1)  = nroots
           ld(2)  = numi
           call nga_get(g_xmy(ip),blo,bhi,dbl_mb(iptr_b),ld)
           ic = 0
           do b = 1, numb
            do i = 1, numi
              imo = nfc(ip)+mini-1+i
              do ir = 1, nroots
                dbl_mb(iptr_b+ic) = 0.5d0*eps(imo,ip)*dbl_mb(iptr_b+ic)
                ic = ic + 1
              enddo
            enddo
           enddo
c
           do icount = 0, nproc-1
            jproc = mod(iproc+icount,nproc)
c
            call nga_distribution(g_xmy(ip),jproc,alo,ahi)
            doitxmy2 = .not.(alo(1).ne.1.or.ahi(1).ne.nroots)
c
            if (doitxmy2) then
            if (alo(2).eq.mini.and.ahi(2).eq.maxi) then
            mina = alo(3)
            maxa = ahi(3)
            numa = maxa-mina+1
            alo(2) = mini
            ahi(2) = maxi
            if (.not.ma_push_get(mt_dbl,nroots*numa*numi,'block a',
     +                           ihdl_a,iptr_a)) 
     +      call errquit(pname//'failed to allocate blockb',0, UERR)
            if (.not.ma_push_get(mt_dbl,nroots*numa*numb,'block r',
     +                          ihdl_r,iptr_r)) 
     +      call errquit(pname//'failed to allocate blockr',0, UERR)
            ld(1)  = nroots
            ld(2)  = numi
            call nga_get(g_xmy(ip),alo,ahi,dbl_mb(iptr_a),ld)
            call dfill(nroots*numa*numb,0.0d0,dbl_mb(iptr_r),1)
            do b = 0, numb-1
              do a = 0, numa-1
                do i = 0, numi-1
                  do ir = 0, nroots-1
                    dbl_mb(iptr_r+ir+nroots*a+nroots*numa*b) 
     +              = dbl_mb(iptr_r+ir+nroots*a+nroots*numa*b) 
     +              + dbl_mb(iptr_a+ir+nroots*i+nroots*numi*a)
     +              * dbl_mb(iptr_b+ir+nroots*i+nroots*numi*b)
                  enddo
                enddo
              enddo
            enddo
            clo(1) = 1
            chi(1) = nroots
            clo(2) = naoc(ip)+mina
            chi(2) = naoc(ip)+maxa
            clo(3) = naoc(ip)+minb
            chi(3) = naoc(ip)+maxb
            ld(1) = nroots
            ld(2) = numa
            call nga_acc(g_w(ip),clo,chi,dbl_mb(iptr_r),ld,1.0d0)
            if (.not.ma_pop_stack(ihdl_r)) 
     +      call errquit(pname//'failed to deallocate blockr',0,MA_ERR)
            if (.not.ma_pop_stack(ihdl_a)) 
     +      call errquit(pname//'failed to deallocate blocka',0,MA_ERR)
            endif
           endif ! doitxmy1 and doitxmy2
          enddo ! icount
c
            if (.not.ma_pop_stack(ihdl_b)) 
     +       call errquit(pname//'failed to deallocate blockb',0,MA_ERR)
          endif !doitxmy1

          call ga_sync()

          endif  ! tda
c
      enddo
c
c     Copy Wia to Wai
c
      do ip = 1, ipol
        alo(1) = 1
        ahi(1) = nroots
        alo(2) = 1
        ahi(2) = naoc(ip)
        alo(3) = naoc(ip)+1
        ahi(3) = naoc(ip)+nav(ip)
        call nga_scale_patch(g_w(ip),alo,ahi,0.5d0)
        idim(1) = naoc(ip)
        idim(2) = nav(ip)
        ichnk(1) = -1
        ichnk(2) = -1
        if (.not.nga_create(mt_dbl,2,idim,"transpose",ichnk,g_t))
     +     call errquit(pname//'failed to create g_t',0,GA_ERR)
        do ir = 1, nroots
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = naoc(ip)
          alo(3) = naoc(ip)+1
          ahi(3) = naoc(ip)+nav(ip)
          blo(1) = 1
          bhi(1) = naoc(ip)
          blo(2) = 1
          bhi(2) = nav(ip)
          call nga_copy_patch('n',g_w(ip),alo,ahi,g_t,blo,bhi)
          alo(2) = naoc(ip)+1
          ahi(2) = naoc(ip)+nav(ip)
          alo(3) = 1
          ahi(3) = naoc(ip)
          bhi(2) = nav(ip)
          bhi(1) = naoc(ip)
          call nga_copy_patch('t',g_t,blo,bhi,g_w(ip),alo,ahi)
        enddo
        if (.not.ga_destroy(g_t))
     +     call errquit(pname//'failed to destroy g_t',0,GA_ERR)
      enddo
      call ga_sync()
c
      if (util_print('tddft grad w',print_debug)) then
        oroot = ga_nodeid().eq.0
        if (oroot) write(LuOut,*)'DEBUG: '//pname//'W'
        call tddft_grad_print_array(ipol,nroots,g_w,dble(ipol))
      endif
c
      end
c $Id$
