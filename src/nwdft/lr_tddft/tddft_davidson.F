      subroutine tddft_davidson(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1  nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     2  g_buffer1,g_buffer2,d_trials,d_prod_p,d_prod_m,
     3  g_movecs,evl,ntrials,nproducts,
     4  tda,oskel,singlet,triplet,tol2e,kfac,lhashf,thresh,
     5  maxiter,target,targetsym,algorithm,mult,l_irs,k_irs)
c
c $Id$
c
c Iteratively solve TDDFT matrix eigenvalue problem using 
c a modified Davidson trial vector subspace diagonalization
c algorithm.
c
c ! Caution ! shell_fock_build must work with antisym=.true. 
c for CIS, TDHF, and hybrid TDDFT.
c
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
c ---------------------------------
c TDDFT matrix eigenvalue equations 
c ---------------------------------
c
c (Real orbitals are assumed):
c
c (A-B)^(1/2)(A+B)(A-B)^(1/2)Z = omega^2 Z, 
c Z = (A-B)^(-1/2)(X+Y),
c
c A_ia(a),jb(a) = delta_ij delta_ab (epsilon_a(a)-epsilon_i(a))
c               + (i(a)a(a)|j(a)b(a)) - kfac(i(a)j(a)|a(a)b(a))
c               + (i(a)a(a)|w(aa)|j(a)b(a)),
c A_ia(a),jb(b) = (i(a)a(a)|j(b)b(b))
c               + (i(a)a(a)|w(ab)|j(b)b(b)),
c B_ia(a),jb(a) = (i(a)a(a)|j(a)b(a)) - kfac(i(a)b(a)|j(a)a(a))
c               + (i(a)a(a)|w(aa)|j(a)b(a)),
c B_ia(a),jb(b) = (i(a)a(a)|j(b)b(b))
c               + (i(a)a(a)|w(ab)|j(b)b(b)),
c
c where (a),(b) mean alpha and beta spins.  When (tda) B=0, so
c
c AZ = omega Z.
c
c When there is no HF exchange contribution,
c
c (A-B) = delta_ij delta_ab (epsilon_a(a)-epsilon_i(a)),
c
c so
c 
c C(A+B)Z = omega^2 Z
c
c C_ia(a),ia(a) = (epsilon_a(a)-epsilon_i(a)) otherwise 0,
c
c also all trial density matrices can be symmetrized for fock_2e,
c since (pq|rs)=(pq|sr), (pq|w|rs)=(pq|w|sr) for real orbitals.
c
c When there is HF exchange contribution, trial density matrices
c that are in general non-symmetric cannot be symmetrized.  One
c must decompose such non-symmetric matrices into sums of symmetric
c and anti-symmetric matrices and rely upon fock_2e with asym=.true.
c to work.
c
c -------------------------------------------
c Subspace representation of A and B matrices
c -------------------------------------------
c
c p-th trial vectors
c 
c T^p_ia(a) and T^p_ia(b),
c
c are independent for unrestricted TDDFT (open-shells)
c and are dependent for restricted TDDFT (closed-shells)
c in such ways as
c
c T(a) = T(b)  : singlet excited states
c T(a) = -T(b) : triplet excited states
c
c [+]_ia,jb = A_ia,jb + B_ia,jb
c [-]_ia,jb = A_ia,jb - B_ia,jb (diagonal if pure DFT)
c
c p-th product vectors
c
c P^p_ia(a) = [+]_ia(a),jb(a) T^p_jb(a) + [+]_ia(a),jb(b) T^p_jb(b)
c Q^p_ia(a) = [-]_ia(a),jb(a) T^p_jb(a) + [-]_ia(a),jb(b) T^p_jb(b)
c
c These are carried out with fock_2e with the option requesting
c a CPHF LHS calculation.
c
c Subspace [+]_pq = T^p_ia(a) P^q_ia(q) + T^p_ia(b) P^q_ia(b)
c Subspace [-]_pq = T^p_ia(a) Q^q_ia(q) + T^p_ia(b) Q^q_ia(b)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
c
      integer rtdb             ! Runtime database handle
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer nroots           ! Number of roots sought
      integer nresiduals       ! Number of residual vectors
      integer maxvecs          ! Maximum number of trial vectors
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for A+B product vectors
      integer g_prod_m(2)      ! GA handle for A-B product vectors
      integer g_buffer1(2)     ! GA handle for buffer vectors
      integer g_buffer2(2)     ! GA handle for buffer vectors
      integer d_trials(2)      ! DRA handle for trial vectors
      integer d_prod_p(2)      ! DRA handle for A+B product vectors
      integer d_prod_m(2)      ! DRA handle for A-B product vectors
      integer g_residuals_1(2) ! GA handle for residual vectors 1
      integer g_residuals_2(2) ! GA handle for residual vectors 2
      integer g_movecs(2)      ! GA handle for MO eigenvectors
      double precision evl(nbf_ao,2)   ! Eigenvalues
      integer ntrials          ! Current number of trial vectors
      integer nproducts        ! Current number of product vectors
      logical tda              ! True if Tamm-Dancoff approximation
      logical singlet          ! True if singlet excited states
      logical triplet          ! True if triplet excited states
      logical oskel            ! True if use symmetry in fock_2e
      double precision tol2e   ! 2-e integral tolerance
      double precision kfac    ! Prefactor for HF exchange contribution
      logical lhashf           ! True if XC has HF contribution
      double precision thresh  ! Convergence threshold for Davidson iteration
      integer maxiter          ! Maximum number of Davidson iterations
      integer target           ! Target excited state
      character*4 targetsym    ! Symmetry of target excited state
      double precision cwall   ! Wall timings
      double precision ccpu    ! CPU timings
      double precision dcpu,dwall ! Davidson    ! CPU timings
      integer algorithm        ! =1 multiple fock build; =2: conventional
      integer mult             ! Ground state spin multiplicity
      integer l_irs(2),k_irs(2)! Irreps
c
      integer nov(2)           ! Number of occupied virtuals pairs
      integer g_td(2)          ! GA handles for work trial densities
      integer g_pd(2)          ! GA handles for work A+B product densities
      integer g_qd(2)          ! GA handles for work A-B product densities
      integer l_residual_1,k_residual_1 ! MA handles for norm of residuals
      integer l_residual_2,k_residual_2 ! MA handles for norm of residuals
      integer iter             ! Iteration counter
      integer l_apb,k_apb      ! MA handles for (A+B) matrix
      integer l_apbvec,k_apbvec! MA handles for (A+B) eigenvector matrix
      integer l_apbval,k_apbval! MA handles for (A+B) eigenvalue vector
      integer l_amb,k_amb      ! MA handles for (A-B) matrix
      integer l_ambvec,k_ambvec! MA handles for (A-B) eigenvector matrix
      integer l_ambval,k_ambval! MA handles for (A-B) eigenvalue vector
      integer l_valsav,k_valsav! MA handles for eigenvalues from previous iteration
      double precision maxresidual       ! Max norm of residual vectors
      double precision maxomega          ! Max change in excitation energies
      logical converged        ! True if convergence criterion is met
      integer nconverged       ! Number of converged roots
      integer nnewtrials       ! Number of trial vectors in the next iteration
      integer memleft          ! GA memory space available
      integer g_ntrl           ! NGA handle for trial vectors
      integer g_nprd_p         ! NGA handle for (A+B) product vectors
      integer g_nprd_m         ! NGA handle for (A-B) product vectors
      integer idim(3)          ! Dimension array for NGA create
      integer ichunk(3)        ! Chunk array for NGA create
      integer alo(2),ahi(2)    ! Regions for patch
      integer blo(3),bhi(3)    ! Regions for patch
      integer nvectors         ! Number of trial vectors to be processed at once
      logical nodezero
      integer i,n
      integer dummy,request
      integer vshift
      logical lecut
      double precision ecut
      double precision diff_max
      logical lewin
      double precision ewinl,ewinh,etmp
      logical lowin
      integer owstart(2)
      integer owend(2)
      logical lcheck
      integer ispintag
      integer iwhich
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      nodezero=(ga_nodeid().eq.0)
      iwhich = 0   ! call to tddft_nga_cont()
c
c     Set spintag
      if (singlet) then 
          ispintag = 1
      else if (triplet) then 
          ispintag = 2
      else 
          ispintag = 3
      end if
c
c -------------------------
c Allocate density matrices
c -------------------------
c
c Determine the length of the trial vector: Occ*Virt
      call tddft_lentrialvec(ipol,nocc,nmo,nfc,nfv,nov)
c
c     nresiduals=min(nroots+min(nroots/5,10),maxvecs)
      nresiduals=nroots
      do i=1,ipol
        g_td(i)=ga_create_atom_blocked
     1    (geom,ao_bas_han,'Trial density')
        g_pd(i)=ga_create_atom_blocked
     1    (geom,ao_bas_han,'A product density')
        if (.not.tda) 
     1    g_qd(i)=ga_create_atom_blocked
     2      (geom,ao_bas_han,'B product density')
        if (.not.ga_create(mt_dbl,max(1,nov(i)),nresiduals,
     1    'residual vectors 1',-1,nresiduals,g_residuals_1(i)))
     2    call errquit
     3    ('tddft_davidson: failed to create g_residuals_1',0, GA_ERR)
        call ga_zero(g_residuals_1(i))
        if (.not.tda) then
          if (.not.ga_create(mt_dbl,max(1,nov(i)),nresiduals,
     1      'residual vectors 2',-1,nresiduals,g_residuals_2(i)))
     2      call errquit
     3      ('tddft_davidson: failed to create g_residuals_2',0, GA_ERR)
          call ga_zero(g_residuals_2(i))
        endif
      enddo
c
      if (.not.ma_push_get(mt_dbl,nroots,'last eigenvalues',
     1  l_valsav,k_valsav)) call errquit
     2  ('tddft_davidson: failed to allocate valsav',0, MA_ERR)
      do i=1,nroots
        dbl_mb(k_valsav+i-1)=1.0d99
      enddo
      if (.not.ma_push_get(mt_dbl,nresiduals,'residuals',
     1  l_residual_1,k_residual_1)) call errquit
     2  ('tddft_davidson: failed to allocate residual_1',0, MA_ERR)
      if (.not.tda) then
        if (.not.ma_push_get(mt_dbl,nresiduals,'residuals',
     1    l_residual_2,k_residual_2)) call errquit
     2    ('tddft_davidson: failed to allocate residual_2',0, MA_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tddft:check',mt_log,1,lcheck))
     &   lcheck=.true.
c
c --------------------------------------------
c Get the energy cutoff on the occupied states
c --------------------------------------------
      if (.not.rtdb_get(rtdb,'tddft:lecut',mt_log,1,lecut))
     &   lecut = .false.
      if (.not.rtdb_get(rtdb,'tddft:ecut',mt_dbl,1,ecut)) then
         ecut  = 0.d0
         lecut = .false.
      end if
      if (.not.rtdb_get(rtdb,'tddft:diff_max',mt_dbl,1,diff_max))
     &   diff_max = 100000.d0
c
c --------------------------------------------
c Get reference virtual state
c --------------------------------------------
      if (.not.rtdb_get(rtdb,'tddft:vshift',mt_int,1,vshift))
     &   vshift = 0
c
c --------------------------------------------
c Get the energy window on the occupied states
c --------------------------------------------
      if (.not.rtdb_get(rtdb,'tddft:lewin',mt_log,1,lewin))
     &   lewin = .false.
      if (.not.rtdb_get(rtdb,'tddft:ewinl',mt_dbl,1,ewinl)) then
         ewinl  = 0.d0
         lewin = .false.
      end if
      if (.not.rtdb_get(rtdb,'tddft:ewinh',mt_dbl,1,ewinh)) then
         ewinh  = 0.d0
         lewin = .false.
      end if
      if (ewinl .gt. ewinh) then
         etmp = ewinl
         ewinl = ewinh
         ewinh = etmp
      end if
c
c ------------------------------------------------
c Get window of occupied orbitals and set defaults
c ------------------------------------------------
      if (.not.rtdb_get(rtdb,'tddft:lowin',mt_log,1,lowin))
     &  lowin = .false.
      if (.not.rtdb_get(rtdb,'tddft:owstart',mt_int,2,owstart))
     &  lowin = .false.
      if (.not.rtdb_get(rtdb,'tddft:owend',mt_int,2,owend))
     &  lowin = .false.
c
c -------------------------
c Top of Davidson iteration
c -------------------------
c
      dcpu=-util_cpusec()
      dwall=-util_wallsec()
      ccpu=-util_cpusec()
      cwall=-util_wallsec()
      if (nodezero.and.util_print('iteration',print_default)) then
        write(LuOut,9000)
        if (ipol.eq.1) then
          if (singlet) write(LuOut,9010)
          if (triplet) write(LuOut,9020)
          if (lowin) write(LuOut,9040) owstart(1), owend(1)
        else
          write(LuOut,9030)
        endif
        write(LuOut,9100)
        call util_flush(LuOut)
      endif
 9000 format(/,2x,'Entering Davidson iterations')
 9010 format(2x,'Restricted singlet excited states')
 9020 format(2x,'Restricted triplet excited states')
 9030 format(2x,'Unrestricted')
 9040 format(2x,'Orbital Window',2x,i4,2x,i4)
 9100 format(/,
     1  '  Iter ','  NTrls ','  NConv ',
     2  '   DeltaV  ','   DeltaE  ','    Time   ',/,
     3  '  ---- ',' ------ ',' ------ ',
     4  ' --------- ',' --------- ',' --------- ')
      iter=0
   10 continue
        iter=iter+1
        if (iter.gt.maxiter) call errquit
     1    ('tddft_davidson: exceeded maxiter',0, CALC_ERR)
c
c ================================================
c Algorithms 1 & 3 (multiple product vector build)
c ================================================
c
        if ((algorithm.eq.1).or.(algorithm.eq.3)) then
c
c -----------------------------------------------
c Allocate NGA trial and product density matrices
c -----------------------------------------------
c
          nvectors=ntrials-nproducts
          idim(1)=ipol*nvectors
          idim(2)=nbf_ao
          idim(3)=nbf_ao
          ichunk(1)=ipol*nvectors
          ichunk(2)=-1
          ichunk(3)=-1
          if (.not.nga_create(mt_dbl,3,idim,'Work trial',
     1      ichunk,g_ntrl)) call errquit
     2      ('tddft_davidson: failed to create g_ntrl',0, GA_ERR)
          if (.not.nga_create(mt_dbl,3,idim,'Work product',
     1      ichunk,g_nprd_p)) call errquit
     2      ('tddft_davidson: failed to create g_nprd_p',0, GA_ERR)
          if (.not.tda) then
            if (.not.nga_create(mt_dbl,3,idim,'Work product',
     1        ichunk,g_nprd_m)) call errquit
     2        ('tddft_davidson: failed to create g_nprd_m',0, GA_ERR)
          endif
c
c -----------------------
c MO to AO transformation
c -----------------------
c
          do n=nproducts+1,ntrials
            if (algorithm.eq.1) then
              call tddft_transfm(n,g_trials,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_td,1)
            else
              do i=1,ipol
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_trials(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
              enddo
              call tddft_transfm(1,g_buffer1,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_td,1)
            endif
            do i=1,ipol
              alo(1)=1
              ahi(1)=nbf_ao
              alo(2)=1
              ahi(2)=nbf_ao
              blo(1)=(i-1)*nvectors+(n-nproducts)
              bhi(1)=(i-1)*nvectors+(n-nproducts)
c             blo(1)=(n-nproducts-1)*ipol+i
c             bhi(1)=(n-nproducts-1)*ipol+i
              blo(2)=1
              bhi(2)=nbf_ao
              blo(3)=1
              bhi(3)=nbf_ao
              call nga_copy_patch('N',g_td(i),alo,ahi,g_ntrl,blo,bhi)
            enddo
          enddo
          ccpu=ccpu+util_cpusec()
          cwall=cwall+util_wallsec()
          if (nodezero.and.util_print('timings',print_high))
     1      write(LuOut,9130) 'Transform',ccpu,cwall
          ccpu=-util_cpusec()
          cwall=-util_wallsec()
          memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
          if (nodezero.and.util_print('trial vectors',print_high))
     1      write(LuOut,9200) memleft
c
c ---------------------------------------------------------
c Contract with integrals and make product density matrices
c ---------------------------------------------------------
c
          call tddft_nga_cont(rtdb,geom,ao_bas_han,
     1      g_ntrl,g_nprd_p,g_nprd_m,nbf_ao,ipol,tol2e,tda,oskel,kfac,
     2      lhashf,triplet,nvectors,iwhich)
          ccpu=ccpu+util_cpusec()
          cwall=cwall+util_wallsec()
          if (nodezero.and.util_print('timings',print_high))
     1      write(LuOut,9130) 'Contract',ccpu,cwall
          ccpu=-util_cpusec()
          cwall=-util_wallsec()
          memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
          if (nodezero.and.util_print('trial vectors',print_high))
     1      write(LuOut,9200) memleft
c
c -----------------------
c AO to MO transformation
c -----------------------
c
          do n=nproducts+1,ntrials
            do i=1,ipol
              alo(1)=1
              ahi(1)=nbf_ao
              alo(2)=1
              ahi(2)=nbf_ao
              blo(1)=(i-1)*nvectors+(n-nproducts)
              bhi(1)=(i-1)*nvectors+(n-nproducts)
c             blo(1)=(n-nproducts-1)*ipol+i
c             bhi(1)=(n-nproducts-1)*ipol+i
              blo(2)=1
              bhi(2)=nbf_ao
              blo(3)=1
              bhi(3)=nbf_ao
              call nga_copy_patch('N',g_nprd_p,blo,bhi,g_pd(i),alo,ahi)
              if (.not.tda) call nga_copy_patch
     1          ('N',g_nprd_m,blo,bhi,g_qd(i),alo,ahi)
            enddo
            if (algorithm.eq.1) then
              call tddft_transfm(n,g_prod_p,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_pd,2)
              if (.not.tda) then
                call tddft_transfm(n,g_prod_m,g_movecs,
     1            nbf_ao,nocc,nmo,nfc,nfv,ipol,g_qd,2)
              endif
            else
              call tddft_transfm(1,g_buffer1,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_pd,2)
              do i=1,ipol
                dummy=dra_write_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_prod_p(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
              enddo
              if (.not.tda) then
                call tddft_transfm(1,g_buffer1,g_movecs,
     1            nbf_ao,nocc,nmo,nfc,nfv,ipol,g_qd,2)
                do i=1,ipol
                  dummy=dra_write_section(.false.,
     1              g_buffer1(i),1,nov(i),1,1,
     2              d_prod_m(i),1,nov(i),n,n,request)
                  dummy=dra_wait(request)
                enddo
              endif
            endif
          enddo
          ccpu=ccpu+util_cpusec()
          cwall=cwall+util_wallsec()
          if (nodezero.and.util_print('timings',print_high))
     1      write(LuOut,9130) 'Transform',ccpu,cwall
          ccpu=-util_cpusec()
          cwall=-util_wallsec()
          memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
          if (nodezero.and.util_print('trial vectors',print_high))
     1      write(LuOut,9200) memleft
c
c --------------------------------------
c Add (epsilon_a-epsilon_i) contribution
c --------------------------------------
c
          do n=nproducts+1,ntrials
            if (algorithm.eq.1) then
              call tddft_epsilon(n,g_trials,g_prod_p,
     1          ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2          lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
              if (.not.tda) then
                call tddft_epsilon(n,g_trials,g_prod_m,
     1            ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2            lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
              endif
            else
              do i=1,ipol
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_trials(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
                dummy=dra_read_section(.false.,
     1            g_buffer2(i),1,nov(i),1,1,
     2            d_prod_p(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
              enddo
              call tddft_epsilon(1,g_buffer1,g_buffer2,
     1          ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2          lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
              do i=1,ipol
                dummy=dra_write_section(.false.,
     1            g_buffer2(i),1,nov(i),1,1,
     2            d_prod_p(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
              enddo
              if (.not.tda) then
                do i=1,ipol
                  dummy=dra_read_section(.false.,
     1              g_buffer1(i),1,nov(i),1,1,
     2              d_trials(i),1,nov(i),n,n,request)
                  dummy=dra_wait(request)
                  dummy=dra_read_section(.false.,
     1              g_buffer2(i),1,nov(i),1,1,
     2              d_prod_m(i),1,nov(i),n,n,request)
                  dummy=dra_wait(request)
                enddo
                call tddft_epsilon(1,g_buffer1,g_buffer2,
     1            ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2            lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
                do i=1,ipol
                  dummy=dra_write_section(.false.,
     1              g_buffer2(i),1,nov(i),1,1,
     2              d_prod_m(i),1,nov(i),n,n,request)
                  dummy=dra_wait(request)
                enddo
              endif
            endif
            if (util_print('trial vectors',print_debug))
     1        call tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     2          g_prod_p,g_prod_m,g_buffer1,d_trials,d_prod_p,
     3          d_prod_m,n,ntrials,n,tda,algorithm)
          enddo
          nproducts=ntrials
          ccpu=ccpu+util_cpusec()
          cwall=cwall+util_wallsec()
          if (nodezero.and.util_print('timings',print_high))
     1      write(LuOut,9130) 'Epsilon',ccpu,cwall
          ccpu=-util_cpusec()
          cwall=-util_wallsec()
          memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
          if (nodezero.and.util_print('trial vectors',print_high))
     1      write(LuOut,9200) memleft
c
c --------------
c Deallocate NGA
c --------------
c
          if (.not.tda) then
            if (.not.ga_destroy(g_nprd_m)) call errquit
     1        ('tddft_davidson: failed to deallocate g_nprd_m',0,
     &       GA_ERR)
          endif
          if (.not.ga_destroy(g_nprd_p)) call errquit
     1      ('tddft_davidson: failed to deallocate g_nprd_p',0, GA_ERR)
          if (.not.ga_destroy(g_ntrl)) call errquit
     1      ('tddft_davidson: failed to deallocate g_ntrl',0, GA_ERR)
c
c ===============================================
c Algorithms 2 & 4 (one product vector at a time)
c ===============================================
c
        else if ((algorithm.eq.2).or.(algorithm.eq.4)) then
c
c ---------------------------
c Loop over all trial vectors
c ---------------------------
c
          do n=nproducts+1,ntrials
c
c ------------------------------------------------
c MO to AO transformation of the i-th trial vector       
c ------------------------------------------------
c
            if (algorithm.eq.2) then
              call tddft_transfm(n,g_trials,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_td,1)
            else
              do i=1,ipol
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_trials(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
              enddo
              call tddft_transfm(1,g_buffer1,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_td,1)
            endif
            ccpu=ccpu+util_cpusec()
            cwall=cwall+util_wallsec()
            if (nodezero.and.util_print('timings',print_high))
     1        write(LuOut,9130) 'Transform',ccpu,cwall
            ccpu=-util_cpusec()
            cwall=-util_wallsec()
            memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
            if (nodezero.and.util_print('trial vectors',print_high))
     1        write(LuOut,9200) memleft
 9200       format(2x,'Memleft=',i10,' doubles')
c
c ---------------------------------------------------------
c Contract with integrals and make product density matrices
c ---------------------------------------------------------
c
            call tddft_contract(rtdb,geom,ao_bas_han,
     1        g_td,g_pd,g_qd,nbf_ao,ipol,tol2e,tda,oskel,kfac,
     2        lhashf,triplet)
            ccpu=ccpu+util_cpusec()
            cwall=cwall+util_wallsec()
            if (nodezero.and.util_print('timings',print_high))
     1        write(LuOut,9130) 'Contract',ccpu,cwall
            ccpu=-util_cpusec()
            cwall=-util_wallsec()
            memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
            if (nodezero.and.util_print('trial vectors',print_high))
     1        write(LuOut,9200) memleft
c
c ------------------------------------------------
c MO to AO transformation of the i-th trial vector       
c ------------------------------------------------
c
            if (algorithm.eq.2) then
              call tddft_transfm(n,g_prod_p,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_pd,2)
              if (.not.tda) then
                call tddft_transfm(n,g_prod_m,g_movecs,
     1            nbf_ao,nocc,nmo,nfc,nfv,ipol,g_qd,2)
              endif
            else
              call tddft_transfm(1,g_buffer1,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_pd,2)
              do i=1,ipol
                dummy=dra_write_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_prod_p(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
              enddo
              if (.not.tda) then
                call tddft_transfm(1,g_buffer1,g_movecs,
     1            nbf_ao,nocc,nmo,nfc,nfv,ipol,g_qd,2)
                do i=1,ipol
                  dummy=dra_write_section(.false.,
     1              g_buffer1(i),1,nov(i),1,1,
     2              d_prod_m(i),1,nov(i),n,n,request)
                  dummy=dra_wait(request)
                enddo
              endif
            endif
            ccpu=ccpu+util_cpusec()
            cwall=cwall+util_wallsec()
            if (nodezero.and.util_print('timings',print_high))
     1        write(LuOut,9130) 'Transform',ccpu,cwall
            ccpu=-util_cpusec()
            cwall=-util_wallsec()
            memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
            if (nodezero.and.util_print('trial vectors',print_high))
     1        write(LuOut,9200) memleft
c
c --------------------------------------
c Add (epsilon_a-epsilon_i) contribution
c --------------------------------------
c
            if (algorithm.eq.2) then
              call tddft_epsilon(n,g_trials,g_prod_p,
     1          ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2          lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
              if (.not.tda) then
                call tddft_epsilon(n,g_trials,g_prod_m,
     1            ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2            lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
              endif
            else
              do i=1,ipol
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_trials(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
                dummy=dra_read_section(.false.,
     1            g_buffer2(i),1,nov(i),1,1,
     2            d_prod_p(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
              enddo
              call tddft_epsilon(1,g_buffer1,g_buffer2,
     1          ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2          lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
              do i=1,ipol
                dummy=dra_write_section(.false.,
     1            g_buffer2(i),1,nov(i),1,1,
     2            d_prod_p(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
              enddo
              if (.not.tda) then
                do i=1,ipol
                  dummy=dra_read_section(.false.,
     1              g_buffer1(i),1,nov(i),1,1,
     2              d_trials(i),1,nov(i),n,n,request)
                  dummy=dra_wait(request)
                  dummy=dra_read_section(.false.,
     1              g_buffer2(i),1,nov(i),1,1,
     2              d_prod_m(i),1,nov(i),n,n,request)
                  dummy=dra_wait(request)
                enddo
                call tddft_epsilon(1,g_buffer1,g_buffer2,
     1            ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2            lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
                do i=1,ipol
                  dummy=dra_write_section(.false.,
     1              g_buffer2(i),1,nov(i),1,1,
     2              d_prod_m(i),1,nov(i),n,n,request)
                  dummy=dra_wait(request)
                enddo
              endif
            endif
            ccpu=ccpu+util_cpusec()
            cwall=cwall+util_wallsec()
            if (nodezero.and.util_print('timings',print_high))
     1        write(LuOut,9130) 'Epsilon',ccpu,cwall
            ccpu=-util_cpusec()
            cwall=-util_wallsec()
            memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
            if (nodezero.and.util_print('trial vectors',print_high))
     1        write(LuOut,9200) memleft
            if (util_print('trial vectors',print_debug))
     1        call tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     2          g_prod_p,g_prod_m,g_buffer1,d_trials,d_prod_p,
     3          d_prod_m,n,ntrials,n,tda,algorithm)
          enddo
          nproducts=ntrials
        else
          call errquit('tddft_davidson: unknown algorithm',algorithm,
     &       INPUT_ERR)
        endif
c
c ------------------------------------
c Form (A+B) and (A-B) in the subspace
c ------------------------------------
c
        if (.not.ma_push_get(mt_dbl,ntrials*ntrials,'(A+B)',
     1    l_apb,k_apb)) call errquit
     2    ('tddft_davidson: failed to allocate apb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,ntrials*ntrials,
     1    '(A+B) eigenvectors',l_apbvec,k_apbvec)) call errquit
     2    ('tddft_davidson: failed to allocate apbvec',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,ntrials,'(A+B) eigenvalues',
     1    l_apbval,k_apbval)) call errquit
     2    ('tddft_davidson: failed to allocate apbval',0, MA_ERR)
        if (.not.tda) then
          if (.not.ma_push_get(mt_dbl,ntrials*ntrials,'(A-B)',
     1      l_amb,k_amb)) call errquit
     2      ('tddft_davidson: failed to allocate amb',0, MA_ERR)
          if (.not.ma_push_get(mt_dbl,ntrials*ntrials,
     1      '(A-B) eigenvectors',l_ambvec,k_ambvec)) call errquit
     2      ('tddft_davidson: failed to allocate ambvec',0, MA_ERR)
          if (.not.ma_push_get(mt_dbl,ntrials,'(A-B) eigenvalues',
     1      l_ambval,k_ambval)) call errquit
     2      ('tddft_davidson: failed to allocate ambval',0, MA_ERR)
        endif
c
c Make (A+B) and (A-B)
c
        call tddft_subspace(ipol,nocc,nmo,nfc,nfv,g_trials,
     1    g_prod_p,g_prod_m,g_buffer1,g_buffer2,d_trials,
     2    d_prod_p,d_prod_m,ntrials,tda,dbl_mb(k_apb),dbl_mb(k_amb),
     3    algorithm)
        ccpu=ccpu+util_cpusec()
        cwall=cwall+util_wallsec()
        if (nodezero.and.util_print('timings',print_high))
     1    write(LuOut,9130) 'Subspace',ccpu,cwall
        ccpu=-util_cpusec()
        cwall=-util_wallsec()
        memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
        if (nodezero.and.util_print('trial vectors',print_high))
     1    write(LuOut,9200) memleft
c
c -----------------------------------------------------------
c Diagonalize (A-B) then (A-B)^(1/2)(A+B)(A-B)^(1/2) matrices
c -----------------------------------------------------------
c
        call tddft_diagon(ntrials,tda,dbl_mb(k_apb),dbl_mb(k_apbvec),
     1    dbl_mb(k_apbval),dbl_mb(k_amb),dbl_mb(k_ambvec),
     2    dbl_mb(k_ambval))
        ccpu=ccpu+util_cpusec()
        cwall=cwall+util_wallsec()
        if (nodezero.and.util_print('timings',print_high))
     1    write(LuOut,9130) 'Diagon',ccpu,cwall
        ccpu=-util_cpusec()
        cwall=-util_wallsec()
        memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
        if (nodezero.and.util_print('trial vectors',print_high))
     1    write(LuOut,9200) memleft
c
c ----------
c Checkpoint
c ----------
c
        if (lcheck) 
     1 call tddft_checkpoint(rtdb,ipol,ispintag,nroots,ntrials,
     2    nov,algorithm,dbl_mb(k_apbvec),dbl_mb(k_ambvec),g_trials,
     3    g_buffer1,d_trials,tda)
c
c ---------------------------------------------------------------
c Form residual vectors, new trial vectors, and check convergence
c ---------------------------------------------------------------
c 
        call tddft_residual(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1    ntrials,tda,g_trials,g_prod_p,g_prod_m,g_buffer1,g_buffer2,
     2    d_trials,d_prod_p,d_prod_m,dbl_mb(k_apbvec),
     3    dbl_mb(k_apbval),dbl_mb(k_ambvec),dbl_mb(k_residual_1),
     4    dbl_mb(k_residual_2),g_residuals_1,g_residuals_2,thresh,
     5    dbl_mb(k_valsav),maxresidual,maxomega,converged,maxvecs,
     6    evl,nconverged,nnewtrials,nresiduals,algorithm,lecut,ecut,
     7    diff_max,lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
        nproducts = ntrials ! ntrials may have changed do due restart
        ccpu=ccpu+util_cpusec()
        cwall=cwall+util_wallsec()
        if (nodezero.and.util_print('timings',print_high))
     1    write(LuOut,9130) 'Residual',ccpu,cwall
        ccpu=-util_cpusec()
        cwall=-util_wallsec()
        memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
        if (nodezero.and.util_print('trial vectors',print_high))
     1    write(LuOut,9200) memleft
        dcpu=dcpu+util_cpusec()
        dwall=dwall+util_wallsec()
        if (nodezero.and.util_print('iteration',print_default)) then
          write(LuOut,9110) iter,ntrials,nconverged,maxresidual,
     1      maxomega,dwall
          call util_flush(LuOut)
        endif
        dcpu=-util_cpusec()
        dwall=-util_wallsec()
 9110   format(i5,2x,i5,3x,i5,4x,e9.2,2x,e9.2,2x,f9.1)
c
c ---------
c Converged
c ---------
c
        if (converged) then
          if (nodezero.and.util_print('iteration',print_default)) then
            write(LuOut,9120)
          endif
 9120     format(
     1      '  ---- ',' ------ ',' ------ ',
     2      ' --------- ',' --------- ',' --------- ',/,
     3      2x,'Convergence criterion met')
 9130 format(2x,a20,1x,f9.1,' CPU seconds ',f9.1,' Wall seconds')
c
c -----------------------
c Analysis of the results
c -----------------------
c
          call tddft_analysis(rtdb,geom,ao_bas_han,nroots,nbf_ao,ipol,
     1      nocc,nmo,nfc,nfv,ntrials,tda,oskel,g_trials,g_buffer1,
     2      d_trials,dbl_mb(k_apbvec),dbl_mb(k_apbval),dbl_mb(k_ambvec),
     3      g_movecs,evl,singlet,triplet,target,targetsym,mult,
     4      algorithm,l_irs,k_irs)
          ccpu=ccpu+util_cpusec()
          cwall=cwall+util_wallsec()
          if (nodezero.and.util_print('timings',print_high))
     1      write(LuOut,9130) 'Analysis',ccpu,cwall
          ccpu=-util_cpusec()
          cwall=-util_wallsec()
          memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
          if (nodezero.and.util_print('trial vectors',print_high))
     1      write(LuOut,9200) memleft
          if (.not.tda) then
            if (.not.ma_pop_stack(l_ambval)) call errquit
     1        ('tddft_davidson: failed to deallocate ambval',0, MA_ERR)
            if (.not.ma_pop_stack(l_ambvec)) call errquit
     1        ('tddft_davidson: failed to deallocate ambvec',0, MA_ERR)
            if (.not.ma_pop_stack(l_amb)) call errquit
     1        ('tddft_davidson: failed to deallocate amb',0, MA_ERR)
          endif
          if (.not.ma_pop_stack(l_apbval)) call errquit
     1      ('tddft_davidson: failed to deallocate apbval',0, MA_ERR)
          if (.not.ma_pop_stack(l_apbvec)) call errquit
     1      ('tddft_davidson: failed to deallocate apbvec',0, MA_ERR)
          if (.not.ma_pop_stack(l_apb)) call errquit
     1      ('tddft_davidson: failed to deallocate apb',0, MA_ERR)
          goto 20
        else
          ntrials=nnewtrials
        endif
c
c Deallocate (A+B) and (A-B) matrices
c
        if (.not.tda) then
          if (.not.ma_pop_stack(l_ambval)) call errquit
     1      ('tddft_davidson: failed to deallocate ambval',0, MA_ERR)
          if (.not.ma_pop_stack(l_ambvec)) call errquit
     1      ('tddft_davidson: failed to deallocate ambvec',0, MA_ERR)
          if (.not.ma_pop_stack(l_amb)) call errquit
     1      ('tddft_davidson: failed to deallocate amb',0, MA_ERR)
        endif
        if (.not.ma_pop_stack(l_apbval)) call errquit
     1    ('tddft_davidson: failed to deallocate apbval',0, MA_ERR)
        if (.not.ma_pop_stack(l_apbvec)) call errquit
     1    ('tddft_davidson: failed to deallocate apbvec',0, MA_ERR)
        if (.not.ma_pop_stack(l_apb)) call errquit
     1    ('tddft_davidson: failed to deallocate apb',0, MA_ERR)
c
c End of Davidson iterations
c
      goto 10
c           
c ------
c Return
c ------
c
   20 continue
      if (.not.tda) then
        if (.not.ma_pop_stack(l_residual_2)) call errquit
     1    ('tddft_davidson: failed to deallocate residual_2',0, MA_ERR)
      endif
      if (.not.ma_pop_stack(l_residual_1)) call errquit
     1  ('tddft_davidson: failed to deallocate residual_1',0, MA_ERR)
      if (.not.ma_pop_stack(l_valsav)) call errquit
     1  ('tddft_davidson: failed to deallocate valsav',0, MA_ERR)
      do i=ipol,1,-1
        if (.not.tda) then
          if (.not.ga_destroy(g_residuals_2(i))) call errquit
     1      ('tddft_davidson: failed to destroy g_residuals_2',0,
     &       GA_ERR)
        endif
        if (.not.ga_destroy(g_residuals_1(i))) call errquit
     1    ('tddft_davidson: failed to destroy g_residuals_1',0, GA_ERR)
        if (.not.tda) then
          if (.not.ga_destroy(g_qd(i))) call errquit
     1      ('tddft_davidson: failed to destroy g_qd',0, GA_ERR)
        endif
        if (.not.ga_destroy(g_pd(i))) call errquit
     1    ('tddft_davidson: failed to destroy g_pd',0, GA_ERR)
        if (.not.ga_destroy(g_td(i))) call errquit
     1    ('tddft_davidson: failed to destroy g_td',0, GA_ERR)
      enddo
      return
      end
