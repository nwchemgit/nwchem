      subroutine tddft_prod_vec(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1  nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     3  g_movecs,evl,ntrials,nproducts,
     4  tda,oskel,singlet,triplet,tol2e,kfac,lhashf,thresh,
     5  target,targetsym,algorithm,mult,l_irs,k_irs,iwhich,lhascd)
c
c
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
      logical lhascd
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
      logical lecut
      double precision ecut
      double precision diff_max
      logical lewin
      double precision ewinl,ewinh,etmp
      logical lowin
      integer owstart(2)
      integer owend(2)
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
cjb
c      double precision, allocatable, dimension(:,:) :: eye
      integer nab,iblk, ierr, nblks, ie, ntrials_last
      integer mylo(2), myhi(2)
      logical llast
      character*255 filename
      integer unitn,l_aux,k_aux
      integer l_t,k_t,j,oblk
      integer iwhich
c
      nodezero=(ga_nodeid().eq.0)
      ccpu=0d0
      cwall=0d0
cjb
c      call ga_print(g_movecs(1))
c
c -------------------------
c Allocate density matrices
c -------------------------
c
      do i=1,ipol
        nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
      enddo
      if (ipol.eq.1) nov(2)=0
cjb
      nab = nov(1)
c        if (.not.ma_push_get(mt_dbl,nov(1),'Work trial vector',
c    1    l_t,k_t)) call errquit
c     2    ('tddft_iguess: failed to allocate t',0, MA_ERR)     
cjb
c
c     nresiduals=min(nroots+min(nroots/5,10),maxvecs)
      nresiduals=nroots
      do i=1,ipol
        g_td(i)=ga_create_atom_blocked
     1    (geom,ao_bas_han,'Trial density')

      if(iwhich.eq.1) then
        g_pd(i)=ga_create_atom_blocked
     1    (geom,ao_bas_han,'A product density')
      else
        if (.not.tda) 
     1    g_qd(i)=ga_create_atom_blocked
     2      (geom,ao_bas_han,'B product density')
      endif

      enddo
c
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
c ================================================
c Algorithms 1 & 3 (multiple product vector build)
c ================================================
c
cc        call ga_zero(g_trials(1))

cc              mylo(1) = 1
cc              mylo(2) = 1
cc              myhi(1) = nab
cc              myhi(2) = 1

cc        do i=1,nab
cc          dbl_mb(k_t+i-1) = 0.0d0
cc        enddo
cc          dbl_mb(k_t) = 1.0d0
cc
cc        call nga_put(g_trials(1), mylo, myhi,
cc     &           dbl_mb(k_t), nab)

c
c -----------------------------------------------
c Allocate NGA trial and product density matrices
c -----------------------------------------------
c
c          write(6,*)'Check 2',ntrials,nproducts
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

         if(iwhich.eq.1) then

          if (.not.nga_create(mt_dbl,3,idim,'Work product',
     1      ichunk,g_nprd_p)) call errquit
     2      ('tddft_davidson: failed to create g_nprd_p',0, GA_ERR)
         else

          if (.not.tda) then
            if (.not.nga_create(mt_dbl,3,idim,'Work product',
     1        ichunk,g_nprd_m)) call errquit
     2        ('tddft_davidson: failed to create g_nprd_m',0, GA_ERR)
         endif
          endif
c          write(6,*)'Check 3'

c
c -----------------------
c MO to AO transformation
c -----------------------
c
cc          if(nodezero)write(6,*)'Before MO->AO'

          do n=nproducts+1,ntrials
              call tddft_transfm(n,g_trials,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_td,1)
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
c
c ---------------------------------------------------------
c Contract with integrals and make product density matrices
c ---------------------------------------------------------
c
c          if(nodezero)write(6,*)'Before contractions'
          ccpu=-util_cpusec()
          cwall=-util_wallsec()

          call tddft_nga_cont(rtdb,geom,ao_bas_han,
     1      g_ntrl,g_nprd_p,g_nprd_m,nbf_ao,ipol,tol2e,tda,oskel,kfac,
     2      lhashf,triplet,nvectors,iwhich,lhascd)
 
          ccpu=ccpu+util_cpusec()
          cwall=cwall+util_wallsec()
cjb         if (nodezero.and.util_print('timings',print_high))
c           if(nodezero)write(LuOut,9130) 'Contract',ccpu,cwall
          ccpu=-util_cpusec()
          cwall=-util_wallsec()
          memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
cjb          if (nodezero.and.util_print('trial vectors',print_high))
cjb     1      write(LuOut,9200) memleft
c
c -----------------------
c AO to MO transformation
c -----------------------
c
cc          if(nodezero)write(6,*)'BEfore AO->MO'

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

           if(iwhich.eq.1) then
              call nga_copy_patch('N',g_nprd_p,blo,bhi,g_pd(i),alo,ahi)

           else 
              if (.not.tda) call nga_copy_patch
     1          ('N',g_nprd_m,blo,bhi,g_qd(i),alo,ahi)
           endif

            enddo
          
           if(iwhich.eq.1) then
 
              call tddft_transfm(n,g_prod_p,g_movecs,
     1          nbf_ao,nocc,nmo,nfc,nfv,ipol,g_pd,2)
            else
              if (.not.tda) then
                call tddft_transfm(n,g_prod_m,g_movecs,
     1            nbf_ao,nocc,nmo,nfc,nfv,ipol,g_qd,2)
              endif
            endif
          enddo
          ccpu=-util_cpusec()
          cwall=-util_wallsec()

c --------------------------------------
c Add (epsilon_a-epsilon_i) contribution
c --------------------------------------
c
cc          if(nodezero)write(6,*)'Before epsilon'

          do n=nproducts+1,ntrials
           if(iwhich.eq.1) then
              call tddft_epsilon(n,g_trials,g_prod_p,
     1          ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2          lowin,owstart,owend,lewin,ewinl,ewinh,0)
           else
              if (.not.tda) then
                call tddft_epsilon(n,g_trials,g_prod_m,
     1            ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2            lowin,owstart,owend,lewin,ewinl,ewinh,0)
              endif
           endif
          enddo
          ccpu=-util_cpusec()
          cwall=-util_wallsec()

c          if(iwhich.eq.1) call ga_print(g_prod_p(1))
c          if(iwhich.eq.2) call ga_print(g_prod_m(1))

c          if(iwhich.eq.1) then
c            call ga_copy(g_prod_p(1),g_trials(1))
c          else
c            call ga_copy(g_prod_m(1),g_trials(1))
c          endif
c
c --------------
c Deallocate NGA
c --------------
c

         if(iwhich.eq.2) then
          if (.not.tda) then
            if (.not.ga_destroy(g_nprd_m)) call errquit
     1        ('tddft_davidson: failed to deallocate g_nprd_m',0,
     &       GA_ERR)
          endif
         endif

           if(iwhich.eq.1) then
          if (.not.ga_destroy(g_nprd_p)) call errquit
     1      ('tddft_davidson: failed to deallocate g_nprd_p',0, GA_ERR)
           endif
          if (.not.ga_destroy(g_ntrl)) call errquit
     1      ('tddft_davidson: failed to deallocate g_ntrl',0, GA_ERR)

cc        if(nodezero)write(LuOut,"('Finished')")

 9130 format(2x,a20,1x,f9.1,' CPU seconds ',f9.1,' Wall seconds')
      do i=ipol,1,-1
        if(iwhich.eq.2) then
        if (.not.tda) then
          if (.not.ga_destroy(g_qd(i))) call errquit
     1      ('tddft_davidson: failed to destroy g_qd',0, GA_ERR)
        endif
        endif

        if(iwhich.eq.1) then
        if (.not.ga_destroy(g_pd(i))) call errquit
     1    ('tddft_davidson: failed to destroy g_pd',0, GA_ERR)
        endif

        if (.not.ga_destroy(g_td(i))) call errquit
     1    ('tddft_davidson: failed to destroy g_td',0, GA_ERR)
      enddo
cjb
c        if (.not.ma_pop_stack(l_t)) call errquit
c     1    ('tddft_iguess: failed to deallocate t',0, MA_ERR)

c      stop
      return
      end
