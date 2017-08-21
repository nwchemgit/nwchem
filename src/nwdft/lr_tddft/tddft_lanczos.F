      subroutine tddft_lanczos(rtdb,ipolar)
c
c  Lanczos symmetric solver
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
#include "xc.fh"
#include "cosmo.fh"
#include "prop.fh"
c
      integer rtdb             ! Runtime database handle
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      logical oskel,oadapt     ! Symmetry control parameters
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      character*255 movecs_out ! MO vector file name
      character*255 civecs_tmp ! CI vector file name
      character*40 xname(numfunc)       ! Names of X functionals
      character*40 cname(numfunc)       ! Names of C functionals
      character*40 xcname(numfunc)      ! Names of XC functionals
      double precision xfac(numfunc)    ! Prefactors for X functionals
      double precision cfac(numfunc)    ! Prefactors for C functionals
      logical lxfac(numfunc)            ! True if local X functional
      logical nlxfac(numfunc)           ! True if gc X functional
      logical lcfac(numfunc)            ! True if local C functional
      logical nlcfac(numfunc)           ! True if gc C functional
      logical xccomb(numfunc)           ! True if composite XC functional
      logical lhf                       ! True if HF (no DFT)
      logical ldft                      ! True if pure DFT (no HF exch)
      logical lhybrid                   ! True if hybrid DFT and HF
      logical lhasdft                   ! True if it has DFT
      logical lhashf                    ! True if it has HF
      logical tda                       ! True if Tamm-Dancoff approximation
      logical singlet                   ! True if singlet excited states
      logical triplet                   ! True if triplet excited states
      integer nroots           ! Number of roots sought
      integer maxvecs          ! Maximum number of trial vectors
      integer maxiter          ! Maximum number of Davidson iterations
      integer target           ! Target excited state
      character*4 targetsym    ! Symmetry of target excited state
      logical symmetry         ! Symmetry restriction
      double precision tol2e   ! 2-e integral tolerance
      integer nov(2)           ! Number of occupied virtual pairs
      integer g_trials(2)      ! GA handle for trial vectors (XY)
      integer g_prod_p(2)      ! GA handle for A+B product vectors (XY)
      integer g_prod_m(2)      ! GA handle for A-B product vectors (XY)
      integer g_movecs(2)      ! GA handle for MO vectors
      integer l_evl,k_evl      ! MA handles for eigenvalues (AB)
      integer l_occ,k_occ      ! MA handles for occupancies (AB)
      integer l_irs(2),k_irs(2)! MA handles for irreps (AB)
      integer ntrials          ! Current number of trial vectors
      integer nproducts        ! Current number of product vectors
      double precision thresh  ! Convergence threshold for Davidson iteration
      double precision energy  ! Target excited state energy
      logical xc_active_save   ! Previous value for xc_active
      integer calc_type_save   ! Previous value for calc_type
      integer algorithm        ! =1 incore-multiple fock build; =2: incore-single fock build
                               ! =3 disk-multiple fock build;   =4:
                               ! disk-single fock build
      integer mult             ! Ground state spin multiplicity
c
      logical dft_energy
      external dft_energy
c
      integer i
      integer memleft
      integer dummy
      logical nodezero
c
      integer vshift
      logical lecut
      double precision ecut
      logical lewin
      double precision ewinl,ewinh
      double precision diff_max
      logical lowin(2)
      integer owstart(2)
      integer owend(2)


cjb  Lanczos specific variables
      integer l_t(2),k_t(2),j
      integer mdeg,ideg,l_h,k_h
      double precision beta
      integer g_bigv(2),g_bigvp(2),g_sh,g_sh2
      integer k, l_hv, k_hv
      integer g_h,g_hvec,g_wgt2
      integer l_gamma2,k_gamma2
      logical llanczos
      integer nvec,isave,s
      integer l_ritz,k_ritz

      double precision ccpu,cwall    ! CPU timings
      double precision ecpu,ewall    ! CPU timings

      double precision cntr(3) ! Center of dipole (arbitrary)
      double precision tmom(4) ! Transition moments
      integer len_moments,g_dens(2),iv
      logical llandip
      integer ipolar, l, unitn, reclen
      integer g_s
      double precision dlanshift, dtheta
      logical ldolshift,lpolar
c
      nodezero=(ga_nodeid().eq.0)

        do i=1,3
           cntr(i)=0.0d0
        enddo
c
c ---------------------
c Initialize parameters
c ---------------------
c
      call tddft_init(rtdb,geom,oskel,oadapt,ao_bas_han,
     1  nbf_ao,movecs_out,ipol,nocc,nmo,nfc,nfv,numfunc,
     2  xname,cname,xcname,xfac,cfac,lxfac,nlxfac,lcfac,nlcfac,xccomb,
     3  lhf,ldft,lhybrid,lhasdft,lhashf,
     4  tda,singlet,triplet,nroots,maxvecs,tol2e,thresh,
     5  maxiter,target,targetsym,symmetry,algorithm,g_movecs,
     6  l_evl,k_evl,l_occ,k_occ,l_irs,k_irs,mult,lecut,ecut,diff_max,
     7  lowin,owstart,owend,lewin,ewinl,ewinh,vshift)

      if (lhasdft) then
        if (.not.rtdb_get(rtdb,'bgj:xc_active',mt_log,1,
     1    xc_active_save)) xc_active_save=.false.
        if (.not.rtdb_get(rtdb,'fock_xc:calc_type',mt_int,1,
     1    calc_type_save)) calc_type_save=0
        if (.not.rtdb_put(rtdb,'bgj:xc_active',mt_log,1,.true.))
     1    call errquit('tddft_energy: failed to set xc_active true',0,
     &       RTDB_ERR)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,2))
     1    call errquit('tddft_energy: failed to set calc_type 2',0,
     &       RTDB_ERR)
      endif

c
c ---------------------
c Initialize parameters - Lanczos
c ---------------------
c
       if (.not.rtdb_get(rtdb,'tddft:mdeg',mt_int,1,mdeg))
     1 call errquit('tddft_lanczos: lsteps is a mandatory parameter',0,
     & GA_ERR)

        lpolar = .true.
        if (.not.rtdb_get(rtdb,'tddft:ipolar',mt_int,1,ipolar))
     1 lpolar = .false.

       if (.not.rtdb_get(rtdb,'tddft:lanshift',mt_dbl,1,dlanshift)) then
         dlanshift = 0.d0
         ldolshift = .false.
       else
         ldolshift = .true.
       endif

       dlanshift = dlanshift**2;

c
c Initial check
c
       if(algorithm.gt.1) then
       call errquit
     1('tddft: Only algorithm 1 is compatible with Lanczos solver',0,
     & GA_ERR)
        endif

       if(ipol.gt.1) then
c        singlet = .false.
c       call errquit
c     1('Lanczos solver implemented for closed shell only',0,
c     & GA_ERR)
       endif


c
c ----------------------------------
c Cosmo initialization, if necessary
c ----------------------------------
c
      cosmo_on = .false.
c
c ----------------------------------------------------
c Allocate global arrays for trial and product vectors
c ----------------------------------------------------
c
        call tddft_lentrialvec(ipol,nocc,nmo,nfc,nfv,nov)
        nproducts = 0

        maxvecs = 1
c       if(nodezero)write(6,*)'CHCK 0'

        do i=1,ipol
          if (.not.ga_create(mt_dbl,max(1,nov(i)),maxvecs,
     1      'trial vectors',-1,maxvecs,g_trials(i)))
     2      call errquit('tddft_energy: failed to create g_trials',0,
     &       GA_ERR)

          if (.not.ga_duplicate(g_trials(i),g_prod_p(i),
     1      'product vectors'))
     2      call errquit('tddft_energy: failed to create g_prod_p',0,
     &       GA_ERR)

          if (.not.tda) then
            if (.not.ga_duplicate(g_trials(i),g_prod_m(i),
     1        'product vectors'))
     2        call errquit('tddft_energy: failed to create g_prod_m',0,
     &       GA_ERR)
          endif
        enddo

       
        do i=1,ipol
          call ga_zero(g_trials(i))
          call ga_zero(g_prod_p(i))
          call ga_zero(g_prod_m(i))
        enddo
c
c -----------------------------------------------------
c Iterative solution of TDDFT matrix eigenvalue problem
c -----------------------------------------------------
c

          if(nodezero)write(LuOut,
     1"(/,'Starting Lanczos for estimating absorption spectra')")

          if(nodezero) then
            write(LuOut,"(/,'Number of Lanczos steps:',I7,/)")mdeg
            write(LuOut,"(/,'Polarization:',I7,/)")ipolar
          endif

      ccpu=-util_cpusec()
      cwall=-util_wallsec()

        if (.not.ma_push_get(mt_dbl,mdeg,'ritz',
     1    l_ritz,k_ritz)) call errquit('tddft_iguess: hv',0, MA_ERR)

        do j =1,mdeg
          dbl_mb(k_ritz+j-1) = 0.0d0
        enddo
c
c INITIATE
c
c       if(nodezero)write(6,*)'CHCK 1'

       do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nov(i),'Work trial vector',
     1    l_t(i),k_t(i))) call errquit
     2    ('tddft_iguess: failed to allocate t',0, MA_ERR)
       enddo

c       if(nodezero)write(6,*)'CHCK 2'

        if (.not.ma_push_get(mt_dbl,mdeg*mdeg,'H matrix',
     1    l_h,k_h)) call errquit
     2    ('tddft_iguess: failed to allocate t',0, MA_ERR)

c       if(nodezero)write(6,*)'CHCK 3'

c
c Lanczos vectors
c
        do i=1,ipol
          if (.not.ga_create(mt_dbl,nov(i),mdeg,
     1      'big V',-1,-1,g_bigv(i)))
     2      call errquit('tddft_energy: failed to create g_trials',0,
     &       GA_ERR)
          if (.not.ga_create(mt_dbl,nov(i),mdeg,
     1      'big prod V',-1,-1,g_bigvp(i)))
     2      call errquit('tddft_energy: failed to create g_trials',0,
     &       GA_ERR)

          call ga_zero(g_bigv(i))
          call ga_zero(g_bigvp(i))

        enddo

c       if(nodezero)write(6,*)'CHCK 4'


        do j=1,mdeg*mdeg
          dbl_mb(k_h+j-1) = 0.0d0
        enddo
        do i=1,ipol
        do j=1,nov(i)
          dbl_mb(k_t(i)+j-1) = 0.0d0
        enddo
        enddo

c       if(nodezero)write(6,*)'CHCK 5'

c
c Get <|mu|> integrals as a trial vectors
c
       do i=1,ipol

         if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'AO ovl',
     &     -1, -1, g_dens(i)))
     &     call errquit('Error creating ga',0,GA_ERR)
         call ga_zero(g_dens(i))

c       if(nodezero)write(6,*)'CHCK 6'

      call tddft_dipole(ao_bas_han, cntr, 1,
     &     g_dens(i), tmom, 4, k_evl, ipolar, g_movecs(i))

c       if(nodezero)write(6,*)'CHCK 6a'

           iv = 1
           call ga_zero(g_trials(i))

       enddo

      call ga_sync()

c       if(nodezero)write(6,*)'CHCK 6b',singlet,triplet

      call tddft_transfm(iv,g_trials,g_movecs,
     1  nbf_ao,nocc,nmo,nfc,nfv,ipol,g_dens,2)

c       if(nodezero)write(6,*)'CHCK 6c',ipol

      do i=1,ipol

      call ga_get(g_trials(i),1,nov(i),1,1,dbl_mb(k_t(i)),nov(i))

c       if(nodezero)write(6,*)'CHCK 6d WWWWWWWWWWWWWWWWWW Warning'

          if (.not.ga_destroy(g_dens(i))) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)

      enddo

      call ga_sync()

c       if(nodezero)write(6,*)'CHCK 7'


c********************************
c Perform Lanczos steps
c********************************

       do ideg = 1,mdeg

       if(nodezero)write(6,"('Working on step = ',I6,'/',I6)")ideg,mdeg

       ecpu=-util_cpusec()
       ewall=-util_wallsec()

       call tddft_prod_vec(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1  nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     3  g_movecs,dbl_mb(k_evl),1,nproducts,
     4  tda,oskel,singlet,triplet,tol2e,xfac(1),lhashf,thresh,
     5  target,targetsym,algorithm,mult,l_irs,k_irs,1)
c
c keep result in memory
c
        beta = 0.0d0

        do i=1,ipol
          beta = beta + ga_ddot(g_trials(i),g_prod_p(i))
        enddo

        beta = dsqrt(beta)

       if(ideg.gt.1) then
        dbl_mb(k_h+(ideg-1)*mdeg+(ideg-1)-1) = beta
       endif

        do i=1,ipol

          call ga_scale(g_trials(i),1.0d0/beta)

        call ga_copy_patch('N',g_trials(i),1,nov(i),1,1,
     1 g_bigv(i),1,nov(i),ideg,ideg)

         call ga_scale(g_prod_p(i),1.0d0/beta)

        call ga_copy_patch('N',g_prod_p(i),1,nov(i),1,1,
     1 g_bigvp(i),1,nov(i),ideg,ideg)

        enddo

              call tddft_prod_vec(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1  nocc,nmo,nfc,nfv,nroots,maxvecs,g_prod_p,g_prod_p,g_prod_m,
     3  g_movecs,dbl_mb(k_evl),1,nproducts,
     4  tda,oskel,singlet,triplet,tol2e,xfac(1),lhashf,thresh,
     5  target,targetsym,algorithm,mult,l_irs,k_irs,2)

cjb spectral transform

      if(ldolshift) then

       do i = 1,ipol

         call ga_scale(g_trials(i),-dlanshift)
         call ga_add(1.0d0,g_trials(i),1.0d0,g_prod_m(i),g_prod_m(i))
         call ga_copy(g_prod_m(i),g_trials(i))

       enddo

       call tddft_prod_vec(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1  nocc,nmo,nfc,nfv,nroots,maxvecs,g_prod_m,g_prod_p,g_prod_m,
     3  g_movecs,dbl_mb(k_evl),1,nproducts,
     4  tda,oskel,singlet,triplet,tol2e,xfac(1),lhashf,thresh,
     5  target,targetsym,algorithm,mult,l_irs,k_irs,1)

       call tddft_prod_vec(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1  nocc,nmo,nfc,nfv,nroots,maxvecs,g_prod_p,g_prod_p,g_prod_m,
     3  g_movecs,dbl_mb(k_evl),1,nproducts,
     4  tda,oskel,singlet,triplet,tol2e,xfac(1),lhashf,thresh,
     5  target,targetsym,algorithm,mult,l_irs,k_irs,2)

       do i = 1,ipol

         call ga_scale(g_trials(i),-dlanshift)
         call ga_add(1.0d0,g_trials(i),1.0d0,g_prod_m(i),g_prod_m(i))

       enddo

      endif

cjb spectral end

       if (.not.ga_create(mt_dbl,ideg,1,'small h',-1,-1,g_sh))
     2   call errquit('tddft_energy: failed',0,GA_ERR)
       call ga_zero(g_sh)

c h    = V'*(B*Av);

       do i=1,ipol

         call ga_dgemm('T','N',ideg,1,nov(i),1.0d0,g_bigvp(i),
     1 g_prod_m(i),1.0d0,g_sh)

       enddo
c V*h
       do i=1,ipol

         call ga_dgemm('N','N',nov(i),1,ideg,1.0d0,g_bigv(i),g_sh,
     1 0.0d0,g_trials(i))

c f    = Av - V*h;
        call ga_add(1.0d0,g_prod_m(i),-1.0d0,g_trials(i),g_trials(i))

       enddo

c =============== Reorthogonalization begin
c    s    = V'*(B*f);
c    h    = h + s;
c    f    = f - V*s;
c       if(0.eq.1) then

       if (.not.ga_create(mt_dbl,ideg,1,'small h',-1,-1,g_s))
     2   call errquit('tddft_energy: failed',0,GA_ERR)
         call ga_zero(g_s)
c B*f
c s = V'*(B*f)
       do i=1,ipol

         call ga_dgemm('T','N',ideg,1,nov(i),1.0d0,g_bigvp(i),
     1 g_trials(i),1.0d0,g_s)

       enddo

c h = h + s
        call ga_add(1.0d0,g_sh,1.0d0,g_s,g_sh)

c V*s;
        do i=1,ipol

         call ga_dgemm('N','N',nov(i),1,ideg,1.0d0,g_bigv(i),g_s,
     1 0.0d0,g_prod_m(i))

c f = f - V*s
        call ga_add(1.0d0,g_trials(i),-1.0d0,g_prod_m(i),g_trials(i))

        enddo

          if (.not.ga_destroy(g_s)) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)

c ======================= Reortho end

        if (.not.ma_push_get(mt_dbl,mdeg*mdeg,'h vec',
     1    l_hv,k_hv)) call errquit('tddft_iguess: hv',0, MA_ERR)

         call ga_get(g_sh,1,ideg,1,1,dbl_mb(k_hv),ideg)

c H(1:j,j) = h;
         do j = 1,ideg
           dbl_mb(k_h+ideg-1+(j-1)*mdeg) = dbl_mb(k_hv+j-1)
         enddo

        if(ideg.eq.1) then
         if(nodezero) write(6,"('New H elements:',2F18.12)")
     1 dbl_mb(k_hv), 0.0d0
        else
         if(nodezero) write(6,"('New H elements:',2F18.12)")
     1 dbl_mb(k_hv+ideg-1),
     1 dbl_mb(k_hv+ideg-2)
        endif

        if (.not.ma_pop_stack(l_hv)) call errquit
     1    ('tddft_iguess: failed to deallocate hv',0, MA_ERR)

          if (.not.ga_destroy(g_sh)) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)


          call ga_sync()

          ecpu=ecpu+util_cpusec()
          ewall=ewall+util_wallsec()

          if(nodezero)write(LuOut,"(2x,a20,1x,f9.1,' CPU seconds ',
     1f9.1,' Wall seconds')") 'Step time: ',ecpu,ewall

        enddo !ideg

c ***************ideg END ***************
c
c   Print H
c
       if(nodezero) then
         write(6,"(/,'Printing tridiagonal H matrix')")
         write(6,"('=============================')")
           do j=1,mdeg
              if(j.ne.mdeg) then
               write(6,"('H element (i,i) (i,i+1)',I6,2F22.12)")
     1 j,dbl_mb(k_h+(j-1)*mdeg+j-1),dbl_mb(k_h+(j-1)*mdeg+j)
              else
               write(6,"('H element (i,i)        ',I6,F22.12)")
     1 j,dbl_mb(k_h+(j-1)*mdeg+j-1)
              endif
           enddo
       endif
c
c   Diagonalize H
c
       if (.not.ga_create(mt_dbl,mdeg,mdeg,'Global H',-1,-1,g_h))
     2   call errquit('tddft_energy: failed',0,GA_ERR)

       if (.not.ga_create(mt_dbl,mdeg,mdeg,'H evec',-1,-1,g_hvec))
     2   call errquit('tddft_energy: failed',0,GA_ERR)

       if(nodezero) then
         call ga_put(g_h,1,mdeg,1,mdeg,dbl_mb(k_h),mdeg)
       endif

       call ga_sync()

         do j=1,mdeg*mdeg
           dbl_mb(k_h+j-1) = 0.0d0
         enddo

         call ga_diag_std(g_h,g_hvec,dbl_mb(k_h))

cccc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CJB Analyse
cccc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       if(nodezero) then
        do i=1,ipol
         call ga_put(g_trials(i),1,nov(i),1,1,dbl_mb(k_t(i)),nov(i))
        enddo
       endif

         call ga_sync()

       if (.not.ga_create(mt_dbl,mdeg,1,'weig h',-1,-1,g_sh))
     2   call errquit('tddft_energy: failed',0,GA_ERR)
c2       call ga_zero(g_sh)

       if (.not.ga_create(mt_dbl,mdeg,1,'weig h',-1,-1,g_sh2))
     2   call errquit('tddft_energy: failed',0,GA_ERR)

         call ga_zero(g_sh2)
         call ga_zero(g_sh)

         do i=1,ipol
          call ga_dgemm('T','N',mdeg,1,nov(i),1.0d0,g_bigv(i),
     1 g_trials(i),1.0d0,g_sh)
         enddo

c         call ga_sync()

c wgt1 = eigvec'*(V'*v0);

         call ga_dgemm('T','N',mdeg,1,mdeg,1.0d0,g_hvec,g_sh,
     1 0.0d0,g_sh2)

c2         call ga_sync()

         call ga_zero(g_sh)

         call ga_elem_multiply(g_sh2,g_sh2,g_sh)

cccc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CJB DOS
cccc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if (.not.ma_push_get(mt_dbl,mdeg,'gamma2',
     1    l_gamma2,k_gamma2)) call errquit('tddft',0, MA_ERR)

        call ga_sync()

        call ga_get(g_sh,1,mdeg,1,1,dbl_mb(k_gamma2),mdeg)

        if(nodezero) then
            write(6,"(/,'Ritz values:',
     1/,'==============E[a.u.]==========E[eV]=============Osc.s.===='
     2 )")
          do j=1,mdeg

            dtheta = dsqrt(dbl_mb(k_h+j-1))

           if(ldolshift) then
            if(dtheta.lt.dlanshift) then
              dtheta = dsqrt(dabs(dtheta-dlanshift)) 
            else
              dtheta = dsqrt(dtheta + dlanshift)
            endif
           endif

            write(6,"('Ritz v.: ',3F18.12)")
     1 dtheta,
     1 dtheta*27.2113851,
     1 dbl_mb(k_gamma2+j-1)*dtheta*
     1 dtheta*(4.0d0/3.0d0)/(DBLE(ipol))
          enddo
        endif

        if (.not.ma_pop_stack(l_gamma2)) call errquit
     1    ('tddft_iguess: failed to deallocate hv',0, MA_ERR)

cjb end DOS
          if (.not.ga_destroy(g_sh)) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)
          if (.not.ga_destroy(g_sh2)) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)
cjb end analyse
       if (.not.ga_destroy(g_hvec)) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)
       if (.not.ga_destroy(g_h)) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)

cccc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (.not.ma_pop_stack(l_h)) call errquit
     1    ('tddft_iguess: failed to deallocate t',0, MA_ERR)

       do i=1,ipol
        if (.not.ma_pop_stack(l_t(ipol+1-i))) call errquit
     1    ('tddft_iguess: failed to deallocate t',0, MA_ERR)
       enddo

       do i=1,ipol
          if (.not.ga_destroy(g_bigv(i))) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)

          if (.not.ga_destroy(g_bigvp(i))) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)
       enddo  
        
c =================================================================
c                Output
c =================================================================


c
c DEALLOC
c
        if (.not.ma_pop_stack(l_ritz)) call errquit
     1    ('tddft_iguess: failed to deallocate t',0, MA_ERR)

       if(nodezero) write(LuOut
     1 ,"('=======================================================')")
       if(nodezero)write(LuOut,"(/,'End of Lanczos subroutine')")

cjb  reinitialize Global Arrays

       do i=1,ipol

       if (.not.ga_destroy(g_trials(i))) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)
       if (.not.ga_destroy(g_prod_m(i))) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)
       if (.not.ga_destroy(g_prod_p(i))) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)

       enddo
c
c And we close the 2-e integral evaluators
c
      call fock_2e_tidy(rtdb)
      call schwarz_tidy()
      call int_terminate()
c
c -----------------
c Print utility end
c -----------------
c
      if(lpolar.or.(ipolar.eq.3))call util_print_pop
c
c ------
c Return
c ------
c
      if (lhasdft) then
        if (.not.rtdb_put(rtdb,'bgj:xc_active',mt_log,1,
     1    xc_active_save))
     2    call errquit('tddft_energy: failed to restore xc_active',0,
     &       RTDB_ERR)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,
     1    calc_type_save))
     2    call errquit('tddft_energy: failed to restore calc_type',0,
     &       RTDB_ERR)
      endif
      if (.not.rtdb_put(rtdb,'tddft:energy',mt_dbl,1,0.0d0))
     1  call errquit('tddft_lanczos: failed to put tddft energy',0,
     &       RTDB_ERR)
      if(.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1  call errquit('tddft_energy: failed to unset triplet',0,
     &       RTDB_ERR)
      call tddft_finalize(geom,ao_bas_han,ipol,g_movecs,l_evl,k_evl,
     1  l_occ,k_occ,l_irs,k_irs)

c       stop
      return
      end
