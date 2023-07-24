      subroutine bse_davidson(pars)

      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bse.fh"
#include "cdft.fh"
      type(bse_params_t) :: pars
      character(*), parameter :: pname = 'bse_analytic: '

      integer lSize,isp
      integer l_omega, k_omega
      integer l_dia,k_dia
      integer lW,kW
      integer maxpoles,maxnpoles,totnpoles,nri

      integer indeces(pars%nroots),ntrials,g_trials,npoles,ndim
      integer g_x,g_y,g_ym,g_s,g_ks,g_mks,g_w,g_kw,g_mkw,off1
      integer g_wia,g_tmp,info,lwork,kwork,idum,mdim
      integer k_x,k_y,k_ym,k_s,k_ks,k_mks,k_w,k_kw,k_mkw,k_r,l_r
      integer myld,vkv,l_vkv,myhi,mylo,k_trials,l_ksw,ksw,jdum
      double precision minimum,norm(pars%nroots),rdum,factor
      double precision,allocatable :: diagonal(:),diagonal2(:)

      integer ipole,jpole,itrial,imo,amo,idavidson
      logical found,dohartree,stat,debug
      double precision,external :: dasum

      debug = pars%me.eq.0
      debug = .false.

      ! synchronize parallel processes
      call ga_sync()

      !initialize useful variables
      nri = pars%nri
      maxpoles = maxval(pars%mynpoles(1:pars%ipol))
      maxnpoles = maxval(pars%npoles(1:pars%ipol))
      totnpoles = maxnpoles*pars%ipol
      npoles = sum(pars%npoles(1:pars%ipol))
      ntrials = pars%nroots


      allocate(diagonal(npoles),diagonal2(npoles))

      ! allocation
      if(.not.ma_push_get(mt_dbl,totnpoles,'dia',l_dia,k_dia))
     &  call errquit(pname//'failed to allocate Delta_ia',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,npoles,'omega',l_omega,k_omega))
     &  call errquit(pname//'failed to allocate omega',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,npoles,'residual',l_r,k_r))
     &  call errquit(pname//'failed to allocate residual',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,ntrials*ntrials,'vkv',l_vkv,vkv))
     &  call errquit(pname//'failed to allocate vkv',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nri**2,'W',lW,kW))
     &  call errquit(pname//'failed to allocate W',0,MA_ERR)

      stat = ga_create(mt_dbl,npoles,ntrials,'x',0,ntrials,g_x).and.
     $  ga_create(mt_dbl,npoles,ntrials,'y',0,ntrials,g_y).and.
     $  ga_create(mt_dbl,npoles,ntrials,'ym',0,ntrials,g_ym).and.
     $  ga_create(mt_dbl,npoles,ntrials,'s',0,ntrials,g_s).and.
     $  ga_create(mt_dbl,npoles,ntrials,'ks',0,ntrials,g_ks).and.
     $  ga_create(mt_dbl,npoles,ntrials,'mks',0,ntrials,g_mks).and.
     $  ga_create(mt_dbl,npoles,ntrials,'w',0,ntrials,g_w).and.
     $  ga_create(mt_dbl,npoles,ntrials,'kw',0,ntrials,g_kw).and.
     $  ga_create(mt_dbl,npoles,ntrials,'mkw',0,ntrials,g_mkw).and.
     $  ga_create(mt_dbl,npoles,ntrials,'t',0,ntrials,g_trials).and.
     $  nga_create(mt_dbl,1,npoles,'wia',0,g_wia)

      if(.not.stat) call errquit('ga create failed',66,GA_ERR)

      ! obtain "old" eigenvalue difference to build W
      do isp=1,pars%ipol
         lSize = (isp-1)*maxpoles
         call bse_get_eia(dbl_mb(pars%k_mf_evals+(isp-1)*nbf_ao),
     &                   dbl_mb(k_dia+lSize),pars%nocc(isp),
     &                   pars%nvir(isp),pars%ovlo(isp),pars%ovhi(isp))
      enddo 

      ! obtain screened Coulomb matrix
      call bse_buildw(pars,dbl_mb(k_dia),dbl_mb(kW),pars%nmo,nri,
     &                maxpoles,pars%ipol)

      ! transform ERIs using screened Coulomb matrix
      call bse_wmn(pars,dbl_mb(kW),pars%nmo,nri,pars%ipol)
      if(.not.ma_chop_stack(lW))
     &  call errquit(pname//'failed to chop stack',0,MA_ERR) 

      ! obtain "GW" eigenvalue difference to build Casida matrix
      do isp=1,pars%ipol
         lSize = (isp-1)*maxpoles
         call bse_get_eia(dbl_mb(pars%k_gw_evals+(isp-1)*nbf_ao),
     &                   dbl_mb(k_dia+lSize),pars%nocc(isp),
     &                   pars%nvir(isp),pars%ovlo(isp),
     &                   pars%ovhi(isp))    
      enddo 

      ! Get diagonal of Casida matrix
      call bse_davidson_diagonal(pars,dbl_mb(k_dia),maxpoles)

      ! Get Eigenvalue differences
      do isp=1,pars%ipol
        off1 = (isp-1)*pars%npoles(1)
        lSize = (isp-1)*maxpoles
        call nga_put(g_wia,pars%ovlo(isp)+off1,pars%ovhi(isp)+off1,
     $           dbl_mb(k_dia+lSize),pars%ovhi(isp)-pars%ovlo(isp)+1)
      enddo

      ! Get guess vectos
      if(pars%me.eq.0) 
     $  call bse_davidson_guess(pars,dbl_mb(k_dia),g_trials,npoles,
     $                          ntrials)
      call ga_brdcst(mt_dbl,diagonal,ma_sizeof(mt_dbl,npoles,mt_byte),0)

      call ga_distribution(g_s,pars%me,mylo,myhi,idum,idum)
      myld = myhi - mylo + 1

      ! K-orthogonalization
      call bse_davidson_kvec(pars,g_trials,g_wia,g_x,npoles,
     $                       ntrials,pars%singlet)
      call ga_access(g_x,mylo,myhi,1,ntrials,k_x,myld)
      call ga_access(g_trials,mylo,myhi,1,ntrials,k_trials,myld)
      call bse_davidson_kortho(pars,dbl_mb(k_trials),dbl_mb(k_x),
     $    dbl_mb(k_x),dbl_mb(vkv),dbl_mb(k_x),myld,ntrials,mdim,.true.)
      call ga_release_update(g_trials,mylo,myhi,1,ntrials)
      call ga_release_update(g_x,mylo,myhi,1,ntrials)

      ! Starting vectors
      call bse_davidson_kvec(pars,g_x,g_wia,g_y,npoles,
     $                       ntrials,pars%singlet)
      call bse_davidson_mvec(pars,g_y,g_wia,g_ym,npoles,
     $                       ntrials)
      call ga_copy(g_x,g_s)
      call ga_copy(g_y,g_ks)
      call ga_copy(g_ym,g_mks)

      if (pars%me.eq.0) then
        write(*,*)
        write(*,*) "--------------------------------"
        write(*,*) " Iter  Ntrials  Nconv     Res   "
        write(*,*) "--------------------------------"
      endif

      do idavidson=1,20

        ! Get access to local chunks
        call ga_access(g_s,mylo,myhi,1,ntrials,k_s,myld)
        call ga_access(g_ks,mylo,myhi,1,ntrials,k_ks,myld)
        call ga_access(g_mks,mylo,myhi,1,ntrials,k_mks,myld)
        call ga_access(g_x,mylo,myhi,1,ntrials,k_x,myld)
        call ga_access(g_y,mylo,myhi,1,ntrials,k_y,myld)
        call ga_access(g_ym,mylo,myhi,1,ntrials,k_ym,myld)

        ! Build projected matrix
        call dgemm('t','n',ntrials,ntrials,myld,1d0,dbl_mb(k_ks),
     $              myld,dbl_mb(k_mks),myld,0d0,dbl_mb(vkv),ntrials)
        call ga_dgop((/194/),dbl_mb(vkv),ntrials**2,'+')
        if (debug) call ma_print(dbl_mb(vkv),ntrials,ntrials,'Proj')

        ! Diagonalize projected matrix
        call dsyev('v','l',ntrials,dbl_mb(vkv),ntrials,dbl_mb(k_omega),
     $              rdum,-1,info)
        if(.not.ma_push_get(mt_dbl,int(rdum),'work',lwork,kwork))
     &    call errquit(pname//'failed to allocate work',199,MA_ERR)
        call dsyev('v','l',ntrials,dbl_mb(vkv),ntrials,dbl_mb(k_omega),
     $             dbl_mb(kwork),int(rdum),info)
        if(.not.ma_chop_stack(lwork))
     &    call errquit(pname//'failed to allocate work',192,MA_ERR)
        if (info.ne.0)
     $    call errquit(pname//'Diagonalization failed',206,0)
        if (debug) call ma_print(dbl_mb(vkv),ntrials,ntrials,'Eigen')


        ! Project back vectors
        call dgemm('n','n',myld,ntrials,ntrials,1d0,dbl_mb(k_s),myld,
     $              dbl_mb(vkv),ntrials,0d0,dbl_mb(k_x),myld)
        call dgemm('n','n',myld,ntrials,ntrials,1d0,dbl_mb(k_ks),myld,
     $              dbl_mb(vkv),ntrials,0d0,dbl_mb(k_y),myld)
        call dgemm('n','n',myld,ntrials,ntrials,1d0,dbl_mb(k_mks),myld,
     $              dbl_mb(vkv),ntrials,0d0,dbl_mb(k_ym),myld)

        ! Get residual
        call ga_access(g_w,mylo,myhi,1,pars%nroots,k_w,myld)

#if 1
        call bse_davidson_residual(pars%nroots,ntrials,dbl_mb(k_omega),
     $         dbl_mb(k_x),dbl_mb(k_ym),dbl_mb(k_r),
     $         dbl_mb(k_dia+mylo-1),dbl_mb(k_w),myld,ndim,norm)

        ! Print info
        if (pars%me.eq.0) then
          write(*,9010) idavidson,ntrials,pars%nroots-ndim,maxval(norm)
        endif
 9010   format(I5,3X,I5,3X,I5,G16.5)
#else        
        minimum = 0d0
        ndim = 0
        do itrial=1,pars%nroots
          call dcopy(myld,dbl_mb(k_ym+(itrial-1)*myld),1,dbl_mb(k_r),1)
          call daxpy(myld,-dbl_mb(k_omega+itrial-1),
     $               dbl_mb(k_x+(itrial-1)*myld),1,dbl_mb(k_r),1)
          norm = maxval(abs(dbl_mb(k_r:k_r+myhi-mylo)))
          call ga_dgop((/229/),norm,1,'maxabs')
          write(*,*) "       ",itrial,norm(1),
     $       dsqrt(dbl_mb(k_omega+itrial-1))*27.2114d0
          if (norm(1).gt.1d-5) then
            ndim = ndim + 1
            do ipole=mylo,myhi
              factor = (dbl_mb(k_dia+ipole-1)-dbl_mb(k_omega+itrial-1))
              dbl_mb(k_r+ipole-mylo) = dbl_mb(k_r+ipole-mylo)/factor
            enddo
            call ga_put(g_w,mylo,myhi,ndim,ndim,dbl_mb(k_r),myld)
          endif
        enddo
#endif        

        call ga_release_update(g_w,mylo,myhi,1,pars%nroots)
        call ga_release_update(g_x,mylo,myhi,1,ntrials)
        call ga_release_update(g_y,mylo,myhi,1,ntrials)
        call ga_release_update(g_ym,mylo,myhi,1,ntrials)
        call ga_release_update(g_s,mylo,myhi,1,ntrials)
        call ga_release_update(g_ks,mylo,myhi,1,ntrials)
        call ga_release_update(g_mks,mylo,myhi,1,ntrials)

        ! All residuals are below the threshold
        if (ndim.eq.0) exit

        ! Restart calculation if dimensions exceeded
        if (ndim+ntrials.gt.pars%nmax) then
          call bse_davidson_restart(ntrials,pars%nroots,npoles,
     $          g_x,g_y,g_ym,g_s,g_ks,g_mks,l_vkv,vkv)
          cycle
        endif

        ! Project-out current vectors
        call ga_access(g_ks,mylo,myhi,1,ntrials,k_ks,myld)
        call ga_access(g_s,mylo,myhi,1,ntrials,k_s,myld)
        call ga_access(g_w,mylo,myhi,1,pars%nroots,k_w,myld)
        if(.not.ma_push_get(mt_dbl,ndim*ntrials,'ksw',l_ksw,ksw))
     $    call errquit(pname//'ksw allocation failed',264,MA_ERR)
        call dgemm('t','n',ntrials,ndim,myld,1d0,dbl_mb(k_ks),myld,
     $              dbl_mb(k_w),myld,0d0,dbl_mb(ksw),ntrials)
        call ga_dgop((/268/),dbl_mb(ksw),ntrials*ndim,'+')
        call dgemm('n','n',myld,ndim,ntrials,-1d0,dbl_mb(k_s),myld,
     $              dbl_mb(ksw),ntrials,1d0,dbl_mb(k_w),myld)
        if(.not.ma_chop_stack(l_ksw))
     $    call errquit(pname//'ksw deallocation failed',272,MA_ERR)
        call ga_release(g_ks,mylo,myhi,1,ntrials)
        call ga_release(g_s,mylo,myhi,1,ntrials)
        call ga_release_update(g_w,mylo,myhi,1,pars%nroots)


        ! K-orthogonalization
        call bse_davidson_kvec(pars,g_w,g_wia,g_kw,npoles,
     $                         ndim,pars%singlet)
        call bse_davidson_mvec(pars,g_kw,g_wia,g_mkw,npoles,
     $                         ndim)
        call ga_access(g_w,mylo,myhi,1,pars%nroots,k_w,myld)
        call ga_access(g_kw,mylo,myhi,1,pars%nroots,k_kw,myld)
        call ga_access(g_mkw,mylo,myhi,1,pars%nroots,k_mkw,myld)
        call ga_access(g_x,mylo,myhi,1,mdim,k_x,myld)
        call bse_davidson_kortho(pars,dbl_mb(k_w),dbl_mb(k_kw),
     $    dbl_mb(k_mkw),dbl_mb(vkv),dbl_mb(k_x),myld,ndim,mdim,
     $    .false.)
        call ga_release_update(g_x,mylo,myhi,1,mdim)
        call ga_release_update(g_w,mylo,myhi,1,pars%nroots)
        call ga_release_update(g_kw,mylo,myhi,1,pars%nroots)
        call ga_release_update(g_mkw,mylo,myhi,1,pars%nroots)

        ! All trial vectors were linearly-dependent
        ! restart calculation with current best guess
        if (mdim.eq.0) then
          call bse_davidson_restart(ntrials,pars%nroots,npoles,
     $          g_x,g_y,g_ym,g_s,g_ks,g_mks,l_vkv,vkv)
          cycle
        endif

        ndim = mdim

        ! Extend the size of the current vectors
        call bse_davidson_newvecs(ntrials,ndim,npoles,g_x,g_y,g_ym,
     $                            g_s,g_ks,g_mks,g_w,g_kw,g_mkw,
     $                            l_vkv,vkv)

      enddo

      stat = ga_destroy(g_x).and.ga_destroy(g_y).and.ga_destroy(g_ym)
      stat = stat.and.ga_destroy(g_s).and.ga_destroy(g_ks)
      stat = stat.and.ga_destroy(g_mks).and.ga_destroy(g_w)
      stat = stat.and.ga_destroy(g_kw).and.ga_destroy(g_mkw)
      stat = stat.and.ga_destroy(g_trials).and.ga_destroy(g_wia)
      stat = stat.and.ga_destroy(pars%g_apb).and.ga_destroy(pars%g_amb)
      do isp=1,pars%ipol
        stat = stat.and.ga_destroy(pars%g_erioo(isp))
        stat = stat.and.ga_destroy(pars%g_eriov(isp))
        stat = stat.and.ga_destroy(pars%g_erivv(isp))
        stat = stat.and.ga_destroy(pars%g_wov(isp))
      enddo
      if(.not.stat) call errquit('deallocation failed',0,GA_ERR) 
      stat = ma_chop_stack(l_dia)
      if(.not.stat) call errquit('deallocation failed',0,MA_ERR) 

      end subroutine bse_davidson




