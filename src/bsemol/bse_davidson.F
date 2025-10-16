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

      integer ntrials,g_trials,npoles,ndim
      integer g_x,g_y,g_ym,g_s,g_ks,g_mks,g_w,g_kw,g_mkw,off1
      integer g_wia,g_tmp,info,lwork,kwork,idum,mdim,k_trials
      integer k_x,k_y,k_ym,k_s,k_ks,k_mks,k_w,k_kw,k_mkw,k_r,l_r
      integer l_map,k_map,l_tmp1,k_tmp1,l_tmp2,k_tmp2,k_x2,k_trials2
      integer k_w2,l_tmp3,k_tmp3,l_tmp4,k_tmp4,k_kw2,k_mkw2
      integer myld,vkv,l_vkv,myhi,mylo,l_ksw,ksw,jdum,myld2
      integer subscript(2)
      double precision minimum,norm(pars%nroots),rdum,factor

      integer inode,nprocs,nrestart
      integer ipole,jpole,itrial,imo,amo,idavidson
      logical found,dohartree,stat,debug
      double precision,external :: dasum

      !debug = pars%me.eq.0
      debug = .false.

      ! synchronize parallel processes
      call ga_sync()

      !initialize useful variables
      nri = pars%nri
      maxpoles = maxval(pars%mynpoles(1:pars%ipol))
      maxnpoles = maxval(pars%npoles(1:pars%ipol))
      totnpoles = maxnpoles*pars%ipol
      npoles = sum(pars%npoles(1:pars%ipol))
      ntrials = pars%nroots*pars%nspace
      nprocs = ga_nnodes()

      ! allocation
      if(.not.ma_push_get(mt_dbl,totnpoles,'dia',l_dia,k_dia))
     &  call errquit(pname//'failed to allocate Delta_ia',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,npoles,'omega',l_omega,k_omega))
     &  call errquit(pname//'failed to allocate omega',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,npoles,'residual',l_r,k_r))
     &  call errquit(pname//'failed to allocate residual',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,maxpoles*pars%nmax,'tmp1',
     &                    l_tmp1,k_tmp1))
     &  call errquit(pname//'failed to allocate tmp1',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,maxpoles*pars%nmax,'tmp2',
     &                    l_tmp2,k_tmp2))
     &  call errquit(pname//'failed to allocate tmp2',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,maxpoles*pars%nmax,'tmp3',
     &                    l_tmp3,k_tmp3))
     &  call errquit(pname//'failed to allocate tmp3',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,maxpoles*pars%nmax,'tmp4',
     &                    l_tmp4,k_tmp4))
     &  call errquit(pname//'failed to allocate tmp4',0,MA_ERR)
      if(.not.ma_push_get(mt_int,2*ga_nnodes()+1,'map',l_map,k_map))
     &  call errquit(pname//'failed to allocate map',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,ntrials*ntrials,'vkv',l_vkv,vkv))
     &  call errquit(pname//'failed to allocate vkv',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nri**2,'W',lW,kW))
     &  call errquit(pname//'failed to allocate W',0,MA_ERR)

      ! Get map for irreg array
      pars%nblocks = ga_nnodes()
      do inode=0,ga_nnodes()-1
        idum = pars%ovlo(1)
        call ga_brdcst(mt_int,idum,ma_sizeof(mt_int,1,mt_byte),inode)
        int_mb(k_map+inode) = idum
      enddo
      if (pars%ipol.gt.1) then
        pars%nblocks = pars%nblocks + ga_nnodes()
        do inode=0,ga_nnodes()-1
          idum = pars%ovlo(2)
          call ga_brdcst(mt_int,idum,ma_sizeof(mt_int,1,mt_byte),inode)
          int_mb(k_map+ga_nnodes()+inode) = idum + pars%npoles(1)
        enddo
      endif
      int_mb(k_map+pars%nblocks) = 1

      ! Allocate GA arrays
      g_x = ga_create_handle()
      call ga_set_data(g_x,2,(/npoles,ntrials/),mt_dbl)
      call ga_set_array_name(g_x,'x')
      call ga_set_tiled_irreg_proc_grid(g_x,int_mb(k_map),
     $   (/pars%nblocks,1/),(/ga_nnodes(),1/))
      stat = ga_allocate(g_x)
      stat = stat.and.ga_duplicate(g_x,g_y,'y').and.
     $                ga_duplicate(g_x,g_ym,'ym').and. 
     $                ga_duplicate(g_x,g_s,'s').and. 
     $                ga_duplicate(g_x,g_ks,'ks').and. 
     $                ga_duplicate(g_x,g_mks,'mks').and. 
     $                ga_duplicate(g_x,g_w,'w').and. 
     $                ga_duplicate(g_x,g_kw,'kw').and. 
     $                ga_duplicate(g_x,g_mkw,'mkw').and. 
     $                ga_duplicate(g_x,g_trials,'trials')
      if(.not.stat) call errquit('ga create failed',103,GA_ERR)

      g_wia = ga_create_handle()
      call ga_set_data(g_wia,1,npoles,mt_dbl)
      call ga_set_array_name(g_wia,'wia')
      call ga_set_tiled_irreg_proc_grid(g_wia,int_mb(k_map),
     $                pars%nblocks,ga_nnodes())
      stat = stat.and.ga_allocate(g_wia)
      if(.not.stat) call errquit('ga create failed',110,GA_ERR)


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
     $                          ntrials,pars%nroots)
      call ga_brdcst(mt_dbl,dbl_mb(k_dia),
     $                           ma_sizeof(mt_dbl,npoles,mt_byte),0)
      do isp=1,pars%ipol
        off1 = (isp-1)*pars%npoles(1)
        call dcopy(pars%mynpoles(isp),
     $             dbl_mb(k_dia+off1+pars%ovlo(isp)-1),
     $             1,dbl_mb(k_dia+(isp-1)*pars%mynpoles(1)),1)
      enddo


      do isp=1,pars%ipol
        if(.not.ga_create(mt_dbl,nri*pars%nvir(isp),pars%nvir(isp),
     $          'g_erivv',0,pars%nvir(isp),g_tmp))
     $  call errquit('could not create GA',0,GA_ERR)
        call ga_copy_patch('n',pars%g_erivv(isp),1,nri,1,
     $                      pars%nvir(isp)*pars%nvir(isp),
     $                      g_tmp,1,nri*pars%nvir(isp),1,
     $                      pars%nvir(isp))
        if(.not.ga_destroy(pars%g_erivv(isp)))
     $      call errquit('could not destroy GA',0,GA_ERR)
        pars%g_erivv(isp) = g_tmp
      enddo

      ! Start K-Davidson Algorithm
      idavidson = 0

      ! Restarts should land here
    1 continue

      ! K-orthogonalization
      call bse_davidson_kvec(pars,g_trials,g_wia,g_x,npoles,
     $                       ntrials,pars%singlet,1d0)

      call nga_access_block_grid(g_x, (/pars%me,0/), k_x, myld)
      call nga_access_block_grid(g_trials, (/pars%me,0/),k_trials,myld)
      if (pars%ipol.eq.1) then
        myld2 = 0
        call dcopy(myld*ntrials,dbl_mb(k_x),1,dbl_mb(k_tmp1),1)
        call dcopy(myld*ntrials,dbl_mb(k_trials),1,dbl_mb(k_tmp2),1)
      else
        call nga_access_block_grid(g_x,(/pars%me+nprocs,0/),k_x2,myld2)
        call nga_access_block_grid(g_trials,(/pars%me+nprocs,0/),
     $                            k_trials2,myld2)
        do itrial=1,ntrials
          call dcopy(myld,dbl_mb(k_x+(itrial-1)*myld),1,
     $                 dbl_mb(k_tmp1+(itrial-1)*(myld+myld2)),1)
          call dcopy(myld2,dbl_mb(k_x2+(itrial-1)*myld2),1,
     $                 dbl_mb(k_tmp1+(itrial-1)*(myld+myld2)+myld),1)
          call dcopy(myld,dbl_mb(k_trials+(itrial-1)*myld),1,
     $                 dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)),1)
          call dcopy(myld2,dbl_mb(k_trials2+(itrial-1)*myld2),1,
     $                 dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)+myld),1)
        enddo
      endif

      call bse_davidson_kortho(pars,dbl_mb(k_tmp2),dbl_mb(k_tmp1),
     $    dbl_mb(k_tmp1),dbl_mb(vkv),dbl_mb(k_tmp1),
     $    myld+myld2,ntrials,mdim,.true.)

      if (pars%ipol.eq.1) then
        call dcopy(myld*ntrials,dbl_mb(k_tmp2),1,dbl_mb(k_trials),1)
      else
        do itrial=1,ntrials
          call dcopy(myld,dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)),1,
     $               dbl_mb(k_trials+(itrial-1)*myld),1)
          call dcopy(myld2,dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)+myld),
     $                1,dbl_mb(k_trials2+(itrial-1)*myld2),1)
        enddo
        call nga_release_update_block_grid(g_trials,
     $          (/pars%me+nprocs,0/))
        call nga_release_block_grid(g_x,(/pars%me+nprocs,0/))
      endif
      call nga_release_update_block_grid(g_trials,(/pars%me,0/))
      call nga_release_block_grid(g_x,(/pars%me,0/))

      ! Starting vectors
      call ga_copy(g_trials,g_x)
      call bse_davidson_kvec(pars,g_x,g_wia,g_y,npoles,
     $                       ntrials,pars%singlet,1d0)
      call bse_davidson_kvec(pars,g_y,g_wia,g_ym,npoles,
     $                       ntrials,.false.,-1d0)
      call ga_copy(g_x,g_s)
      call ga_copy(g_y,g_ks)
      call ga_copy(g_ym,g_mks)

      if (pars%me.eq.0 .and. idavidson.eq.0) then
        write(*,*)
        write(*,*) "--------------------------------"
        write(*,*) " Iter  Ntrials  Nconv     Res   "
        write(*,*) "--------------------------------"
      endif

    2 idavidson = idavidson + 1
      if(idavidson.gt.pars%maxiter) goto 3

        ! Build projected matrix
        call dcopy(ntrials*ntrials,0d0,0,dbl_mb(vkv),1)
        do isp=1,pars%ipol
          subscript = (/pars%me+(isp-1)*ga_nnodes(),0/)
          call nga_access_block_grid(g_ks,subscript,k_ks,myld)
          call nga_access_block_grid(g_mks,subscript,k_mks,myld)
          call dgemm('t','n',ntrials,ntrials,myld,1d0,dbl_mb(k_ks),
     $                myld,dbl_mb(k_mks),myld,1d0,dbl_mb(vkv),ntrials)
          call nga_release_block_grid(g_ks,subscript)
          call nga_release_block_grid(g_mks,subscript)
        enddo
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
        if (debug) call ma_print(dbl_mb(k_omega),ntrials,1,'Omega')


        ! Project back vectors
        do isp=1,pars%ipol
          subscript = (/pars%me+(isp-1)*nprocs,0/)
          call nga_access_block_grid(g_s,subscript,k_s,myld)
          call nga_access_block_grid(g_ks,subscript,k_ks,myld)
          call nga_access_block_grid(g_mks,subscript,k_mks,myld)
          call nga_access_block_grid(g_x,subscript,k_x,myld)
          call nga_access_block_grid(g_y,subscript,k_y,myld)
          call nga_access_block_grid(g_ym,subscript,k_ym,myld)

          call dgemm('n','n',myld,ntrials,ntrials,1d0,dbl_mb(k_s),myld,
     $              dbl_mb(vkv),ntrials,0d0,dbl_mb(k_x),myld)
          call dgemm('n','n',myld,ntrials,ntrials,1d0,dbl_mb(k_ks),myld,
     $              dbl_mb(vkv),ntrials,0d0,dbl_mb(k_y),myld)
          call dgemm('n','n',myld,ntrials,ntrials,1d0,dbl_mb(k_mks),
     $              myld,dbl_mb(vkv),ntrials,0d0,dbl_mb(k_ym),myld)

          call nga_release_block_grid(g_s,subscript)
          call nga_release_block_grid(g_ks,subscript)
          call nga_release_block_grid(g_mks,subscript)
          call nga_release_update_block_grid(g_x,subscript)
          call nga_release_update_block_grid(g_y,subscript)
          call nga_release_update_block_grid(g_ym,subscript)
        enddo

        ! Get residual
        subscript = (/pars%me,0/)
        call nga_access_block_grid(g_x,subscript,k_x,myld)
        call nga_access_block_grid(g_ym,subscript,k_trials,myld)
        if (pars%ipol.eq.1) then
          myld2 = 0
          call dcopy(myld*ntrials,dbl_mb(k_x),1,dbl_mb(k_tmp1),1)
          call dcopy(myld*ntrials,dbl_mb(k_trials),1,dbl_mb(k_tmp2),1)
        else
          subscript = (/pars%me+nprocs,0/)
          call nga_access_block_grid(g_x,subscript,k_x2,myld2)
          call nga_access_block_grid(g_ym,subscript,k_trials2,myld2)
          do itrial=1,pars%nroots
            call dcopy(myld,dbl_mb(k_x+(itrial-1)*myld),1,
     $                  dbl_mb(k_tmp1+(itrial-1)*(myld+myld2)),1)
            call dcopy(myld2,dbl_mb(k_x2+(itrial-1)*myld2),1,
     $                  dbl_mb(k_tmp1+(itrial-1)*(myld+myld2)+myld),1)
            call dcopy(myld,dbl_mb(k_trials+(itrial-1)*myld),1,
     $                  dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)),1)
            call dcopy(myld2,dbl_mb(k_trials2+(itrial-1)*myld2),1,
     $                  dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)+myld),1)
          enddo
        endif

        call bse_davidson_residual(pars%nroots,ntrials,dbl_mb(k_omega),
     $         dbl_mb(k_tmp1),dbl_mb(k_tmp2),dbl_mb(k_r),
     $         dbl_mb(k_dia),dbl_mb(k_tmp3),myld+myld2,ndim,norm)

        call nga_release_block_grid(g_x,(/pars%me,0/))
        call nga_release_block_grid(g_ym,(/pars%me,0/))
        call nga_access_block_grid(g_w,(/pars%me,0/),k_w,myld)
        if(pars%ipol.eq.1) then
          call dcopy(pars%nroots*myld,dbl_mb(k_tmp3),1,dbl_mb(k_w),1)
        else
          call nga_release_block_grid(g_x,subscript)
          call nga_release_block_grid(g_ym,subscript)
          call nga_access_block_grid(g_w,subscript,k_w2,myld2)
          do itrial=1,pars%nroots
            call dcopy(myld,dbl_mb(k_tmp3+(itrial-1)*(myld+myld2)),1,
     $                      dbl_mb(k_w+(itrial-1)*myld),1)
            call dcopy(myld2,dbl_mb(k_tmp3+
     $                     (itrial-1)*(myld+myld2)+myld),1,
     $                     dbl_mb(k_w2+(itrial-1)*myld2),1)
          enddo
          call nga_release_update_block_grid(g_w,subscript)
        endif
        call nga_release_update_block_grid(g_w,(/pars%me,0/))

        ! Print info
        if (pars%me.eq.0) then
          write(*,9010) idavidson,ntrials,pars%nroots-ndim,maxval(norm)
        endif
 9010   format(I5,3X,I5,3X,I5,G16.5)

        ! All residuals are below the threshold
        if (ndim.eq.0) goto 3

        ! Restart calculation if dimensions exceeded
        if (ndim+ntrials.gt.pars%nmax) then
          nrestart = min(ntrials, pars%nroots*pars%nspace)
          call bse_davidson_restart(ntrials,nrestart,npoles,
     $          g_x,g_y,g_ym,g_s,g_ks,g_mks,g_trials,l_vkv,vkv,
     $          int_mb(k_map),pars%nblocks)
          goto 1
        endif

        ! Project-out current vectors
        if(.not.ma_push_get(mt_dbl,ndim*ntrials,'ksw',l_ksw,ksw))
     $    call errquit(pname//'ksw allocation failed',264,MA_ERR)
        call dcopy(ntrials*ndim,0d0,0,dbl_mb(ksw),1)
        do isp=1,pars%ipol
          subscript = (/pars%me+(isp-1)*nprocs,0/)
          call nga_access_block_grid(g_ks,subscript,k_ks,myld)
          call nga_access_block_grid(g_w,subscript,k_w,myld)
          call dgemm('t','n',ntrials,ndim,myld,1d0,dbl_mb(k_ks),myld,
     $              dbl_mb(k_w),myld,1d0,dbl_mb(ksw),ntrials)
          call nga_release_block_grid(g_ks,subscript)
          call nga_release_block_grid(g_w,subscript)
        enddo
        call ga_dgop((/268/),dbl_mb(ksw),ntrials*ndim,'+')

        do isp=1,pars%ipol
          subscript = (/pars%me+(isp-1)*nprocs,0/)
          call nga_access_block_grid(g_s,subscript,k_s,myld)
          call nga_access_block_grid(g_w,subscript,k_w,myld)
          call dgemm('n','n',myld,ndim,ntrials,-1d0,dbl_mb(k_s),myld,
     $                dbl_mb(ksw),ntrials,1d0,dbl_mb(k_w),myld)
          call nga_release_update_block_grid(k_w,subscript)
          call nga_release_block_grid(k_s,subscript)
        enddo
        if(.not.ma_chop_stack(l_ksw))
     $    call errquit(pname//'ksw deallocation failed',272,MA_ERR)


        ! K-orthogonalization
        call bse_davidson_kvec(pars,g_w,g_wia,g_kw,npoles,
     $                         ndim,pars%singlet,1d0)
        call bse_davidson_kvec(pars,g_kw,g_wia,g_mkw,npoles,
     $                         ndim,.false.,-1d0)

        subscript = (/pars%me,0/)
        call nga_access_block_grid(g_w,subscript,k_w,myld)
        call nga_access_block_grid(g_kw,subscript,k_kw,myld)
        call nga_access_block_grid(g_mkw,subscript,k_mkw,myld)
        if (pars%ipol.eq.1) then
          myld2 = 0
          call dcopy(ndim*myld,dbl_mb(k_w),1,dbl_mb(k_tmp1),1)
          call dcopy(ndim*myld,dbl_mb(k_kw),1,dbl_mb(k_tmp2),1)
          call dcopy(ndim*myld,dbl_mb(k_mkw),1,dbl_mb(k_tmp3),1)
        else
          subscript = (/pars%me+nprocs,0/)
          call nga_access_block_grid(g_w,subscript,k_w2,myld2)
          call nga_access_block_grid(g_kw,subscript,k_kw2,myld2)
          call nga_access_block_grid(g_mkw,subscript,k_mkw2,myld2)
          do itrial=1,ndim
            call dcopy(myld,dbl_mb(k_w+(itrial-1)*myld),1,
     $         dbl_mb(k_tmp1+(itrial-1)*(myld+myld2)),1)
            call dcopy(myld2,dbl_mb(k_w2+(itrial-1)*myld2),1,
     $         dbl_mb(k_tmp1+(itrial-1)*(myld+myld2)+myld),1)
            call dcopy(myld,dbl_mb(k_kw+(itrial-1)*myld),1,
     $         dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)),1)
            call dcopy(myld2,dbl_mb(k_kw2+(itrial-1)*myld2),1,
     $         dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)+myld),1)
            call dcopy(myld,dbl_mb(k_mkw+(itrial-1)*myld),1,
     $         dbl_mb(k_tmp3+(itrial-1)*(myld+myld2)),1)
            call dcopy(myld2,dbl_mb(k_mkw2+(itrial-1)*myld2),1,
     $         dbl_mb(k_tmp3+(itrial-1)*(myld+myld2)+myld),1)
          enddo
        endif

        call bse_davidson_kortho(pars,dbl_mb(k_tmp1),dbl_mb(k_tmp2),
     $    dbl_mb(k_tmp3),dbl_mb(vkv),dbl_mb(k_tmp4),myld+myld2,
     $    ndim,mdim,.false.)

        if(pars%ipol.eq.1) then
          call dcopy(mdim*myld,dbl_mb(k_tmp1),1,dbl_mb(k_w),1)
          call dcopy(mdim*myld,dbl_mb(k_tmp2),1,dbl_mb(k_kw),1)
          call dcopy(mdim*myld,dbl_mb(k_tmp3),1,dbl_mb(k_mkw),1)
        else
          do itrial=1,mdim
            call dcopy(myld,dbl_mb(k_tmp1+(itrial-1)*(myld+myld2)),1,
     $         dbl_mb(k_w+(itrial-1)*myld),1)
            call dcopy(myld2,
     $         dbl_mb(k_tmp1+(itrial-1)*(myld+myld2)+myld),1,
     $         dbl_mb(k_w2+(itrial-1)*myld2),1)
            call dcopy(myld,dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)),1,
     $         dbl_mb(k_kw+(itrial-1)*myld),1)
            call dcopy(myld2,
     $         dbl_mb(k_tmp2+(itrial-1)*(myld+myld2)+myld),1,
     $         dbl_mb(k_kw2+(itrial-1)*myld2),1)
            call dcopy(myld,dbl_mb(k_tmp3+(itrial-1)*(myld+myld2)),1,
     $         dbl_mb(k_mkw+(itrial-1)*myld),1)
            call dcopy(myld2,
     $         dbl_mb(k_tmp3+(itrial-1)*(myld+myld2)+myld),1,
     $         dbl_mb(k_mkw2+(itrial-1)*myld2),1)
          enddo
          call nga_release_update_block_grid(g_w,subscript)
          call nga_release_update_block_grid(g_kw,subscript)
          call nga_release_update_block_grid(g_mkw,subscript)
        endif
        call nga_release_update_block_grid(g_w,(/pars%me,0/))
        call nga_release_update_block_grid(g_kw,(/pars%me,0/))
        call nga_release_update_block_grid(g_mkw,(/pars%me,0/))

        ! All trial vectors were linearly-dependent
        ! restart calculation with current best guess
        if (mdim.lt.3) then
          nrestart = min(ntrials, pars%nroots*pars%nspace)
          call bse_davidson_restart(ntrials,nrestart,npoles,
     $          g_x,g_y,g_ym,g_s,g_ks,g_mks,g_trials,l_vkv,vkv,
     $          int_mb(k_map),pars%nblocks)
          goto 1
        endif

        ndim = mdim

        ! Extend the size of the current vectors
        call bse_davidson_newvecs(ntrials,ndim,npoles,g_x,g_y,g_ym,
     $                            g_s,g_ks,g_mks,g_w,g_kw,g_mkw,
     $                            l_vkv,vkv,int_mb(k_map),pars%nblocks)
        goto 2

    3 continue

      if(pars%me.eq.0) write(*,*)

      stat = ga_destroy(pars%g_apb) .and. ga_destroy(pars%g_amb)
      if(.not.stat) call errquit('deallocation failed',481,GA_ERR) 

      stat = ga_create(mt_dbl,npoles,pars%nroots,'X',0,pars%nroots,
     $                 pars%g_apb)
      stat = ga_create(mt_dbl,npoles,pars%nroots,'Y',0,pars%nroots,
     $                 pars%g_amb)
      if(.not.stat) call errquit('failed to creare GA',484,GA_ERR)
      call nga_copy_patch('n',g_x,(/1,1/),(/npoles,pars%nroots/),
     $                 pars%g_apb,(/1,1/),(/npoles,pars%nroots/))
      call nga_copy_patch('n',g_y,(/1,1/),(/npoles,pars%nroots/),
     $                 pars%g_amb,(/1,1/),(/npoles,pars%nroots/))

      stat = ga_destroy(g_x).and.ga_destroy(g_y).and.ga_destroy(g_ym)
      stat = stat.and.ga_destroy(g_s).and.ga_destroy(g_ks)
      stat = stat.and.ga_destroy(g_mks).and.ga_destroy(g_w)
      stat = stat.and.ga_destroy(g_kw).and.ga_destroy(g_mkw)
      stat = stat.and.ga_destroy(g_trials).and.ga_destroy(g_wia)
      do isp=1,pars%ipol
        stat = stat.and.ga_destroy(pars%g_erioo(isp))
        stat = stat.and.ga_destroy(pars%g_eriov(isp))
        stat = stat.and.ga_destroy(pars%g_erivv(isp))
        stat = stat.and.ga_destroy(pars%g_wov(isp))
      enddo
      if(.not.stat) call errquit('deallocation failed',496,GA_ERR) 

      if (.not.pars%tda) then
        do ipole=0,pars%nroots-1
          dbl_mb(k_omega+ipole) = dsqrt(dbl_mb(k_omega+ipole))
          call ga_scale_patch(pars%g_apb,1,npoles,ipole+1,ipole+1,
     $              dsqrt(dbl_mb(k_omega+ipole)))
          call ga_scale_patch(pars%g_amb,1,npoles,ipole+1,ipole+1,
     $              1.0/dsqrt(dbl_mb(k_omega+ipole))) 
        enddo
      endif

      call ga_add(0.5d0,pars%g_apb,0.5d0,pars%g_amb,pars%g_apb)
      call ga_add(1d0,pars%g_apb,-1d0,pars%g_amb,pars%g_amb)
      call bse_oscstr(pars,dbl_mb(k_omega),npoles)

      stat = ma_chop_stack(l_dia)
      if(.not.stat) call errquit('deallocation failed',501,MA_ERR) 

      ! Destroy GAs
      if (.not.ga_destroy(pars%g_apb))
     &    call errquit(pname//'could not destroy temp GA',519,GA_ERR)
      if (.not.ga_destroy(pars%g_amb))
     &    call errquit(pname//'could not destroy temp GA',519,GA_ERR)

      end subroutine bse_davidson




