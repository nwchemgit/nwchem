
      subroutine bse_davidson_mvec(pars,vec,wia,mv,npoles,ntrials)
      implicit none 
#include "bse.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      type(bse_params_t) :: pars

      integer npoles,ntrials
      integer vec,wia,kv,mv

      integer ipol,npoles1,npoles2,uli,lli,ipole
      integer nocc(2),nvir(2),amo,imo,itrial,ula,lla,nri

      integer g_tmp1(2),g_tmp2(2)
      double precision temp(npoles),factor

      ipol = pars%ipol
      nri = pars%nri
      nocc = pars%nocc
      nvir = pars%nvir
      npoles1 = nocc(1)*nvir(1)
      npoles2 = nocc(2)*nvir(2)

      ! Create temporary arrays
      if (.not.ga_create(mt_dbl,nvir(1),nocc(1),'g_tmp1',0,nocc(1),
     $                   g_tmp1))
     $  call errquit('bse_davidson_kvec: ga create failed',34,GA_ERR)
      if (.not.ga_create(mt_dbl,nri,nocc(1),'g_tmp2',0,nocc(1),
     $                   g_tmp2))
     $  call errquit('bse_davidson_kvec: ga create failed',37,GA_ERR)

      if (pars%ipol.gt.1) then
        if (.not.ga_create(mt_dbl,nvir(2),nocc(2),'g_tmp1',0,nocc(2),
     $                   g_tmp1(2)))
     $    call errquit('bse_davidson_kvec: ga create failed',48,GA_ERR)
        if (.not.ga_create(mt_dbl,nri,nocc(2),'g_tmp2',0,nocc(2),
     $                   g_tmp2(2)))
     $    call errquit('bse_davidson_kvec: ga create failed',51,GA_ERR)
      endif


      call ga_copy(vec,mv)
      call ga_scale_rows(mv,wia)

      do itrial=1,ntrials
        call ga_copy_patch('n',vec,1,npoles1,itrial,itrial,
     $                        g_tmp1,1,nvir(1),1,nocc(1))
        ! Screened Coulomb A
        do amo=1,nvir(1)
          lla = (amo-1)*nvir(1) + 1
          ula = amo*nvir(1)
          call ga_matmul_patch('n','n',1d0,0d0,
     $            pars%g_erivv(1),1,nri,lla,ula,
     $            g_tmp1,1,nvir(1),1,nocc(1),
     $            g_tmp2,1,nri,1,nocc(1))
          do imo=1,nocc(1)
            lli = (imo-1)*nocc(1) + 1
            uli = imo*nocc(1)
            ipole = (imo-1)*nvir(1) + amo
            temp(ipole) = ga_ddot_patch(g_tmp2,'n',1,nri,1,nocc(1),
     $                      pars%g_erioo(1),'n',1,nri,lli,uli)
          enddo
        enddo
        if (ipol.gt.1) then
          call ga_copy_patch('n',vec,1+npoles1,npoles,itrial,itrial,
     $                        g_tmp1,1,nvir(2),1,nocc(2))
          do amo=1,nvir(2)
            lla = (amo-1)*nvir(2) + 1
            ula = amo*nvir(2)
            call ga_matmul_patch('n','n',1d0,0d0,
     $              pars%g_erivv(2),1,nri,lla,ula,
     $              g_tmp1,1,nvir(2),1,nocc(2),
     $              g_tmp2,1,nri,1,nocc(2))
            do imo=1,nocc(2)
              lli = (imo-1)*nocc(2) + 1
              uli = imo*nocc(2)
              ipole = (imo-1)*nvir(2) + amo + npoles1
              temp(ipole) = ga_ddot_patch(g_tmp2,'n',1,nri,1,nocc(2),
     $                      pars%g_erioo(2),'n',1,nri,lli,uli)
            enddo
          enddo
        endif
        if (pars%me.eq.0)
     $    call ga_acc(mv,1,npoles,itrial,itrial,temp,npoles,-1d0)
      enddo

      if(.not. (ga_destroy(g_tmp1(1)) .and. ga_destroy(g_tmp2(1))))
     $  call errquit('bse_davidson_kvec: ga_destroy_failed',93,GA_ERR)
      if(pars%ipol.gt.1) then
        if(.not. (ga_destroy(g_tmp1(2)) .and. ga_destroy(g_tmp2(2))))
     $  call errquit('bse_davidson_kvec: ga_destroy_failed',96,GA_ERR)
      endif

      if (.not.ga_create(mt_dbl,nri,ntrials,'g_tmp1',0,ntrials,
     $                   g_tmp1))
     $  call errquit('bse_davidson_kvec: ga create failed',102,GA_ERR)
      if (.not.ga_create(mt_dbl,npoles,ntrials,'g_tmp2',0,ntrials,
     $                   g_tmp2))
     $  call errquit('bse_davidson_kvec: ga create failed',105,GA_ERR)

      ! Screened Coulomb B
      call ga_matmul_patch('n','n',1d0,0d0,
     $     pars%g_wov(1),1,nri,1,npoles1,
     $     vec,1,npoles1,1,ntrials,
     $     g_tmp1,1,nri,1,ntrials)
      call ga_matmul_patch('t','n',-1d0,0d0,
     $     pars%g_wov(1),1,npoles1,1,nri,
     $     g_tmp1,1,nri,1,ntrials,
     $     g_tmp2,1,npoles1,1,ntrials)
      if (pars%ipol.gt.1) then
        call ga_matmul_patch('n','n',1d0,0d0,
     $     pars%g_wov(2),1,nri,1,npoles2,
     $     vec,1+npoles1,npoles,1,ntrials,
     $     g_tmp1,1,nri,1,ntrials)
        call ga_matmul_patch('t','n',-1d0,0d0,
     $     pars%g_wov(2),1,npoles2,1,nri,
     $     g_tmp1,1,nri,1,ntrials,
     $     g_tmp2,1+npoles1,npoles,1,ntrials)
      endif
      call ga_add_patch(-1d0,g_tmp2,1,npoles,1,ntrials,
     $                  1d0,mv,1,npoles,1,ntrials,
     $                      mv,1,npoles,1,ntrials)

      if(.not.(ga_destroy(g_tmp1) .and. ga_destroy(g_tmp2)))
     $  call errquit('bse_davidson_kvec: ga_destroy failed',156,GA_ERR)

      end

