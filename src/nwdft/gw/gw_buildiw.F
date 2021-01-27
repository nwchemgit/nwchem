      subroutine gw_buildiw(pars,eia,x,iw,iwidx,nmo,nri,ngl,mynmo,
     $                      maxpoles,maxqp,ipol)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: nmo,nri,ngl,mynmo,maxpoles,maxqp,ipol

      double precision,intent(in)    :: eia(maxpoles,ipol),x(ngl)
      integer,intent(in)             :: iwidx(0:mynmo,maxqp,ipol)
      double precision,intent(out)   :: iw(mynmo,ngl,maxqp,ipol)

      character(len=12),parameter :: pname = 'gw_buildiw: '

      logical ok
      logical main
      integer molo,mohi,ldpi,sizepi,l_sol,k_sol
      integer ilo,ihi,me,ld,mynpoles,igl,ipole,info
      integer pilo,pihi,k_pi,l_pi,idx,idy,l_tmp,tmp
      integer imo,jmo,kmo,lmo,nprocs,kOV1,kOV2,isp
      double precision w2
      double precision factor(maxpoles,ipol)

      double precision,external :: ddot
      integer,external :: ga_nodeid, ga_nnodes

      me = ga_nodeid()
      nprocs = ga_nnodes()
      main = me.eq.0

      if (mod(nri,2).eq.0) then
        ldpi = nri + 1
      else
        ldpi = nri
      endif
      sizepi = (nri*(nri+1))/2

#ifdef PARALLEL_DIAG
#endif      

      if (.not.ma_push_get(mt_dbl,sizepi,'Pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate Pi',0,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,nri*mynmo,'temp',l_tmp,tmp))
     $  call errquit(pname//'can''t allocate tmp',0,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,nri*mynmo,'sol',l_sol,k_sol))
     $  call errquit(pname//'can''t allocate sol',0,MA_ERR) 


      call ga_access(pars%g_eriov(1),1,nri,pars%ovlo(1),pars%ovhi(1),
     $               kOV1,ld)
      if (ipol.gt.1) then
        call ga_access(pars%g_eriov(2),1,nri,pars%ovlo(2),pars%ovhi(2),
     $                 kOV2,ld)
      else
        kOV2 = 1
      endif

      do igl=1,ngl
c
        w2 = x(igl)
        do isp=1,ipol
          do ipole=1,pars%mynpoles(isp)
            factor(ipole,isp) = 
     $          dsqrt(eia(ipole,isp)/(eia(ipole,isp)**2 + w2))
          enddo
        enddo
c
        call gw_buildpi('w',dbl_mb(kOV1),dbl_mb(kOV2),dbl_mb(k_pi),
     $                  factor,ldpi,nri,pars%mynpoles,main,ipol,
     $                  maxpoles)
        call ga_dgop(mt_dbl,dbl_mb(k_pi),sizepi,'+')
        !
        ! All cores invert the same matrix in order to avoid
        ! communication. This is good for small matrices but
        ! will be very slow for large ones.
        !

        call dpftrf('n','l',nri,dbl_mb(k_pi),info)
        if (info.ne.0) then
          call errquit(pname//'Cholesky decomposition failed',0,0)
          !call dft_invdiag
        endif
        !
        do isp=1,ipol

          do imo=pars%llmo(isp),pars%llmo(isp)+pars%nqp(isp)-1
            jmo = imo - pars%llmo(isp) + 1
            if (imo.le.pars%nocc(isp)) then

              do kmo=1,iwidx(0,jmo,isp)
                lmo = iwidx(kmo,jmo,isp)
                idx = tmp + (kmo-1)*nri
                if (lmo.le.pars%nocc(isp)) then
                  idy = ((imo-1)*pars%nocc(isp)+lmo-1)/nprocs + 
     $                  pars%oolo(isp)
                  call ga_get(pars%g_erioo(isp),1,nri,idy,idy,
     $                 dbl_mb(idx),nri)
                else
                  idy = ((imo-1)*pars%nvir(isp)+lmo-
     $                   pars%nocc(isp)-1)/nprocs + pars%ovlo(isp)
                  call ga_get(pars%g_eriov(isp),1,nri,idy,idy,
     $                  dbl_mb(idx),nri)
                endif
              enddo


            else

              do kmo=1,iwidx(0,jmo,isp)
                lmo = iwidx(kmo,jmo,isp)
                idx = tmp + (kmo-1)*nri
                if (lmo.le.pars%nocc(isp)) then
                  idy = ((lmo-1)*pars%nvir(isp)+imo-
     $                   pars%nocc(isp)-1)/nprocs + pars%ovlo(isp)
                  call ga_get(pars%g_eriov(isp),1,nri,idy,idy,
     $                   dbl_mb(idx),nri)
                else
                  idy = ((imo-pars%nocc(isp)-1)*pars%nvir(isp)+lmo-
     $                    pars%nocc(isp)-1)/nprocs + pars%vvlo(isp)
                  call ga_get(pars%g_erivv(isp),1,nri,idy,idy,
     $                    dbl_mb(idx),nri)
                endif
              enddo

            endif

            call dcopy(nri*iwidx(0,jmo,isp),dbl_mb(tmp),1,
     $                 dbl_mb(k_sol),1)

            call dpftrs('n','l',nri,iwidx(0,jmo,isp),dbl_mb(k_pi),
     $                  dbl_mb(k_sol),nri,info)
c
            if (info.ne.0) then
              call errquit(pname//' triangular solve failed',0,0)
            endif
c
            do kmo=1,iwidx(0,jmo,isp)
              idx = tmp + (kmo-1)*nri
              idy = k_sol + (kmo-1)*nri
              lmo = iwidx(kmo,jmo,isp)
              iw(kmo,igl,jmo,isp) = 
     $           ddot(nri,dbl_mb(idx),1,dbl_mb(idy),1) -
     $           ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
          enddo
c
        enddo
      enddo

      if (.not.ma_chop_stack(l_pi))
     $  call errquit(pname//'can''t pop stack',0,MA_ERR) 
c
      end subroutine
