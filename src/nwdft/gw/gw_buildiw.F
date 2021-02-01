      subroutine gw_buildiw(pars,eia,x,iw,iwidx,nmo,nri,ngl,mynmo,
     $                      maxpoles,maxqp,ipol)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"      
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: nmo,nri,ngl,mynmo,maxpoles,maxqp,ipol

      double precision,intent(in)    :: eia(maxpoles,*),x(ngl)
      integer,intent(in)             :: iwidx(0:mynmo,maxqp,*)
      double precision,intent(out)   :: iw(mynmo,ngl,maxqp,*)

      character(len=12),parameter :: pname = 'gw_buildiw: '

      logical ok
      logical main
      integer molo,mohi,ldpi,sizepi,l_sol,k_sol,k_work,l_work,lwork
      integer ilo,ihi,me,ld,mynpoles,igl,ipole,info,k_ipiv,l_ipiv
      integer pilo,pihi,k_pi,l_pi,idx,idy,l_tmp,tmp
      integer imo,jmo,kmo,lmo,nprocs,kOV1,kOV2,isp,addr,idz
      double precision w2
      double precision factor(maxpoles,ipol)

      logical cholesky
      double precision,external :: ddot
      integer,external :: ga_nodeid, ga_nnodes,ilaenv

      call ga_sync()

      me = ga_nodeid()
      nprocs = ga_nnodes()
      main = me.eq.0

      if (mod(nri,2).eq.0) then
        ldpi = nri + 1
      else
        ldpi = nri
      endif
      sizepi = (nri*(nri+1))/2

      if (.not.ma_alloc_get(mt_dbl,nri**2,'Pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate Pi',0,MA_ERR) 
      if (.not.ma_alloc_get(mt_dbl,nri*mynmo,'temp',l_tmp,tmp))
     $  call errquit(pname//'can''t allocate tmp',0,MA_ERR) 
      if (.not.ma_alloc_get(mt_dbl,nri*mynmo,'sol',l_sol,k_sol))
     $  call errquit(pname//'can''t allocate sol',0,MA_ERR) 
      if (.not.ma_alloc_get(mt_int,nri,'IPIV',l_ipiv,k_ipiv))
     $  call errquit(pname//'can''t allocate IPIV',0,MA_ERR)
      lwork = ilaenv(1,'DSYTRF','L',nri,-1,-1,-1)*nri
      if (.not.ma_alloc_get(mt_dbl,lwork,'Work',l_work,k_work))
     $  call errquit(pname//'can''t allocate Work',0,MA_ERR) 

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
     $            dsqrt(eia(ipole,isp)/(eia(ipole,isp)**2 + w2))
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

        cholesky = .true.
        call dpftrf('n','l',nri,dbl_mb(k_pi),info)
        call ga_brdcst(mt_dbl,info,mitob(1),0)
        if (info.ne.0) then
          cholesky = .false.
          call gw_buildpi('r',dbl_mb(kOV1),dbl_mb(kOV2),dbl_mb(k_pi),
     $                    factor,nri,nri,pars%mynpoles,main,ipol,
     $                    maxpoles)
          call ga_dgop(mt_dbl,dbl_mb(k_pi),nri**2,'+')
          call dsytrf('l',nri,dbl_mb(k_pi),nri,int_mb(k_ipiv),
     $                 dbl_mb(k_work),lwork,info)
          if (info.gt.0) then
            call gw_corrldl(dbl_mb(k_pi),int_mb(k_ipiv),nri,main)
          else
            call errquit(pname//'LDL factorization failed',info,0)
          endif
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

            if (cholesky) then
              call dtfsm('n','l','l','n','n',nri,iwidx(0,jmo,isp),1d0,
     $                  dbl_mb(k_pi),dbl_mb(k_sol),nri)       
              addr = k_sol
            else
              call dsytrs('L',nri,iwidx(0,jmo,isp),dbl_mb(k_pi),nri,
     $                     int_mb(k_ipiv),dbl_mb(k_sol),nri,info)
              if (info.ne.0) then
                call errquit(pname//'LS solution failed',info,0)
              endif
              addr = tmp
            endif
c
            do kmo=1,iwidx(0,jmo,isp)
              idx = tmp + (kmo-1)*nri
              idy = k_sol + (kmo-1)*nri
              idz = addr + (kmo-1)*nri
              lmo = iwidx(kmo,jmo,isp)
              iw(kmo,igl,jmo,isp) = 
     $            ddot(nri,dbl_mb(idy),1,dbl_mb(idz),1) -
     $            ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
          enddo
c
        enddo
      enddo

      if (.not.(ma_free_heap(l_sol)))
     $  call errquit(pname//'can''t free SOL heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_tmp)))
     $  call errquit(pname//'can''t free TMP heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_pi)))
     $  call errquit(pname//'can''t free PI heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_ipiv)))
     $  call errquit(pname//'can''t free IPIV heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_work)))
     $  call errquit(pname//'can''t free WORK heap',me,MA_ERR) 
c
      end subroutine

      subroutine gw_corrldl(matrix,ipiv,n,main)
      implicit none
      logical,intent(in) :: main
      integer,intent(in) :: n
      double precision,dimension(n,n) :: matrix
      integer,dimension(n)            :: ipiv

      integer i,s
      double precision eig1,eig2,a,b,c

      i = 1
      do while (i.le.n)

      if (ipiv(i).lt.0) then
        s = 2
      else
        s = 1
      endif

      if ( (s.eq.1) .and. (matrix(i,i).lt.1d-6) ) then
        matrix(i,i) = 1d-6
      else if (s.eq.2) then
        a = matrix(i,i)
        b = matrix(i+1,i+1)
        c = matrix(i+1,i)
        eig1 = 0.5d0*(a+b+dsqrt((a-b)**2 + 4d0*c**2))
        eig2 = 0.5d0*(a+b-dsqrt((a-b)**2 + 4d0*c**2))
        if (main) write(*,*) eig1,eig2
        if (eig2.lt.1d-6) then
          matrix(i,i) = matrix(i,i) + 1d-6 - eig2
          matrix(i+1,i+1) = matrix(i+1,i+1) + 1d-6 - eig2
        endif
      endif

      i = i + s
      end do

      end subroutine


