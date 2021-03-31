#define USE_ALLOC 1
      subroutine gw_buildiw(pars,eia,x,wgt,iw,iwidx,nmo,nri,ngl,mynmo,
     $                      maxpoles,maxqp,ipol)
#ifdef USE_ALLOC
      USE ISO_C_BINDING
#endif
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: nmo,nri,ngl,mynmo,maxpoles,maxqp,ipol

      double precision,intent(in)    :: eia(maxpoles,*),x(ngl),wgt(ngl)
      integer,intent(in)             :: iwidx(0:mynmo,maxqp,*)
      double precision,intent(out)   :: iw(mynmo,ngl,maxqp,*)

      character(len=12),parameter :: pname = 'gw_buildiw: '

      logical ok,main,diag
      integer molo,mohi,ldpi,sizepi,l_sol,k_sol
      integer ilo,ihi,me,ld,mynpoles,igl,ipole,info
      integer pilo,pihi,k_pi,l_pi,idx,idy,l_tmp,tmp,l_vx,k_vx
      integer imo,jmo,kmo,lmo,nprocs,kOV1,kOV2,isp,addr,idz,kOO(2)
      integer kVV(2),kOV(2),mynvir,mynocc,nlocal,nmod
      double precision w2,ecrpa,ecgm
#ifdef USE_ALLOC
      real(kind=C_DOUBLE), allocatable    :: factor(:,:)
      integer(kind=c_int)         :: stat_alloc
#else
      double precision factor(maxpoles,ipol)
#endif
      double precision t0, t1

#ifdef USE_OPENMP
      integer iMaxThreads,ithread
      integer,external :: omp_get_max_threads,omp_get_thread_num
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif
#ifdef USE_ALLOC
      allocate(factor(maxpoles,ipol),STAT = stat_alloc)
      if(stat_alloc.ne.0)
     &     call errquit(' gw_buildiw: allocate failed for size',
     M     maxpoles*ipol, MA_ERR)
#endif
      call ga_sync()

      
      me = ga_nodeid()
      nprocs = ga_nnodes()
      main = me.eq.0
      pars%ecrpa = 0d0
      pars%ecgm = 0d0

      if (mod(nri,2).eq.0) then
        ldpi = nri + 1
      else
        ldpi = nri
      endif
      sizepi = (nri*(nri+1))/2

      diag = .false.
      if (nri.gt.4000) diag = .true.

      if (.not.ma_alloc_get(mt_dbl,sizepi,'Pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate Pi',0,MA_ERR) 
      if (.not.ma_alloc_get(mt_dbl,sizepi,'vx',l_vx,k_vx))
     $  call errquit(pname//'can''t allocate vx',0,MA_ERR) 
      if (.not.ma_alloc_get(mt_dbl,nri*mynmo,'temp',l_tmp,tmp))
     $  call errquit(pname//'can''t allocate tmp',0,MA_ERR) 
      if (.not.ma_alloc_get(mt_dbl,nri*mynmo,'sol',l_sol,k_sol))
     $  call errquit(pname//'can''t allocate sol',0,MA_ERR) 


      call ga_access(pars%g_eriov(1),1,nri,pars%ovlo(1),pars%ovhi(1),
     $               kOV(1),ld)
      call ga_access(pars%g_erioo(1),1,nri,pars%oolo(1),pars%oohi(1),
     $               kOO(1),ld) 
      if (pars%nvqp(1).gt.0) then
        call ga_access(pars%g_erivv(1),1,nri,pars%vvlo(1),pars%vvhi(1),
     $                 kVV(1),ld) 
      end if

      if (ipol.gt.1 .and. pars%mynpoles(2).gt.0 ) then
        call ga_access(pars%g_eriov(2),1,nri,pars%ovlo(2),pars%ovhi(2),
     $                 kOV(2),ld)
        call ga_access(pars%g_erioo(2),1,nri,pars%oolo(2),pars%oohi(2),
     $                 kOO(2),ld) 
        if (pars%nvqp(2).gt.0) then
          call ga_access(pars%g_erivv(2),1,nri,pars%vvlo(2),
     $                   pars%vvhi(2),kVV(2),ld) 
        endif
      else
        kOV(2) = 1
        kOO(2) = 1
      endif
      t0=util_wallsec()


      do igl=ngl,1,-1
c
        w2 = x(igl)

        do isp=1,ipol
!$omp   parallel do
           do ipole=1,pars%mynpoles(isp)
              if(eia(ipole,isp).le.0d0) then
                 write(6,123) ' negative eia',ipole,eia(ipole,isp),w2,
     P                (eia(ipole,isp)/(eia(ipole,isp)**2 + w2))
              endif
 123          format(i5,3e12.2)
              factor(ipole,isp) = 
     $            dsqrt(eia(ipole,isp)/(eia(ipole,isp)**2 + w2))
          enddo
!$omp   end parallel do        
        enddo

c
        if (igl.eq.1) call dfill(sizepi,0d0,dbl_mb(k_pi),1)

        call gw_buildpi('w',dbl_mb(kOV(1)),dbl_mb(kOV(2)),dbl_mb(k_pi),
     $                  factor,ldpi,nri,pars%mynpoles,main,ipol,
     $                  maxpoles)

#ifdef GWDEBUG
        if(me.eq.0) then
 1234      format(a,i3,a,i3,a,f20.3)
           write(luout,1234) ' ##  igl ',igl,' out of ',ngl,
     C          ' calling dgop at ',
     T          util_wallsec()-t0
        endif
#endif
        call ga_dgop(mt_dbl,dbl_mb(k_pi),sizepi,'+')
#ifdef GWDEBUG
        if(me.eq.0) then
           write(luout,1234) ' ##  igl ',igl,' out of ',ngl,
     C          ' called dgop at ',
     T          util_wallsec()-t0
        endif
#endif

        call dcopy(sizepi,dbl_mb(k_pi),1,dbl_mb(k_vx),1)
        !
        ! All cores invert the same matrix in order to avoid
        ! communication. This is good for small matrices but
        ! will be very slow for large ones.
        !

        if (diag) then
          call gw_invdiag(dbl_mb(k_vx),ldpi,nri,ecrpa,.true.)
          if (me.eq.0) pars%ecrpa = pars%ecrpa + wgt(igl)*ecrpa
        else
          call dpftrf('n','l',nri,dbl_mb(k_pi),info)
          if (me.eq.0) then
            call gw_cdgw_ecrpa(dbl_mb(k_pi),ldpi,nri,ecrpa)
            pars%ecrpa = pars%ecrpa + wgt(igl)*ecrpa
          endif
        endif
        !
        !
        do isp=1,ipol
          if (pars%nstates(isp).eq.0) cycle

          do imo=pars%llmo(isp),pars%llmo(isp)+pars%nqp(isp)-1
            jmo = imo - pars%llmo(isp) + 1
            

#ifdef USE_OPENMP
            call util_blas_set_num_threads(1)
#endif
            if (imo.le.pars%nocc(isp)) then

!$omp         parallel do private(kmo,lmo,idx,idy)
              do kmo=1,iwidx(0,jmo,isp)
                lmo = iwidx(kmo,jmo,isp)
                idx = tmp + (kmo-1)*nri
                if (lmo.le.pars%nocc(isp)) then
                  idy = (imo-1)*pars%nocc(isp)+lmo-1
                  idy = idy/nprocs*nri + kOO(isp)
                else
                  idy = ((imo-1)*pars%nvir(isp)+lmo-pars%nocc(isp)-1)
                  idy = idy/nprocs*nri + kOV(isp)
                endif
                call dcopy(nri,dbl_mb(idy),1,dbl_mb(idx),1)
              enddo
!$omp         end parallel do

            else

!$omp         parallel do private(kmo,lmo,idx,idy)
              do kmo=1,iwidx(0,jmo,isp)
                lmo = iwidx(kmo,jmo,isp)
                idx = tmp + (kmo-1)*nri
                if (lmo.le.pars%nocc(isp)) then
                  idy = ((lmo-1)*pars%nvir(isp)+imo-pars%nocc(isp)-1)
                  idy = idy/nprocs*nri + kOV(isp)
                else
                  idy = ((imo-pars%nocc(isp)-1)*pars%nvir(isp)+lmo-
     $                    pars%nocc(isp)-1)/nprocs*nri + kVV(isp)
                endif
                call dcopy(nri,dbl_mb(idy),1,dbl_mb(idx),1)
              enddo
!$omp         end parallel do
            endif
#ifdef USE_OPENMP
            call util_blas_set_num_threads(iMaxThreads)
#endif


            if (diag) then
              call gw_rfpmm(dbl_mb(k_vx),dbl_mb(tmp),dbl_mb(k_sol),
     $                      ldpi,nri,iwidx(0,jmo,isp),.true.)
            else
              call dcopy(nri*iwidx(0,jmo,isp),dbl_mb(tmp),1,
     $                   dbl_mb(k_sol),1)

              call dtfsm('n','l','l','n','n',nri,iwidx(0,jmo,isp),1d0,
     $                   dbl_mb(k_pi),dbl_mb(k_sol),nri)       
            endif
c
#ifdef USE_OPENMP
            call util_blas_set_num_threads(1)
#endif


!$omp       parallel do private(kmo,idx,idy)
            do kmo=1,iwidx(0,jmo,isp)
              idx = tmp + (kmo-1)*nri
              idy = k_sol + (kmo-1)*nri
              if (diag) then
                iw(kmo,igl,jmo,isp) = 
     $            ddot(nri,dbl_mb(idx),1,dbl_mb(idy),1) -
     $            ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
              else
                iw(kmo,igl,jmo,isp) = 
     $            ddot(nri,dbl_mb(idy),1,dbl_mb(idy),1) -
     $            ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
              endif
            enddo
!$omp       end parallel do


#ifdef USE_OPENMP
            call util_blas_set_num_threads(iMaxThreads)
#endif            
          enddo
c
        enddo

        if (me.eq.0) then
          call gw_cdgw_ecgm(dbl_mb(k_vx),dbl_mb(k_pi),nri,ldpi,ecgm,
     $                      ecrpa,diag)
          pars%ecgm = pars%ecgm + wgt(igl)*ecgm
          pars%ecrpa = pars%ecrpa + wgt(igl)*ecrpa
        endif

      enddo

      deallocate(factor,STAT = stat_alloc)
      if(stat_alloc.ne.0)
     &     call errquit(' gw_buildiw: deallocate failed for size',
     M     maxpoles*ipol, MA_ERR)
      if (.not.(ma_free_heap(l_sol)))
     $  call errquit(pname//'can''t free SOL heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_tmp)))
     $  call errquit(pname//'can''t free TMP heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_pi)))
     $  call errquit(pname//'can''t free PI heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_vx)))
     $  call errquit(pname//'can''t free vX heap',me,MA_ERR) 

      pars%ecrpa = pars%ecrpa*0.5d0
      pars%ecgm  = pars%ecgm*0.5d0

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif      
c
      end subroutine

