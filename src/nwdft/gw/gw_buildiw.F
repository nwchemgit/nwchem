#define USE_ALLOC 1
      subroutine gw_buildiw(pars,eia,x,iw,iwidx,nmo,nri,ngl,mynmo,
     $     maxpoles,maxqp,ipol)
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

      double precision,intent(in)    :: eia(maxpoles,*),x(ngl)
      integer,intent(in)             :: iwidx(0:mynmo,maxqp,*)
      double precision,intent(out)   :: iw(mynmo,ngl,maxqp,*)

      character(len=12),parameter :: pname = 'gw_buildiw: '

      logical ok
      logical main
      integer molo,mohi,ldpi,sizepi,l_sol,k_sol,k_work,l_work,lwork
      integer ilo,ihi,me,ld,mynpoles,igl,ipole,info,k_ipiv,l_ipiv
      integer pilo,pihi,k_pi,l_pi,idx,idy,l_tmp,tmp
      integer imo,jmo,kmo,lmo,nprocs,kOV1,kOV2,isp,addr,idz,kOO(2)
      integer kVV(2),kOV(2),mynvir,mynocc,nlocal,nmod
      double precision w2
#ifdef USE_ALLOC
      real(kind=C_DOUBLE), allocatable    :: factor(:,:)
      integer(kind=c_int)         :: stat_alloc
#else
      double precision factor(maxpoles,ipol)
#endif
      logical cholesky
      double precision t0, t1
      integer,external :: ilaenv

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
      write(6,*)  ' stat_alloc ',stat_alloc
#endif
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
     $               kOV(1),ld)
      call ga_access(pars%g_erioo(1),1,nri,pars%oolo(1),pars%oohi(1),
     $               kOO(1),ld) 
      if (pars%nvqp(1).gt.0) then
        call ga_access(pars%g_erivv(1),1,nri,pars%vvlo(1),pars%vvhi(1),
     $                 kVV(1),ld) 
      end if

      if (ipol.gt.1) then
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
        call gw_buildpi('w',dbl_mb(kOV(1)),dbl_mb(kOV(2)),dbl_mb(k_pi),
     $                  factor,ldpi,nri,pars%mynpoles,main,ipol,
     $                  maxpoles)
        if(me.eq.0) then
 1234      format(a,i3,a,i3,a,f20.3)
           write(luout,1234) ' ##  igl ',igl,' out of ',ngl,
     C          ' calling dgop at ',
     T          util_wallsec()-t0
        endif
        call ga_dgop(mt_dbl,dbl_mb(k_pi),sizepi,'+')
        if(me.eq.0) then
           write(luout,1234) ' ##  igl ',igl,' out of ',ngl,
     C          ' called dgop at ',
     T          util_wallsec()-t0
        endif
        !
        ! All cores invert the same matrix in order to avoid
        ! communication. This is good for small matrices but
        ! will be very slow for large ones.
        !

        cholesky = .true.
        call dpftrf('n','l',nri,dbl_mb(k_pi),info)
        !
        do isp=1,ipol

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


            call dcopy(nri*iwidx(0,jmo,isp),dbl_mb(tmp),1,
     $                 dbl_mb(k_sol),1)

            call dtfsm('n','l','l','n','n',nri,iwidx(0,jmo,isp),1d0,
     $                 dbl_mb(k_pi),dbl_mb(k_sol),nri)       
c
#ifdef USE_OPENMP
            call util_blas_set_num_threads(1)
#endif


!$omp       parallel do private(kmo,idx,idy)
            do kmo=1,iwidx(0,jmo,isp)
              idx = tmp + (kmo-1)*nri
              idy = k_sol + (kmo-1)*nri
              iw(kmo,igl,jmo,isp) = 
     $            ddot(nri,dbl_mb(idy),1,dbl_mb(idy),1) -
     $            ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
!$omp       end parallel do


#ifdef USE_OPENMP
            call util_blas_set_num_threads(iMaxThreads)
#endif            
          enddo
c
        enddo
      enddo
      deallocate(factor,STAT = stat_alloc)
      if(stat_alloc.ne.0)
     &     call errquit(' gw_buildiw: deallocate failed for size',
     M     maxpoles*ipol, MA_ERR)
      write(6,*)  ' stat_alloc ',stat_alloc 
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

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif      
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


