      subroutine gw_buildiw_minres(pars,eia,x,iw,iwidx,nmo,nri,ngl,
     $                      mynmo,maxpoles,maxqp,ipol)
      use iso_c_binding
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"      
#include "util.fh"
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
      integer kVV(2),kOV(2),ga,root,stat_alloc
      double precision w2,timer
      real(kind=c_double) xvec(nri),bvec(nri)

      real(kind=c_double), allocatable, dimension(:,:,:) :: factor

      logical cholesky
      double precision,external :: ddot
      integer,external :: ga_nodeid, ga_nnodes,ilaenv
      logical,external :: ga_locate

#ifdef USE_OPENMP
      integer iMaxThreads
      integer,external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif

      call ga_sync()

      allocate(factor(maxpoles,ipol,ngl),stat=stat_alloc)
      if (stat_alloc.ne.0) then
        call errquit('gw_buildiw_minres: allocation failed',0,0)
      endif

      me = ga_nodeid()
      nprocs = ga_nnodes()
      main = me.eq.0

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
        kOV(2) = kOV(1)
        kOO(2) = kOO(1)
      endif

      do igl=1,ngl
        do isp=1,ipol
          do ipole=1,pars%mynpoles(isp)
            factor(ipole,isp,igl) = 
     $        eia(ipole,isp)/(eia(ipole,isp)**2 + x(igl))
          enddo
        enddo
      enddo

      do isp=1,ipol
        do imo=pars%llmo(isp),pars%llmo(isp)+pars%nqp(isp)-1
          jmo = imo - pars%llmo(isp) + 1
          lmo = 0

          do kmo=1,nmo

            if (imo.le.pars%nocc(isp)) then
              if (kmo.le.pars%nocc(isp)) then
                root = mod((imo-1)*pars%nocc(isp)+kmo-1,nprocs) 
                idx  = ((imo-1)*pars%nocc(isp)+kmo-1)/nprocs+
     $                 pars%oolo(isp)
                ga = pars%g_erioo(isp)
              else
                root = mod((imo-1)*pars%nvir(isp)+kmo-pars%nocc(isp)-1,
     $                     nprocs)           
                idx  = ((imo-1)*pars%nvir(isp)+kmo-pars%nocc(isp)-1)/
     $                 nprocs + pars%ovlo(isp)           
                ga = pars%g_eriov(isp)
              endif
            else
              if (kmo.le.pars%nocc(isp)) then
                root = mod((kmo-1)*pars%nvir(isp)+imo-pars%nocc(isp)-1,
     $                 nprocs)           
                idx  = ((kmo-1)*pars%nvir(isp)+imo-pars%nocc(isp)-1)/
     $                 nprocs + pars%ovlo(isp)           
                ga = pars%g_eriov(isp)
              else
                root = mod((imo-pars%nocc(isp)-1)*pars%nvir(isp)+
     $                 kmo-pars%nocc(isp)-1,nprocs)
                idx  =((imo-pars%nocc(isp)-1)*pars%nvir(isp)+
     $                 kmo-pars%nocc(isp)-1)/nprocs + pars%vvlo(isp)
                ga = pars%g_erivv(isp)
              endif
            endif


            if (me.eq.root) call ga_get(ga,1,nri,idx,idx,bvec,nri)
            call ga_sync()
            call ga_brdcst(mt_dbl,bvec,ma_sizeof(mt_dbl,nri,mt_byte),
     $                     root)

            xvec(:) = bvec(:)
            timer = util_wallsec()
            call gw_minres(dbl_mb(kOV(1)),dbl_mb(kOV(2)),
     $                     factor(1,1,ngl),xvec,bvec,pars%mynpoles,
     $                     maxpoles,nri,ipol)

#ifdef GWDEBUG
            if (main) write(luout,9000) 200,jmo,kmo,util_wallsec()-timer
#endif
            if (root.eq.me) then
              lmo = lmo + 1
              iw(lmo,ngl,jmo,isp) = 
     $                 dot_product(xvec,bvec) - dot_product(bvec,bvec)
            endif
            

            do igl=ngl-1,1,-1
              timer = util_wallsec()
              call gw_minres(dbl_mb(kOV(1)),dbl_mb(kOV(2)),
     $                       factor(1,1,igl),xvec,bvec,pars%mynpoles,
     $                       maxpoles,nri,ipol)

#ifdef GWDEBUG
              if (main) write(luout,9000) igl,jmo,kmo,
     $                                    util_wallsec()-timer
#endif
              if (root.eq.me) then
                iw(lmo,igl,jmo,isp) = 
     $                 dot_product(xvec,bvec) - dot_product(bvec,bvec)
              endif

            enddo
          enddo
        enddo
      enddo

 9000 format(' ## igl ',I3,' nqp ',I4,' imo ',I4,' time ',F8.3)


#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif      
c
      deallocate(factor)

      end subroutine

