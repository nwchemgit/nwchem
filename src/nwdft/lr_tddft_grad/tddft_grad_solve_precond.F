      subroutine tddft_grad_solve_precond(g_x,shift)
c
c     $Id$
c
c     This subroutine applies a preconditioner to the new trial vectors
c     stored in the global array g_x. The purpose of this is to 
c     enhance the convergence of the linear system solver in the 
c     optimisation of Ria.
c
c     Because of the way the linear system solver works the data in
c     the global array is stored in a vector format:
c
c        g_x(1:nocca*nvirta+noccb*nvirtb,1:nroots)
c
c     The preconditioner used is simply to divide by the orbital 
c     energy difference:
c
c        g_x(ia,iroot) = g_x(ia,iroot)/(e(a)-e(i))
c
c     This is simply a guess for a reasonable preconditioner. So
c     perhaps there is a better option.
c
c     The parallelisation is driven off the distribution of the
c     global array.
c
c     Written by Huub van Dam, March 2007.
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tddft_grad_solve_common.fh"
c
c     Input:
c
      double precision shift ! a level shift
c
c     Input/Output:
c
      integer g_x ! the new trial vectors
c
c     Local:
c
      integer ir  ! counter running over roots
      integer ioc ! counter over active occupied orbitals
      integer ivt ! counter over active virtual orbitals
      integer iov ! counter over occupied-virtual orbital pairs
c
      integer ilo,  ihi  ! limits on CPHF vector block
      integer irlo, irhi ! limits on roots block
      integer k_x        ! index for the local CPHF block
      integer ld_x       ! leading dimension for the local CPHF block
c
      logical oalpha ! .true. if there is an alpha orbital part
      logical obeta  ! .true. if there is an beta orbital part
c
      integer ialo, iahi ! limits on alpha CPHF vector block
      integer iblo, ibhi ! limits on beta CPHF vector block
c
      integer iav_start ! the first local alpha virtual orbital
      integer iao_start ! the first local alpha occupied orbital
      integer ibv_start ! the first local beta virtual orbital
      integer ibo_start ! the first local beta occupied orbital
c
      integer iproc     ! my processor rank
      logical doit
c
      call ga_sync()
      iproc = ga_nodeid()
      call ga_distribution(g_x,iproc,ilo,ihi,irlo,irhi)
      doit = .not. (ilo.eq.0.and.ihi.eq.-1.and.
     &     irlo.eq.0.and.irhi.eq.-1)
      if (doit) then
       call ga_access(g_x,ilo,ihi,irlo,irhi,k_x,ld_x)
c
c      Do we have alpha or beta orbitals or both?
c
       oalpha = ilo.le.tdg_naoc(1)*tdg_nav(1)
       obeta  = ihi.gt.tdg_naoc(1)*tdg_nav(1)
c
c      If we have alpha orbitals then work on those elements
c
       ialo =  0
       iahi = -1
       if (oalpha) then
c
c       First we work out the limits on the loops
c
        ialo      = ilo
        iahi      = min(tdg_naoc(1)*tdg_nav(1),ihi)
        iav_start = (ialo+tdg_naoc(1)-1)/tdg_naoc(1)
        iao_start = mod(ialo-1,tdg_naoc(1))+1
c
c       Loop over the roots
c
        do ir = irlo, irhi
          ioc = iao_start
          ivt = iav_start
          iov = 0
c
c         Loop over occupied (fastest) and virtual orbitals
c
          do while (iov.le.iahi-ialo)
            do while (iov.le.iahi-ialo.and.ioc.le.tdg_naoc(1))
              dbl_mb(k_x+(ir-1)*ld_x+iov) 
     &        = dbl_mb(k_x+(ir-1)*ld_x+iov)/
     &         (dbl_mb(tdg_e_mo+tdg_nocc(1)+ivt-1)-
     &          dbl_mb(tdg_e_mo+tdg_nfc(1)+ioc-1))
              ioc = ioc + 1
              iov = iov + 1
            enddo
            ioc = 1
            ivt = ivt + 1
          enddo
        enddo
      endif
      if (obeta) then
        iblo      = max(ilo-tdg_naoc(1)*tdg_nav(1),1)
        ibhi      = ihi-tdg_naoc(1)*tdg_nav(1)
        ibv_start = (iblo+tdg_naoc(2)-1)/tdg_naoc(2)
        ibo_start = mod(iblo-1,tdg_naoc(2))+1
        do ir = irlo, irhi
          ioc = ibo_start
          ivt = ibv_start
          iov = iahi-ialo+1  ! the first after the alpha spin elements
          do while (iov.le.iahi-ialo+1+ibhi-iblo)
            do while (iov.le.iahi-ialo+1+ibhi-iblo.and.
     &                ioc.le.tdg_naoc(2))
              dbl_mb(k_x+(ir-1)*ld_x+iov) 
     &        = dbl_mb(k_x+(ir-1)*ld_x+iov)/
     &         (dbl_mb(tdg_e_mo+tdg_nmo(1)+tdg_nocc(2)+ivt-1)-
     &          dbl_mb(tdg_e_mo+tdg_nmo(1)+tdg_nfc(2)+ioc-1))
              ioc = ioc + 1
              iov = iov + 1
            enddo
            ioc = 1
            ivt = ivt + 1
          enddo
        enddo
      endif
      call ga_release(g_x,ilo,ihi,irlo,irhi)
      endif ! doit
c
      call ga_sync
c
      return
      end
