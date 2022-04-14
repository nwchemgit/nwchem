c
c     == calculate scalar x2c contribution ==
      subroutine calc_x2c_w_scalar  (
     &   ao_bas_han,  ! in: AO basis handle 
     &   geom,        ! in: geometry handle 
     &   ipol,        ! in: nr. of spin components in g_dens
     &   g_dens,      ! in: atomic superposition density
     &   skipnuc,     ! in: whether to exclude nuclear potential
     &   skipcoul,    ! in: whether to skip atomic Hartree pot.
     &   delchi_ao,   ! in: derivatives of basis fcts.
     &   qxyz,        ! in: grid points 
     &   qwght,       ! in: weight coefficients
     &   nbf,         ! in: nr. basis functions
     &   npts,        ! in: nr. grid points
     &   natoms,      ! in: nr. atoms
     &   x2c0)        ! out: W0

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
cinclude "x2c.fh"
#include "global.fh"
#include "bas.fh"
c
      integer nbf,npts,ao_bas_han,natoms,geom
      logical skipnuc, skipcoul
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision pot(npts)
      double precision qxyz(3,npts)
      double precision delchi_ao(npts,3,nbf)

      double precision x2c0(nbf,nbf)

      integer i,j,k
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot
      double precision denom
      double precision fac1_arr(npts)
      double precision ac_sf 
      external get_ints_x2c_scalar

      integer dbg

c     =================================================================

      dbg = 0

c     debug:
      if (ga_nodeid().eq.0 .and. dbg>0) then
        write (luout,*) 'hello from calc_x2c_w_scalar'
        write (luout,*) 'skipnuc, skipcoul = ',skipnuc, skipcoul
      end if

c     == preliminaries ==
      do ipt = 1,npts
        amat_coul(ipt,1) = 0.d0            
        if (ipol.gt.1) amat_coul(ipt,2) = 0.d0            
        amat_nucl(ipt) = 0.d0
        closegridpts(ipt) = 0
      end do

      tol = 1d-8


c     calculate sum of atomic Hartree potentials, if
c     the skipcoul option is not set:

      if (.not.skipcoul) then
        do i=1,ipol
          call potential_list(ao_bas_han, g_dens(i), npts, qxyz, 
     &       amat_coul(1,i), tol)
        enddo
        if (ipol.gt.1) then     ! enter-ONLY-if ipol>1
          do k = 1,npts
            amat_coul(k,1)=amat_coul(k,1)+amat_coul(k,2)
          enddo
        endif
        
c       write (luout,*) 'x2chscale =', x2chscale
        
c       scale Hartree potential globally if requested by input
c$$$        if (scale_x2c_hartree) then
c$$$          do k = 1,npts
c$$$            amat_coul(k,1)=amat_coul(k,1) * x2chscale
c$$$          end do
c$$$        end if

      else
        if (ga_nodeid().eq.0 .and. dbg>0) 
     &     write (luout,*) 'calc_x2c_w_sc: skipped Coul. pot'
      end if                    ! .not.skipcoul


c     calculate the total point nuclear potential on the grid
c     if the option skipnuc is not set

      if (.not.skipnuc) then
        call gridNuclearPotentialPoint(geom,natoms,npts,qxyz,qwght,
     &     closegridpts,amat_nucl)
        do k = 1,npts
          if (k.eq.closegridpts(k)) qwght(k) = 0.d0
        end do
      else
        if (ga_nodeid().eq.0 .and. dbg>0) 
     &     write (luout,*) 'calc_x2c_w_sc: skipped nuc. pot'
      end if ! .not.skipnuc

c     === define fac1_arr
      do k = 1,npts
c      == assemble hartree and nuclear contributions     ==
c      == and pre-multiply with grid point weight factor ==
       amat_tot =-amat_coul(k,1)+amat_nucl(k)
       fac1_arr(k)=amat_tot * qwght(k)
      end do

c     scale Hartree potential if requested by input
c     (caution, it may have been scaled already once, see above)
c$$$      if (scale_x2c_so_hartree) then
c$$$        write (luout,*) 'x2csoscale =', x2csoscale
c$$$        do k = 1,npts
c$$$          amat_coul(k,1)=amat_coul(k,1) * x2csoscale
c$$$        end do
c$$$      end if


c     == assemble x2c correction ==
c ----- main diagonal --- START
      do i = 1, nbf
          j=i
          call get_ints_x2c_scalar(nbf,npts,delchi_ao,i,j,      
     &                          fac1_arr,   
     &                          ac_sf) ! out 
          x2c0(i,j) = x2c0(i,j) + ac_sf
      enddo ! end-loop-i
c ----- main diagonal --- END
c ----- off diagonal --- START
      do i = 1, nbf
        do j = i+1, nbf
          call get_ints_x2c_scalar(nbf,npts,delchi_ao,i,j,      
     &                          fac1_arr,
     &                          ac_sf) ! out 
          x2c0(i,j) = x2c0(i,j) + 2.0d0*ac_sf
        enddo ! end-loop-j
      enddo ! end-loop-i
c ----- off diagonal --- END

      if (ga_nodeid().eq.0 .and. dbg>0) then
        write (luout,*) 'leaving calc_x2c_w_scalar'
      end if

      return
      end

      subroutine get_ints_x2c_scalar(
     &   nbf,                   ! in: # basis fncts.
     &   npts,           ! in: grid points 
     &   delchi_ao,      ! in: deriv. of basis fcts.
     &   i,j,            ! in: (i,j) indices for delchi_ao
     &   fac1_arr,       ! in
     &   ac_sf)          ! out
c
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts)
      double precision ac_sf
      double precision prod0

      ac_sf    = 0.0d0
      
      do k = 1, npts
        prod0 = delchi_ao(k,1,i)*delchi_ao(k,1,j)
     &     +delchi_ao(k,2,i)*delchi_ao(k,2,j) 
     &     +delchi_ao(k,3,i)*delchi_ao(k,3,j)  
        
        ac_sf    = ac_sf    + fac1_arr(k)*prod0
        
      end do                    ! k
      return
      end

