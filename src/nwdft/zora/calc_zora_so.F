c
c     == calculate spin-orbit zora contribution ==
      subroutine calc_zora_so(ao_bas_han,  ! in: AO basis handle 
     &                           geom,        ! in: geometry handle 
     &                           ipol,        ! in: nr. of polarizations 
     &                           g_dens,      ! in: superposit. atomic densities 
     &                           delchi_ao,   ! in: deriv. of basis functions
     &                           qxyz,        ! in: grid points 
     &                           qwght,       ! in: weighting coeffs.
     &                           nbf,         ! in: nr. basis functions
     &                           npts,        ! in: nr. grid points
     &                           natoms,      ! in: nr. atoms 
     &                           zora0, scal0,! out  
     &                           zorax, scalx,! out 
     &                           zoray, scaly,! out 
     &                           zoraz, scalz,! out
     &                           ofinite,    ! in: =.true. if finite Nuclear Model (Gaussian)
     &                           zetanuc_arr,! in: zetanuc(i) i=1,natoms
     &                           Knucl,      ! in: = .true. if including ONLY nuclear part in K ZORA
     &                           use_modelpotential,gexpo,gcoef)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "zora.fh"
#include "global.fh"
#include "bas.fh"
c
      integer nbf,npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision pot(npts)
      double precision qxyz(3,npts)
      double precision delchi_ao(npts,3,nbf)
      double precision zora0(nbf,nbf)  ! scalar component
      double precision scal0(nbf,nbf)  ! scalar component
      double precision zorax(nbf,nbf)
      double precision scalx(nbf,nbf)
      double precision zoray(nbf,nbf)
      double precision scaly(nbf,nbf)
      double precision zoraz(nbf,nbf)
      double precision scalz(nbf,nbf)
      integer i,j,k
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot
      double precision denom
      double precision fac1_arr(npts),fac2_arr(npts)
      logical use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
      logical ofinite,Knucl
      double precision zetanuc_arr(natoms)
      double precision ac_sf,ac_scl,      ! scalar part
     &                 ac_sx,ac_scl_x,    ! x term
     &                 ac_sy,ac_scl_y,    ! y term
     &                 ac_sz,ac_scl_z     ! z term
      logical do_so
      external get_ints_zora
c
c     == calculate zora factors ==
      call calc_zora_facs(ao_bas_han,geom,
     &                     ipol,g_dens,natoms,npts,qxyz,qwght,
     &                     use_modelpotential,gexpo,gcoef,
     &                     ofinite,
     &                     zetanuc_arr,
     &                     Knucl,
     &                     fac1_arr,
     &                     fac2_arr)
c
c     == assemble zora correction ==
c     diagonal
      do_so = .true.  ! for spin-orbit zora
      do i = 1, nbf
          j=i
          call get_ints_zora(do_so,nbf,npts,delchi_ao,i,j,      
     &                          fac1_arr,fac2_arr,       
     &                          ac_sf,ac_scl,   ! out
     &                          ac_sx,ac_scl_x, ! out 
     &                          ac_sy,ac_scl_y, ! out 
     &                          ac_sz,ac_scl_z) ! out 
c
          zora0(i,j) = zora0(i,j) + ac_sf
          scal0(i,j) = scal0(i,j) + ac_scl
          zorax(i,j) = zorax(i,j) + ac_sx
          scalx(i,j) = scalx(i,j) + ac_scl_x
          zoray(i,j) = zoray(i,j) + ac_sy
          scaly(i,j) = scaly(i,j) + ac_scl_y
          zoraz(i,j) = zoraz(i,j) + ac_sz
          scalz(i,j) = scalz(i,j) + ac_scl_z
      enddo ! end-loop-i
c
c     off diagonal
      do i = 1, nbf
        do j = i+1, nbf
          call get_ints_zora(do_so,nbf,npts,delchi_ao,i,j,      
     &                          fac1_arr,fac2_arr,       
     &                          ac_sf,ac_scl,   ! out
     &                          ac_sx,ac_scl_x, ! out 
     &                          ac_sy,ac_scl_y, ! out 
     &                          ac_sz,ac_scl_z) ! out 
c
          zora0(i,j) = zora0(i,j) + 2.0d0*ac_sf
          scal0(i,j) = scal0(i,j) + 2.0d0*ac_scl
          zorax(i,j) = zorax(i,j) + 2.0d0*ac_sx
          scalx(i,j) = scalx(i,j) + 2.0d0*ac_scl_x
          zoray(i,j) = zoray(i,j) + 2.0d0*ac_sy
          scaly(i,j) = scaly(i,j) + 2.0d0*ac_scl_y
          zoraz(i,j) = zoraz(i,j) + 2.0d0*ac_sz
          scalz(i,j) = scalz(i,j) + 2.0d0*ac_scl_z
        enddo ! end-loop-j
      enddo ! end-loop-i
c
      return
      end
