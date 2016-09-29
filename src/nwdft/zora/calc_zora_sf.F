c
c     == calculate spin-free zora contribution ==
      subroutine calc_zora_sf(ao_bas_han,            ! in: AO basis handle 
     &                        geom,                  ! in: geometry handle 
     &                        ipol,                  ! in: nr. of polarizations 
     &                        g_dens,                ! in: superposit. atomic densities 
     &                        chi_ao,                ! in: basis functions
     &                        delchi_ao,             ! in: deriv. of basis functions
     &                        qxyz,                  ! in: grid points 
     &                        qwght,                 ! in: weighting coeffs.
     &                        nbf,                   ! in: nr. basis functions
     &                        npts,                  ! in: nr. grid points
     &                        natoms,                ! in: nr. atoms 
     &                        zorasf,                ! out
     &                        zorascal,              ! out
     &                        ofinite,               ! in: =.true. if finite Nuclear Model (Gaussian)
     &                        zetanuc_arr,           ! in: zetanuc(i) i=1,natoms
     &                        Knucl,                 ! in: = .true. if including ONLY nuclear part in K ZORA  
     &                        use_modelpotential,    ! in
     &                        gexpo,                 ! in
     &                        gcoef)                 ! in
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
      integer nbf, npts, ao_bas_han, natoms, geom
      integer g_dens(2),ipol
      double precision chi_ao(npts,nbf),qwght(npts),pot(npts)
      double precision qxyz(3,npts)
      double precision delchi_ao(npts,3,nbf)
      double precision zorasf(nbf,nbf)
      double precision zorascal(nbf,nbf)
      integer i,j,k
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot
      double precision denom
      double precision fac1_arr(npts),fac2_arr(npts)
      double precision ac_sf,ac_scl,prod
      logical ofinite,Knucl,use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
      double precision tmp
      logical do_so
      double precision zetanuc_arr(natoms)
c
c     == calculate zora factors ==
      call calc_zora_facs(ao_bas_han,geom,
     &                    ipol,g_dens,natoms,npts,qxyz,qwght,
     &                    use_modelpotential,gexpo,gcoef,
     &                    ofinite,zetanuc_arr,Knucl,
     &                    fac1_arr,fac2_arr)
c
c     == assemble zora correction ==
c     main diagonal
      do_so = .false.  ! for scalar zora
      tmp = 0.d0
      do i = 1, nbf
         j=i
         call get_ints_zora(do_so,nbf,npts,delchi_ao,i,j,          
     &                   fac1_arr,fac2_arr,      
     &                   ac_sf,ac_scl, !out
     &                   tmp, tmp, tmp, tmp, tmp, tmp) 
c
         zorasf(i,j)   = zorasf(i,j)   + ac_sf
         zorascal(i,j) = zorascal(i,j) + ac_scl
c
      enddo ! end-loop-i
c
c     off diagonal
      do_so = .false.  ! for scalar zora
      tmp = 0.d0
      do i = 1, nbf
        do j = i+1, nbf
         call get_ints_zora(do_so,nbf,npts,delchi_ao,i,j,          
     &                         fac1_arr,fac2_arr,      
     &                         ac_sf,ac_scl,  ! out
     &                   tmp, tmp, tmp, tmp, tmp, tmp) 
c
         zorasf(i,j)   = zorasf(i,j)   + 2.0d0*ac_sf
         zorascal(i,j) = zorascal(i,j) + 2.0d0*ac_scl
c
        enddo ! end-loop-j
      enddo ! end-loop-i
c
      return
      end
