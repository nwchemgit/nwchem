c     +++++++++++++++++++++++++
c     + calculate all EPR AOs +
c     +++++++++++++++++++++++++
c     1. (scalx1,scaly1,scalz1):
c     H^{01,EPR}_{mu nu,t}=< dchi_mu/dr_m| K | dchi_nu/dr_n >-
c                          < dchi_mu/dr_n| K | dchi_nu/dr_m >
c
c     2. (scalxx2,scalxy2,scalxz2,
c         scalyx2,scalyy2,scalyz2,
c         scalzx2,scalzy2,scalzz2):
c     A_{mu nu,at}        =< dchi_mu/dr_m| K r_a| dchi_nu/dr_n >-
c                          < dchi_mu/dr_n| K r_a| dchi_nu/dr_m >
c     a=1,2,3 (tmn)=123,231,312
c
c     3. (scalx3,scaly3,scalz3):
c     C_{mu nu,k}         =< chi_mu| K | dchi_nu/dr_k >
c     k=1,2,3
c
c     4. (scalxx4,scalxy4,scalxz4,
c         scalyx4,scalyy4,scalyz4,
c         scalzx4,scalzy4,scalzz4):
c     B_{mu nu,kt}        =< chi_mu| K r_t| dchi_nu/dr_k >
c     t=1,2,3 k=1,2,3
c     
c     (tmn)=123,231,312=xyz,yzx,zxy
      subroutine calc_zora_EPR(ao_bas_han,   ! in: AO basis handle 
     &                         geom,         ! in: geometry handle 
     &                         ipol,         ! in: nr. of polarizations 
     &                         g_dens,       ! in: superposit. atomic densities 
     &                         chi_ao,       ! in:           basis functions
     &                         delchi_ao,    ! in: deriv. of basis functions
     &                         qxyz,         ! in: grid points 
     &                         qwght,        ! in: weighting coeffs.
     &                         nbf,          ! in: nr. basis functions
     &                         npts,         ! in: nr. grid points
     &                         natoms,       ! in: nr. atoms 
     &                         chi_cntr,     ! in: basis centers
     &                         Knucl,        ! in: flag =T do Kzora with Nucl pot ONLY
     C     mbf,ibf,
c     do we need atoms definition from  grid_nbfm ??
     &                         use_modelpotential,
     &                         gexpo,
     &                         gcoef,
     &                         scalx1,       ! out 
     &                         scaly1,       ! out 
     &                         scalz1,       ! out
     &                         scalxx2,      ! out 
     &                         scalxy2,      ! out 
     &                         scalxz2,      ! out 
     &                         scalyx2,      ! out 
     &                         scalyy2,      ! out 
     &                         scalyz2,      ! out 
     &                         scalzx2,      ! out 
     &                         scalzy2,      ! out 
     &                         scalzz2,      ! out 
     &                         scalx3,       ! out 
     &                         scaly3,       ! out 
     &                         scalz3,       ! out
     &                         scalxx4,      ! out 
     &                         scalxy4,      ! out 
     &                         scalxz4,      ! out 
     &                         scalyx4,      ! out 
     &                         scalyy4,      ! out 
     &                         scalyz4,      ! out 
     &                         scalzx4,      ! out 
     &                         scalzy4,      ! out 
     &                         scalzz4)       ! out          
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "bas.fh"
#include "zora.fh"
      logical Knucl
      integer nbf,npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision qxyz(3,npts)
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision scalx1(nbf,nbf),
     &                 scaly1(nbf,nbf),
     &                 scalz1(nbf,nbf)
      double precision scalxx2(nbf,nbf),
     &                 scalxy2(nbf,nbf),
     &                 scalxz2(nbf,nbf),
     &                 scalyx2(nbf,nbf),
     &                 scalyy2(nbf,nbf),
     &                 scalyz2(nbf,nbf),
     &                 scalzx2(nbf,nbf),
     &                 scalzy2(nbf,nbf),
     &                 scalzz2(nbf,nbf)
      double precision scalx3(nbf,nbf),
     &                 scaly3(nbf,nbf),
     &                 scalz3(nbf,nbf)
      double precision scalxx4(nbf,nbf),
     &                 scalxy4(nbf,nbf),
     &                 scalxz4(nbf,nbf),
     &                 scalyx4(nbf,nbf),
     &                 scalyy4(nbf,nbf),
     &                 scalyz4(nbf,nbf),
     &                 scalzx4(nbf,nbf),
     &                 scalzy4(nbf,nbf),
     &                 scalzz4(nbf,nbf)
      integer i,j,k,n
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot,Kzora
      double precision fac1_arr(npts),fac2_arr(3,npts)
      double precision fac3_arr(3,npts)
      double precision ac_epr1(3),ac_epr2(3,3),
     &                 ac_epr3(3),ac_epr4(3,3)
      double precision chi_cntr(3,nbf)
      logical use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
      integer mbf,ibf(*)
      integer i0,j0
c
      external get_ints_zora_epr_12,get_ints_zora_epr_34
      clight_au2 = clight_au*clight_au
c     == preliminaries ==
      do ipt = 1,npts
        do i=1,ipol
         amat_coul(ipt,i) = 0.d0            
        end do        
        amat_nucl(ipt) = 0.d0
        closegridpts(ipt) = 0
      end do
c
c     == calculate the total hartree potential on the grid ==
      call gridHartreePotential(use_modelpotential,
     &    ao_bas_han, geom, natoms, ipol, g_dens, npts, qxyz, qwght,
     &    closegridpts, gexpo, gcoef, amat_coul)
c
c     == calculate the total nuclear potential on the grid ==
      call gridNuclearPotentialPoint(geom,natoms,npts,qxyz,qwght,
     &                          closegridpts,amat_nucl)
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c     === define fac_arr
      if (.not.(Knucl)) then
       do k = 1,npts
c      == assemble hartree and nuclear contributions ==
         amat_tot = amat_nucl(k) + amat_coul(k,1)
c        Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) ! for g-factors
         Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2)-1.0d0 ! for g-shifts
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (do_NonRel) then             ! remove it after TEST
           Kzora=1.0d0                   ! remove it after TEST
         endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        fac1_arr(k)=Kzora*qwght(k)
        do n=1,3
         fac2_arr(n,k)=fac1_arr(k)*qxyz(n,k)
         fac3_arr(n,k)=qwght(k)*qxyz(n,k)
        enddo ! end-loop-n
       enddo ! end-loop-k
      else if (Knucl) then
       do k = 1,npts
c      == assemble hartree and nuclear contributions ==
          amat_tot = amat_nucl(k)
c         Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) ! for g-factors
          Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2)-1.0d0 ! for g-shifts
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (do_NonRel) then             ! remove it after TEST
          Kzora=1.0d0                   ! remove it after TEST
         endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        fac1_arr(k)=Kzora*qwght(k)
        do n=1,3
         fac2_arr(n,k)=fac1_arr(k)*qxyz(n,k)
         fac3_arr(n,k)=qwght(k)*qxyz(n,k)
        enddo ! end-loop-n
       enddo ! end-loop-k
      else
       write(*,*) 'Error in calc_zora_EPR: Knucl not T of F ...'
       stop
      endif
c     == assemble zora correction ==
c ----- main diagonal --- START
c      do i = 1, nbf
      do i0 = 1, mbf
         i=ibf(i0)
          j=i
          call get_ints_zora_epr_12(nbf,npts,
     &                              delchi_ao,i,j,i0,i0,      
     &                              fac1_arr,
     &                              fac2_arr,  
     &                              ac_epr1, ! out 
     &                              ac_epr2) ! out
          scalx1(i,j)  = scalx1(i,j)  + ac_epr1(1)
          scaly1(i,j)  = scaly1(i,j)  + ac_epr1(2)
          scalz1(i,j)  = scalz1(i,j)  + ac_epr1(3)
          scalxx2(i,j) = scalxx2(i,j) + ac_epr2(1,1)
          scalxy2(i,j) = scalxy2(i,j) + ac_epr2(1,2)
          scalxz2(i,j) = scalxz2(i,j) + ac_epr2(1,3)
          scalyx2(i,j) = scalyx2(i,j) + ac_epr2(2,1)
          scalyy2(i,j) = scalyy2(i,j) + ac_epr2(2,2)
          scalyz2(i,j) = scalyz2(i,j) + ac_epr2(2,3)
          scalzx2(i,j) = scalzx2(i,j) + ac_epr2(3,1)
          scalzy2(i,j) = scalzy2(i,j) + ac_epr2(3,2)
          scalzz2(i,j) = scalzz2(i,j) + ac_epr2(3,3)
      enddo ! end-loop-i
c ----- main diagonal --- END
c ----- off diagonal --- START
      do i0 = 1, mbf
         i=ibf(i0)
        do j0 = i0+1, mbf
           j=ibf(j0)
          call get_ints_zora_epr_12(nbf,npts,
     &                              delchi_ao,i,j,i0,j0,      
     &                              fac1_arr,
     &                              fac2_arr,  
     &                              ac_epr1, ! out 
     &                              ac_epr2) ! out
          scalx1(i,j)  = scalx1(i,j)  + 2.0d0*ac_epr1(1)
          scaly1(i,j)  = scaly1(i,j)  + 2.0d0*ac_epr1(2)
          scalz1(i,j)  = scalz1(i,j)  + 2.0d0*ac_epr1(3)
          scalxx2(i,j) = scalxx2(i,j) + 2.0d0*ac_epr2(1,1)
          scalxy2(i,j) = scalxy2(i,j) + 2.0d0*ac_epr2(1,2)
          scalxz2(i,j) = scalxz2(i,j) + 2.0d0*ac_epr2(1,3)
          scalyx2(i,j) = scalyx2(i,j) + 2.0d0*ac_epr2(2,1)
          scalyy2(i,j) = scalyy2(i,j) + 2.0d0*ac_epr2(2,2)
          scalyz2(i,j) = scalyz2(i,j) + 2.0d0*ac_epr2(2,3)
          scalzx2(i,j) = scalzx2(i,j) + 2.0d0*ac_epr2(3,1)
          scalzy2(i,j) = scalzy2(i,j) + 2.0d0*ac_epr2(3,2)
          scalzz2(i,j) = scalzz2(i,j) + 2.0d0*ac_epr2(3,3)
        enddo ! end-loop-j
      enddo ! end-loop-i
c ----- off diagonal --- END
c ---- full matrix calc -------- START
      do i0 = 1, mbf
         i=ibf(i0)
        do j0 = 1, mbf
           j=ibf(j0)
          call get_ints_zora_epr_34(nbf,npts,chi_ao,delchi_ao,i,j,i0,j0,      
     &                              fac1_arr,
     &                              fac2_arr,  
     &                              fac3_arr,  
     &                              chi_cntr, ! in     
     &                              ac_epr3,  ! out 
     &                              ac_epr4)  ! out
          scalx3(i,j)  = scalx3(i,j)  + ac_epr3(1)
          scaly3(i,j)  = scaly3(i,j)  + ac_epr3(2)
          scalz3(i,j)  = scalz3(i,j)  + ac_epr3(3)
          scalxx4(i,j) = scalxx4(i,j) + ac_epr4(1,1)
          scalxy4(i,j) = scalxy4(i,j) + ac_epr4(1,2)
          scalxz4(i,j) = scalxz4(i,j) + ac_epr4(1,3)
          scalyx4(i,j) = scalyx4(i,j) + ac_epr4(2,1)
          scalyy4(i,j) = scalyy4(i,j) + ac_epr4(2,2)
          scalyz4(i,j) = scalyz4(i,j) + ac_epr4(2,3)
          scalzx4(i,j) = scalzx4(i,j) + ac_epr4(3,1)
          scalzy4(i,j) = scalzy4(i,j) + ac_epr4(3,2)
          scalzz4(i,j) = scalzz4(i,j) + ac_epr4(3,3)
        enddo ! end-loop-j
      enddo ! end-loop-i
c ---- full matrix calc -------- END
      return
      end
      subroutine get_ints_zora_epr(nbf,       ! in: nr. basis functions
     &                             npts,      ! in: grid points 
     &                             chi_ao,    ! in:           basis functions
     &                             delchi_ao, ! in: deriv. of basis functions
     &                             i,j,       ! in: (i,j) indices for delchi_ao
     &                                i0,j0,     ! restricted indices after grid_nbfm
     &                             fac1_arr,  ! in
     &                             fac2_arr,  ! in
     &                             ac_epr1,   ! out 
     &                             ac_epr2,   ! out 
     &                             ac_epr3,   ! out 
     &                             ac_epr4)   ! out
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k,m,n,i0,j0
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts),fac2_arr(3,npts)
      double precision ac_epr1(3),ac_epr2(3,3),
     &                 ac_epr3(3),ac_epr4(3,3) 
      double precision prod(3),prod1(3)
      do n=1,3 ! reset
       ac_epr1(n) = 0.0d0
       ac_epr3(n) = 0.0d0
       do m=1,3
        ac_epr2(m,n) = 0.0d0
        ac_epr4(m,n) = 0.0d0
       enddo
      enddo    
      do k = 1, npts
       prod(1)= delchi_ao(k,2,i0)*delchi_ao(k,3,j0)
     &         -delchi_ao(k,3,i0)*delchi_ao(k,2,j0) 
       prod(2)= delchi_ao(k,3,i0)*delchi_ao(k,1,j0)
     &         -delchi_ao(k,1,i0)*delchi_ao(k,3,j0) 
       prod(3)= delchi_ao(k,1,i0)*delchi_ao(k,2,j0)
     &         -delchi_ao(k,2,i0)*delchi_ao(k,1,j0) 
       do n=1,3
        prod1(n) = chi_ao(k,i0)*delchi_ao(k,n,j0)
       enddo
       do n=1,3
        ac_epr1(n) = ac_epr1(n) + fac1_arr(k)*prod(n)
        ac_epr3(n) = ac_epr3(n) + fac1_arr(k)*prod1(n)
        do m=1,3
         ac_epr2(m,n) = ac_epr2(m,n) + fac2_arr(m,k)*prod(n)
         ac_epr4(m,n) = ac_epr4(m,n) + fac2_arr(m,k)*prod1(n)
        enddo ! end-loop-m
       enddo ! end-loop-n
      enddo ! end-loo-k
      return
      end
      subroutine get_ints_zora_epr_12(nbf,       ! in: nr. basis functions
     &                                npts,      ! in: grid points 
     &                                delchi_ao, ! in: deriv. of basis functions
     &                                i,j,       ! in: (i,j) indices for delchi_ao
     &                                i0,j0,     ! restricted indices after grid_nbfm
     &                                fac1_arr,  ! in
     &                                fac2_arr,  ! in
     &                                ac_epr1,   ! out 
     &                                ac_epr2)   ! out
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k,m,n,i0,j0
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts),fac2_arr(3,npts)
      double precision fac3_arr(3,npts)
      double precision ac_epr1(3),ac_epr2(3,3)
      double precision prod(3)

      do n=1,3 ! reset
      ac_epr1(n) = 0.0d0
       do m=1,3
        ac_epr2(m,n) = 0.0d0
       enddo
      enddo    
      do k = 1, npts
       prod(1)= delchi_ao(k,2,i0)*delchi_ao(k,3,j0)
     &         -delchi_ao(k,3,i0)*delchi_ao(k,2,j0) 
       prod(2)= delchi_ao(k,3,i0)*delchi_ao(k,1,j0)
     &         -delchi_ao(k,1,i0)*delchi_ao(k,3,j0) 
       prod(3)= delchi_ao(k,1,i0)*delchi_ao(k,2,j0)
     &         -delchi_ao(k,2,i0)*delchi_ao(k,1,j0) 
       do n=1,3
        ac_epr1(n) = ac_epr1(n) + fac1_arr(k)*prod(n)
        do m=1,3
          ac_epr2(m,n) = ac_epr2(m,n) + fac2_arr(m,k)*prod(n)
        enddo ! end-loop-m
       enddo ! end-loop-n
      enddo ! end-loo-k
      return
      end

      subroutine get_ints_zora_epr_34(nbf,       ! in: nr. basis functions
     &                                npts,      ! in: grid points 
     &                                chi_ao,    ! in:           basis functions
     &                                delchi_ao, ! in: deriv. of basis functions
     &                                i,j,       ! in: (i,j) indices for delchi_ao
     &                                i0,j0,     ! restricted indices after grid_nbfm
     &                                fac1_arr,  ! in
     &                                fac2_arr,  ! in
     &                                fac3_arr,  ! in 
     &                                chi_cntr,  ! in : basis funct. centers
     &                                ac_epr3,   ! out 
     &                                ac_epr4)   ! out
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k,m,n,i0,j0
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts),fac2_arr(3,npts)
      double precision fac3_arr(3,npts)
      double precision ac_epr3(3),ac_epr4(3,3) 
      double precision prod(3),val
      double precision prod1
      double precision chi_cntr(3,nbf)

      do n=1,3 ! reset
       ac_epr3(n) = 0.0d0
       do m=1,3
       ac_epr4(m,n) = 0.0d0
       enddo
      enddo    
      do k = 1, npts
       prod1 = chi_ao(k,i0)*chi_ao(k,j0)
       do n=1,3
        prod(n) = chi_ao(k,i0)*delchi_ao(k,n,j0)
        ac_epr3(n) = ac_epr3(n) + fac1_arr(k)*prod(n)
        prod(n)=prod(n)+chi_ao(k,j0)*delchi_ao(k,n,i0)
        do m=1,3
         ac_epr4(m,n)=ac_epr4(m,n)+(fac2_arr(m,k)-
     &                              fac1_arr(k)*chi_cntr(m,j))*
     &                              prod(n)
        enddo ! end-loop-m
       enddo ! end-loop-n
      enddo ! end-loo-k
      return
      end

c $Id$
