c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     + calculate all NMR Chemical Shieldings matricial AOs +
c     + for scalar relativistic ZORA                        +
c     + FA-11-09-10                                         +
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine calc_zora_NMRCS_SR(ao_bas_han,    ! in: AO basis handle 
     &                              geom,          ! in: geometry handle 
     &                              ipol,          ! in: nr. of polarizations 
     &                              g_dens,        ! in: superposit. atomic densities 
     &                              chi_ao,        ! in:           basis functions
     &                              delchi_ao,     ! in: deriv. of basis functions
     &                              qxyz,          ! in: grid points 
     &                              qwght,         ! in: weighting coeffs.
     &                              nbf,           ! in: nr. basis functions
     &                              npts,          ! in: nr. grid points
     &                              natoms,        ! in: nr. atoms 
     &                              xyz_NMRcoords, ! in : nuclear coordinates
     &                              use_modelpotential,
     &                              gexpo,
     &                              gcoef,
     &                              dia1,          ! out
     &                              dia2x,         ! out
     &                              dia2y,         ! out
     &                              dia2z,         ! out
     &                              dia3xx,        ! out
     &                              dia3xy,        ! out
     &                              dia3xz,        ! out
     &                              dia3yx,        ! out
     &                              dia3yy,        ! out
     &                              dia3yz,        ! out
     &                              dia3zx,        ! out
     &                              dia3zy,        ! out
     &                              dia3zz,        ! out
     &                              scalxy1,       ! out 
     &                              scalxz1,       ! out 
     &                              scalyx1,       ! out 
     &                              scalyz1,       ! out 
     &                              scalzx1,       ! out 
     &                              scalzy1,       ! out 
     &                              scalxy21,      ! out 
     &                              scalxz21,      ! out 
     &                              scalyx21,      ! out 
     &                              scalyz21,      ! out 
     &                              scalzx21,      ! out 
     &                              scalzy21,      ! out 
     &                              scalxy22,      ! out 
     &                              scalxz22,      ! out 
     &                              scalyx22,      ! out 
     &                              scalyz22,      ! out 
     &                              scalzx22,      ! out 
     &                              scalzy22,      ! out 
     &                              scalxy23,      ! out 
     &                              scalxz23,      ! out 
     &                              scalyx23,      ! out 
     &                              scalyz23,      ! out 
     &                              scalzx23,      ! out 
     &                              scalzy23)      ! out 
c     Calculates:
c     scalyz1,scalzx1,scalxy1,    --> B1_{munu,mn}
c     scalzy1,scalxz1,scalyx1
c     scalyz21,scalzx21,scalxy21, --> B2_{munu,mn,a=1}
c     scalzy21,scalxz21,scalyx21
c     scalyz22,scalzx22,scalxy22, --> B2_{munu,mn,a=2}
c     scalzy22,scalxz22,scalyx22
c     scalyz23,scalzx23,scalxy23, --> B2_{munu,mn,a=3}
c     scalzy23,scalxz23,scalyx23
c
c     B1_{munu,mn}  =<chi_{mu}|    K r_{Q,m}/r_Q^3| dchi_{nu}/dr_n> -> ac_nmr1(m,n)
c     B2_{munu,mn,a}=<chi_{mu}|r_a K r_{Q,m}/r_Q^3| dchi_{nu}/dr_n> -> ac_nmr2(a,m,n)
c       a=1,2,3=x,y,z
c     tmn=123,231,312,
c         132,213,321
c     FA-11-07-10
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "zora.fh"
#include "global.fh"
#include "bas.fh"
      integer nbf,npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision qxyz(3,npts)
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision dia1(nbf,nbf),
     &                 dia2x(nbf,nbf),
     &                 dia2y(nbf,nbf),
     &                 dia2z(nbf,nbf),
     &                 dia3xx(nbf,nbf),
     &                 dia3xy(nbf,nbf),
     &                 dia3xz(nbf,nbf),
     &                 dia3yx(nbf,nbf),
     &                 dia3yy(nbf,nbf),
     &                 dia3yz(nbf,nbf),
     &                 dia3zx(nbf,nbf),
     &                 dia3zy(nbf,nbf),
     &                 dia3zz(nbf,nbf)
      double precision scalxy1(nbf,nbf),
     &                 scalxz1(nbf,nbf),
     &                 scalyx1(nbf,nbf),
     &                 scalyz1(nbf,nbf),
     &                 scalzx1(nbf,nbf),
     &                 scalzy1(nbf,nbf),
     &                 scalxy21(nbf,nbf),
     &                 scalxz21(nbf,nbf),
     &                 scalyx21(nbf,nbf),
     &                 scalyz21(nbf,nbf),
     &                 scalzx21(nbf,nbf),
     &                 scalzy21(nbf,nbf),
     &                 scalxy22(nbf,nbf),
     &                 scalxz22(nbf,nbf),
     &                 scalyx22(nbf,nbf),
     &                 scalyz22(nbf,nbf),
     &                 scalzx22(nbf,nbf),
     &                 scalzy22(nbf,nbf),
     &                 scalxy23(nbf,nbf),
     &                 scalxz23(nbf,nbf),
     &                 scalyx23(nbf,nbf),
     &                 scalyz23(nbf,nbf),
     &                 scalzx23(nbf,nbf),
     &                 scalzy23(nbf,nbf)
      integer i,j,k,m,n
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts),amat_NMRnucl(3,npts),
     &                 amat_NMRnucl4dia(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot,Kzora
      double precision fac1,fac2,
     &                 fac3(3,npts),fac4(3,3,npts),
     &                 fac5(npts),fac6(3,npts),fac7(3,3,npts)
      double precision ac_nmr1(3,3),ac_nmr2(3,3,3)
      double precision ac_dia1,ac_dia2(3),ac_dia3(3,3) 
      double precision xyz_NMRcoords(3)   
      external get_ints_zorasr_nmr12,
     &         get_ints_zora_nmr_dia,
     &         gridNMRPotential,gridNMRPotential_dia
c
      logical use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
c
      clight_au2 = clight_au*clight_au
c     == preliminaries ==
      do ipt = 1,npts
        do i=1,ipol
         amat_coul(ipt,i) = 0.d0            
        end do        
        amat_nucl(ipt)        = 0.d0
        amat_NMRnucl4dia(ipt) = 0.d0
        closegridpts(ipt) = 0
        do i=1,3
         amat_NMRnucl(i,ipt) = 0.d0
        enddo
      enddo
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
      call gridNMRPotential(amat_NMRnucl,  ! out: NMR potential
     &                      xyz_NMRcoords, 
     &                      npts,qxyz,closegridpts)
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
      call gridNMRPotential_dia(amat_NMRnucl4dia,  ! out: NMR potential
     &                          xyz_NMRcoords, 
     &                          npts,qxyz,closegridpts)
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c     === define fac_arr
      do k = 1,npts
c      == assemble hartree and nuclear contributions ==
        amat_tot = amat_coul(k,1) + amat_nucl(k)
        Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) 
c       Kzora=1.0d0+0.5d0*amat_tot/clight_au2 ! Using Pauli approx
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c WARNING: When asking for NON REL will use K ne 1
       if (do_NonRel) then             ! remove it after TEST
         Kzora=1.0d0                   ! remove it after TEST
       endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       fac1=Kzora*qwght(k)
       fac5(k)=fac1*amat_NMRnucl4dia(k)
       do n=1,3
        fac2=fac1*qxyz(n,k)
        fac3(n,k)=fac1*amat_NMRnucl(n,k)
        fac6(n,k)=fac2*amat_NMRnucl4dia(k)
        do m=1,3
         fac4(m,n,k)=fac2*amat_NMRnucl(m,k)
         fac7(m,n,k)=fac6(n,k)*qxyz(m,k)
        enddo ! end-loop-m
       enddo ! end-loop-n
      enddo ! end-loop-k 
c ----- main diagonal --- START
      do i = 1, nbf
          j=i
          call get_ints_zora_nmr_dia(nbf,npts,chi_ao,i,j,      
     &                              fac5,fac6,fac7,       
     &                              ac_dia1,ac_dia2,ac_dia3) ! out 
         dia1(i,j)   = dia1(i,j)   + ac_dia1
         dia2x(i,j)  = dia2x(i,j)  + ac_dia2(1)
         dia2y(i,j)  = dia2y(i,j)  + ac_dia2(2)
         dia2z(i,j)  = dia2z(i,j)  + ac_dia2(3)
         dia3xx(i,j) = dia3xx(i,j) + ac_dia3(1,1)
         dia3xy(i,j) = dia3xy(i,j) + ac_dia3(1,2)
         dia3xz(i,j) = dia3xz(i,j) + ac_dia3(1,3)
         dia3yx(i,j) = dia3yx(i,j) + ac_dia3(2,1)
         dia3yy(i,j) = dia3yy(i,j) + ac_dia3(2,2)
         dia3yz(i,j) = dia3yz(i,j) + ac_dia3(2,3)
         dia3zx(i,j) = dia3zx(i,j) + ac_dia3(3,1)
         dia3zy(i,j) = dia3zy(i,j) + ac_dia3(3,2)
         dia3zz(i,j) = dia3zz(i,j) + ac_dia3(3,3)
      enddo ! end-loop-i
c ----- main diagonal --- END
c ----- off diagonal --- START
      do i = 1, nbf
        do j = i+1, nbf
          call get_ints_zora_nmr_dia(nbf,npts,chi_ao,i,j,      
     &                              fac5,fac6,fac7,       
     &                              ac_dia1,ac_dia2,ac_dia3) ! out 
         dia1(i,j)   = dia1(i,j)   + 2.0d0*ac_dia1
         dia2x(i,j)  = dia2x(i,j)  + 2.0d0*ac_dia2(1)
         dia2y(i,j)  = dia2y(i,j)  + 2.0d0*ac_dia2(2)
         dia2z(i,j)  = dia2z(i,j)  + 2.0d0*ac_dia2(3)
         dia3xx(i,j) = dia3xx(i,j) + 2.0d0*ac_dia3(1,1)
         dia3xy(i,j) = dia3xy(i,j) + 2.0d0*ac_dia3(1,2)
         dia3xz(i,j) = dia3xz(i,j) + 2.0d0*ac_dia3(1,3)
         dia3yx(i,j) = dia3yx(i,j) + 2.0d0*ac_dia3(2,1)
         dia3yy(i,j) = dia3yy(i,j) + 2.0d0*ac_dia3(2,2)
         dia3yz(i,j) = dia3yz(i,j) + 2.0d0*ac_dia3(2,3)
         dia3zx(i,j) = dia3zx(i,j) + 2.0d0*ac_dia3(3,1)
         dia3zy(i,j) = dia3zy(i,j) + 2.0d0*ac_dia3(3,2)
         dia3zz(i,j) = dia3zz(i,j) + 2.0d0*ac_dia3(3,3)
        enddo ! end-loop-j
      enddo ! end-loop-i
c ----- off diagonal --- END
c ---- full matrix calc -------- START
      do i = 1, nbf
        do j = 1, nbf
          call get_ints_zorasr_nmr12(nbf,npts,chi_ao,delchi_ao,i,j,      
     &                               fac3,fac4,
     &                               ac_nmr1,ac_nmr2) ! out  
          scalxy1(i,j)  = scalxy1(i,j)  + ac_nmr1(1,2)
          scalxz1(i,j)  = scalxz1(i,j)  + ac_nmr1(1,3)
          scalyx1(i,j)  = scalyx1(i,j)  + ac_nmr1(2,1)
          scalyz1(i,j)  = scalyz1(i,j)  + ac_nmr1(2,3)
          scalzx1(i,j)  = scalzx1(i,j)  + ac_nmr1(3,1)
          scalzy1(i,j)  = scalzy1(i,j)  + ac_nmr1(3,2)   
          scalxy21(i,j) = scalxy21(i,j) + ac_nmr2(1,1,2)
          scalxz21(i,j) = scalxz21(i,j) + ac_nmr2(1,1,3)
          scalyx21(i,j) = scalyx21(i,j) + ac_nmr2(1,2,1)
          scalyz21(i,j) = scalyz21(i,j) + ac_nmr2(1,2,3)
          scalzx21(i,j) = scalzx21(i,j) + ac_nmr2(1,3,1)
          scalzy21(i,j) = scalzy21(i,j) + ac_nmr2(1,3,2)  
          scalxy22(i,j) = scalxy22(i,j) + ac_nmr2(2,1,2)
          scalxz22(i,j) = scalxz22(i,j) + ac_nmr2(2,1,3)
          scalyx22(i,j) = scalyx22(i,j) + ac_nmr2(2,2,1)
          scalyz22(i,j) = scalyz22(i,j) + ac_nmr2(2,2,3)
          scalzx22(i,j) = scalzx22(i,j) + ac_nmr2(2,3,1)
          scalzy22(i,j) = scalzy22(i,j) + ac_nmr2(2,3,2)
          scalxy23(i,j) = scalxy23(i,j) + ac_nmr2(3,1,2)
          scalxz23(i,j) = scalxz23(i,j) + ac_nmr2(3,1,3)
          scalyx23(i,j) = scalyx23(i,j) + ac_nmr2(3,2,1)
          scalyz23(i,j) = scalyz23(i,j) + ac_nmr2(3,2,3)
          scalzx23(i,j) = scalzx23(i,j) + ac_nmr2(3,3,1)
          scalzy23(i,j) = scalzy23(i,j) + ac_nmr2(3,3,2) 
        enddo ! end-loop-j
      enddo ! end-loop-i
c ---- full matrix calc -------- END
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine calc_NMRCS_SR_F1ij(ao_bas_han,    ! in: AO basis handle 
     &                              geom,          ! in: geometry handle 
     &                              ipol,          ! in: nr. of polarizations 
     &                              g_dens,        ! in: superposit. atomic densities 
     &                              chi_ao,        ! in:           basis functions
     &                              delchi_ao,     ! in: deriv. of basis functions
     &                              qxyz,          ! in: grid points 
     &                              qwght,         ! in: weighting coeffs.
     &                              nbf,           ! in: nr. basis functions
     &                              npts,          ! in: nr. grid points
     &                              natoms,        ! in: nr. atoms 
     &                              use_modelpotential,
     &                              gexpo,
     &                              gcoef,
     &                              scalx,         ! out 
     &                              scaly,         ! out 
     &                              scalz,         ! out
     &                              scalx3,        ! out 
     &                              scaly3,        ! out 
     &                              scalz3,        ! out
     &                              scalxy4,       ! out 
     &                              scalxz4,       ! out 
     &                              scalyx4,       ! out 
     &                              scalyz4,       ! out 
     &                              scalzx4,       ! out 
     &                              scalzy4)       ! out 
c     Calculates:
c     scalx,scaly,scalz           --> A_{munu,a}
c     scalx3,scaly3,scalz3        --> B3_{munu,a} 
c     scalyz4,scalzx4,scalxy4,    --> B4_{munu,mn}
c     scalzy4,scalxz4,scalyx4
c     A_{munu,a}    = \sum_{s=1,2,3} <dchi_{mu}/dr_s|r_a K | dchi_{nu}/dr_s>
c     B3_{munu,a}   =<chi_{mu}|    K              | dchi_{nu}/dr_a> -> ac_nmr3(a)
c     B4_{munu,mn}  =<chi_{mu}|r_m K              | dchi_{nu}/dr_n> -> ac_nmr4(m,n)
c       a=1,2,3=x,y,z
c     tmn=123,231,312,
c         132,213,321
c     --> Those matricial AOs are used in evaluation of 1st order derivative
c         of Fock operator (F_{ij}^{1k})
c     FA-11-07-10
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "zora.fh"
#include "global.fh"
#include "bas.fh"
      integer nbf,npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision qxyz(3,npts)
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision scalx(nbf,nbf),
     &                 scaly(nbf,nbf),
     &                 scalz(nbf,nbf),
     &                 scalx3(nbf,nbf),
     &                 scaly3(nbf,nbf),
     &                 scalz3(nbf,nbf),
     &                 scalxy4(nbf,nbf),
     &                 scalxz4(nbf,nbf),
     &                 scalyx4(nbf,nbf),
     &                 scalyz4(nbf,nbf),
     &                 scalzx4(nbf,nbf),
     &                 scalzy4(nbf,nbf)
      integer i,j,k,m,n
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot,Kzora
      double precision fac1(npts),fac2(3,npts)
      double precision ac_nmr3(3),ac_nmr4(3,3)
      double precision ac_nmrsf(3) 
      external get_ints_zora_nmr_sf,get_ints_zorasr_nmr34
c
      logical use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
c
      clight_au2 = clight_au*clight_au
c     == preliminaries ==
      do ipt = 1,npts
        do i=1,ipol
         amat_coul(ipt,i) = 0.d0            
        end do        
        amat_nucl(ipt)        = 0.d0
        closegridpts(ipt) = 0
      enddo
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
      do k = 1,npts
c      == assemble hartree and nuclear contributions ==
        amat_tot = amat_coul(k,1) + amat_nucl(k)
        Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) 
c        Kzora=1.0d0+0.5d0*amat_tot/clight_au2 ! Using Pauli
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c WARNING: When asking for NON REL will use K ne 1
       if (do_NonRel) then ! remove it after TEST
         Kzora=1.0d0                   ! remove it after TEST
       endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       fac1(k)=Kzora*qwght(k)
       do n=1,3
        fac2(n,k)=fac1(k)*qxyz(n,k)
       enddo ! end-loop-n
      enddo ! end-loop-k
  
c ----- main diagonal --- START
      do i = 1, nbf
          j=i
          call get_ints_zora_nmr_sf(nbf,npts,delchi_ao,i,j,      
     &                              fac2,       
     &                              ac_nmrsf) ! out 
          scalx(i,j) = scalx(i,j) + ac_nmrsf(1)
          scaly(i,j) = scaly(i,j) + ac_nmrsf(2)
          scalz(i,j) = scalz(i,j) + ac_nmrsf(3)
      enddo ! end-loop-i
c ----- main diagonal --- END
c ----- off diagonal --- START
      do i = 1, nbf
        do j = i+1, nbf
          call get_ints_zora_nmr_sf(nbf,npts,delchi_ao,i,j,      
     &                              fac2,       
     &                              ac_nmrsf) ! out 
          scalx(i,j) = scalx(i,j) + 2.0d0*ac_nmrsf(1)
          scaly(i,j) = scaly(i,j) + 2.0d0*ac_nmrsf(2)
          scalz(i,j) = scalz(i,j) + 2.0d0*ac_nmrsf(3)
        enddo ! end-loop-j
      enddo ! end-loop-i
c ----- off diagonal --- END
c ---- full matrix calc -------- START
      do i = 1, nbf
        do j = 1, nbf
          call get_ints_zorasr_nmr34(nbf,npts,chi_ao,delchi_ao,i,j,      
     &                               fac1,fac2,
     &                               ac_nmr3,ac_nmr4) ! out   
          scalx3(i,j)   = scalx3(i,j)   + ac_nmr3(1)
          scaly3(i,j)   = scaly3(i,j)   + ac_nmr3(2)
          scalz3(i,j)   = scalz3(i,j)   + ac_nmr3(3)
          scalxy4(i,j)  = scalxy4(i,j)  + ac_nmr4(1,2)
          scalxz4(i,j)  = scalxz4(i,j)  + ac_nmr4(1,3)
          scalyx4(i,j)  = scalyx4(i,j)  + ac_nmr4(2,1)
          scalyz4(i,j)  = scalyz4(i,j)  + ac_nmr4(2,3)
          scalzx4(i,j)  = scalzx4(i,j)  + ac_nmr4(3,1)
          scalzy4(i,j)  = scalzy4(i,j)  + ac_nmr4(3,2)
        enddo ! end-loop-j
      enddo ! end-loop-i
c ---- full matrix calc -------- END
      return
      end

       subroutine get_ints_zorasr_nmr12(nbf,       ! in: nr. basis functions
     &                                 npts,      ! in: grid points 
     &                                 chi_ao,    ! in:           basis functions
     &                                 delchi_ao, ! in: deriv. of basis functions
     &                                 i,j,       ! in: (i,j) indices for delchi_ao
     &                                 fac3,      ! in
     &                                 fac4,      ! in
     &                                 ac_nmr1,   ! out 
     &                                 ac_nmr2)   ! out
c     Calculates:
c     B1_{munu,mn}  =<chi_{mu}|    K r_{Q,m}/r_Q^3| dchi_{nu}/dr_n> -> ac_nmr1(m,n)
c     B2_{munu,mn,a}=<chi_{mu}|r_a K r_{Q,m}/r_Q^3| dchi_{nu}/dr_n> -> ac_nmr2(a,m,n)
c       a=1,2,3
c     tmn=123,231,312,
c         132,213,321
c     FA-11-04-10
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k,m,n,a,t
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision ac_nmr1(3,3),ac_nmr2(3,3,3)
      double precision fac3(3,npts),fac4(3,3,npts)
      double precision prod1(3)
      integer ind_tmn(2,3)  
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312 
       do n=1,3 ! reset
        do m=1,3
         ac_nmr1(m,n) = 0.0d0
         do a=1,3
          ac_nmr2(a,m,n) = 0.0d0
         enddo ! end-loop-a
        enddo ! end-loop-m
       enddo ! end-loop-n
       do k = 1, npts
        do n=1,3
         prod1(n)   = chi_ao(k,i)*delchi_ao(k,n,j)
        enddo    
        do t=1,3
         m=ind_tmn(1,t)
         n=ind_tmn(2,t)
         ac_nmr1(m,n) = ac_nmr1(m,n) + fac3(m,k)*prod1(n)
         ac_nmr1(n,m) = ac_nmr1(n,m) + fac3(n,k)*prod1(m)
         do a=1,3 ! = r_a = x,y,z
          ac_nmr2(a,m,n) = ac_nmr2(a,m,n) + fac4(m,a,k)*prod1(n)
          ac_nmr2(a,n,m) = ac_nmr2(a,n,m) + fac4(n,a,k)*prod1(m)
         enddo ! end-loop-a
        enddo ! end-loop-t
       enddo ! end-loo-k
      return
      end

      subroutine get_ints_zorasr_nmr34(nbf,       ! in: nr. basis functions
     &                                 npts,      ! in: grid points 
     &                                 chi_ao,    ! in:           basis functions
     &                                 delchi_ao, ! in: deriv. of basis functions
     &                                 i,j,       ! in: (i,j) indices for delchi_ao
     &                                 fac1,      ! in
     &                                 fac2,      ! in
     &                                 ac_nmr3,   ! out 
     &                                 ac_nmr4)   ! out
c     Calculates:
c     B3_{munu,a}   =<chi_{mu}|    K              | dchi_{nu}/dr_a> -> ac_nmr3(a)
c     B4_{munu,m,n} =<chi_{mu}|r_m K              | dchi_{nu}/dr_n> -> ac_nmr4(m,n)
c       a=1,2,3
c     tmn=123,231,312,
c         132,213,321
c     FA-11-07-10
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,k,t,m,n,i,j
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision ac_nmr3(3),ac_nmr4(3,3)
      double precision fac1(npts),fac2(3,npts)
      double precision prod1(3)
      integer ind_tmn(2,3)  
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312 
          do n=1,3 ! reset
           ac_nmr3(n) = 0.0d0
           do m=1,3
            ac_nmr4(m,n) = 0.0d0
           enddo ! end-loop-m
          enddo ! end-loop-n
          do k = 1, npts
            do n=1,3
             prod1(n)   = chi_ao(k,i)*delchi_ao(k,n,j)
             ac_nmr3(n) = ac_nmr3(n) + fac1(k)*prod1(n)
            enddo    
            do t=1,3
             m=ind_tmn(1,t)
             n=ind_tmn(2,t)
             ac_nmr4(m,n) = ac_nmr4(m,n) + fac2(m,k)*prod1(n)
             ac_nmr4(n,m) = ac_nmr4(n,m) + fac2(n,k)*prod1(m)
            enddo ! end-loop-t
          enddo ! end-loo-k
      return
      end

      subroutine get_ints_zora_nmr_sf(nbf,      ! in: nr. basis functions
     &                                npts,     ! in: grid points 
     &                                delchi_ao,! in: deriv. of basis functions
     &                                i,j,      ! in: (i,j) indices for delchi_ao
     &                                fac,      ! in
     &                                ac_nmrsf) ! out
c     Calculates:
c     A_{munu,a} = \sum_{s=1,2,3} Y_{munu,ss,a}
c     Y_{munu,ss,a}=<dchi_{mu}/dr_s|r_a K | dchi_{nu}/dr_s>
c     a=1,2,3=x,y,z
c     FA-11-04-10
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k,n
      double precision delchi_ao(npts,3,nbf)
      double precision fac(3,npts)
      double precision ac_nmrsf(3) 
      double precision prod
      do n=1,3
        ac_nmrsf(n) = 0.0d0
      enddo
      do k = 1, npts
       prod = delchi_ao(k,1,i)*delchi_ao(k,1,j)
     &       +delchi_ao(k,2,i)*delchi_ao(k,2,j) 
     &       +delchi_ao(k,3,i)*delchi_ao(k,3,j)  
       do n=1,3
        ac_nmrsf(n) = ac_nmrsf(n) + fac(n,k)*prod
       enddo ! end-loo-n
      enddo ! end-loo-k
      return
      end

      subroutine get_ints_zora_nmr_dia(nbf,      ! in: nr. basis functions
     &                                 npts,     ! in: grid points 
     &                                 chi_ao,   ! in: deriv. of basis functions
     &                                 i,j,      ! in: (i,j) indices for delchi_ao
     &                                 fac5,     ! in
     &                                 fac6,     ! in
     &                                 fac7,     ! in
     &                                 ac_dia1,  ! out
     &                                 ac_dia2,  ! out
     &                                 ac_dia3)  ! out
c     Calculates:
c   1. < chi_mu | K        /r_Q^3 | chi_nu>            (1 matricial AOs) -> ac_dia1
c   2. < chi_mu | K r_k    /r_Q^3 | chi_nu>  k=1,2,3   (3 matricial AOs) -> ac_dia2
c   3. < chi_mu | K r_k r_t/r_Q^3 | chi_nu>  k,t=1,2,3 (9 matricial AOs) -> ac_dia3
c     FA-11-05-10
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k,n,m
      double precision chi_ao(npts,nbf)
      double precision ac_dia1,ac_dia2(3),ac_dia3(3,3) 
      double precision prod
      double precision fac5(npts),fac6(3,npts),fac7(3,3,npts)
        ac_dia1 = 0.0d0
      do n=1,3
        ac_dia2(n) = 0.0d0
        do k=1,3
         ac_dia3(n,k) = 0.0d0
        enddo
      enddo
      do k = 1, npts
       prod = chi_ao(k,i)*chi_ao(k,j) 
       ac_dia1 = ac_dia1 + fac5(k)*prod
       do n=1,3
        ac_dia2(n) = ac_dia2(n) + fac6(n,k)*prod
        do m=1,3
         ac_dia3(m,n) = ac_dia3(m,n) + fac7(m,n,k)*prod        
        enddo ! end-loo-m        
       enddo ! end-loo-n
      enddo ! end-loo-k
      return
      end

      subroutine gridNMRPotential(amat_NMRnucl,  ! out: EFG potential
     &                            xyz_NMRcoords, ! in : EFG-nuclear coord.
     &                            nqpts,         ! in : nr. grid points 
     &                            qxyz,          ! in : grid points                                   
     &                            closegridpts)
      implicit none
#include "global.fh" 
#include "stdio.fh"
#include "zora.fh"
      integer i,igrid,nqpts
      integer closegridpts(*)
      double precision xyz_NMRcoords(3)
      double precision qxyz(3,nqpts)
      double precision rxyz(3),dist,dist3,ac_prod
      double precision amat_NMRnucl(3,nqpts)
c     == loop over the grid points ==      
       do igrid = 1,nqpts
        ac_prod=0.0d0
        do i=1,3
           rxyz(i) = qxyz(i,igrid)-xyz_NMRcoords(i) 
           ac_prod=ac_prod+rxyz(i)*rxyz(i)
        enddo
        dist = dsqrt(ac_prod)
        if (dist.gt.zoracutoff_NMR) then ! check-cutoff
          dist3=dist*dist*dist
          do i=1,3
           amat_NMRnucl(i,igrid) = rxyz(i)/dist3
          enddo 
         goto 3
         if (dist.lt.zoracutoff_NMR*100.5d0) then
          if (ga_nodeid().eq.0) then
          write(*,2) igrid,nqpts,
     &     xyz_NMRcoords(1),xyz_NMRcoords(2),xyz_NMRcoords(3),
     &     qxyz(1,igrid),qxyz(2,igrid),qxyz(3,igrid),
     &     dist,zoracutoff_NMR,
     &     amat_NMRnucl(1,igrid),amat_NMRnucl(2,igrid),
     &     amat_NMRnucl(3,igrid)
 2       format('IN:(igrid,nqpts)=(',i8,',',i8,') ',
     &          'xyz_NMRcoords=(',f15.8,',',f15.8,',',f15.8,') ',
     &          'qxyz=(',f15.8,',',f15.8,',',f15.8,') ',
     &          '(dist,zoracutoff_NMR,100.5)=(',f15.8,',',f15.8,') ',
     &          ' NMRpot=(',f25.8,',',f25.8,',',f25.8,')')
          endif 
         endif
 3      continue
        else

         goto 4
         if (ga_nodeid().eq.0) then
          write(*,1) igrid,nqpts,
     &     xyz_NMRcoords(1),xyz_NMRcoords(2),xyz_NMRcoords(3),
     &     qxyz(1,igrid),qxyz(2,igrid),qxyz(3,igrid),
     &     dist,zoracutoff_NMR,
     &     amat_NMRnucl(1,igrid),amat_NMRnucl(2,igrid),
     &     amat_NMRnucl(3,igrid)
 1       format('OFF:(igrid,nqpts)=(',i8,',',i8,') ',
     &          'xyz_NMRcoords=(',f15.8,',',f15.8,',',f15.8,') ',
     &          'qxyz=(',f15.8,',',f15.8,',',f15.8,') ',
     &          '(dist,zoracutoff_NMR)=(',f15.8,',',f15.8,
     &          ' NMRpot=(',f25.8,',',f25.8,',',f25.8,')')
         endif
 4       continue

          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
       return
       end
      subroutine gridNMRPotential_dia(amat_NMRnucl4dia, ! out: EFG potential
     &                                xyz_NMRcoords,    ! in : EFG-nuclear coord.
     &                                nqpts,            ! in : nr. grid points 
     &                                qxyz,             ! in : grid points                                   
     &                                closegridpts)
      implicit none
#include "global.fh" 
#include "stdio.fh"
#include "zora.fh"
      integer i,igrid,nqpts
      integer closegridpts(*)
      double precision xyz_NMRcoords(3)
      double precision qxyz(3,nqpts)
      double precision rxyz(3),dist,dist3,ac_prod
      double precision amat_NMRnucl4dia(nqpts)
c     == loop over the grid points ==      
       do igrid = 1,nqpts
        amat_NMRnucl4dia(igrid) = 0.0d0
        ac_prod=0.0d0
        do i=1,3
           rxyz(i) = qxyz(i,igrid)-xyz_NMRcoords(i) 
           ac_prod=ac_prod+rxyz(i)*rxyz(i)
        enddo
        dist = dsqrt(ac_prod)
        if (dist.gt.zoracutoff_NMR) then ! check-cutoff
          dist3=dist*dist*dist
          amat_NMRnucl4dia(igrid) = 1.0d0/dist3
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
       return
       end
c $Id$
