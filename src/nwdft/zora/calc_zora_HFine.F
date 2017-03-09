c     +++++++++++++++++++++++++++++++
c     + calculate all Hyperfine AOs +
c     +++++++++++++++++++++++++++++++
c     1. (zpsox,zpsoy,zpsoz):
c     H^{ZPSO}_{mu nu,Aj}= \int dr K/r_A^3 
c                     \vec{r}_A x [chi_{mu}^* \nabla chi_{nu} -
c                                  chi_{nu}^* \nabla chi_{mu}^* ]_j
c     (Eq. 56 in J. Autschbach's write-up of 
c      ZORA-NMR spin-spin coupling constants
c      Sept. 17, 2007's write-up)

c     2. (fcsdxx,fcsdxy,fcsdxz,
c         fcsdyx,fcsdyy,fcsdyz,
c         fcsdzx,fcsdzy,fcsdzz):
c     H^{FC+SD}_{uv}=\int dr K U_{N,v} \nabla_{u} (chi_{mu}^* chi_{nu})
c     where U_{N,v} is,
c     U_{N,v} c^{-2} r_{N,v}/r_N^3  (Eq. 10 in JA's draft of
c     'Calculation of hyperfine tensor using zeroth-order regular approximation
c      and density functional theory: Expectation value versus linear response
c      approaches')

      subroutine calc_zora_HFine_slow(
     &                           ao_bas_han,   ! in: AO basis handle 
     &                           geom,         ! in: geometry handle 
     &                           ipol,         ! in: nr. of polarizations 
     &                           g_dens,       ! in: superposit. atomic densities 
     &                           chi_ao,       ! in:           basis functions
     &                           delchi_ao,    ! in: deriv. of basis functions
     &                           qxyz,         ! in: grid points 
     &                           qwght,        ! in: weighting coeffs.
     &                           nbf,          ! in: nr. basis functions
     &                           npts,         ! in: nr. grid points
     &                           natoms,       ! in: nr. atoms 
     &                           ofinite,      ! in: = .true. if Gaussian Nucl. Model of charges requested
     &                           zetanuc_arr,  ! in: zetanuc(i) i=1,natoms for Gaussian Nuclear Model
     &                           atmass,       ! in: atomic mass
     &                           xyz_NMRcoords,! in : nuclear coordinates
     &                           use_modelpotential,
     &                           gexpo,
     &                           gcoef,
     &                           zpsox,        ! out 
     &                           zpsoy,        ! out 
     &                           zpsoz,        ! out          
     &                           fcsdxx,       ! out
     &                           fcsdxy,       ! out
     &                           fcsdxz,       ! out
c     &                           fcsdyx,       ! out
     &                           fcsdyy,       ! out
     &                           fcsdyz,       ! out
c     &                           fcsdzx,       ! out
c     &                           fcsdzy,       ! out
     &                           fcsdzz)       ! out
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "bas.fh"
#include "zora.fh"
      integer nbf,npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision qxyz(3,npts)
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision zpsox(nbf,nbf),
     &                 zpsoy(nbf,nbf),
     &                 zpsoz(nbf,nbf)
      double precision fcsdxx(nbf,nbf),
     &                 fcsdxy(nbf,nbf),
     &                 fcsdxz(nbf,nbf),
     &                 fcsdyy(nbf,nbf),
     &                 fcsdyz(nbf,nbf),
     &                 fcsdzz(nbf,nbf)
      double precision ac_fcsd(3,3)
      integer i,j,k,n
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts),amat_NMRnucl(3,npts),
     &                 amat_Pnucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot,Kzora
      double precision fac1_arr(npts),fac2_arr(3,npts)
      double precision ac_zpso(3)
      double precision xyz_NMRcoords(3),atmass  
      double precision chi_cntr(3,nbf),threehalf
      data threehalf /1.5d0/
      logical ofinite
c ------- for Gaussian Nuclear Model --- START
      double precision zetanuc_arr(natoms)
c ------- for Gaussian Nuclear Model --- START
c
      logical use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
c
c     == preliminaries ==
      clight_au2 = clight_au*clight_au
      do ipt = 1,npts
        do i=1,ipol
         amat_coul(ipt,i) = 0.d0            
        end do        
        amat_nucl(ipt)  = 0.d0
        amat_Pnucl(ipt) = 0.0d0
        closegridpts(ipt) = 0
        do i=1,3
         amat_NMRnucl(i,ipt) = 0.d0
        enddo
      end do
c
c     == calculate the total hartree potential on the grid ==
      call gridHartreePotential(use_modelpotential,
     &    ao_bas_han, geom, natoms, ipol, g_dens, npts, qxyz, qwght,
     &    closegridpts, gexpo, gcoef, amat_coul)
c
c     == calculate the total nuclear potential on the grid ==
      if (ofinite) then
c ------ Choosing Nuclear Model: erf(zetanuc^0.5 r_L)
        call gridNuclearPotentialFinite(geom,natoms,npts,qxyz,qwght,
     &                             zetanuc_arr,
     &                             closegridpts,
     &                             amat_nucl)
c ------ Choosing Nuclear Model: P(1/2,zetanuc r_L^2)
c        call gridNuclearPotentialFinite2(geom,natoms,npts,qxyz,qwght,
c     &                             closegridpts,amat_nucl)
      else ! default : point charge model for nuclei
        call gridNuclearPotentialPoint(geom,natoms,npts,qxyz,qwght,
     &                            closegridpts,amat_nucl)
      endif
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
      call gridNMRPotential(amat_NMRnucl,  ! out: NMR potential
     &                      xyz_NMRcoords, 
     &                      npts,qxyz,closegridpts)
      if (ofinite) then ! ====> GAUSSIAN charge nuclear model 
       call get_Pnucl(amat_Pnucl,    ! out: P(3/2,r_N^2)
     &                atmass,          ! in : atomic mass
     &                xyz_NMRcoords, ! in : EFG-nuclear coord.
     &                threehalf,
     &                npts,qxyz)
c     === define fac_arr
       do k = 1,npts
c      == assemble hartree and nuclear contributions ==
        amat_tot = amat_nucl(k) + amat_coul(k,1)
        Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (do_NonRel) then             ! remove it after TEST
          Kzora=1.0d0                   ! remove it after TEST
        endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        fac1_arr(k)=Kzora*qwght(k)*amat_Pnucl(k)
        do n=1,3
         fac2_arr(n,k)=fac1_arr(k)*amat_NMRnucl(n,k) ! for ZPSO
        enddo ! end-loop-n
       enddo ! end-loop-k
      else             ! ====> POINT charge nuclear model (default)---START
c     === define fac_arr
       do k = 1,npts
c      == assemble hartree and nuclear contributions ==
        amat_tot = amat_nucl(k) + amat_coul(k,1)
        Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (do_NonRel) then             ! remove it after TEST
          Kzora=1.0d0                   ! remove it after TEST
        endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        fac1_arr(k)=Kzora*qwght(k)
        do n=1,3
         fac2_arr(n,k)=fac1_arr(k)*amat_NMRnucl(n,k) ! for ZPSO
        enddo ! end-loop-n
       enddo ! end-loop-k
      endif             ! ====> POINT charge nuclear model (default)---END
c     == assemble zora correction ==
c ---- full matrix calc -------- START
      do i = 1, nbf
        do j = 1, nbf
          call get_ints_zora_hfine_slow(
     &                      nbf,npts,chi_ao,delchi_ao,i,j,      
     &                      fac2_arr,  
     &                      ac_zpso,  ! out
     &                      ac_fcsd)  ! out
          zpsox(i,j)  = zpsox(i,j)  + ac_zpso(1)
          zpsoy(i,j)  = zpsoy(i,j)  + ac_zpso(2)
          zpsoz(i,j)  = zpsoz(i,j)  + ac_zpso(3)
          fcsdxx(i,j) = fcsdxx(i,j) + ac_fcsd(1,1)
          fcsdxy(i,j) = fcsdxy(i,j) + ac_fcsd(1,2)
          fcsdxz(i,j) = fcsdxz(i,j) + ac_fcsd(1,3)
c         fcsdyx(i,j) = fcsdyx(i,j) + ac_fcsd(2,1)
          fcsdyy(i,j) = fcsdyy(i,j) + ac_fcsd(2,2)
          fcsdyz(i,j) = fcsdyz(i,j) + ac_fcsd(2,3)
c         fcsdzx(i,j) = fcsdzx(i,j) + ac_fcsd(3,1)
c         fcsdzy(i,j) = fcsdzy(i,j) + ac_fcsd(3,2)
          fcsdzz(i,j) = fcsdzz(i,j) + ac_fcsd(3,3)
        enddo ! end-loop-j
      enddo ! end-loop-i
c ---- full matrix calc -------- END
      return
      end

      subroutine get_ints_zora_hfine_slow(
     &                               nbf,       ! in: nr. basis functions
     &                               npts,      ! in: grid points 
     &                               chi_ao,    ! in:           basis functions
     &                               delchi_ao, ! in: deriv. of basis functions
     &                               i,j,       ! in: (i,j) indices for delchi_ao
     &                               fac2_arr,  ! in
     &                               ac_zpso,   ! out : ZPSO  term 
     &                               ac_fcsd)   ! out : FC+SD term (n,m) component
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k,m,n,a,b
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision fac2_arr(3,npts)
      double precision ac_zpso(3),
     &                 ac_fcsd(3,3)
      double precision prod(3),prod1(3)
      integer ind_nab(2,3)
      data ind_nab / 2, 3,  ! nab=123
     &               3, 1,  ! nab=231
     &               1, 2 / ! nab=312
      do n=1,3 ! reset
       ac_zpso(n) = 0.0d0
        do m=1,3
         ac_fcsd(n,m) = 0.0d0
        enddo
      enddo    
      do k = 1, npts
       do n=1,3
        prod1(n) = chi_ao(k,i)*delchi_ao(k,n,j)
     &            +chi_ao(k,j)*delchi_ao(k,n,i)
       enddo ! end-loop-n
       do m=n,3
          do n=1,3
         ac_fcsd(n,m) = ac_fcsd(n,m) +
     &                  fac2_arr(n,k)*prod1(m)         
        enddo ! end-loop-m
       enddo ! end-loop-n
      enddo ! end-loo-k
      do k = 1, npts
       do n=1,3
        prod(n)  = chi_ao(k,i)*delchi_ao(k,n,j)
     &            -chi_ao(k,j)*delchi_ao(k,n,i)
       enddo ! end-loop-n
       do n=1,3
        a=ind_nab(1,n)
        b=ind_nab(2,n)      
        ac_zpso(n) = ac_zpso(n) + 
     &               fac2_arr(a,k)*prod(b)-
     &               fac2_arr(b,k)*prod(a) 
       enddo ! end-loop-n
      enddo ! end-loo-k
      return
      end
c +++++++++++++++++++++++++++++++++++++++
c +++++++++++++++++++++++++++++++++++++++
      subroutine calc_zora_HFine_fast(
     &                           ao_bas_han,   ! in: AO basis handle 
     &                           geom,         ! in: geometry handle 
     &                           ipol,         ! in: nr. of polarizations 
     &                           g_dens,       ! in: superposit. atomic densities 
     &                           chi_ao,       ! in:           basis functions
     &                           delchi_ao,    ! in: deriv. of basis functions
     &                           qxyz,         ! in: grid points 
     &                           qwght,        ! in: weighting coeffs.
     &                           nbf, mbf,ibf,          ! in: nr. basis functions
     &                           npts,         ! in: nr. grid points
     &                           natoms,       ! in: nr. atoms 
     &                           ofinite,      ! in: = .true. if Gaussian Nucl. Model of charges requested
     &                           zetanuc_arr,  ! in: sqrt(zetanuc(i)) i=1,natoms for Gaussian Nuclear Model
     &                           zetanuc_slc,  ! in: zetanuc(i) 
     &                           Knucl,
     &                           xyz_NMRcoords,! in: nuclear coordinates
     &                           use_modelpotential,
     &                           gexpo,
     &                           gcoef,
     &                           zpsoxy,       ! out 
     &                           zpsoyz,       ! out 
     &                           zpsozx,       ! out 
     &                           fcsdxx,       ! out
     &                           fcsdxy,       ! out
     &                           fcsdxz,       ! out
     &                           fcsdyy,       ! out
     &                           fcsdyz,       ! out
     &                           fcsdzz)       ! out
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "bas.fh"
#include "zora.fh"
      integer nbf,npts,ao_bas_han,natoms,geom
      integer mbf,ibf(*)
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision qxyz(3,npts)
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision zpsoxy(nbf,nbf),
     &                 zpsoxz(nbf,nbf),
     &                 zpsoyx(nbf,nbf),
     &                 zpsoyz(nbf,nbf),
     &                 zpsozx(nbf,nbf),
     &                 zpsozy(nbf,nbf)
      double precision fcsdxx(nbf,nbf),
     &                 fcsdxy(nbf,nbf),
     &                 fcsdxz(nbf,nbf),
c     &                 fcsdyx(nbf,nbf),
     &                 fcsdyy(nbf,nbf),
     &                 fcsdyz(nbf,nbf),
c     &                 fcsdzx(nbf,nbf),
c     &                 fcsdzy(nbf,nbf),
     &                 fcsdzz(nbf,nbf)
      double precision ac_fcsd(3,3)
      integer i,j,k,n
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts),amat_NMRnucl(3,npts),
     &                 amat_Pnucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot,Kzora
      double precision fac1_arr(npts),fac2_arr(3,npts)
      double precision ac_zpso(3,3)
      double precision xyz_NMRcoords(3)
      double precision chi_cntr(3,nbf),qxyz1(3)
      double precision threehalf
      data threehalf /1.5/
      logical ofinite,Knucl
c
      double precision zetanuc_arr(natoms),Pnucl
      double precision zetanuc_slc
      integer count_pt ! ONLY for checking get_Pnucl
c
      integer i0,j0
      logical use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
c
c     == preliminaries ==
      clight_au2 = clight_au*clight_au
      do ipt = 1,npts
        do i=1,ipol
         amat_coul(ipt,i) = 0.d0            
        end do        
        amat_nucl(ipt) = 0.d0
        closegridpts(ipt) = 0
        do i=1,3
         amat_NMRnucl(i,ipt) = 0.d0
        enddo
      end do
c
c     == calculate the total hartree potential on the grid ==
      call gridHartreePotential(use_modelpotential,
     &    ao_bas_han, geom, natoms, ipol, g_dens, npts, qxyz, qwght,
     &    closegridpts, gexpo, gcoef, amat_coul)
c
c     == calculate the total nuclear potential on the grid ==
      if (ofinite) then
c
c ------ Choosing Nuclear Model: erf(zetanuc^0.5 r_L)
        call gridNuclearPotentialFinite(geom,natoms,npts,qxyz,qwght,
     &                         zetanuc_arr,
     &                         closegridpts,amat_nucl)
c ------ Choosing Nuclear Model: P(1/2,zetanuc r_L^2)
c        call gridNuclearPotentialFinite2(geom,natoms,npts,qxyz,qwght,
c     &                             closegridpts,amat_nucl)
      else ! default : point charge model for nuclei
c
        call gridNuclearPotentialPoint(geom,natoms,npts,qxyz,qwght,
     &                            closegridpts,amat_nucl)
      endif
c
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c
      call gridNMRPotential(amat_NMRnucl,  ! out: NMR potential
     &                      xyz_NMRcoords, 
     &                      npts,qxyz,closegridpts)
      if (ofinite) then ! ====> GAUSSIAN charge nuclear model 
       if (Knucl) then !-- V=Vnucl     (amat_tot) 
        count_pt=1 ! ONLY for checking get_Pnucl
        do k = 1,npts
         qxyz1(1)=qxyz(1,k)
         qxyz1(2)=qxyz(2,k)
         qxyz1(3)=qxyz(3,k)
         call get_Pnucl1(Pnucl,         ! out: P(3/2,r_N^2)
     &                   zetanuc_slc,   ! in : atomic mass
     &                   xyz_NMRcoords, ! in : EFG-nuclear coord.
     &                   threehalf,
     &                   npts,qxyz1,count_pt)
c      == assemble hartree and nuclear contributions ==
         amat_tot = amat_nucl(k) ! V = Vnucl (ONLY)
         Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (do_NonRel) then             ! remove it after TEST
           Kzora=1.0d0                   ! remove it after TEST
         endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         fac1_arr(k)=Kzora*qwght(k)*Pnucl
         do n=1,3
          fac2_arr(n,k)=fac1_arr(k)*amat_NMRnucl(n,k) ! for ZPSO
         enddo ! end-loop-n
        enddo ! end-loop-k
       else ! ------------ V=Vnucl+Vee (amat_tot) (default)
c     === define fac_arr
        count_pt=1 ! ONLY for checking get_Pnucl
        do k = 1,npts
         qxyz1(1)=qxyz(1,k)
         qxyz1(2)=qxyz(2,k)
         qxyz1(3)=qxyz(3,k)
         call get_Pnucl1(Pnucl,         ! out: P(3/2,r_N^2)
     &                   zetanuc_slc,   ! in : atomic mass
     &                   xyz_NMRcoords, ! in : EFG-nuclear coord.
     &                   threehalf,
     &                   npts,qxyz1,count_pt)
c      == assemble hartree and nuclear contributions ==
         amat_tot = amat_nucl(k) + amat_coul(k,1) ! V=Vnucl+Vee
         Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (do_NonRel) then             ! remove it after TEST
           Kzora=1.0d0                   ! remove it after TEST
         endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         fac1_arr(k)=Kzora*qwght(k)*Pnucl
         do n=1,3
          fac2_arr(n,k)=fac1_arr(k)*amat_NMRnucl(n,k) ! for ZPSO
         enddo ! end-loop-n
        enddo ! end-loop-k
       endif ! end-if-Knucl
      else              ! ====> POINT charge nuclear model (default)---START
       if (Knucl) then !-- V=Vnucl     (amat_tot) 
c     === define fac_arr
         do k = 1,npts
c      == assemble hartree and nuclear contributions ==
          amat_tot = amat_nucl(k) ! V=Vnucl (ONLY)
          Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if (do_NonRel) then             ! remove it after TEST
            Kzora=1.0d0                   ! remove it after TEST
          endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
          fac1_arr(k)=Kzora*qwght(k)
          do n=1,3
           fac2_arr(n,k)=fac1_arr(k)*amat_NMRnucl(n,k) ! for ZPSO
          enddo ! end-loop-n
         enddo ! end-loop-k
       else ! ------------ V=Vnucl+Vee (amat_tot) (default)
c     === define fac_arr
         do k = 1,npts
c      == assemble hartree and nuclear contributions ==
          amat_tot = amat_nucl(k) + amat_coul(k,1) ! V=Vnucl+Vee
          Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2) 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if (do_NonRel) then             ! remove it after TEST
            Kzora=1.0d0                   ! remove it after TEST
          endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
          fac1_arr(k)=Kzora*qwght(k)
          do n=1,3
           fac2_arr(n,k)=fac1_arr(k)*amat_NMRnucl(n,k) ! for ZPSO
          enddo ! end-loop-n
         enddo ! end-loop-k
       endif ! end-if-Knucl
      endif              ! ====> POINT charge nuclear model (default)---END
c     == assemble zora correction ==
c ---- full matrix calc -------- START
      do i0 = 1, mbf
         i=ibf(i0)
        do j0 = i0, mbf
           j=ibf(j0)
          call get_ints_zora_hfine_fast(
     &                  nbf,npts,chi_ao,delchi_ao,i0,j0,      
     &                  fac2_arr,  
     &                  ac_zpso,  ! out
     &                  ac_fcsd)  ! out
          zpsoxy(i,j) = zpsoxy(i,j) + ac_zpso(1,2)
          zpsoxz(i,j) = zpsoxz(i,j) + ac_zpso(1,3)
          zpsoyx(i,j) = zpsoyx(i,j) + ac_zpso(2,1)
          zpsoyz(i,j) = zpsoyz(i,j) + ac_zpso(2,3)
          zpsozx(i,j) = zpsozx(i,j) + ac_zpso(3,1)
          zpsozy(i,j) = zpsozy(i,j) + ac_zpso(3,2)   
          fcsdxx(i,j) = fcsdxx(i,j) + ac_fcsd(1,1)
          fcsdxy(i,j) = fcsdxy(i,j) + ac_fcsd(1,2)
          fcsdxz(i,j) = fcsdxz(i,j) + ac_fcsd(1,3)
c          fcsdyx(i,j) = fcsdyx(i,j) + ac_fcsd(2,1)
          fcsdyy(i,j) = fcsdyy(i,j) + ac_fcsd(2,2)
          fcsdyz(i,j) = fcsdyz(i,j) + ac_fcsd(2,3)
c          fcsdzx(i,j) = fcsdzx(i,j) + ac_fcsd(3,1)
c          fcsdzy(i,j) = fcsdzy(i,j) + ac_fcsd(3,2)
          fcsdzz(i,j) = fcsdzz(i,j) + ac_fcsd(3,3)
        enddo ! end-loop-j
      enddo ! end-loop-i
crecover upper triangle
      do i0 = 1, mbf
         i=ibf(i0)
        do j0 = i0+1, mbf
           j=ibf(j0)
          zpsoxy(j,i) = zpsoxy(i,j)
          zpsoxz(j,i) = zpsoxz(i,j)
          zpsoyx(j,i) = zpsoyx(i,j)
          zpsoyz(j,i) = zpsoyz(i,j)
          zpsozx(j,i) = zpsozx(i,j)
          zpsozy(j,i) = zpsozy(i,j)
          fcsdxx(j,i) = fcsdxx(i,j)
          fcsdxy(j,i) = fcsdxy(i,j)
          fcsdxz(j,i) = fcsdxz(i,j)
          fcsdyy(j,i) = fcsdyy(i,j)
          fcsdyz(j,i) = fcsdyz(i,j)
          fcsdzz(j,i) = fcsdzz(i,j)
        enddo ! end-loop-j
      enddo ! end-loop-i
c ---- full matrix calc -------- END
      return
      end

      subroutine get_ints_zora_hfine_fast(
     &                               nbf,       ! in: nr. basis functions
     &                               npts,      ! in: grid points 
     &                               chi_ao,    ! in:           basis functions
     &                               delchi_ao, ! in: deriv. of basis functions
     &                               i,j,       ! in: (i,j) indices for delchi_ao
     &                               fac2_arr,  ! in
     &                               ac_zpso,   ! out : ZPSO  term 
     &                               ac_fcsd)   ! out : FC+SD term (n,m) component
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k,m,n,a,b
      double precision chi_ao(npts,nbf)
      double precision delchi_ao(npts,3,nbf)
      double precision fac2_arr(3,npts)
      double precision ac_zpso(3,3),
     &                 ac_fcsd(3,3)
      double precision prod(3),prod1(3),prod13
      integer ind_nab(2,3)
      data ind_nab / 2, 3,  ! nab=123
     &               3, 1,  ! nab=231
     &               1, 2 / ! nab=312
      do n=1,3 ! reset
        do m=1,3
         ac_zpso(n,m) = 0.0d0
         ac_fcsd(n,m) = 0.0d0
        enddo
      enddo    
      do k = 1, npts
            prod13 = chi_ao(k,i)*delchi_ao(k,1,j)
     &           +chi_ao(k,j)*delchi_ao(k,1,i)
               ac_fcsd(1,1) = ac_fcsd(1,1)+fac2_arr(1,k)*prod13
               ac_fcsd(1,2) = ac_fcsd(1,2)+fac2_arr(2,k)*prod13
               ac_fcsd(1,3) = ac_fcsd(1,3)+fac2_arr(3,k)*prod13
            prod13 = chi_ao(k,i)*delchi_ao(k,2,j)
     &           +chi_ao(k,j)*delchi_ao(k,2,i)
               ac_fcsd(2,2) = ac_fcsd(2,2)+fac2_arr(2,k)*prod13
               ac_fcsd(2,3) = ac_fcsd(2,3)+fac2_arr(3,k)*prod13

            prod13 = chi_ao(k,i)*delchi_ao(k,3,j)
     &           +chi_ao(k,j)*delchi_ao(k,3,i)
               ac_fcsd(3,3) = ac_fcsd(3,3)+fac2_arr(3,k)*prod13
      enddo ! end-loo-k
#if 0
      do k = 1, npts
       do n=1,3
        prod(n)  = chi_ao(k,i)*delchi_ao(k,n,j)
       enddo ! end-loop-n
       do n=1,3
         a=ind_nab(1,n)
         b=ind_nab(2,n)      
         ac_zpso(a,b) = ac_zpso(a,b)+fac2_arr(a,k)*prod(b)
         ac_zpso(b,a) = ac_zpso(b,a)+fac2_arr(b,k)*prod(a)
       enddo ! end-loop-n
      enddo ! end-loo-k
#else
      do n=1,3
         a=ind_nab(1,n)
         b=ind_nab(2,n)      
         do k = 1, npts
            ac_zpso(a,b) = ac_zpso(a,b)+fac2_arr(a,k)*
     B           chi_ao(k,i)*delchi_ao(k,b,j)
            ac_zpso(b,a) = ac_zpso(b,a)+fac2_arr(b,k)*
     A           chi_ao(k,i)*delchi_ao(k,a,j)
         enddo                  ! end-loop-n
      enddo                     ! end-loo-k
#endif
      return
      end

      subroutine calc_NMRHFine_F1ij(
     &                           ao_bas_han,   ! in: AO basis handle 
     &                           geom,         ! in: geometry handle 
     &                           ipol,         ! in: nr. of polarizations 
     &                           g_dens,       ! in: superposit. atomic densities 
     &                           delchi_ao,    ! in: deriv. of basis functions
     &                           qxyz,         ! in: grid points 
     &                           qwght,        ! in: weighting coeffs.
     &                           nbf,          ! in: nr. basis functions
     &                           npts,         ! in: nr. grid points
     &                           natoms,       ! in: nr. atoms 
     &                           ofinite,      ! in: = .true. if Gaussian Nucl. Model of charges requested
     &                           zetanuc_arr,  ! in: sqrt(zetanuc(i)) i=1,natoms for Gaussian Nuclear Model
     &                           Knucl,
     &                           use_modelpotential,
     &                           gexpo,
     &                           gcoef,
     &                           zsrkineticx,  ! out 
     &                           zsrkineticy,  ! out 
     &                           zsrkineticz)  ! out
c Purpose : Evaluates AO matrix for operator
c           [\vec{p} K x \vec{p}]_v  v=1,2,3=x,y,z
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "bas.fh"
#include "zora.fh"
      integer nbf,npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision qxyz(3,npts)
      double precision delchi_ao(npts,3,nbf)
      double precision zsrkineticx(nbf,nbf),
     &                 zsrkineticy(nbf,nbf),
     &                 zsrkineticz(nbf,nbf)
      integer i,j,k,n
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot,Kzora
      double precision fac1_arr(npts)
      double precision ac_hfineF1ji(3)
      double precision chi_cntr(3,nbf)
      logical ofinite,Knucl
      double precision zetanuc_arr(natoms)
      external get_ints_zora_hfine_F1ji, 
     &         gridNuclearPotentialFinite,
     &         gridNuclearPotentialPoint
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
      if (ofinite) then
c
c ------ Choosing Nuclear Model: erf(zetanuc^0.5 r_L)
        call gridNuclearPotentialFinite(geom,natoms,npts,qxyz,qwght,
     &                             zetanuc_arr,
     &                             closegridpts,amat_nucl)
c ------ Choosing Nuclear Model: P(1/2,zetanuc r_L^2)
c        call gridNuclearPotentialFinite2(geom,natoms,npts,qxyz,qwght,
c     &                             closegridpts,amat_nucl)
      else ! default : point charge model for nuclei
        call gridNuclearPotentialPoint(geom,natoms,npts,qxyz,qwght,
     &                            closegridpts,amat_nucl)
      endif
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c     === define fac_arr
      if (Knucl) then !-- V=Vnucl     (amat_tot) 
       do k = 1,npts
c      == assemble hartree and nuclear contributions ==
        amat_tot = amat_nucl(k)
        Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2)-1.0d0 ! Alternative expression
                                                            ! gives same value as K for this AO
                                                            ! but it is suppose to cancel noise
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (do_NonRel) then             ! remove it after TEST
           Kzora=1.0d0                   ! remove it after TEST
         endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         fac1_arr(k)=Kzora*qwght(k)
       enddo ! end-loop-k
      else ! default  V=Vnucl+Vhartee
       do k = 1,npts
c      == assemble hartree and nuclear contributions ==
        amat_tot = amat_nucl(k) + amat_coul(k,1)
        Kzora=1.0d0/(1.0d0-0.5d0*amat_tot/clight_au2)-1.0d0 ! Alternative expression
                                                            ! gives same value as K for this AO
                                                            ! but it is suppose to cancel noise
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (do_NonRel) then             ! remove it after TEST
           Kzora=1.0d0                   ! remove it after TEST
         endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
         fac1_arr(k)=Kzora*qwght(k)
       enddo ! end-loop-k
      endif
c     == assemble zora correction ==
c ---- main diagonal -------- START
      do i = 1, nbf
          j=i
          call get_ints_zora_hfine_F1ji(nbf,npts,delchi_ao,i,j,      
     &                             fac1_arr,  
     &                             ac_hfineF1ji)  ! out
          zsrkineticx(i,j)  = zsrkineticx(i,j)  + ac_hfineF1ji(1)
          zsrkineticy(i,j)  = zsrkineticy(i,j)  + ac_hfineF1ji(2)
          zsrkineticz(i,j)  = zsrkineticz(i,j)  + ac_hfineF1ji(3)
      enddo ! end-loop-i
c ---- main diagonal -------- END
c ---- off diagonal -------- START
      do i = 1, nbf
        do j = i+1, nbf
          call get_ints_zora_hfine_F1ji(nbf,npts,delchi_ao,i,j,      
     &                                  fac1_arr,  
     &                                  ac_hfineF1ji)  ! out
          zsrkineticx(i,j)  = zsrkineticx(i,j)  + 2.0d0*ac_hfineF1ji(1)
          zsrkineticy(i,j)  = zsrkineticy(i,j)  + 2.0d0*ac_hfineF1ji(2)
          zsrkineticz(i,j)  = zsrkineticz(i,j)  + 2.0d0*ac_hfineF1ji(3)
        enddo ! end-loop-j
      enddo ! end-loop-i
c ---- off diagonal -------- END
      return
      end

      subroutine get_ints_zora_hfine_F1ji(
     &                          nbf,         ! in: nr. basis functions
     &                          npts,        ! in: grid points 
     &                          delchi_ao,   ! in: deriv. of basis functions
     &                          i,j,         ! in: (i,j) indices for delchi_ao
     &                          fac1_arr,    ! in
     &                          ac_hfineF1ji)! out
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,k,n,i,j
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts)
      double precision ac_hfineF1ji(3)
      double precision prod(3)

      do n=1,3 ! reset
      ac_hfineF1ji(n) = 0.0d0
      enddo    
      do k = 1, npts
       prod(1)= delchi_ao(k,2,i)*delchi_ao(k,3,j)
     &         -delchi_ao(k,3,i)*delchi_ao(k,2,j) 
       prod(2)= delchi_ao(k,3,i)*delchi_ao(k,1,j)
     &         -delchi_ao(k,1,i)*delchi_ao(k,3,j) 
       prod(3)= delchi_ao(k,1,i)*delchi_ao(k,2,j)
     &         -delchi_ao(k,2,i)*delchi_ao(k,1,j) 
       do n=1,3
        ac_hfineF1ji(n) = ac_hfineF1ji(n) + fac1_arr(k)*prod(n)
       enddo ! end-loop-n
      enddo ! end-loo-k
      return
      end
c
      subroutine get_Pnucl(amat_Pnucl,    ! out: potential
     &                     atmass,        ! in : atomic mass
     &                     xyz_NMRcoords, ! in : EFG-nuclear coord.
     &                     a_coeff,       ! in : =3/2 for AOs =1/2 for Vnucl
     &                     nqpts,         ! in : nr. grid points 
     &                     qxyz)
c    About: a_coeff,       
c    =3/2 when adding finite size charge Gaussian model in evaluation
c         of hyperfine AOs (calc_zora_HFine)
c    =1/2 for Vnucl (in gridNuclearPotential)
      implicit none
#include "geom.fh"
#include "global.fh" 
#include "msgids.fh"
#include "stdio.fh"
      integer i,igrid,nqpts
      double precision xyz_NMRcoords(3)
      double precision qxyz(3,nqpts)
      double precision rxyz(3),dist,dist2,ac_prod
      double precision amat_Pnucl(nqpts)
      character*16 element
      character*2  symbol
      character*16 tags
      logical is_atom
      double precision atmass,zetanuc
      double precision rtemp,a_coeff
c
      double precision dgami
      external dgami, 
     &         get_znuc
c
      call get_znuc(atmass,zetanuc)
      do igrid = 1,nqpts
        ac_prod=0.0d0
        do i=1,3
         rxyz(i) = qxyz(i,igrid)-xyz_NMRcoords(i) 
         ac_prod=ac_prod+rxyz(i)*rxyz(i)
        enddo
        rtemp = zetanuc*ac_prod  ! dist*dist
        amat_Pnucl(igrid) = dgami(a_coeff,rtemp)  ! P(3/2,\tilde{r}_N^2)
      end do ! igrid
c
      return
      end
c
c------- get_Pnucl1() ------------ START
c Purpose : Get one single value of Pnucl
      subroutine get_Pnucl1(Pnucl,        ! out: potential
     &                      zetanuc_slc,  ! in : zetanuc
     &                      xyz_NMRcoords,! in : EFG-nuclear coord.
     &                      a_coeff,      ! in : =3/2 for AOs =1/2 for Vnucl
     &                      nqpts,        ! in : nr. grid points 
     &                      qxyz,         ! in : one single grid point
     &                      count_pt)     ! TO CHECK
c    About: a_coeff,       
c    =3/2 when adding finite size charge Gaussian model in evaluation
c         of hyperfine AOs (calc_zora_HFine)
c    =1/2 for Vnucl (in gridNuclearPotential)
      implicit none
#include "msgids.fh"
#include "stdio.fh"
#include "global.fh"
   
      integer count_pt ! to check

      integer i,igrid,nqpts
      double precision xyz_NMRcoords(3)
      double precision qxyz(3),Pnucl
      double precision rxyz(3),dist,dist2,ac_prod
      double precision zetanuc_slc
      double precision rtemp,a_coeff
      double precision dgami
      external dgami ! Incomplete Gamma
c ---------- Defining values at hand to check --------- START    

c      xyz_NMRcoords(1)=  0.07090063d0
c      xyz_NMRcoords(2)= -0.12532286d0
c      xyz_NMRcoords(3)=  0.00000000d0
c      qxyz(1)= 0.07090403d0
c      qxyz(2)=-0.12532425d0
c      qxyz(3)=-0.00000339d0
c      zetanuc_slc=140130060.38598028d0
cNWChem: Grid coord=(0.07090403,    -0.12532425,    -0.00000339)
c        Nuclpos   =(0.07090063,    -0.12532286,     0.00000000)
c        (zetanuc,rtemp,gammap)=(140130060.38598028,0.00350104,0.00015551)
cADF:    Grid coord=(0.07090403,    -0.12532425,    -0.00000339) 
c        Nuclpos   =(0.07090063,    -0.12532286,     0.00000000)
c        (zetanuc,rtemp,gammap)=(140130060.38598028,0.00349682,0.00015523)
c ---------- Defining values at hand to check --------- ENd
        ac_prod=0.0d0
        do i=1,3
         rxyz(i) = qxyz(i)-xyz_NMRcoords(i) 
         ac_prod=ac_prod+rxyz(i)*rxyz(i)
        enddo
        rtemp = zetanuc_slc*ac_prod  ! dist*dist
        Pnucl = dgami(a_coeff,rtemp)  ! P(3/2,\tilde{r}_N^2)
       return
       end
c
c------- get_Pnucl1() ------------ END
c     Purpose: Evaluation of zetanuc
c              to be used in evaluation of Incomplete
c              Gamma Function [gratio(...)]
c              rtemp = zetanuc*ac_prod  ! dist*dist
c              call dgratio(threehalf,rtemp,gammap,gammaq,0)
c
       subroutine get_znuc(atmass, zetanuc)
c
       implicit none
#include "msgids.fh"
#include "stdio.fh"
c
       double precision atmass
       double precision parnuc1,parnuc2,
     &      one,two,three,threehalf,
     &      fm2bohr,zetanuc,rtemp
       data parnuc1   /0.836d0/
       data parnuc2   /0.570d0/
       data one       /1.0d0/
       data two       /2.0d0/
       data three     /3.0d0/
       data threehalf /1.5d0/
       data fm2bohr   /52917.7249d0/
       rtemp = (parnuc1 * atmass**(one/three) + parnuc2) / fm2bohr
       zetanuc = three / ( two * (rtemp**2))
       return
       end
c   
c     Purpose: Evaluation of zetanuc arr
c              to be used in evaluation of Incomplete
c              Gamma Function [gratio(...)]
c              rtemp = zetanuc*ac_prod  ! dist*dist
c              call dgratio(threehalf,rtemp,gammap,gammaq,0)
c              This routine is used in gridNuclearPotential()
c
       subroutine get_zetanuc_arr(geom, natoms, zetanuc_arr) 
c
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "global.fh" 
#include "msgids.fh"
#include "stdio.fh"
        integer geom ! handle for geometry
        integer i,natoms
        double precision zetanuc_arr(natoms),
     &                   atmass,zetanuc
        external get_znuc
        do i=1,natoms
          if(.not.geom_mass_get(geom,i,atmass)) call
     &       errquit(' mass_get  failed ',i,GEOM_ERR)   
          call get_znuc(atmass,zetanuc)
          zetanuc_arr(i)=zetanuc      
        enddo ! end-loop-i
       return
       end   
c $Id$
