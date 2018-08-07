C>
C> \ingroup cosmo
C> @{
C>
C> \file cosmo_charges.F
C> The COSMO surface generator
C>
C> \brief Compute the COSMO charges
C>
      subroutine cosmo_charges(rtdb,basis,geom,ecos,some,
     &         itype_wfn,vectors,nclosed,nopen,nbf,nmo)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "prop.fh"
#include "cosmo_params.fh"
c
      integer rtdb          !< [Input] The RTDB handle
      integer basis         !< [Input] The basis set handle
      integer geom          !< [Input] The geometry handle
      integer nbf           !< [Input] The number of basis functions
      integer nmo           !< [Input] The number of molecular orbitals
      integer itype_wfn     !< [Input] The wavefunction type
                            !< - 1: RHF type wavefunction
                            !< - 2: UHF type wavefunction
      integer vectors(2)    !< [Input] The GA handles for the SCF vectors
      integer nclosed(2)    !< [Input] The number of closed shell orbitals
                            !< (only used with RHF)
      integer nopen(2)      !< [Input] The number of open shell orbitals
c
      double precision ecos !< [Output] The intra COSMO charges
                            !< interaction energy
      logical some          !< [Input] Do you want "some" additional output?
      logical dbug, out, status
      character*8 scfruntyp
c
      character*255 cosmo_file
c
c
      integer l_efcc, k_efcc, l_efcs, k_efcs, l_efcz, k_efcz
      integer l_efciat, k_efciat
      integer l_rad,  k_rad,  nrad
      integer l_occ,  k_occ
      integer nefc, ief, i, nat
      integer g_dens(3)  ! 1: up, 2: down, 3: total
      integer ndens,nvirt(2)
      logical stat
c
      double precision bohr
      parameter (bohr=0.529177249d+00)
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      dbug=.false..and.ga_nodeid().eq.0
      out =.false..and.ga_nodeid().eq.0
      out =out.or.dbug
      some=some.or.out.and.ga_nodeid().eq.0
c
      if(some) then
         write(luout,9999)
      endif
c
      if(out) then
         write(luout,*) 'in cosmo_charges ... geom = ',geom
      endif
c
c     ----- retrieve the number of atoms from -geom- -----
c
      if (.not.geom_ncent(geom,nat))
     &   call errquit("cosmo_charges: geom_ncent failed",0,UERR)
c
c     ----- read -efc- coordinates from -rtdb- -----
c
      if(.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1     ,nefc))
     &   call errquit('cosmo_charges: rtdb get failed for nefc  ',911,
     &       RTDB_ERR)
c
c     ----- allocate memory for efc's -----
c
      if(.not.ma_push_get(mt_dbl,nefc*3,'cosmo efcc',l_efcc,k_efcc))
     & call errquit('cosmo_charges malloc k_efcc failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcs',l_efcs,k_efcs))
     & call errquit('cosmo_charges malloc k_efcs failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcz',l_efcz,k_efcz))
     & call errquit('cosmo_charges malloc k_efcz failed',911,MA_ERR)
      if(.not.ma_push_get(mt_int,nefc,'cosmo efciat',l_efciat,k_efciat))
     & call errquit('cosmo_charges malloc k_efciat failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nat,'cosmo rad',l_rad,k_rad))
     & call errquit('cosmo_charges malloc k_rad failed',nat,MA_ERR)
c
      if(.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,dbl_mb(k_efcc)))
     &   call errquit('cosmo_charges: rtdb get failed for efcc  ',912,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,  nefc,dbl_mb(k_efcz)))
     &   call errquit('cosmo_charges: rtdb get failed for efcz  ',913,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcs',mt_dbl,  nefc,dbl_mb(k_efcs)))
     &   call errquit('cosmo_charges: rtdb get failed for efcs  ',914,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efciat',
     &                 mt_int,nefc,int_mb(k_efciat)))
     &   call errquit('cosmo_charges: rtdb get failed for efciat',914,
     &       RTDB_ERR)
      call cosmo_def_radii(rtdb,geom,nat,dbl_mb(k_rad))
      status = rtdb_get(rtdb,'cosmo:radius',mt_dbl, nat,dbl_mb(k_rad))
      do i = 0, nat-1
        dbl_mb(k_rad+i) = dbl_mb(k_rad+i)/bohr
      enddo
c
      if(out) then
         write(luout,*) 'in cosmo_charges, nefc = ',nefc
         do ief=1,nefc
            write(luout,*) dbl_mb(k_efcc+(ief-1)*3), 
     &            dbl_mb(k_efcc+(ief-1)*3+1),dbl_mb(k_efcc+(ief-1)*3+2)
         enddo
      endif
c
c      ----- calculate the density matrices from molecular orbitals -----
cc
cc     ----- get density matrix -----
cc
c      scfruntyp='RHF'
c      nocc=nclosed(1)+nopen(1)
c      if (itype_wfn.eq.2) then
c         scfruntyp='UHF'
c         nocc=max(nocc,nclosed(2)+nopen(2))
c      endif
c      if(.not.ma_push_get(mt_dbl,nocc*2,'cosmo occ',l_occ,k_occ))
c     &   call errquit('cosmo_charges malloc k_occ failed',911,MA_ERR)
c      do i=1,nocc*2
c         dbl_mb(k_occ+i-1)=1.0d0
c      enddo
c
c     changes to make cosmo work for open shell DFT (MV)
c
      if(.not.ma_push_get(mt_dbl,nbf*2,'cosmo occ',l_occ,k_occ))
     &   call errquit('cosmo_charges malloc k_occ failed',911,MA_ERR)
         call dfill(2*nbf, 0.0d0, dbl_mb(k_occ), 1)
c
c     Only need to set occupation numbers for UHF
c     occupation numbers for RHF are done inside hnd_prop_dens_make
c
      if (itype_wfn.eq.2) then
         scfruntyp='UHF'
         do i = 1, nopen(1)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
         do i = nbf+1, nbf+nopen(2)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
      else if (itype_wfn.eq.1) then
         scfruntyp='RHF'
         do i = 1, nclosed(1)
            dbl_mb(i-1+k_occ) = 2.0d0
         enddo
         do i = nclosed(1)+1, nclosed(1)+nopen(1)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
      else
         call errquit("unknown function type",0,0)
      endif
c
c     end of changes (MV)
c
      call hnd_prop_dens_make(rtdb,geom,basis,nbf,nmo,nclosed,nopen,
     &                        nvirt,scfruntyp,vectors,dbl_mb(k_occ),
     &                        g_dens,ndens)
c
      if(out) then
         write(luout,*) 'in -cosmo_charges, wfntyp = ',itype_wfn
         write(luout,*) 'in -cosmo_charges, scftyp = ',scfruntyp
         if(dbug) then
            write(luout,*) 'in -cosmo_charges, -da- ...'
            call ga_print(g_dens(1))
            if(scfruntyp.eq.'UHF     ') then
               write(luout,*) 'in -cosmo_charges, -db- ...'
               call ga_print(g_dens(2))
            endif
         endif ! dbug
      endif !out
c
c     ----- get -cosmo- charges -----
c
      cosmo_file = "cosmo.xyz"  ! default name
      call hnd_coschg(g_dens,ndens,rtdb,geom,basis,nat,nefc,
     &                dbl_mb(k_efcc),dbl_mb(k_efcs),dbl_mb(k_efcz),
     &                int_mb(k_efciat),dbl_mb(k_rad),ecos,cosmo_file)
c
c     ----- release memory block -----
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('cosmo_charges: ga_destroy failed g_dens',0,GA_ERR)
      enddo
      if(.not.ma_chop_stack(l_efcc))
     & call errquit('cosmo_charges, ma_chop_stack of l_efcc failed',911,
     &       ma_err)
c
      return
 9999 format(/,10X,15(1H-),
     1       /,10X,'-cosmo- charges',
     2       /,10X,15(1H-))
      end
c
C> \brief Compute the COSMO charges from an input density matrix
c
      subroutine cosmo_charges_from_dmat(rtdb,basis,geom,ecos,some,
     &         ipol,g_dens_in,cosmo_file)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "prop.fh"
c
      integer rtdb          !< [Input] The RTDB handle
      integer basis         !< [Input] The basis set handle
      integer geom          !< [Input] The geometry handle
      double precision ecos !< [Output] The intra COSMO charges, interaction energy
      logical some          !< [Input] Do you want "some" additional output?
      integer ipol          !< [Input] 1: RHF type wavefunction, 2: UHF type wavefunction
      integer g_dens_in(ipol)  !< [Input] 1: up, 2: down  (input density matrix to get the charges)
c
      logical dbug, out, status
c
!      double precision dielec,dielecinf,screen,rsolv,adiag,dsurf,dvol
!      double precision srfmol,volmol,zeta,ptspatm,gammas,swtol
!      integer do_cosmo_model,lineq,minbem,maxbem,ificos,cosmo_sccor
!      common/hnd_cospar/dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
!      common/hnd_cosmod/do_cosmo_model,lineq,minbem,maxbem,ificos,
!     &                  cosmo_sccor
!      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol,ptspatm
#include "cosmo_params.fh"
c
      integer l_efcc, k_efcc, l_efcs, k_efcs, l_efcz, k_efcz
      integer l_efciat, k_efciat
      integer l_rad,  k_rad,  nrad
      integer l_occ,  k_occ
      integer nefc, ief, i, nat
      integer g_dens(3)  ! 1: up, 2: down, 3: total
      integer ndens,nvirt(2)
      logical stat
c
      character*255 cosmo_file
c
      double precision bohr
      parameter (bohr=0.529177249d+00)
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      dbug=.false..and.ga_nodeid().eq.0
      out =.false..and.ga_nodeid().eq.0
      out =out.or.dbug
      some=some.or.out.and.ga_nodeid().eq.0
c
      if(some) then
         write(luout,9999)
      endif
c
      if(out) then
         write(luout,*) 'in cosmo_charges ... geom = ',geom
      endif
c
c     ----- retrieve the number of atoms from -geom- -----
c
      if (.not.geom_ncent(geom,nat))
     &   call errquit("cosmo_charges: geom_ncent failed",0,UERR)
c
c     ----- read -efc- coordinates from -rtdb- -----
c
      if(.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1     ,nefc))
     &   call errquit('cosmo_charges: rtdb get failed for nefc  ',911,
     &       RTDB_ERR)
c
c     ----- allocate memory for efc's -----
c
      if(.not.ma_push_get(mt_dbl,nefc*3,'cosmo efcc',l_efcc,k_efcc))
     & call errquit('cosmo_charges malloc k_efcc failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcs',l_efcs,k_efcs))
     & call errquit('cosmo_charges malloc k_efcs failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcz',l_efcz,k_efcz))
     & call errquit('cosmo_charges malloc k_efcz failed',911,MA_ERR)
      if(.not.ma_push_get(mt_int,nefc,'cosmo efciat',l_efciat,k_efciat))
     & call errquit('cosmo_charges malloc k_efciat failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nat,'cosmo rad',l_rad,k_rad))
     & call errquit('cosmo_charges malloc k_rad failed',nat,MA_ERR)
c
      if(.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,dbl_mb(k_efcc)))
     &   call errquit('cosmo_charges: rtdb get failed for efcc  ',912,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,  nefc,dbl_mb(k_efcz)))
     &   call errquit('cosmo_charges: rtdb get failed for efcz  ',913,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcs',mt_dbl,  nefc,dbl_mb(k_efcs)))
     &   call errquit('cosmo_charges: rtdb get failed for efcs  ',914,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efciat',
     &                 mt_int,nefc,int_mb(k_efciat)))
     &   call errquit('cosmo_charges: rtdb get failed for efciat',914,
     &       RTDB_ERR)
      call cosmo_def_radii(rtdb,geom,nat,dbl_mb(k_rad))
      status = rtdb_get(rtdb,'cosmo:radius',mt_dbl, nat,dbl_mb(k_rad))
      do i = 0, nat-1
        dbl_mb(k_rad+i) = dbl_mb(k_rad+i)/bohr
      enddo
c
      if(out) then
         write(luout,*) 'in cosmo_charges, nefc = ',nefc
         do ief=1,nefc
            write(luout,*) dbl_mb(k_efcc+(ief-1)*3), 
     &            dbl_mb(k_efcc+(ief-1)*3+1),dbl_mb(k_efcc+(ief-1)*3+2)
         enddo
      endif
c
c     ----- create density matrix handles -----
c
      if (ipol.eq.1) ndens = 1
      if (ipol.eq.2) ndens = 3
c
      do i = 1, ndens
        g_dens(i) = ga_create_atom_blocked(geom,basis,'density matrix')
        call ga_zero(g_dens(i))
      end do
c
      if (ipol.eq.2) then
           call ga_copy(g_dens_in(1),g_dens(1))
           call ga_copy(g_dens_in(2),g_dens(2))
           call ga_dadd (1.D0,g_dens(1),1.D0,g_dens(2),g_dens(3)) ! sum up components
      else if (ipol.eq.1) then
           call ga_copy(g_dens_in(1),g_dens(1))
      else
           call errquit("unknown wavefunction function type",0,0)
      endif
c
c     ----- get -cosmo- charges -----
c
      call hnd_coschg(g_dens,ndens,rtdb,geom,basis,nat,nefc,
     &                dbl_mb(k_efcc),dbl_mb(k_efcs),dbl_mb(k_efcz),
     &                int_mb(k_efciat),dbl_mb(k_rad),ecos,cosmo_file)
c
c     ----- release memory block -----
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('cosmo_charges: ga_destroy failed g_dens',0,GA_ERR)
      enddo
      if(.not.ma_chop_stack(l_efcc))
     & call errquit('cosmo_charges, ma_chop_stack of l_efcc failed',911,
     &       ma_err)
c
      return
 9999 format(/,10X,15(1H-),
     1       /,10X,'-cosmo- charges',
     2       /,10X,15(1H-))
      end
C>
C> @}
c $Id$
