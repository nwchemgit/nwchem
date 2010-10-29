      subroutine hnd_spind(rtdb,basis,geom)
c
c $Id$
c
c     This routine calculates the Fermi Contact Term and
c     Spin-Dipolar terms for a given (spin) density at the
c     atomic positions
c
      implicit none 
#include "errquit.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
c
      integer rtdb      ! [Input] rtdb         
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
c
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, g_spin, i, j, ij, isonr
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt, l_efgs, k_efgs
      integer l_elds, k_elds
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      double precision xp, yp, zp, fac, fac1, fac2
      double precision pi, betae, betan, convf, con, aiso, spnden
      double precision ge, emf, vl, auev, evmhz, gmhz, gnu
      double precision spndxx, spndyy, spndzz, spndxy, spndxz, spndyz
      double precision a(6),vec(3,3),vec2(3,3),eig(3),hbar,gnuc
c
      logical atom_gfac
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,atom_gfac
c
      data gnuc   /5.05078343d-27/      ! Nuclear magneton
      data hbar   /1.05457168d-34/      ! Planck constant over 2 pi
      data ge     /2.002319304386d+00/  ! Electron g-factor
      data emf    /1836.152701d+00/     ! Proton-electron mass ratio
      data vl     /137.0359895d+00/     ! Speed of light in au
      data auev   /27.2113961d+00/      ! Conversion from au to eV
      data evmhz  /2.41798836d+08/      ! Conversion from eV to MHz
      data gmhz   /2.8025d+00/          ! Conversion from Gauss to MHz
c
c     Get density matrix (needs to be uhf or udft to be useful)
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     Return if ndens equals 1 (RHF)
c
      if (ndens.eq.1) then
          if (ga_nodeid().eq.0) write(luout,7000)
          goto 911
      endif
c
c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
      if (ga_nodeid().eq.0) write(luout,9999)
c
c     --- calculate constants and conversion terms ---
c
      pi   = acos(-1.0d0)
      fac  =(4.0d0*pi/3.0d0)
      betae=1.0d0/(2.0d0*vl)
      betan=betae/emf
      convf=auev*evmhz
      con  =ge*betae*betan*convf
c
c     ----- define points for calculation -----
c
c           1. nuclei
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_spind: geom_ncent',911,GEOM_ERR)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_spind: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat*2,'eld pnt',l_elds,k_elds)) !two time 
     &    call errquit('hnd_spind: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,6*nat,'efg pnt',l_efgs,k_efgs))
     &    call errquit('hnd_spind: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_spind: ma failed',911,MA_ERR)
c
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_eldmap: geom_cent_get',911,GEOM_ERR)
   30 continue
c
c     ----- get spin density matrix -----
c
      g_spin = ga_create_atom_blocked(geom, basis, 'hnd_spind:spindens')
      call ga_add(1.0d0,g_dens(1),-1.0d0,g_dens(2),g_spin)
c
c     ----- calculate spin density at the nuclei    -----
c           this is for the isotropic interaction
c           ( fermi contact term ) of the hyperfine
c           coupling tensor.
c
      call hnd_elfcon(basis,geom,g_spin,dbl_mb(k_xyzpt),nat,
     &                dbl_mb(k_elds),-1)
c
c     Let node 0 do the printing
c
      if (ga_nodeid().gt.0) goto 400
c
c     ----- print spin density and -aiso- at the nuclei -----
c
      write(luout,9994)
      write(luout,9997)
      call ecce_print_module_entry('fermi contact')
c
      do 370 iat=1,nat
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_efgmap: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) call
     &      errquit('hnd_efgmap: geom_tag_to_element failed',0,GEOM_ERR)
         if (.not. atom_gfac(dbl_mb(k_zanpt+iat-1),gnu,isonr)) call
     &       errquit('hnd_efgmap: atom_gfac failed',0, UERR)
         xp  = dbl_mb(k_xyzpt  +3*(iat-1))
         yp  = dbl_mb(k_xyzpt+1+3*(iat-1))
         zp  = dbl_mb(k_xyzpt+2+3*(iat-1))
         spnden=dbl_mb(k_elds+iat-1)
         aiso  =2.0d0*fac*con*gnu*spnden
         write(luout,9995) iat,isonr,symbol,xp,yp,zp,spnden,aiso,
     &                     aiso/gmhz
         call ecce_print1_char('atom name',symbol,1)
         call ecce_print1('fermi contact',MT_DBL,aiso,1)
  370 continue
c
      call ecce_print_module_exit('fermi contact','ok')
c
  400 call ga_sync()
c
c     ----- calculate spin-dipolar term at the nuclei -----
c           this is for the anisotropic interaction
c           of the hyperfine coupling tensor.
c
c     ----- calculate spin-dipolar term at all points -----
c
      call hnd_elfcon(basis,geom,g_spin,dbl_mb(k_xyzpt),nat,
     &                dbl_mb(k_efgs),2)
c
c     Let node 0 do the printing
c
      if (ga_nodeid().gt.0) goto 500
c
      write(luout,9996)
      call ecce_print_module_entry('spin-dipole')
c
      do 430 iat=1,nat
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_efgmap: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) call
     &      errquit('hnd_efgmap: geom_tag_to_element failed',0,GEOM_ERR)
         if (.not. atom_gfac(dbl_mb(k_zanpt+iat-1),gnu,isonr)) call
     &       errquit('hnd_efgmap: atom_gfac failed',0, UERR)
         xp  = dbl_mb(k_xyzpt  +3*(iat-1))
         yp  = dbl_mb(k_xyzpt+1+3*(iat-1))
         zp  = dbl_mb(k_xyzpt+2+3*(iat-1))
c
c        Do remember to take out the isotropic Fermi Contact Term
c        (which you get for free from the EFG)
c
         spndxx  = dbl_mb(k_efgs  +6*(iat-1)) + fac*dbl_mb(k_elds+iat-1)
         spndyy  = dbl_mb(k_efgs+1+6*(iat-1)) + fac*dbl_mb(k_elds+iat-1)
         spndzz  = dbl_mb(k_efgs+2+6*(iat-1)) + fac*dbl_mb(k_elds+iat-1)
         spndxy  = dbl_mb(k_efgs+3+6*(iat-1))
         spndxz  = dbl_mb(k_efgs+4+6*(iat-1))
         spndyz  = dbl_mb(k_efgs+5+6*(iat-1)) 
c
         a(1) = spndxx
         a(2) = spndxy
         a(3) = spndyy
         a(4) = spndxz
         a(5) = spndyz
         a(6) = spndzz
         ij=0
         do 434 i = 1, 3
         do 434 j = 1, i
            ij = ij + 1
            vec(i,j) = a(ij)
            vec(j,i) = a(ij)
            vec2(i,j) = a(ij)*con*gnu
            vec2(j,i) = a(ij)*con*gnu
  434    continue
         call ecce_print1_char('atom name',symbol,1)
         call ecce_print2('spin-dipole tensor',MT_DBL,vec2,3,3,3)
         call hnd_diag(vec,eig,3,.false.,.false.)
         fac1 = con*gnu
         fac2 = con*gnu/gmhz
         write(luout,9992)
         write(luout,9991) iat,isonr,symbol,xp,yp,zp
         write(luout,9984)
         write(luout,9983) spndxx,spndyy,spndzz,spndxy,spndxz,spndyz
         write(luout,9990)
         write(luout,9989) eig(1),eig(2),eig(3),eig(1)*fac1,eig(2)*fac1,
     1                  eig(3)*fac1,eig(1)*fac2,eig(2)*fac2,eig(3)*fac2
         write(luout,9988)
         write(luout,9986) ((vec(i,j),j=1,3),i=1,3)
         write(luout,*)
         write(luout,*)
         do i = 1, 3
            eig(i) = eig(i)*con*gnu
         enddo
         call ecce_print1('spin-dipole eigenvalues',MT_DBL,eig,3)
         call ecce_print2('spin-dipole eigenvectors',MT_DBL,vec,3,3,3)
  430 continue
c
      call ecce_print_module_exit('spin-dipole','ok')
      call util_flush(luout)
c
c     ----- release memory block -----
c
  500 call ga_sync()
c
c     ------- Deallocate MA and GA memory ------
c
      if (.not.ga_destroy(g_spin)) call errquit
     &   ('hnd_spind, ga_destroy of g_spin failed',911,GA_ERR)
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_spind, ma_pop_stack of l_dens failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_efgs)) call errquit
     &   ('hnd_spind, ma_pop_stack of l_dens failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_elds)) call errquit
     &   ('hnd_spind, ma_pop_stack of l_elds failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_spind, ma_pop_stack of l_dens failed',911,MA_ERR)
  911 do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call errquit
     &      ('hnd_spind, ga_destroy of g_dens failed',911,GA_ERR)
      enddo
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
 7000 format(/,10x,'Hyperfine Interactions not calculated for RHF',
     1      ' wave functions')
 9999 format(/,10x,22(1h-),/,10x,'Hyperfine Interactions',
     1       /,10x,22(1h-),/)
 9994 format(/,' Total Spin Density (Fermi Contact Term)',/,
     1         ' ------------------ --------------------')
 9997 format(5x,'Atom ',11x,'x',9x,'y',9x,'z',7x,
     1          'Density (a.u.)',6x,'Aiso(MHz)',5x,'Aiso(Gauss)')
 9995 format(1x,i5,1x,i3,'-',a2,2x,3f10.5,2x,f13.6,5x,f13.6,2x,f13.6)
 9996 format(/,/,10x,' -------------------------------------------',
     1       /,10x,' Anisotropic Interaction (Spin-Dipolar Term)',
     2       /,10x,' -------------------------------------------',/)
 9992 format(1x,60(1h-),/,5x,'Atom ',11x,'x',9x,'y',9x,'z',/,
     1       1x,60(1h-))
 9991 format(1x,i5,1x,i3,'-',a2,2x,3f10.5)
 9984 format(/,3x,'Spin-Dipolar term (a.u.)',/,8x,'xx',13x,'yy',
     1       13x,'zz',13x,'xy',13x,'xz',13x,'yz',/,1x,90(1h-))
 9983 format(6f15.6)
 9990 format(/,3x,'Principal Comp. of -hf- tensor (au) ',3x,
     1          8x,'Anisotropic -hfcc- (MHz)',3x,4x,
     2          8x,'Anisotropic -hfcc- (Gauss)',/,1x,122(1h-))
 9989 format(1x,3f13.6,2x,3f13.6,2x,3f13.6)
 9988 format(/,3x,'Orientation of the principal axis',
     1         ' of Hyperfine Tensor w.r.t absolute (molecular) frame',
     2         /,1x,122(1h-))
 9986 format(1x,3f13.6)
      end
