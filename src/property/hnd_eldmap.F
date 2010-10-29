      subroutine hnd_eldmap(rtdb,basis,geom)
c
c $Id$
c
c     This routine calculates the electron density and spin
c     density for a given density at the atomic positions.
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
c
      integer rtdb          ! [Input] rtdb     
      integer basis         ! [Input] Basis set
      integer geom          ! [Input] Geometry
c
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, i, g_spin
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt, l_eden, k_eden
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      logical spindens
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Get density matrix
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     ----- calculate electron (and spin density) -----
c
      spindens = ndens .gt. 1
c
      if (ga_nodeid().eq.0) write(luout,9999)
      if (ga_nodeid().eq.0) write(luout,9996)
c
       call ecce_print_module_entry('Eldensities')
c
c     ----- define points for calculation -----
c           1. grid points    (not active)
c           2. nuclei
c           3. center of mass (not active)
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_eldmap: geom_ncent',911,GEOM_ERR)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_eldmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nat,'eden pnt',l_eden,k_eden))
     &    call errquit('hnd_eldmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_eldmap: ma failed',911,MA_ERR)
c
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_eldmap: geom_cent_get',911,GEOM_ERR)
   30 continue
c
c     ----- calculate electronic contribution at all points -----
c           do for both total dens and spin dens if ndens.gt.1
c
      if (spindens) then
         call hnd_elfcon(basis,geom,g_dens(3),dbl_mb(k_xyzpt),nat,
     &                   dbl_mb(k_eden),-1)
         g_spin=ga_create_atom_blocked(geom,basis,'hnd_spind:spindens')
         call ga_add(1.0d0,g_dens(1),-1.0d0,g_dens(2),g_spin)
         call hnd_elfcon(basis,geom,g_spin,dbl_mb(k_xyzpt),nat,
     &                   dbl_mb(k_eden+nat),-1)
         if (.not.ga_destroy(g_spin)) call errquit
     &      ('hnd_hnd_eldmap, ga_destroy of g_spin failed',911,GA_ERR)
      else
         call hnd_elfcon(basis,geom,g_dens,dbl_mb(k_xyzpt),nat,
     &                   dbl_mb(k_eden),-1)
      endif
c
c     ----- collect and output results of all points -----
c
      if (ga_nodeid().gt.0) goto 300
c
c     ----- total electron density now -----
c
      write(luout,9993)
      write(luout,9997)
      do 130  iat=1,nat
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_elfmap: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) call
     &      errquit('hnd_elfmap: geom_tag_to_element failed',0,GEOM_ERR)
         write(luout,9995) iat,symbol,dbl_mb(k_xyzpt+3*(iat-1)),
     &                     dbl_mb(k_xyzpt+1+3*(iat-1)),
     &                     dbl_mb(k_xyzpt+2+3*(iat-1)),
     &                     dbl_mb(k_eden+iat-1)
c
c        ----- store ecce data -----
c
         call ecce_print1_char('atom name',symbol,1)
         call ecce_print1('Total electron density',MT_DBL,
     &                    dbl_mb(k_eden+iat-1),1)
  130 continue
c
c     ----- spin density now -----
c
      if (spindens) then
         write(luout,9992)
         write(luout,9997)
         do 230  iat=1,nat
            if(.not. geom_cent_tag(geom,iat,at_tag)) call
     &         errquit('hnd_elfmap: geom_cent_tag failed',0,GEOM_ERR)
            if(.not.geom_tag_to_element(at_tag,symbol,element,atn)) call
     &         errquit('hnd_elfmap: geom_tag_to_element',0,GEOM_ERR)
            write(luout,9995) iat,symbol,dbl_mb(k_xyzpt+3*(iat-1)),
     &                        dbl_mb(k_xyzpt+1+3*(iat-1)),
     &                        dbl_mb(k_xyzpt+2+3*(iat-1)),
     &                        dbl_mb(k_eden+nat+iat-1)
c
c           ----- store ecce data -----
c
            call ecce_print1_char('atom name',symbol,1)
            call ecce_print1('Total spin-density',MT_DBL,
     &                       dbl_mb(k_eden+nat+iat-1),1)
  230    continue
      endif
c
      call ecce_print_module_exit('Eldensities','ok')
      call util_flush(luout)
c
c     ----- release memory block -----
c
  300 call ga_sync()
c
c     ------- Deallocate MA memory ------
c
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_eldmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_eden)) call errquit
     &   ('hnd_eldmap, ma_pop_stack of l_eden failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_eldmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('eldmap: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
 9999 format(/,10x,21(1h-),/,10x,'Electron/spin density',
     1       /,10x,21(1h-),/)
 9996 format(' 1 a.u. = 4.80286 10**(-10) esu ')
 9997 format(3x,'Point',6X,'X',9X,'Y',9X,'Z',6X,'Density (a.u.)')
 9995 format(i5,1x,a2,3f10.5,f15.6)
 9993 format(/,' Total electron density ',/,' ---------------------- ')
 9992 format(/,' Total spin density     ',/,' ------------------     ')
      end
