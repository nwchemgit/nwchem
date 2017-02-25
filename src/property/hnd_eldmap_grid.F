      subroutine hnd_eldmap_grid(rtdb,basis,geom)
c
c $Id: hnd_elpmap.F 27042 2015-04-28 21:14:03Z edo $
c
c     This routine calculates the electrostatic potential
c     for a given density at the atomic positions.
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cosmo.fh"
#include "bas.fh"
c
      integer rtdb      ! [Input] rtdb        
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
c
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, i, g_spin
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt
      integer l_eden, k_eden
      integer nefc, l_efcc, k_efcc, l_efcz, k_efcz
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision elpotn
      double precision rr
      double precision enuc
c     bq variables (MV)
      logical dobq
      integer bq_ncent
      integer i_cbq
      integer i_qbq
      double precision elpotbq
c     property grid points variables (MV)
      integer h_prp_c,i_prp_c
      integer ma_prp_type          
      integer nprp
      character*26 prp_date
      logical do_points
      logical do_output
c
      character*30 theory
      integer nbf
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical ao_1prdm_read
      external ao_1prdm_read
      logical ocube, spindens
      integer nsp(3)
      integer iptr,nprp_pass,ngrid(3)
      character*255 cube_name
      character*255 cube_name_total
      character*255 cube_name_spin
      character*255 cube_name_alpha
      character*255 cube_name_beta
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
      spindens = ndens .gt. 1
c
c     ----- define points for calculation -----
c           1. grid points    (not active)
c           2. nuclei
c           3. center of mass (not active)
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elpmap_gr: geom_ncent',911,GEOM_ERR)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_elpmap_gr: ma failed',1,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_elpmap_gr: ma failed',2,MA_ERR)
c
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_elpmap_gr: geom_cent_get',911,GEOM_ERR)
   30 continue
c
      if(.not.rtdb_get(rtdb,'prop:cubefile',mt_log,1,ocube)) 
     +     ocube = .false.
c     
c     define points for the calculation now 
c     either custom grid or (default) nuclei positions (M.V.)
c     -------------------------------------------------
      if(ocube) then
         call prop_grid_initialize(rtdb,nat,dbl_mb(k_xyzpt))
         call prop_grid_get_r_ptr(nprp,i_prp_c,ngrid)
         do_points = .false.
         do_output = .false.
      else if(rtdb_get_info(rtdb, "prop:xyz", ma_prp_type, 
     >                 nprp, prp_date)) then
        nprp = nprp/3
        if (.not. ma_push_get(mt_dbl,3*nprp,'prop:xyz',h_prp_c,i_prp_c))
     &    call errquit('hnd_elfmap_gr: prop:xyz',911,MA_ERR)
        if (.not. rtdb_get(rtdb,'prop:xyz',mt_dbl,
     >                      3*nprp,dbl_mb(i_prp_c)))
     &    call errquit('hnd_elfmap_gr: prop:xyz failed',911,RTDB_ERR)
        do_points = .true.
        do_output = .true.
      else
        nprp = nat
        if (.not. ma_push_get(mt_dbl,3*nat,'prop:xyz',h_prp_c,i_prp_c))
     &    call errquit('hnd_eldmap_gr: ma failed',3,MA_ERR)
        call dcopy(3*nat,dbl_mb(k_xyzpt),1,dbl_mb(i_prp_c),1)
        do_points = .false.
        do_output = .true.
      end if
c
c     density array
      if (.not. ma_push_get(mt_dbl,4*nprp,'eden',l_eden,k_eden))
     &    call errquit('hnd_elpmap_gr: eden ma failed',4,MA_ERR)
c
      if (spindens) then
         call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(i_prp_c),nprp,
     &                dbl_mb(k_eden),-1)
c
         g_spin=ga_create_atom_blocked(geom,basis,'hnd_spind:spindens')

c        spin density  
         call ga_zero(g_spin)
         call ga_add(1.0d0,g_dens(1),-1.0d0,g_dens(2),g_spin)
         call hnd_elfcon(basis,geom,g_spin,dbl_mb(i_prp_c),nprp,
     &                dbl_mb(k_eden+nprp),-1)
c
c        alpha density  
         call ga_zero(g_spin)
         call ga_add(1.0d0,g_dens(1),0.0d0,g_spin,g_spin)
         call hnd_elfcon(basis,geom,g_spin,dbl_mb(i_prp_c),nprp,
     &                dbl_mb(k_eden+2*nprp),-1)
c
c        beta density  
         call ga_zero(g_spin)
         call ga_add(1.0d0,g_dens(2),0.0d0,g_spin,g_spin)
         call hnd_elfcon(basis,geom,g_spin,dbl_mb(i_prp_c),nprp,
     &                dbl_mb(k_eden+3*nprp),-1)
c
         if (.not.ga_destroy(g_spin)) call errquit
     &      ('hnd_hnd_eldmap, ga_destroy of g_spin failed',911,GA_ERR)
      else
         call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(i_prp_c),nprp,
     &                dbl_mb(k_eden),-1)
      endif
c
c     ----- collect and output results of all points -----
c
      if(ocube) then
       call util_file_prefix("eld.total.cube",cube_name_total)
       call util_file_prefix("eld.spin.cube",cube_name_spin)
       call util_file_prefix("eld.alpha.cube",cube_name_alpha)
       call util_file_prefix("eld.beta.cube",cube_name_beta)
      end if
c
      if (ga_nodeid().eq.0.and.ocube) then
         call prop_grid_write_cube(geom,nprp,dbl_mb(k_eden),
     &        cube_name_total)
         call prop_grid_write_cube(geom,nprp,dbl_mb(k_eden+nprp),
     &        cube_name_spin)
         call prop_grid_write_cube(geom,nprp,dbl_mb(k_eden+2*nprp),
     &        cube_name_alpha)
         call prop_grid_write_cube(geom,nprp,dbl_mb(k_eden+3*nprp),
     &        cube_name_beta)
      end if 
c
c     ----- release memory block -----
c
  300 call ga_sync()
c
c     ------- Deallocate MA memory ------
c
      if (.not.ma_pop_stack(l_eden)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_eden failed',911,MA_ERR)
      if(.not.ocube)  then
        if (.not.ma_pop_stack(h_prp_c)) call errquit
     &      ('hnd_elpmap, ma_pop_stack of h_prp_c failed',911,MA_ERR)
      else
         call prop_grid_destroy()
      endif
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('elpmap_gr: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
 9999 format(/,10x,45(1H-),
     1       /,10x,'Electrostatic potential/diamagnetic shielding',
     2       /,10x,45(1H-),/)
 9998 format(' Not enough core in -elpmap-')
 9997 format(3x,'Point',6x,'X',9x,'Y',9x,'Z',5x,'Total Potential(a.u.)',
     1       3x,'Diamagnetic shielding(a.u.)')
 9996 format(' --- Warning - electrostatic potential at ',
     1 3f10.5,' . contribution from nucleus ',i3,' ignored')
 9995 format(i5,1x,a2,3F10.5,f15.6,6x,f15.6)
 9994 format(' 1 a.u. = 9.07618 esu/cm ( or statvolts ) ')
 9993 format(' --- Warning - electrostatic potential at ',
     1 3f10.5,' . contribution from  -efc-  ',i3,' ignored')
c   
 9992 format(3x,'Point',6x,'X',9x,'Y',9x,'Z',5x,'Total Potential(a.u.)',
     1       3x,'Diamagnetic shielding(a.u.)',
     2       3x,'External Bq potential')
 9991 format(i5,1x,a2,3F10.5,f15.6,6x,f15.6,12x,f15.6)
 
      END
