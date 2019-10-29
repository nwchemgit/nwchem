      subroutine hnd_elpmap(rtdb,basis,geom)
c
c $Id$
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
#include "msgids.fh"
c
      integer rtdb      ! [Input] rtdb        
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
c
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, i
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt, l_epot, k_epot
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
      integer l_tepot,k_tepot
c
      character*30 theory
      integer nbf
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical ao_1prdm_read
      external ao_1prdm_read
      logical ocube
      integer nsp(3)
      integer un_grid,avail_ma
      integer npasses,i_pass
      integer iptr,nprp_pass,ngrid(3)
      character*255 cube_name
      character*255 dmat_file
c
c
c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Get density matrix
c
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     + call errquit('task: no task input for theory?',0, RTDB_ERR)
      if (.not. rtdb_cget(rtdb,'scf:scftype',1,scftyp)) scftyp = "RHF"
      ndens = 1
      if (scftyp.eq.'UHF') ndens = 3
c
c     if (ga_nodeid().eq.0) write(luout,*) "scftyp:",scftyp
c     if (ga_nodeid().eq.0) write(luout,*) "theory:",theory
c
c     Read density matrix from a file
c
      if (theory.eq.'tce'.or.theory.eq.'tddft') then
c
        do i = 1, ndens
         g_dens(i) = ga_create_atom_blocked(geom,basis,'density matrix')
         call ga_zero(g_dens(i))
        enddo
c
        if (.not. bas_numbf(basis,nbf))
     &    call errquit('hnd_elfmap: could not get nbf',0, BASIS_ERR)
c
        call util_file_name('dmat',.false.,.false.,dmat_file) ! get filename
        if (ga_nodeid().eq.0)
     &        write(luout,*) "Reading file from:",dmat_file
        if(.not.ao_1prdm_read(nbf,g_dens(ndens),dmat_file))
     &      call errquit('hnd_elfmap: ao_1prdm_read failed',0,0)
c
      else  ! generate from movecs for scf or dft
         call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
      endif
c
      call ecce_print_module_entry('Elpoten')
c
c     ----- define points for calculation -----
c           1. grid points    (not active)
c           2. nuclei
c           3. center of mass (not active)
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elpmap: geom_ncent',911,GEOM_ERR)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_elpmap: ma failed',1,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_elpmap: ma failed',2,MA_ERR)
c
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_elpmap: geom_cent_get',911,GEOM_ERR)
   30 continue
c
      if(.not.rtdb_get(rtdb,'prop:cubefile',mt_log,1,ocube)) 
     +     ocube = .false.

      if(ocube) then
        if (.not. rtdb_cget(rtdb, "prop:grid:output",1,cube_name))
     >     call util_file_prefix("elp.cube",cube_name)
      end if
cc     
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
     &    call errquit('hnd_elpmap: prop:xyz',911,MA_ERR)
        if (.not. rtdb_get(rtdb,'prop:xyz',mt_dbl,
     >                      3*nprp,dbl_mb(i_prp_c)))
     &    call errquit('hnd_elpmap: prop:xyz failed',911,RTDB_ERR)
        do_points = .true.
        do_output = .true.
      else
        nprp = nat
        if (.not. ma_push_get(mt_dbl,3*nat,'prop:xyz',h_prp_c,i_prp_c))
     &    call errquit('hnd_elpmap: ma failed',3,MA_ERR)
        call dcopy(3*nat,dbl_mb(k_xyzpt),1,dbl_mb(i_prp_c),1)
        do_points = .false.
        do_output = .true.
      end if
c
#ifndef OLDCUBE
      if(ocube)
     c     call prop_grid_writecubehead(geom,un_grid,cube_name)
#endif
c
c     check if enough mem is available
c
      avail_ma = ma_inquire_avail(mt_dbl)*0.9d0
      npasses=nprp*2d0/avail_ma + 1d0
c
      call ga_igop(msg_efgs_col-4,npasses,1,'max')
#ifdef OLDCUBE
      if(npasses.gt.1) call errquit('hnd_elpmap: ma failed',0,MA_ERR)
#endif
      if(npasses.gt.1.and.do_points) 
     C     call errquit('hnd_elpmap: ma failed',0,MA_ERR)
      nprp_pass=nprp/npasses
      if(npasses.gt.1) then
c
c     enforce nprp_pass to be a multiple of ngrid(3) to keep cube happy
c
         nprp_pass=(nprp_pass/ngrid(3))*ngrid(3)
         npasses=nprp/nprp_pass
         if(ga_nodeid().eq.0) then
            write(6,'(a,i5,a,i10,a,i19)') 
     H           ' hnd_elpmap: npasses:',npasses,
     W           ' nprp_pass:',nprp_pass,' nprp:',nprp
         endif
      endif
      iptr=0
      do i_pass=1,npasses
         iptr=(i_pass-1)*nprp_pass
         if(i_pass.eq.npasses) nprp_pass=nprp-iptr

      if (.not. ma_push_get(mt_dbl,nprp_pass,'epot pnt',l_epot,k_epot))
     &    call errquit('hnd_elpmap: ma failed',4,MA_ERR)
c
c     total electric field array (M.V.)
c     --------------------------------
      if (.not. ma_push_get(mt_dbl,nprp_pass,'totepot',l_tepot,k_tepot))
     &    call errquit('hnd_elpmap: ma failed',5,MA_ERR)
c
cc
c     ----- calculate electronic contribution at all points -----
c
      call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(i_prp_c+iptr*3),
     N     nprp_pass,   dbl_mb(k_epot),0)
c
      dobq = .false.
      bq_ncent = 0
      i_cbq = 0
      i_qbq = 0
      if(geom_extbq_on()) then
        dobq = .true.
        bq_ncent = geom_extbq_ncenter()
        i_cbq = geom_extbq_coord()
        i_qbq = geom_extbq_charge()
      end if

c
c     ----- collect and output results of all points -----
c
      if (ga_nodeid().gt.0) goto 300
c
c     ----- calculate electrostatic potential -----
c
      if(do_output) then
        if (ga_nodeid().eq.0) write(luout,9999)
        if (ga_nodeid().eq.0) write(luout,9994)
c
c
        if(dobq) then
          write(luout,9992)
        else
          write(luout,9997)
        end if
      end if
      do 230  iat=1,nprp_pass
         xp = dbl_mb(iptr*3+i_prp_c  +3*(iat-1))
         yp = dbl_mb(iptr*3+i_prp_c+1+3*(iat-1))
         zp = dbl_mb(iptr*3+i_prp_c+2+3*(iat-1))
c
c     ----- add nuclear contribution -----
c
         elpotn = -dbl_mb(k_epot+iat-1)
         elpotbq = 0.0d0
         enuc = 0.0
         do 210 i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if (rr.lt.1.0d-3) go to 210
            elpotn = elpotn + zan/rr
            enuc = enuc + zan/rr
  210    continue
c
c     ----- form -efc- contribution -----
c           from cosmo point charges !!!!
c
         if (cosmo_last) then
            if (.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc))
     &         call errquit('hnd_elpmap: rtdb get failed for nefc ',911,
     &         RTDB_ERR)
            if (.not.ma_push_get(mt_dbl,nefc*3,'efcc',l_efcc,k_efcc))
     &         call errquit('hnd_elpmap: malloc k_efcc fail',911,ma_err)
            if (.not.ma_push_get(mt_dbl,nefc,'efcz',l_efcz,k_efcz))
     &         call errquit('hnd_elpmap: malloc k_efcz fail',911,ma_err)
            if (.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,
     &         dbl_mb(k_efcc))) call
     &         errquit('hnd_elpmap: rtdb get failed efcc',912,rtdb_err)
            if (.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,nefc,
     &         dbl_mb(k_efcz))) call
     &         errquit('hnd_elpmap: rtdb get failed efcz',913,rtdb_err)
            do 220 i = 1,nefc
               xn = dbl_mb(k_efcc+3*(i-1)  ) - xp
               yn = dbl_mb(k_efcc+3*(i-1)+1) - yp
               zn = dbl_mb(k_efcc+3*(i-1)+2) - zp
               rr =  sqrt(xn*xn + yn*yn + zn*zn)
               if (rr.lt.1.0d-3) then
                  write(luout,9993) xp,yp,zp,i
                  go to 220
               endif
               elpotn = elpotn + dbl_mb(k_efcz+i-1)/rr
  220       continue
            if (.not.ma_chop_stack(l_efcc)) call
     &         errquit('hnd_elpmap: chop stack l_efcc',913,ma_err)
         endif
c        adding external bq contributions(MV)
c        ----------------------------------
         if (dobq) then
            do i = 1,bq_ncent
               xn = dbl_mb(i_cbq+3*(i-1)  ) - xp
               yn = dbl_mb(i_cbq+3*(i-1)+1) - yp
               zn = dbl_mb(i_cbq+3*(i-1)+2) - zp
               rr =  sqrt(xn*xn + yn*yn + zn*zn)
               elpotbq = elpotbq+dbl_mb(i_qbq+i-1)/rr
            end do
         end if
         elpotn = elpotn + elpotbq
         dbl_mb(k_tepot+iat-1) = elpotn
c        end of external bq contributions
c        -------------------------------
         if(do_output) then
         if(.not.do_points) then
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_elpmap: geom_cent_tag failed',0,GEOM_ERR)
c        geom_tag_to_element returns false for Bq elements (MV)
c        -----------------------------------------------------
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) then 
            if(symbol.ne."bq") call
     &      errquit('hnd_elpmap: geom_tag_to_element failed',0,GEOM_ERR)
         end if
         end if
c
         if(do_points) then
           symbol = "pt"
           if(dobq) then
             write(luout,9991) iat,symbol,xp,yp,zp,elpotn,
     &                         dbl_mb(k_epot+iat-1),elpotbq
           else
             write(luout,9995) iat,symbol,xp,yp,zp,enuc,
     &                         dbl_mb(k_epot+iat-1)
           end if
         else
           if(dobq) then
             write(luout,9991) iat,symbol,xp,yp,zp,elpotn,
     &                         dbl_mb(k_epot+iat-1),elpotbq
           else
             write(luout,9995) iat,symbol,xp,yp,zp,elpotn,
     &                         dbl_mb(k_epot+iat-1)
           end if
         end if
         end if
c
c        ----- store ecce data -----
c
         call ecce_print1_char('atom name',symbol,1)
         call ecce_print1('Electrostatic pot',MT_DBL,elpotn,1)
         call ecce_print1('Diamagnetic shield',MT_DBL,
     &                    dbl_mb(k_epot+iat-1),1)
c
  230 continue ! Assembling and printing next atom
c
      call ecce_print_module_exit('Elpoten','ok')
      call util_flush(luout)
c
c     ----- release memory block -----
c
  300 call ga_sync()
      
      if(ga_nodeid().eq.0.and.ocube) then
        call util_print_centered(luout, 
     >    "writing esp potential to "//cube_name,.true.,.true.)
#ifdef OLDCUBE
        call prop_grid_write_cube(geom,nprp,dbl_mb(k_tepot),cube_name)
#else
c        call prop_grid_writecubehead(geom,un_grid,cube_name)
        call prop_grid_writecubegrid(nprp_pass,dbl_mb(k_tepot),
     F       un_grid)
#endif
      end if
c
c     if custom grid is requested save final electric 
c     into rtdb(M.V.)
c     -----------------------------------------------
      if(do_points) then
        if (.not. rtdb_put(rtdb,'prop:epot_xyz',mt_dbl,
     >                        nprp,dbl_mb(k_tepot)))
     &      call errquit('hnd_elpmap: epot_xyz failed',0,RTDB_ERR)
      end if
cc
c     ------- Deallocate MA memory ------
c
      if (.not.ma_pop_stack(l_tepot)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_tepot failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_epot)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_epot failed',911,MA_ERR)
      enddo ! i_pass
      if(.not.ocube)  then
        if (.not.ma_pop_stack(h_prp_c)) call errquit
     &      ('hnd_elpmap, ma_pop_stack of h_prp_c failed',911,MA_ERR)
      else
         close(un_grid)
         call prop_grid_destroy()
      endif
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('elpmap: ga_destroy failed g_dens',0,GA_ERR)
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
