      subroutine hnd_elpmap(rtdb,basis,geom)
c
c $Id: hnd_elpmap.F,v 1.6 2005-04-13 23:39:38 marat Exp $
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
c     bq variables (MV)
      logical dobq
      integer bq_ncent
      integer i_cbq
      integer i_qbq
      double precision elpotbq
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
c     ----- calculate electrostatic potential -----
c
      if (ga_nodeid().eq.0) write(luout,9999)
      if (ga_nodeid().eq.0) write(luout,9994)
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
     &    call errquit('hnd_elpmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'epot pnt',l_epot,k_epot))
     &    call errquit('hnd_elpmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_elpmap: ma failed',911,MA_ERR)
c
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_elpmap: geom_cent_get',911,GEOM_ERR)
   30 continue
c
c     ----- calculate electronic contribution at all points -----
c
      call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(k_xyzpt),nat,
     &                dbl_mb(k_epot),0)

      dobq = .false.
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
      if(dobq) then
        write(luout,9992)
      else
        write(luout,9997)
      end if
      do 230  iat=1,nat    
         xp = dbl_mb(k_xyzpt  +3*(iat-1))
         yp = dbl_mb(k_xyzpt+1+3*(iat-1))
         zp = dbl_mb(k_xyzpt+2+3*(iat-1))
c
c     ----- add nuclear contribution -----
c
         elpotn = -dbl_mb(k_epot+iat-1)
         elpotbq = 0.0d0
         do 210 i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if (rr.lt.1.0d-3) go to 210
            elpotn = elpotn + zan/rr
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
c        end of external bq contributions
c        -------------------------------
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_elfmap: geom_cent_tag failed',0,GEOM_ERR)
c        geom_tag_to_element returns false for Bq elements (MV)
c        -----------------------------------------------------
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) then 
            if(symbol.ne."bq") call
     &      errquit('hnd_elfmap: geom_tag_to_element failed',0,GEOM_ERR)
         end if
c
         if(dobq) then
           write(luout,9991) iat,symbol,xp,yp,zp,elpotn,
     &                       dbl_mb(k_epot+iat-1),elpotbq
         else
           write(luout,9995) iat,symbol,xp,yp,zp,elpotn,
     &                       dbl_mb(k_epot+iat-1)
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
c
c     ------- Deallocate MA memory ------
c
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_epot)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_epot failed',911,MA_ERR)
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
