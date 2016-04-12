*
* $Id$
*

*     *******************************************
*     *						*
*     *	 	   expand_cell			*
*     *						*
*     *******************************************

      logical function expand_cell(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"
#include "sym.fh"
#include "errquit.fh"

      logical value,default,setamatrix,getsym
      integer geom
      integer geom1,nion1
      integer geom2,nion2
      integer cell_expand(3),i,j,k,l,ind,ngrid(3)
      integer rt1(2),vt1(2),qt1(2),mt1(2),tt1(2)
      integer rt2(2),vt2(2),qt2(2),mt2(2),tt2(2)
      real*8  amatrix(3,3),amatrix2(3,3)

      character*50 cell_name,filename
      character*50 rtdb_unita,rtdb_ngrid

*     **** external functions ****
      logical  geom_amatrix_set,wvfnc_expand_cell,geom_P1_set
      external geom_amatrix_set,wvfnc_expand_cell,geom_P1_set


      value = btdb_get(rtdb,'nwpw:cell_expand',mt_int,3,cell_expand)
      if (.not.value) then
         value = .true.
         goto 201
      end if

*     **** generate new unit cell ****
      setamatrix = .false.
      if (.not.btdb_cget(rtdb,'cgsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
      if (btdb_get(rtdb,rtdb_unita,mt_dbl,9,amatrix)) then
        amatrix2(1,1) = amatrix(1,1)*cell_expand(1)
        amatrix2(2,1) = amatrix(2,1)*cell_expand(1)
        amatrix2(3,1) = amatrix(3,1)*cell_expand(1)
        amatrix2(1,2) = amatrix(1,2)*cell_expand(2)
        amatrix2(2,2) = amatrix(2,2)*cell_expand(2)
        amatrix2(3,2) = amatrix(3,2)*cell_expand(2)
        amatrix2(1,3) = amatrix(1,3)*cell_expand(3)
        amatrix2(2,3) = amatrix(2,3)*cell_expand(3)
        amatrix2(3,3) = amatrix(3,3)*cell_expand(3)
        value = value.and.btdb_put(rtdb,rtdb_unita,mt_dbl,9,amatrix2)

      else
        call dcopy(9,0.0d0,0,amatrix,1)
        

         value = value.and.geom_create(geom,'geometry')
         value = value.and.geom_rtdb_load(rtdb,geom,'geometry')
         value = value.and.geom_amatrix_get(geom,amatrix)

        default =  (amatrix(1,1) .eq. 1.0d0).and.
     >             (amatrix(2,1) .eq. 0.0d0).and.
     >             (amatrix(3,1) .eq. 0.0d0).and.
     >             (amatrix(1,2) .eq. 0.0d0).and.
     >             (amatrix(2,2) .eq. 1.0d0).and.
     >             (amatrix(3,2) .eq. 0.0d0).and.
     >             (amatrix(1,3) .eq. 0.0d0).and.
     >             (amatrix(2,3) .eq. 0.0d0).and.
     >             (amatrix(3,3) .eq. 1.0d0)
        if (default) then
          amatrix(1,1) = 20.0d0
          amatrix(2,1) =  0.0d0
          amatrix(3,1) =  0.0d0
          amatrix(1,2) =  0.0d0
          amatrix(2,2) = 20.0d0
          amatrix(3,2) =  0.0d0
          amatrix(1,3) =  0.0d0
          amatrix(2,3) =  0.0d0
          amatrix(3,3) = 20.0d0
        end if
        amatrix2(1,1) = amatrix(1,1)*cell_expand(1)
        amatrix2(2,1) = amatrix(2,1)*cell_expand(1)
        amatrix2(3,1) = amatrix(3,1)*cell_expand(1)
        amatrix2(1,2) = amatrix(1,2)*cell_expand(2)
        amatrix2(2,2) = amatrix(2,2)*cell_expand(2)
        amatrix2(3,2) = amatrix(3,2)*cell_expand(2)
        amatrix2(1,3) = amatrix(1,3)*cell_expand(3)
        amatrix2(2,3) = amatrix(2,3)*cell_expand(3)
        amatrix2(3,3) = amatrix(3,3)*cell_expand(3)


        if (default) then
          value = value.and.btdb_put(rtdb,rtdb_unita,mt_dbl,9,amatrix2)
        else
          value = value.and.geom_amatrix_set(geom,amatrix2)
          value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
        end if
        value = value.and.geom_destroy(geom)

      end if


c*     **** generate new grid ***
c      l = index(cell_name,' ') - 1
c      rtdb_ngrid = cell_name(1:l)//':ngrid'
c
c      if (.not.btdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)) then
c        ngrid(1) = 32
c        ngrid(2) = 32
c        ngrid(3) = 32
c      end if
c      ngrid(1) = ngrid(1)*cell_expand(1)
c      ngrid(2) = ngrid(2)*cell_expand(2)
c      ngrid(3) = ngrid(3)*cell_expand(3)
ccc      value = value.and.btdb_put(rtdb,rtdb_ngrid,mt_int,3,ngrid)



*     **** generate new geometry ****
      value = value.and.geom_create(geom1,'geometry')
      value = value.and.geom_rtdb_load(rtdb,geom1,'geometry')
      value = value.and.geom_ncent(geom1,nion1)
      if (.not. value) call errquit('opening geometry',0, GEOM_ERR)

*     *** set allocate memory ****
      nion2 = nion1*cell_expand(1)*cell_expand(2)*cell_expand(3)


      value = value.and.
     >        BA_push_get(mt_dbl, (3*nion1), 'rt1',rt1(2),rt1(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (3*nion2), 'rt2',rt2(2),rt2(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (3*nion1), 'vt1',vt1(2),vt1(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (3*nion2), 'vt2',vt2(2),vt2(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (nion1),   'qt1',qt1(2),qt1(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (nion2),   'qt2',qt2(2),qt2(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (nion1),   'mt1',mt1(2),mt1(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (nion2),   'mt2',mt2(2),mt2(1))
      value = value.and.
     >        BA_push_get(mt_byte,(16*nion1),'tt1',tt1(2),tt1(1))
      value = value.and.
     >        BA_push_get(mt_byte,(16*nion2),'tt2',tt2(2),tt2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     *** get geometry from rtdb ****
      value = value.and.
     >        geom_cart_get(geom1,nion1,byte_mb(tt1(1)),
     >                                  dbl_mb(rt1(1)),
     >                                  dbl_mb(qt1(1)))
      value = value.and.geom_masses_get(geom1,nion1,dbl_mb(mt1(1)))
      value = value.and.geom_vel_get(geom1,dbl_mb(vt1(1)))

      call ion_expand_cell(cell_expand,amatrix,
     >            nion1,
     >            dbl_mb(rt1(1)),dbl_mb(vt1(1)),
     >            dbl_mb(mt1(1)),dbl_mb(qt1(1)),byte_mb(tt1(1)),
     >            dbl_mb(rt2(1)),dbl_mb(vt2(1)),
     >            dbl_mb(mt2(1)),dbl_mb(qt2(1)),byte_mb(tt2(1)))


      value = value.and.geom_rtdb_delete(rtdb,'geometry')
      value = value.and.
     >        geom_cart_set(geom1,nion2,byte_mb(tt2(1)),
     >                                  dbl_mb(rt2(1)),
     >                                  dbl_mb(qt2(1)))
      value = value.and.geom_masses_set(geom1,nion2,dbl_mb(mt2(1)))
      value = value.and.geom_vel_set(geom1,dbl_mb(vt2(1)))

*     **** set unitcell to have P1 symmetry - may want to use more intricate logic in the future ***
      value = value.and.geom_P1_set(geom1)

      value = value.and.geom_rtdb_store(rtdb,geom1,'geometry')
      value = value.and.geom_destroy(geom1)

      value = value.and.BA_pop_stack(tt2(2))
      value = value.and.BA_pop_stack(tt1(2))
      value = value.and.BA_pop_stack(mt2(2))
      value = value.and.BA_pop_stack(mt1(2))
      value = value.and.BA_pop_stack(qt2(2))
      value = value.and.BA_pop_stack(qt1(2))
      value = value.and.BA_pop_stack(vt2(2))
      value = value.and.BA_pop_stack(vt1(2))
      value = value.and.BA_pop_stack(rt2(2))
      value = value.and.BA_pop_stack(rt1(2))

*     **** IO ****
      write(*,110) 1,1,1,cell_expand
      write(*,111)  amatrix(1,1), amatrix(2,1), amatrix(3,1),
     >             amatrix2(1,1),amatrix2(2,1),amatrix2(3,1)
      write(*,112)  amatrix(1,2), amatrix(2,2), amatrix(3,2),
     >             amatrix2(1,2),amatrix2(2,2),amatrix2(3,2)
      write(*,113)  amatrix(1,3), amatrix(2,3), amatrix(3,3),
     >             amatrix2(1,3),amatrix2(2,3),amatrix2(3,3)
      write(*,114) nion1,nion2


*     **** generate new wavefunctions ****
      if (.not.btdb_cget(rtdb,'pspw:input vectors',1,filename)) then
         call util_file_name('movecs', .false.,.false.,filename)
      endif
      ind = index(filename,' ') - 1
      if (.not. btdb_cput(rtdb,'xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

      if (.not. btdb_cput(rtdb,'xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)
      value = value.and.wvfnc_expand_cell(rtdb)



  201 continue
      expand_cell = value

      return 
  110 format(' cell conversion : ',I3,'x',I3,'x',I3,' --> ',
     >                     I3,'x',I3,'x',I3)
  111 format(' lattice vectors : a1=<',3f8.3,'>  --> a1new=<',3f8.3,'>')
  112 format('                   a2=<',3f8.3,'>  --> a2new=<',3f8.3,'>')
  113 format('                   a3=<',3f8.3,'>  --> a3new=<',3f8.3,'>')
  114 format(' number of ions  : ',I4,' --> ',I4)
      end
 



*     ***************************************************
*     *							*
*     *	 	   ion_expand_cell	    		*
*     *							*
*     ***************************************************

      subroutine ion_expand_cell(cell_expand,amatrix,
     >                           nion,
     >            coords,    voords,    mass,    charge,    name,
     >            coords_out,voords_out,mass_out,charge_out,name_out)
      implicit none
      integer cell_expand(3)
      real*8  amatrix(3,3)
      integer nion
      real*8 coords(3,*),voords(3,*),mass(*),charge(*)
      character*16 name(*)
      real*8 coords_out(3,*),voords_out(3,*),mass_out(*),charge_out(*)
      character*16 name_out(*)

*     **** local variables ****
      integer indx,i,j,k,ii

      indx = 0
      do k=0,cell_expand(3)-1
      do j=0,cell_expand(2)-1
      do i=0,cell_expand(1)-1
        do ii=1,nion
           coords_out(1,ii+indx) = coords(1,ii)
     >                           + i*amatrix(1,1)
     >                           + j*amatrix(1,2)
     >                           + k*amatrix(1,3)
     >                           - 0.5d0*amatrix(1,1)*(cell_expand(1)-1)
     >                           - 0.5d0*amatrix(1,2)*(cell_expand(2)-1)
     >                           - 0.5d0*amatrix(1,3)*(cell_expand(3)-1)
           coords_out(2,ii+indx) = coords(2,ii)
     >                           + i*amatrix(2,1)
     >                           + j*amatrix(2,2)
     >                           + k*amatrix(2,3)
     >                           - 0.5d0*amatrix(2,1)*(cell_expand(1)-1)
     >                           - 0.5d0*amatrix(2,2)*(cell_expand(2)-1)
     >                           - 0.5d0*amatrix(2,3)*(cell_expand(3)-1)
           coords_out(3,ii+indx) = coords(3,ii)
     >                           + i*amatrix(3,1)
     >                           + j*amatrix(3,2)
     >                           + k*amatrix(3,3)
     >                           - 0.5d0*amatrix(3,1)*(cell_expand(1)-1)
     >                           - 0.5d0*amatrix(3,2)*(cell_expand(2)-1)
     >                           - 0.5d0*amatrix(3,3)*(cell_expand(3)-1)
           voords_out(1,ii+indx) = voords(1,ii)
           voords_out(2,ii+indx) = voords(2,ii)
           voords_out(3,ii+indx) = voords(3,ii)
           mass_out(ii+indx)     = mass(ii)
           charge_out(ii+indx)   = charge(ii)
           name_out(ii+indx)     = name(ii)
        end do
        indx = indx + nion
      end do
      end do
      end do
      return 
      end
 

