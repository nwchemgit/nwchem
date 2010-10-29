      subroutine argos_diana_pltgrd(lprint,isel)
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      logical lprint
      integer isel(nsa)
c
      ngrd(1)=ngrid
      ngrd(2)=ngrid
      ngrd(3)=ngrid
c
c
c     open the charge file
c
      if(me.eq.0) then
      open(unit=lfnchg,file=filchg(1:index(filchg,' ')-1),
     + form='formatted',status='old',err=9999)
c
      read(lfnchg,1000) nsa,nq,box
 1000 format(i7,i3,3f12.6)
      endif
      if(np.gt.1) then
      call ga_brdcst(mag_d01,nsa,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d01,nq,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d01,box,3*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
c     allocate memory for coordinates and charges
c
      if(.not.ma_push_get(mt_dbl,3*nsa,'x',l_x,i_x))
     + call md_abort('Failed to allocate memory for x',me)
      if(.not.ma_push_get(mt_dbl,nq*nsa,'q',l_q,i_q))
     + call md_abort('Failed to allocate memory for q',me)
      if(.not.ma_push_get(mt_byte,2*nsa,'t',l_t,i_t))
     + call md_abort('Failed to allocate memory for t',me)
c
c     read coordinates and charges
c
      call argos_diana_rdchg(dbl_mb(i_x),dbl_mb(i_q),byte_mb(i_t))
c
c     allocate memory for the grid
c
      ngx=ngrd(1)
      ngy=ngrd(2)
      ngz=ngrd(3)
      mgrd=ngx*ngy*ngz
c
      if(me.eq.0) then
      if(lprint) write(*,2000) ngrd,rcut,
     + filchg(1:index(filchg,' ')-1),
     + filplt(1:index(filplt,' ')-1)
 2000 format(' Electrostatic potential on ',i5,'x',i5,'x',i5,
     + ' grid extending',
     + f12.6,' nm',/,' Charges from ',a,' written to ',a)
      if(lprint.and.iesppb.gt.0) write(*,2001)
 2001 format(' Periodic boundary conditions')
c
      endif
c
      if(.not.ma_push_get(mt_dbl,mgrd,'grid',l_grid,i_grid)) then
      call md_abort('Failed to allocate memory for grid',me)
      endif
c
c     calculate the electrostatic potential on the grid
c
      call argos_diana_espgrd(dbl_mb(i_x),dbl_mb(i_q),byte_mb(i_t),
     + dbl_mb(i_grid),lprint,isel)
c
c     deallocate memory
c
      if(.not.ma_pop_stack(l_grid))
     + call md_abort('Failed to deallocate grid',0)
      if(.not.ma_pop_stack(l_t))
     + call md_abort('Failed to deallocate t',0)
      if(.not.ma_pop_stack(l_q))
     + call md_abort('Failed to deallocate q',0)
      if(.not.ma_pop_stack(l_x))
     + call md_abort('Failed to deallocate x',0)
c
      close(unit=lfnchg)
      return
 9999 continue
      call md_abort('argos_diana_pltgrd: error',0)
      return
      end
