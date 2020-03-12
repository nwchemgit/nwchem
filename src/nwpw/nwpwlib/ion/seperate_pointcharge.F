*     
* $Id$
*     

*     ***************************
*     *                         *
*     *   seperate_pointcharge  *
*     *                         *
*     ***************************
      subroutine seperate_pointcharge(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"

      logical mmexist
      common / ion_charge / mmexist

*     *** local variables ***
      integer     i,geom1,geom2,geom3,nion1,nion2,nion3
      integer rt(2),tt(2),qt(2),mt(2)
      double precision q,rxyz(3)
      character*16     t
      logical value

*     **** external functions ****
      logical          parsepointcharge
      external         parsepointcharge


*     *******************************************************************
*     **** seperate ions and charges from molecular pseudopotentials ****
*     *******************************************************************

*     **** generate pspwgeometry from chargepspwgeometry ****
*     **** generate chargegeometry from chargepspwgeometry ****
      value = beom_create(geom1,'chargepspwgeometry')
      value = value.and.beom_create(geom2,'pspwgeometry')
      value = value.and.beom_create(geom3,'chargegeometry')
      value = value.and.beom_rtdb_load(rtdb,geom1,'geometry')
      value = value.and.beom_rtdb_load(rtdb,geom2,'geometry')
      value = value.and.beom_ncent(geom1,nion1)

      !value = value.and.beom_rtdb_store(rtdb,geom1,'pspwgeometry')
      !value = value.and.beom_load(geom2,'pspwgeometry')


      if (.not. value) call errquit('opening geometry',0,
     &       GEOM_ERR)

      value = BA_push_get(mt_dbl, (3*nion1), 'rt',rt(2),rt(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (nion1),   'qt',qt(2),qt(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (nion1),   'mt',mt(2),mt(1))
      value = value.and.
     >        BA_push_get(mt_byte,(16*nion1),'tt',tt(2),tt(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      nion2 = 0
      nion3 = 0
      do i=1,nion1
         value = beom_cent_get(geom1,i,t,rxyz,q)
         if (.not.parsepointcharge(t)) then
            nion2 = nion2 + 1
         else
            nion3 = nion3 + 1
         end if
      end do
      value = value.and.
     >        beom_cart_get(geom1,nion1,byte_mb(tt(1)),
     >                                  dbl_mb(rt(1)),
     >                                  dbl_mb(qt(1)))
      value = value.and.
     >        beom_cart_set(geom2,nion2,byte_mb(tt(1)),
     >                                  dbl_mb(rt(1)),
     >                                  dbl_mb(qt(1)))
      value = value.and.beom_masses_get(geom1,nion1,dbl_mb(mt(1)))
      value = value.and.beom_masses_set(geom2,nion2,dbl_mb(mt(1)))

      if (nion3.gt.0) then
        value = value.and.
     >          beom_cart_set(geom3,nion3,byte_mb(tt(1)+16*nion2),
     >                                   dbl_mb(rt(1) + 3*nion2),
     >                                   dbl_mb(qt(1) +   nion2))
        value = value.and.beom_masses_set(geom3,nion3,
     >                                   dbl_mb(mt(1) +   nion2))
      end if


      call dcopy(nion1,0.0d0,0,dbl_mb(rt(1)),1)
      value = value.and.beom_vel_get(geom1,dbl_mb(rt(1)))
      value = value.and.beom_vel_set(geom2,dbl_mb(rt(1)))
      if (nion3.gt.0)
     > value = value.and.beom_vel_set(geom3,dbl_mb(rt(1)+3*nion2))

      value = value.and.beom_rtdb_store(rtdb,geom2,'pspwgeometry')
      if(nion3.gt.0) then
         value = value.and.beom_rtdb_store(rtdb,geom3,'chargegeometry')
         mmexist = .true.
      else
         mmexist = .false.
      end if
      value = value.and.beom_destroy(geom3)
      value = value.and.beom_destroy(geom2)
      value = value.and.beom_destroy(geom1)
      if (.not. value) 
     >    call errquit('chargepspwgeometry->pspwgeometry write',0,
     &       GEOM_ERR)
      value = BA_pop_stack(tt(2))
      value = value.and.BA_pop_stack(mt(2))
      value = value.and.BA_pop_stack(qt(2))
      value = value.and.BA_pop_stack(rt(2))
      if (.not. value) call errquit('popping stack',0, MA_ERR)

      return
      end


*     ***************************
*     *                         *
*     *    parsepointcharge     *
*     *                         *
*     ***************************
      logical function parsepointcharge(string)
      implicit none
      character*16 string

      logical qmmm

      qmmm = .false.
      if (index(string,'x').eq.1) then
         qmmm = .true.
         if (index(string,'e').eq.2)  qmmm = .false.
         if (index(string,'E').eq.2)  qmmm = .false.
      end if
      if (index(string,'X').eq.1) then
         qmmm = .true.
         if (index(string,'e').eq.2)  qmmm = .false.
         if (index(string,'E').eq.2)  qmmm = .false.
      end if
      if (index(string,'bq').eq.1)  qmmm = .true.
      if (index(string,'Bq').eq.1)  qmmm = .true.
      if (index(string,'bQ').eq.1)  qmmm = .true.
      if (index(string,'BQ').eq.1)  qmmm = .true.

      parsepointcharge = qmmm
      return
      end





*     ***************************
*     *                         *
*     *   combine_pointcharge   *
*     *                         *
*     ***************************
      subroutine combine_pointcharge(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"

      logical mmexist
      common / ion_charge / mmexist

*     **** local variables ****
      integer i,geom1,geom2,geom3,nion1,nion2,nion3,rt(2)
      logical value
      double precision rxyz(3),q
      character*16     t


*       **********************************************************
*       **** put together ions and molecular pseudopotentials ****
*       **********************************************************
        value = beom_create(geom1,'chargepspwgeometry')
        value = value.and.beom_create(geom2,'pspwgeometry')
        if (mmexist)
     >    value = value.and.beom_create(geom3,'chargegeometry')
        value = value.and.
     >          beom_rtdb_load(rtdb,geom1,'geometry')
        value = value.and.beom_ncent(geom1,nion1)
        value = value.and.beom_rtdb_load(rtdb,geom2,'pspwgeometry')
        value = value.and.beom_ncent(geom2,nion2)
        if (mmexist) then
           value = value.and.beom_rtdb_load(rtdb,geom3,'chargegeometry')
           value = value.and.beom_ncent(geom3,nion3)
        else
           nion3 = 0
        end if
        if (.not. value)
     >   call errquit('pspwgeometry->geometry write 1',0, GEOM_ERR)

        value = BA_push_get(mt_dbl,(3*nion1),'rt',rt(2),rt(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)

        do i=1,nion2
           value = value.and.beom_cent_get(geom2,i,t,rxyz,q)
           value = value.and.beom_cent_set(geom1,i,t,rxyz,q)
        end do
        do i=1,nion3
           value = value.and.beom_cent_get(geom3,i,t,rxyz,q)
           value = value.and.beom_cent_set(geom1,(i+nion2),t,rxyz,q)
        end do
        value = value.and.beom_vel_get(geom2,dbl_mb(rt(1)))
        if (nion3.gt.0)
     >     value = value.and.beom_vel_get(geom3,dbl_mb(rt(1+3*nion2)))
        value = value.and.beom_vel_set(geom1, dbl_mb(rt(1)))

        value = value.and.beom_rtdb_delete(rtdb,'geometry')
        value = value.and.beom_rtdb_delete(rtdb,'pspwgeometry')
        if (nion3.gt.0)
     >     value = value.and.beom_rtdb_delete(rtdb,'chargegeometry')
        value = value.and.
     >          beom_rtdb_store(rtdb,geom1,'geometry')
        if (mmexist)
     >     value = value.and.beom_destroy(geom3)
        value = value.and.beom_destroy(geom2)
        value = value.and.beom_destroy(geom1)
        value = value.and.BA_pop_stack(rt(2))
        if (.not. value)
     >   call errquit('pspwgeometry->geometry write 2',0,
     &       GEOM_ERR)

      return
      end






