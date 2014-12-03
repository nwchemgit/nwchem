*     
* $Id$
*     

*     ***************************
*     *                         *
*     *    seperate_molpsp      *
*     *                         *
*     ***************************
      subroutine seperate_molpsp(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"



      logical bqbq
      common / nwpw_bqbq_block / bqbq

      logical mmexist
      common / ion_qmmm/ mmexist

      integer nionall
      common / ion_nionall_block / nionall

*     *** local variables ***
      integer     i,geom1,geom2,geom3,nion1,nion2,nion3
      integer rt(2),tt(2),qt(2),mt(2)
      double precision q,rxyz(3)
      character*16     t
      logical value

*     **** external functions ****
      logical          parseqmmm
      external         parseqmmm


*     *******************************************************************
*     **** seperate ions and charges from molecular pseudopotentials ****
*     *******************************************************************

*     **** generate chargepspwgeometry from geometry ****
*     **** generate qmmmgeometry from geometry ****
      value = geom_create(geom1,'geometry')
      value = value.and.geom_create(geom2,'chargepspwgeometry')
      value = value.and.geom_create(geom3,'qmmmgeometry')
      value = value.and.geom_rtdb_load(rtdb,geom1,'geometry')
      value = value.and.geom_ncent(geom1,nion1)
      if (.not. value) call errquit('opening geometry',0, GEOM_ERR)

*     *** set nionall and bqbq ****
      nionall = nion1
      bqbq    = geom_include_bqbq(geom1)

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
         value = geom_cent_get(geom1,i,t,rxyz,q)
         if (.not.parseqmmm(t)) then
            nion2 = nion2 + 1
         else
            nion3 = nion3 + 1
         end if
      end do
      value = value.and.
     >        geom_cart_get(geom1,nion1,byte_mb(tt(1)),
     >                                  dbl_mb(rt(1)),
     >                                  dbl_mb(qt(1)))
      value = value.and.
     >        geom_cart_set(geom2,nion2,byte_mb(tt(1)),
     >                                  dbl_mb(rt(1)),
     >                                  dbl_mb(qt(1)))
      value = value.and.geom_masses_get(geom1,nion1,dbl_mb(mt(1)))
      value = value.and.geom_masses_set(geom2,nion2,dbl_mb(mt(1)))

      if (nion3.gt.0) then
        value = value.and.
     >          geom_cart_set(geom3,nion3,byte_mb(tt(1)+16*nion2),
     >                                   dbl_mb(rt(1) + 3*nion2),
     >                                   dbl_mb(qt(1) +   nion2))
        value = value.and.geom_masses_set(geom3,nion3,
     >                                   dbl_mb(mt(1) +   nion2))
      end if


      call dcopy(nion1,0.0d0,0,dbl_mb(rt(1)),1)
      value = value.and.geom_vel_get(geom1,dbl_mb(rt(1)))
      value = value.and.geom_vel_set(geom2,dbl_mb(rt(1)))
      if (nion3.gt.0)
     > value = value.and.geom_vel_set(geom3,dbl_mb(rt(1)+3*nion2))

      value = value.and.geom_rtdb_store(rtdb,geom2,'chargepspwgeometry')
      if(nion3.gt.0) then
         value = value.and.geom_rtdb_store(rtdb,geom3,'qmmmgeometry')
         mmexist = .true.
      else
         mmexist = .false.
      end if
      value = value.and.geom_destroy(geom3)
      value = value.and.geom_destroy(geom2)
      value = value.and.geom_destroy(geom1)
      if (.not. value) 
     >   call errquit('geometry->chargepspwgeometry write',0, GEOM_ERR)
      value = BA_pop_stack(tt(2))
      value = value.and.BA_pop_stack(mt(2))
      value = value.and.BA_pop_stack(qt(2))
      value = value.and.BA_pop_stack(rt(2))
      if (.not. value) call errquit('popping stack',0, MA_ERR)

      return
      end

*     ***************************
*     *                         *
*     *      parseqmmm          *
*     *                         *
*     ***************************
      logical function parseqmmm(string)
      implicit none
      character*16 string

      logical qmmm

      qmmm = .false.
      if (index(string,'^').gt.0)   qmmm = .true.

      parseqmmm = qmmm
      return
      end


*     ***************************
*     *                         *
*     *     combine_molpsp      *
*     *                         *
*     ***************************
      subroutine combine_molpsp(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"

      logical mmexist
      common / ion_qmmm/ mmexist

*     **** local variables ****
      integer i,geom1,geom2,geom3,nion1,nion2,nion3,rt(2)
      logical value
      double precision rxyz(3),q
      character*16     t



*       **********************************************************
*       **** put together ions and molecular pseudopotentials ****
*       **********************************************************
        value = geom_create(geom1,'geometry')
        value = value.and.geom_create(geom2,'chargepspwgeometry')
        if (mmexist)
     >    value = value.and.geom_create(geom3,'qmmmgeometry')
        value = value.and.geom_rtdb_load(rtdb,geom1,'geometry')
        value = value.and.geom_ncent(geom1,nion1)
        value = value.and.
     >          geom_rtdb_load(rtdb,geom2,'chargepspwgeometry')
        value = value.and.geom_ncent(geom2,nion2)
        if (mmexist) then
           value = value.and.geom_rtdb_load(rtdb,geom3,'qmmmgeometry')
           value = value.and.geom_ncent(geom3,nion3)
        else
           nion3 = 0
        end if
        if (.not. value)
     >   call errquit('chargepspwgeometry->geometry write 1',0,
     &       GEOM_ERR)

        value = BA_push_get(mt_dbl,(3*nion1),'rt',rt(2),rt(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)

        do i=1,nion2
           value = value.and.geom_cent_get(geom2,i,t,rxyz,q)
           value = value.and.geom_cent_set(geom1,i,t,rxyz,q)
        end do
        do i=1,nion3
           value = value.and.geom_cent_get(geom3,i,t,rxyz,q)
           value = value.and.geom_cent_set(geom1,(i+nion2),t,rxyz,q)
        end do
        value = value.and.geom_vel_get(geom2,dbl_mb(rt(1)))
        if (nion3.gt.0)
     >     value = value.and.geom_vel_get(geom3,dbl_mb(rt(1+3*nion2)))
        value = value.and.geom_vel_set(geom1, dbl_mb(rt(1)))

        value = value.and.geom_rtdb_delete(rtdb,'geometry')
        value = value.and.geom_rtdb_delete(rtdb,'chargepspwgeometry')
        if (nion3.gt.0)
     >     value = value.and.geom_rtdb_delete(rtdb,'qmmmgeometry')
        value = value.and.geom_rtdb_store(rtdb,geom1,'geometry')
        if (mmexist)
     >     value = value.and.geom_destroy(geom3)
        value = value.and.geom_destroy(geom2)
        value = value.and.geom_destroy(geom1)
        value = value.and.BA_pop_stack(rt(2))
        if (.not. value)
     >   call errquit('chargepspwgeometry->geometry write 2',0,
     &       GEOM_ERR)

      return
      end






