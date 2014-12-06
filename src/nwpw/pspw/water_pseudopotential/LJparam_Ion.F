*     *************************
*     *                       *
*     *     LJparam_init	  *
*     *                       *
*     *************************
      subroutine LJparam_init(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"


*     **** LJparam common block ****
      integer ljcount
      integer ion_number(2)
      integer ei(2),si(2)
      real*8 ew,sw
      common / LJparam_Block / ew,sw,ei,si,ion_number,ljcount

*     **** local variables
      logical value
      integer taskid
      integer MASTER
      parameter(MASTER=0)
      integer i,j

      value = btdb_parallel(.true.)
      value = value.and.
     >        btdb_get(rtdb,'water_pseudopotential:lj_count',
     >                 mt_int,1,ljcount)
      if (.not. value) call errquit('error reading rtdb',0, RTDB_ERR)


      if (ljcount.gt.0) then

        value = BA_alloc_get(mt_int,(ljcount),'ion_number',
     >                       ion_number(2),ion_number(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,(ljcount),'ei',ei(2),ei(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,(ljcount),'si',si(2),si(1))
        if (.not. value) call errquit('out of heap memory',0, MA_ERR)

        value = btdb_get(rtdb,'water_pseudopotential:lj_ion_number',
     >                   mt_int,ljcount,int_mb(ion_number(1)))
        value = value.and.
     >          btdb_get(rtdb,'water_pseudopotential:lj_ei',
     >                   mt_dbl,ljcount,dbl_mb(ei(1)))
        value = value.and.
     >          btdb_get(rtdb,'water_pseudopotential:lj_si',
     >                   mt_dbl,ljcount,dbl_mb(si(1)))
        value = value.and.
     >          btdb_get(rtdb,'water_pseudopotential:lj_ew',
     >                   mt_dbl,1,ew)
        value = value.and.
     >          btdb_get(rtdb,'water_pseudopotential:lj_sw',
     >                   mt_dbl,1,sw)
        if (.not. value) call errquit('error reading rtdb',0, RTDB_ERR)
      
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
          write(*,*)
          write(*,*) "LJ Parameters:"
          do i=1,ljcount
            j = int_mb(ion_number(1)+i-1)
            write(*,'(A11,I4,A12,E10.3,A12,E10.3)') 
     >                  "ion #",j,
     >                  "epsilon:",dbl_mb(ei(1)+j-1), 
     >                  "sigma:",  dbl_mb(si(1)+j-1)
          end do
          write(*,'(A11,A16,E10.3,A12,E10.3)') "water", "epsilon:",ew,
     >                                               "sigma:",sw
          write(*,*)

         end if
      end if

      return
      end

*     *************************
*     *                       *
*     *     LJparam_end		  *
*     *                       *
*     *************************
      subroutine LJparam_end()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"

*     **** LJparam common block ****
      integer ljcount
      integer ion_number(2)
      integer ei(2),si(2)
      real*8 ew,sw
      common / LJparam_Block / ew,sw,ei,si,ion_number,ljcount

      logical value

      if (ljcount.gt.0) then
         value = BA_free_heap(ion_number(2))
         value = value.and.BA_free_heap(ei(2))
         value = value.and.BA_free_heap(si(2))
         if (.not.value) call errquit('error BA_free_heap',0, MA_ERR)
      end if

      return
      end


*     ***************************
*     *							*
*     *         LJparam_Ion		*
*     *							*
*     ***************************

      subroutine LJparam_Ion(i,e,s)
      implicit none
      integer i
      real*8  e,s

#include "bafdecls.fh"

*     **** LJparam common block ****
      integer ljcount
      integer ion_number(2)
      integer ei(2),si(2)
      real*8 ew,sw
      common / LJparam_Block / ew,sw,ei,si,ion_number,ljcount

*     **** local variables ****
      logical found
      integer ii,j

      found = .false.
      do ii=1,ljcount
         if (int_mb(ion_number(1)+ii-1) .eq. i) then
             j = ii
             found = .true.
         end if
      end do

      if (found) then
         e = dbl_mb(ei(1)+j-1)
         s = dbl_mb(si(i)+j-1)
      else
         e = 0.0d0
         s = 0.0d0
      end if

      return
      end
      
*     ***************************
*     *							*
*     *         LJparam_water	*
*     *							*
*     ***************************

      subroutine LJparam_water(e,s)
      implicit none
      real*8  e,s

*     **** LJparam common block ****
      integer ljcount
      integer ion_number(2)
      integer ei(2),si(2)
      real*8 ew,sw
      common / LJparam_Block / ew,sw,ei,si,ion_number,ljcount

      e = ew
      s = sw
      return
      end
c $Id$
