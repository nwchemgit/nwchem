*     *************************
*     *                       *
*     *     pspw_init_SIC     *
*     *                       *
*     *************************
      subroutine pspw_init_SIC(ne,rtdb)
      implicit none
      integer ne(2)
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      integer ma_type

      norbs = 0
*     **** get the number of SIC orbitals ****
      if (rtdb_ma_get(rtdb, 'pspw:SIC_orbitals', ma_type,
     >                norbs, orbital_list(2))) then

         if (.not.MA_get_index(orbital_list(2),orbital_list(1))) 
     >      call errquit(
     >       'pspw_init_SIC: ma_get_index failed for actlist',911)

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_screening_radius',
     >                      mt_dbl,1,rcut)) 
     >       rcut = -1.0d0

      else
         norbs = 0
      end if

   
      return
      end


*     *************************
*     *                       *
*     *     pspw_end_SIC      *
*     *                       *
*     *************************
      subroutine pspw_end_SIC()
      implicit none

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      logical value

      if (norbs.gt.0) then
        value =  MA_free_heap(orbital_list(2)) 
        if (.not. value)
     >  call errquit('pspw_end_SIC:error freeing heap memory',0)
      end if


      return
      end

*     *************************
*     *                       *
*     *     pspw_print_SIC    *
*     *                       *
*     *************************
      subroutine pspw_print_SIC(unit)
      implicit none
      integer unit

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      integer i

      if (norbs.ne.0) then
        write(unit,*) "      - SIC orbitals:",
     >      (int_mb(orbital_list(1)+i-1),i=1,norbs)
        write(unit,*)
      end if

      return
      end


