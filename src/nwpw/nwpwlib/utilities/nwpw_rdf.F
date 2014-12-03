c
c     $Id$
c

*     *********************************
*     *                               *
*     *          nwpw_rdf             *
*     *                               *
*     *********************************
      logical function nwpw_rdf(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"
#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer indx1(2),n1,ma1
      integer indx2(2),n2,ma2
      character*50 filename
      character*255 mfilename,gfilename

      value = .true.

*     **** set mfilename and gfilename ****
      if (.not.btdb_cget(rtdb,'cpmd:ion_motion_filename',1,filename))
     >  call util_file_prefix('ion_motion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             mfilename)

      if (.not.btdb_cget(rtdb,'nwpw:gr_filename',1,filename))
     >  call util_file_prefix('rdf',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             gfilename)

      if (
     >  rtdb_ma_get(rtdb,'nwpw:rdf_index1_atoms',ma1,n1,indx1(2)).and.
     >  rtdb_ma_get(rtdb,'nwpw:rdf_index2_atoms',ma2,n2,indx2(2))
     >   ) then
         
        value =           BA_get_index(indx1(2),indx1(1))
        value = value.and.BA_get_index(indx2(2),indx2(1))
        if (.not.value) 
     >  call errquit('nwpw_rdf:error allocatin heap',0,MA_ERR)

*        **** calculate rdf ****
         call pair_dist2(mfilename,gfilename,
     >                   int_mb(indx1(1)),n1,
     >                   int_mb(indx2(1)),n2)

         value =           BA_free_heap(indx1(2))
         value = value.and.BA_free_heap(indx2(2))
         if (.not.value) 
     >   call errquit('nwpw_rdf:error freeing heap',0,MA_ERR)
  
      end if

      nwpw_rdf = value
      return
      end
