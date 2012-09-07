C> \ingroup geom
C> @{
C>
C> \brief Compute the checksum for a geometry instance
C> 
C> Computes a checksum for a geometry instance using MD5 technology.
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function geom_checksum(geom, sum)
*
* $Id$
*
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
      integer geom     !< [Input] The geometry handle
      character*32 sum !< [Output] The checksum
c
      integer dtob, itob, natom, i
      logical geom_check_handle
      external geom_check_handle
c
      geom_checksum = geom_check_handle(geom, 'geom_checksum')
      if (.not. geom_checksum) return
      dtob = ma_sizeof(mt_dbl, 1, mt_byte)
      itob = ma_sizeof(mt_int, 1, mt_byte)
      natom = ncenter(geom)
c
      call checksum_init()
      call checksum_update(dtob*3*natom, coords(1,1,geom))
      call checksum_update(dtob*natom,   charge(1,geom))
      call checksum_update(dtob*natom,   geom_mass(1,geom))
      call checksum_update(dtob*3,       efield(1,geom))
      do i = 1, natom
         call checksum_char_update(tags(i,geom))
      enddo
      call checksum_char_update(group_name(geom))
      call checksum_char_update(user_units(geom))
c
      call checksum_final(sum)
c
      end
C> @}
