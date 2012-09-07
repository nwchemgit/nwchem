C> \ingroup basis
C> @{
C>
C> \brief Compute the checksum of a basis set instance
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function bas_checksum(basisin, sum)
      implicit none
*
* $Id$
*
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
      integer basisin  !< [Input] The basis set handle
      character*32 sum !< [Output] The checksum
c
      integer dtob, itob, ltob, basis, ntags, i
      logical bas_check_handle
      external bas_check_handle
#include "bas_exndcf.fh"
c
      bas_checksum = bas_check_handle(basisin, 'bas_checksum')
      if (.not. bas_checksum) return
c
      dtob = ma_sizeof(mt_dbl, 1, mt_byte)
      itob = ma_sizeof(mt_int, 1, mt_byte)
      ltob = ma_sizeof(mt_log, 1, mt_byte)
c
      basis = basisin + BASIS_HANDLE_OFFSET
      ntags = infbs_head(HEAD_NTAGS, basis)
c
c     Checksum just the unique basis info
c
c     Spherical
c     Header
c     Tags
c     Coefficients and exponents
c
      call checksum_init()
      call checksum_update(ltob, bas_spherical(basis))
      call checksum_update(itob*ndbs_head, infbs_head(1,basis))
      do i = 1, ntags
         call checksum_char_update(bs_tags(i,basis))
      enddo
      call checksum_update(dtob*exndcf(SZ_exndcf,basis),
     $     dbl_mb(exndcf(K_exndcf,basis)))
      call checksum_final(sum)
c
      end
C> @}
