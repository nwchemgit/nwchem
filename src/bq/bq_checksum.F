C> \ingroup bq
C> @{
C>
C> \brief Compute the checksum for a Bq instance
C>
C> Computes the checksum for a Bq instance using MD5 technology. This
C> checksum is used for example to compare the current Bq instance
C> against the one for which movecs were calculated (to test whether
C> we can mark a wavefunction as converged already during a previous
C> calculation).
C>
C> \return .TRUE. if successfull, and .FALSE. otherwise.
C>
      logical function bq_checksum(bq_handle,sum)
      implicit none
C
C     $Id$
C
#include "bq_params.fh"
#include "mafdecls.fh"
C
      integer bq_handle !< [Input] The handle of the Bq instance
      character*32 sum  !< [Output] The checksum
c
      character*32    pname     !< Procedure name
      character*(255) namespace !< Name of Bq instance
      integer ncent             !< The number of centers
      integer i_c               !< The index of the coordinates
      integer i_q               !< The index of the charges
      integer dtob              !< Doubles to bytes factor
c
      logical  bq_check_handle, bq_ncenter
      external bq_check_handle, bq_ncenter
c
      logical  bq_namespace, bq_index_coord, bq_index_charge
      external bq_namespace, bq_index_coord, bq_index_charge
c
      pname = "bq_checksum"
c
      bq_checksum = bq_check_handle(bq_handle)
      if (.not.bq_checksum) return
c
      bq_checksum = bq_ncenter(bq_handle,ncent)
      if (.not.bq_checksum) return
c
      bq_checksum = bq_namespace(bq_handle,namespace)
      if (.not.bq_checksum) return
c
      bq_checksum = bq_index_coord(bq_handle,i_c)
      if (.not.bq_checksum) return
c
      bq_checksum = bq_index_charge(bq_handle,i_q)
      if (.not.bq_checksum) return
c
      dtob = ma_sizeof(mt_dbl, 1, mt_byte)
c
      call checksum_init()
      call checksum_update(dtob*ncent,   dbl_mb(i_q))
      call checksum_update(dtob*3*ncent, dbl_mb(i_c))
      call checksum_char_update(namespace)
      call checksum_final(sum)
c
      bq_checksum = .true.
c
      return
      end
C>
C> @}
