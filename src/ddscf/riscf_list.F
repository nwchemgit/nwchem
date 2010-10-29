      block data riscf_three_eri_list
C$Id$
      implicit none
#include "riscf.fh"
      data riscf_info / dim*0 /
      data nbases / 0 /
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function riscf_put_3c_ints (ao_basis, ri_basis, handle, 
     $     index, npair, odisk)
C     put another set of 3-center ERIs into list

C     parameters:
C     ao_basis, ri_basis: basis set handles [input]
C     handle:             ERI handle (GA)   [input]
C     index:            index handle (GA)   [input]
C     npair:   # of non-zero (ij)-pairs     [input]
C     odisk: disk array (I/O based method) or global array (in-core) [output]

      implicit none

#include "riscf.fh"

      integer ao_basis, ri_basis, handle, npair, index, i
      logical odisk
      external riscf_three_eri_list ! For T3D linker

C     is it already there?
      do i = 1, maxbases
        if (riscf_info(1,i).eq.ao_basis) then
          if (riscf_info(2,i).eq.ri_basis) then
C           nothing to do
            goto 30
          endif
        endif
      enddo

      do i = 1, maxbases
        if (riscf_info(1,i).eq.0) then
          riscf_info(1,i) = ao_basis
          riscf_info(2,i) = ri_basis
          riscf_info(3,i) = handle
          riscf_info(4,i) = index
          riscf_info(5,i) = npair
          if ( odisk ) then
            riscf_info(6,i) = 1
          else
            riscf_info(6,i) = 0
          endif
          nbases = nbases + 1
          goto 30
        endif
      enddo

C     looks like we did not find a free slot
      write (6,*)'put_3c_ints: too many basis sets'
      riscf_put_3c_ints = .FALSE.
      return

 30   continue
      riscf_put_3c_ints = .TRUE.
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function riscf_get_3c_ints (ao_basis, ri_basis, handle, 
     $     index, npair, odisk)
C     find handle to three_center ERIs for a pair of basis sets

C     parameters:
C     ao_basis, ri_basis: basis set handles [input]
C     handle:             ERI handle (GA)   [output]
C     index:            index handle (GA)   [output]
C     npair:   # of non-zero (ij)-pairs     [output]
C     odisk: disk array (I/O based method) or global array (in-core) [output]

      implicit none

#include "riscf.fh"

      integer ao_basis, ri_basis, handle, index, npair, i
      logical odisk

C     search entry
      do i = 1, maxbases
        if (riscf_info(1,i).eq.ao_basis) then
          if (riscf_info(2,i).eq.ri_basis) then
            handle = riscf_info(3,i)
            index = riscf_info(4,i)
            npair = riscf_info(5,i)
            odisk = ( riscf_info(6,i) .gt. 0 )
            riscf_get_3c_ints = .TRUE.
            return
          endif
        endif
      enddo

      riscf_get_3c_ints = .FALSE.
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function riscf_del_3c_ints (handle)
C     delete one or all sets of three_center ERIs

C     parameter:
C     handle >= 0: array to be deleted
C             < 0: delete all

      implicit none

#include "riscf.fh"
#ifdef RISCF_USE_DRA
#include "dra.fh"
#endif

      integer handle, i
      logical status, ga_destroy
      logical riscf_del_index
      external riscf_del_index

      if ( handle.ge.0 ) then
C       search entry
        do i = 1, maxbases
          if (riscf_info(3,i).eq.handle .and. riscf_info(1,i).gt.0) then
            if ( riscf_info(6,i) .eq. 0 ) then
              riscf_del_3c_ints = ga_destroy (handle)
            else
#ifdef RISCF_USE_DRA
              riscf_del_3c_ints = ( dra_delete(handle) .eq. 0 )
#endif
            endif
            riscf_info(1,i) = 0
            riscf_info(2,i) = 0
            riscf_info(3,i) = 0
            riscf_info(5,i) = 0
            status = riscf_del_index (riscf_info(4,i))
            riscf_info(4,i) = 0
            nbases = nbases - 1
            return
          endif
        enddo
C       we did not find it
        write(6,*)'riscf_del_3c_int: handle not in list'
        riscf_del_3c_ints = .FALSE.
      else
        do i = 1, maxbases
          if (riscf_info(1,i).gt.0) then
            if ( riscf_info(6,i) .eq. 0 ) then
              status = ga_destroy (riscf_info(3,i))
            else
#ifdef RISCF_USE_DRA
              status = ( dra_delete (riscf_info(3,i)) .eq. 0 )
#endif
            endif
            if (.not. status ) then
              write(6,*)
     $           'riscf_del_3c_in: could not destroy array'
              riscf_del_3c_ints = .FALSE.
            else
              riscf_del_3c_ints = .TRUE.
            endif
            riscf_info(1,i) = 0
            riscf_info(2,i) = 0
            riscf_info(3,i) = 0
            status = riscf_del_index (riscf_info(4,i))
            riscf_info(4,i) = 0
          endif
        enddo
        nbases = 0
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function riscf_del_index (handle)
C     delete one or all sets of three_center ERIs

C     parameter:
C     handle >= 0: array to be deleted
C             < 0: delete all

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "riscf.fh"

      integer handle, i, ic
C     how often does it occur?
      ic = 0

      do i = 1, maxbases
        if (riscf_info(4,i).ne.0) then
          ic = ic + 1
        endif
      enddo

C     delete if it is only once
      riscf_del_index = .true.
      if ( ic.eq.1 )
     $     riscf_del_index = ga_destroy (handle)
      
      return
      end
