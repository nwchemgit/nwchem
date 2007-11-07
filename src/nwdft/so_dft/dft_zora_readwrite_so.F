c
czora...Inquire if the zora correction file is present
c
      logical function dft_zora_inquire_file_so(filename)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c
      character*(*) filename
      logical found
c
      call ga_sync()
c
c     Inquire if file is present
      inquire(file=filename,exist=found)
      dft_zora_inquire_file_so = found
c
      call ga_sync()
c
      return
      end
c
czora...Read in the zora atomic corrections from disk
c
      logical function dft_zora_read_so(filename, nbf, nsets, nmo, 
     &               mult, g_zora_sf, g_zora_scale_sf, 
     &               g_zora_so, g_zora_scale_so) 
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c
      character*(*) filename
      integer iset              ! restricted or unrestricted
c
      integer g_zora_sf(2)
      integer g_zora_scale_sf(2)
      integer g_zora_so(3)
      integer g_zora_scale_so(3)
c
      integer nsets             ! restricted or unrestricted
      integer nbf               ! No. of functions in basis
      integer nmo(nsets)       
      integer mult
      integer ok, jset, i, j
      integer l_zora, k_zora
c
      integer unitno
      parameter (unitno = 777)
      integer inntsize,ddblsize
c     
      integer nsets_read
      integer nbf_read
      integer mult_read
c
      l_zora = -1               ! An invalid MA handle
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
c
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
c
c     Print a message indicating the file being read
      write(6,22) filename(1:inp_strlen(filename))
 22   format(/' Read atomic ZORA corrections from ',a/)
      call util_flush(luout)
c
c     Open the file
      open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nsets_read
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) mult_read
c
c      Error checks
       if ((nsets_read .ne. nsets) 
     &  .or. (nbf_read .ne. nbf) 
     &  .or. (mult_read .ne. mult) ) goto 1003
c
c     Allocate the temporary buffer
      if (.not.ma_push_get(mt_dbl,nbf,'dft_zora_read_so',l_zora,k_zora))
     $        call errquit('dft_zora_read_so: ma failed', nbf, MA_ERR)
c
c     Read in g_zora_sf
      do iset = 1,2  ! 2 components
       do i = 1,nbf
        call sread(unitno, dbl_mb(k_zora), nbf)
        call ga_put(g_zora_sf(iset), 1, nbf, i, i, dbl_mb(k_zora), 1)
       end do
      end do
c
c     Read in g_zora_scale_sf
      do iset = 1,2  ! 2 components
       do i = 1,nbf
        call sread(unitno, dbl_mb(k_zora), nbf)
        call ga_put(g_zora_scale_sf(iset), 1, nbf, i, i, 
     &   dbl_mb(k_zora), 1)
       end do
      end do
c
c     Read in g_zora_so
      do iset = 1,3  ! 3 components: x,y,z
       do i = 1,nbf
        call sread(unitno, dbl_mb(k_zora), nbf)
        call ga_put(g_zora_so(iset), 1, nbf, i, i, dbl_mb(k_zora), 1)
       end do
      end do
c
c     Read in g_zora_scale_so
      do iset = 1,3  ! 3 components: x,y,z
       do i = 1,nbf
        call sread(unitno, dbl_mb(k_zora), nbf)
        call ga_put(g_zora_scale_so(iset), 1, nbf, i, i, 
     &   dbl_mb(k_zora), 1)
       end do
      end do
c
c     Close the file
      close(unitno,err=1002)
      ok = 1
c
c     Deallocate the temporary buffer
      if (.not. ma_pop_stack(l_zora)) call errquit
     $      ('dft_zora_read_so: pop failed', l_zora, MA_ERR)
c
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
      dft_zora_read_so = ok .eq. 1
      return
c
 1000 write(6,*) 'dft_zora_read_so: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) 'dft_zora_read_so: failed to read ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1003 write(6,*) 'dft_zora_read_so: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) 'dft_zora_read_so: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
c
czora...Write out the zora atomic corrections to disk
c
      logical function dft_zora_write_so(rtdb, basis, filename,
     &     nbf, nsets, nmo, mult, 
     &     g_zora_sf, g_zora_scale_sf,
     &     g_zora_so, g_zora_scale_so)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
c     Temporary routine
c
      integer rtdb              ! [input] RTDB handle (-1 if not accessible)
      integer basis             ! [input] Basis handle(-1 if not accessible)
      character*(*) filename    ! [input] File to write to
      integer nbf               ! [input] No. of functions in basis
      integer nsets             ! [input] No. of sets of matrices
      integer nmo(nsets)        ! [input] No. of mos in each set
      integer mult
c
      integer g_zora_sf(2)    
      integer g_zora_scale_sf(2)  
      integer g_zora_so(3)    
      integer g_zora_scale_so(3)  
c
      integer unitno
      parameter (unitno = 777)
c
      integer lentit
      integer lenbas
      integer l_zora, k_zora
      integer ok, iset, i, j
      integer geom, ma_type, nelem
      character*26 date
      character*32 geomsum, basissum, key
      character*20 scftype20 
      character*128 basis_name, trans_name
      double precision energy, enrep
      integer inntsize
c
      l_zora = -1               ! An invalid MA handle
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
      ok = 0
c
c     Read routines should be consistent with this
c
c     Write out the atomic zora corrections
c
      if (ga_nodeid() .eq. 0) then
c
c     Open the file
      open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c
c     Write out the number of sets and basis functions
      write(unitno, err=1001) nsets
      write(unitno, err=1001) nbf
      write(unitno, err=1001) mult
c
c     Allocate the temporary buffer
      if (.not.ma_push_get(mt_dbl,nbf,'dft_zora_write_so',
     & l_zora,k_zora))
     &   call errquit('dft_zora_write_so: ma failed', nbf, MA_ERR)
c
c     Write out g_zora_sf
      do iset = 1,2  ! 2 components
       do i = 1, nbf
        call ga_get(g_zora_sf(iset), 1, nbf, i, i, dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Write out g_zora_scale_sf
      do iset = 1,2  ! 2 components
       do i = 1, nbf
        call ga_get(g_zora_scale_sf(iset), 1, nbf, i, i, 
     &              dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Write out g_zora_so
      do iset = 1,3  ! 3 components: x,y,z
       do i = 1, nbf
        call ga_get(g_zora_so(iset), 1, nbf, i, i, dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Write out g_zora_scale_so
      do iset = 1,3  ! 3 components: x,y,z
       do i = 1, nbf
        call ga_get(g_zora_scale_so(iset), 1, nbf, i, i, 
     &              dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Deallocate the temporary buffer
      if (.not. ma_pop_stack(l_zora))
     $  call errquit('dft_zora_write_so: ma pop failed', l_zora, MA_ERR)
c
c     Close the file
      close(unitno,err=1002)
c
      ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
      dft_zora_write_so = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote atomic ZORA corrections to ',a/)
         call util_flush(luout)
      endif
c
      return
c
 1000 write(6,*) 'dft_zora_write_so: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) 'dft_zora_write_so: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) 'dft_zora_write_so: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
