      subroutine tce_create_movecs(rtdb)
c
c $Id$
c
c     Initialize for many-electron calculations.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "errquit.fh"
c#include "tce.fh"
c#include "tce_main.fh"
c#include "tce_diis.fh"
c#include "tce_prop.fh"
#include "dra.fh"
c
      integer rtdb             ! RTDB handle
      integer nbf_verify       ! nbf read from movecs_read_header
      integer ipol_verify      ! ipol read from movecs_read_header
      character*255 movecs_out ! MO vector file name
      character*255 title      ! Title
      character*255 basisname  ! Basis set name
      character*255 scftype    ! Type of SCF
      character*4 hftype
      character*3 onoff1,onoff2
      integer i,j,n
      integer nclosed, nopen
      double precision a
      character*60 ctype
      character*60 cperturbative
      integer itol2e
      logical nodezero
      logical dummy
      logical rohf
      logical int_normalize
      external int_normalize
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical xc_gotxc
      external xc_gotxc
      character*10 module
      character*8 multchar(8)
       character*26 date
       integer ma_type
c
c
      nodezero=(ga_nodeid().eq.0)
c
c     ----------
c     MO vectors
c     ----------
c
      if (reference.eq.0) then
        if (.not.rtdb_cget(rtdb,'dft:output vectors',1,movecs_out))
     1   call errquit('tce_create_movecs: no final MO found',0,RTDB_ERR)
      else if (reference.eq.1) then
        if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_out))
     1   call errquit('tce_create_movecs: no final MO found',0,RTDB_ERR)
      endif
c
c     -------------------------------------
c     Allocate global arrays for MO vectors
c     -------------------------------------
c
      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nbf,'eigenvalues',
     1    l_evl(i),k_evl(i))) call errquit
     2    ('tce_create_movecs: failed to allocate evl',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,nbf,'occupancy',
     1    l_occ(i),k_occ(i))) call errquit
     2    ('tce_create_movecs: failed to allocate occ',0,MA_ERR)
        if (.not.ma_push_get(mt_int,nbf,'irreps',
     1    l_irs(i),k_irs(i))) call errquit
     2    ('tce_create_movecs: failed to allocate irs',0,MA_ERR)
      enddo
      if (.not.movecs_read_header(movecs_out,title,basisname,
     1  scftype,nbf_verify,ipol_verify,nmo,2)) call errquit
     2  ('tce_create_movecs: failed to read MO file header',0,RTDB_ERR)
      if (nbf.ne.nbf_verify) call errquit
     1  ('tce_create_movecs: corrupted MO vectors',0,RTDB_ERR)
      if (rohf) nmo(2) = nmo(1)
      do i=1,ipol
        g_movecs(i)=ga_create_atom_blocked
     1    (geom,ao_bas_han,'MO eigenvectors')
        if (rohf.and.(i.eq.2)) then
          call ga_copy(g_movecs(1),g_movecs(2))
          do j = 1,nbf
            dbl_mb(k_evl(2)+j-1) = dbl_mb(k_evl(1)+j-1)
            dbl_mb(k_occ(2)+j-1) = dbl_mb(k_occ(1)+j-1)
          enddo
        else
          if (.not.movecs_read(movecs_out,i,
     1      dbl_mb(k_occ(i)),dbl_mb(k_evl(i)),g_movecs(i)))
     2      call errquit('tce_create_movecs: failed to read MO vectors',
     3      0,RTDB_ERR)
        endif
        if (util_print('MO vectors',print_debug)) then
          if (nodezero) call ma_print(dbl_mb(k_evl(i)),
     1      nbf,1,'Eigenvalues')
          if (nodezero) call ma_print(dbl_mb(k_occ(i)),
     1      nbf,1,'Occupancies')
          call ga_print(g_movecs(i))
        endif
      enddo
      if (reference.eq.0) then
        if (.not.rtdb_get(rtdb,'dft:alpha irreps',mt_int,
     1    nbf,int_mb(k_irs(1))))
     2    call errquit('tce_create_movecs: no alpha irreps found',0,
     3    RTDB_ERR)
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:alpha irreps',mt_int,
     1    nbf,int_mb(k_irs(1))))
     2    call errquit('tce_create_movecs: no alpha irreps found',1,
     3    RTDB_ERR)
      endif
      do i=1,nbf
        if (oadapt) then
          int_mb(k_irs(1)+i-1) = int_mb(k_irs(1)+i-1) - 1
        else
          int_mb(k_irs(1)+i-1) = 0
        endif
      enddo
      if (ipol.eq.2) then
        if (reference.eq.0) then
          if (.not.rtdb_get(rtdb,'dft:beta irreps',mt_int,
     1      nbf,int_mb(k_irs(2))))
     2      call errquit('tce_create_movecs: no beta irreps found',0,
     3      RTDB_ERR)
        else if (reference.eq.1) then
          if (rohf) then
            do j = 1,nbf
              int_mb(k_irs(2)+j-1) = int_mb(k_irs(1)+j-1) + 1
            enddo
          else
            if (.not.rtdb_get(rtdb,'scf:beta irreps',mt_int,
     1        nbf,int_mb(k_irs(2))))
     2        call errquit('tce_create_movecs: no beta irreps found',1,
     3        RTDB_ERR)
          endif
        endif
        do i=1,nbf
          if (oadapt) then
            int_mb(k_irs(2)+i-1) = int_mb(k_irs(2)+i-1) - 1
          else
            int_mb(k_irs(2)+i-1) = 0
          endif
        enddo
      endif
      if (nodezero) call util_flush(LuOut)
c
c     ------
c     Return
c     ------
c
      return
      end
