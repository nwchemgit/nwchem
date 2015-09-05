      subroutine dimqm_printAtomicDipoles(rtdb, lpol, liscmplx, label)
      implicit none
#include "stdio.fh"
#include "util.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
#include "global.fh"
c
c   Input variables
      integer rtdb
      logical lpol
      logical liscmplx
      character(len=*) label
c
c   Local variables
      double precision diptot(3)
      character elements(2,nDIMTypes)
      integer tx(nDIM)
      integer l_muind, k_muind
      character*256 f
      integer i, i3, j, itype
      integer icmplx
      logical stat
      character*50 dd
      character*50 d
      integer l_name, k_name
      character*2 ele
      character aele(2)
      double precision mu(3)
      double precision mu_c(3)
      character*16 mu_str
c      character*16 q_str
c
c   Return without execution if this table is not requested
      if(.not.latmdip) return
      if(ga_nodeid().eq..0) then
        stat = rtdb_parallel(.false.)
c
c   Check if real or complex
      icmplx = 1
      if(liscmplx) icmplx = 2
      i3 = nDIM * 3
c
c   Determine which dipoles to grab
      mu_str = 'dimqm:muind'
c
c   ===============
c   Allocate Arrays
c   ===============
c
c   Allocate memory for induced dipoles 
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimqm muind',
     $                                  l_muind,k_muind))
     $  call errquit('printAtomicDipoles malloc muind failed',1,MA_ERR)
c   Allocate memory for atom names
c   Read in atom names
      if(.not.rtdb_get(rtdb,'dimpar:typeindex',mt_int,nDIM,tx))
     $  call errquit('printAtomicDipoles get tx failed',1,RTDB_ERR)
      if(.not.rtdb_get(rtdb, 'dimpar:name', mt_byte, nDIMTypes*2,
     $                 elements))
     $  call errquit('printAtomicDipoles get names failed',1,RTDB_ERR)

      dd = 
     $ '=================================================='
      d =
     $ '--------------------------------------------------'
c
c -------------------------------------------------
c Print out the dipoles for each atom
c -------------------------------------------------
c
      if(liscmplx) then

c   Read in induced dipoles
      if(.not.rtdb_get(rtdb,mu_str,mt_dbl,i3*icmplx,dbl_mb(k_muind)))
     $  call errquit('dimqm_print get muind(c) failed',1, RTDB_ERR)

        f = '(1X,A)'
        write (luout,f) dd
        write (luout,f) 'Induced dipoles for each DIM atom :'
        if (label .ne. '') write (luout,f) label
        write (luout,f) d
        f = '(5X,A,17X,A)'
        write (luout,f) 'ATOM', 'Dipole'
        f = '(17X,A1,11X,A1,11X,A1)'
        write (luout,f) 'X', 'Y', 'Z'
        f = '(1X,A)'
        write (luout,f) d
        f = '(1X,I5,1X,2A1,3F12.5,1X,A3)'
        do i = 1, nDIM
          itype = tx(i)
          do j = 1, 3
            mu(j)   = dbl_mb(k_muind+3*(i-1)+j-1)
            mu_c(j) = dbl_mb(k_muind+3*(i-1)+j-1+i3)
          end do
          write (luout,f) i, elements(1:2,itype),  mu(:), '(R)'
          write (luout,'(9X,3F12.5,1X,A3)')  mu_c(:), '(I)'
        end do
        f = '(1X,A)'
        write (luout,f) dd
        write (luout,*)

      else

c   Read in induced dipoles
      if(.not.rtdb_get(rtdb,mu_str,mt_dbl,i3,dbl_mb(k_muind)))
     $  call errquit('dimqm_print get muind(r) failed',1, RTDB_ERR)

        f = '(1X,A)'
        write (luout,f) dd
        write (luout,f) 'Induced dipoles for each DIM atom :'
        if (label .ne. '') write (luout,f) label
        write (luout,f) d
        f = '(5X,A,17X,A)'
        write (luout,f) 'ATOM', 'Dipole'
        f = '(17X,A1,11X,A1,11X,A1)'
        write (luout,f) 'X', 'Y', 'Z'
        f = '(1X,A)'
        write (luout,f) d
        f = '(1X,I5,1X,2A1,3F12.5)'
        do i =1, nDIM
          itype = tx(i)
          do j = 1, 3
            mu(j) = dbl_mb(k_muind+3*(i-1)+j-1)
          end do
          write (luout,f) i, elements(1:2,itype), mu(:)
        end do
        f = '(1X,A)'
        write (luout,f) dd
        write (luout,*)

      end if
      call util_flush(LuOut)
      stat = ma_chop_stack(l_muind)
      stat = rtdb_parallel(.true.)
      end if

      end subroutine dimqm_printAtomicDipoles
