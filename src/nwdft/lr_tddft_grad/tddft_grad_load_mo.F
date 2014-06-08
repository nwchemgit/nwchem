c     Load the MO vectors from the file specified by fn_movecs.
c     Apart from the coefficients we also need the eigenvalues.
c     The orbital occupations can be discarded.
c     HvD 5/2006, NG 11/2012
c
      subroutine tddft_grad_load_mo(ipol,nao,fn_movecs,e,g_mo)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
c
c     Input:
c
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nao              ! the number of AOs
      character*255 fn_movecs  ! filename for MO vectors
c
c     Output:
c
      double precision e(nao,ipol) ! the MO energies
      integer          g_mo(2)     ! the MO coefficients
c
c     Functions:
c
      logical  movecs_read
      external movecs_read
c
c     Local:
c
      integer i     ! counter
      integer k_occ ! index for orbital occupations
      integer l_occ ! handle for orbital occupations
      character*32 pname
c
      pname ="tddft_grad_load_mo: "
c
      if (.not.ma_push_get(mt_dbl,nao,'occ',l_occ,k_occ)) 
     +  call errquit(pname//'failed to allocate occ',0,MA_ERR)
c
      do i = 1, ipol
       if (.not.movecs_read(fn_movecs,i,dbl_mb(k_occ),e(1,i),g_mo(i))) 
     +  call errquit(pname//'failed to read MO vectors',0, DISK_ERR)
      enddo
c
      if (.not.ma_pop_stack(l_occ)) 
     + call errquit(pname//'failed to release occ',0, MA_ERR)
c
      end
c $Id$
