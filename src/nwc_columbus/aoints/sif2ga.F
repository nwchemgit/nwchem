      subroutine sif2ga(g,locar,basis,nbft)

#include "global.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "errquit.fh"
#include "mafdecls.fh"


      integer g              ! atom blocked global array
      double precision locar(*) ! local array from sifs with all values
      integer basis
      integer nbft

      integer nshell
      integer ishell, jshell
      integer ijshell, ilo, ihi, jlo, jhi, blkdim
      integer max1e, mem1
      integer l_buf
      integer k_buf
      integer iloc, i, j, ijmap, il, jl
      logical shells_ok
      logical int_chk_sh
      external int_chk_sh

cgk debug
*     write(*,*)'gk: entered sif2ga'
cgk end
      if (.not. bas_numcont(basis, nshell)) call errquit
     &  ('sif2ga: bas_numcont failed for basis', basis,
     &     BASIS_ERR)

      call int_mem_1e(max1e, mem1)
cgk debug
*     write(*,*)'gk: max1e = ', max1e
cgk end

      if (.not. MA_push_get(MT_DBL,max1e,'sif2ga:buf',l_buf,k_buf))
     $  call errquit('sif2ga: ma failed', max1e, MA_ERR)
   
      call ga_zero(g)

cgk debug
*     write(*,*)'gk: nshell = ', nshell
cgk end
      ijmap = 0
      do jshell = 1, nshell
*     write(*,*)'gk: jshell = ', jshell
        if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     $    call errquit('sif2ga: bas_cn2bfr jshell', jshell, 
     &      BASIS_ERR)
*     write(*,*)'gk: jlo = ', jlo, ' jhi = ', jhi
        do ishell = 1, nshell
*     write(*,*)'gk: ishell = ', ishell
          if (.not. bas_cn2bfr(basis, ishell, ilo, ihi))
     $      call errquit('sif2ga: bas_cn2bfr ishell', ishell, 
     &        BASIS_ERR)
*     write(*,*)'gk: ilo = ', ilo, ' ihi = ', ihi
          call dfill(max1e, 0.0d0, dbl_mb(k_buf), 1)
          blkdim=ihi-ilo+1
          iloc=0
          do j = jlo, jhi
            do i = ilo, ihi
              jl=max(i,j)
              il=min(i,j)
              ijmap=(jl-1)*jl/2+il
*             write(*,2000) 'gk: j=',j,' i=',i,' ijmap=',ijmap,
*    &          ' iloc=',iloc
2000   format(4(a,i4))
*             write(*,2010) 'gk:              ',locar(ijmap)
2010   format(a,f20.12)
              dbl_mb(k_buf+iloc)=locar(ijmap)
              iloc=iloc+1
            enddo
          enddo
          call ga_put(g,ilo,ihi,jlo,jhi,dbl_mb(k_buf),blkdim)
        enddo
      enddo

      if (.not. MA_pop_stack(l_buf)) call errquit
     $  ('sif2ga: pop failed', 0, GA_ERR)

      return
      end
