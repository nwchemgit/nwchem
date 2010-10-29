       subroutine moints_op_print(ostart,olo,ohi,nbf,g_a)
*
* $Id$
*
       implicit none
#include "global.fh"
#include "mafdecls.fh"       
       integer ostart,olo,ohi,nbf
       integer g_a
       integer j,k,jk,jj,kk,ofroz
       integer my_id, ilo, ihi, jlo, jhi, ooff
       integer k_local, ld_local

       ofroz = ostart - 1
       my_id = ga_nodeid()
       call ga_distribution(g_a,my_id,ilo,ihi,jlo,jhi)
       ooff = ((olo-ofroz)*(olo-ofroz-1))/2
       if (ga_nodeid().eq.0) then
         do j=olo,ohi
           jj = j - ofroz
           do k=ostart,j
             kk = k - ofroz 
             jk = (jj*(jj-1))/2 + kk - ooff
             if ((jk.ge.jlo).and.(jk.le.jhi)) then
               write(6,901) j,k,jk
 901           format(//,'Operator: [',i2,',',i2,']','  index:',i5/)
               call ga_access(g_a,ilo,ihi,jk,jk,k_local,ld_local)
               call moints_matprint(nbf,nbf,dbl_mb(k_local))
             endif
           enddo
         enddo
       endif
       call ga_sync()

       return
       end




      subroutine moints_matprint(n,m,x)
      implicit none
      integer n,m
      double precision x(n,m)
      integer i, j, mxrow

      mxrow = min(m,10)
      do i=1,min(n,20)
        write(6,901) (x(i,j),j=1,mxrow)
 901    format(10f12.6)
      enddo
      write(6,*)
      return
      end






      subroutine moints_matprintx(n,m,x)
      implicit none
      integer n,m
      double precision x(n,m)
      integer i, j, mxrow

      mxrow = min(m,6)
      do i=1,min(n,20)
        write(6,901) (x(i,j),j=1,mxrow)
 901    format(6f20.14)
      enddo
      write(6,*)
      return
      end



