      subroutine hv_write(lun,ix,len,x)
      implicit none
      integer lun, ix, len
      double precision x(len)
      integer i, k, khi
      
      write(lun,901) ix,(x(i),i=1,8)
 901  format(i5,8f12.6)
      do k=9,len,8
        khi = min((k+7),len)
        write(lun,902) (x(i),i=k,khi)
 902    format(5x,8f12.6)
      enddo
      return
      end
      



      subroutine hv_reorder( nbf, nclosed, nact, g_x, x, y )
      implicit none
      integer nbf, nclosed, nact
      integer g_x
      double precision x(*), y(*)
c     
      integer nvir, vlen, voff, aoff, aend, xoff, xend, ii
      integer nca, i, j
c     
c     
      nvir = nbf - nclosed - nact
      vlen = (nclosed+nact)*nvir + nclosed*nact
      voff = nclosed + nact + 1
      aoff = nclosed + 1
      aend = nclosed + nact
      nca = nclosed + nact
c     
c     
      xoff = nvir*(nclosed+nact) + 1
      xend = vlen
      call ga_get(g_x,xoff,xend,1,1,x,(nclosed*nact))
      call dcopy((nclosed*nact),x,1,y,1)

      xoff = 1
      xend = nvir*nclosed
      call ga_get(g_x,xoff,xend,1,1,x,(nclosed*nvir))
      ii = nact*nclosed
      do i=1,nvir
        do j=1,nclosed
          y(ii+(i-1)*nca+j) = x((j-1)*nvir+i)
        enddo
      enddo

      xoff = xend + 1
      xend = nvir*(nclosed+nact)
      call ga_get(g_x,xoff,xend,1,1,x,(nact*nvir))
      do i=1,nvir
        do j=1,nact
          y(ii+(i-1)*nca+nclosed+j) = x((j-1)*nvir+i)
        enddo
      enddo

      call dscal(vlen,0.5d0,y,1)
c     
c     
      return
      end



c$$$
c$$$       subroutine mattr(n,m,x,y)
c$$$       implicit none
c$$$       integer n,m
c$$$       double precision x(n,m),y(m,n)
c$$$       integer i,j
c$$$
c$$$       do i=1,n
c$$$         do j=1,m
c$$$           y(j,i) = x(i,j)
c$$$         enddo
c$$$       enddo
c$$$       return
c$$$       end
c$$$




       subroutine hv_writev(lun,len,x)
       implicit none
       integer lun, len
       double precision x(len)
       
       write(lun) x
       return
       end



       subroutine hv_readv(lun,len,x)
       implicit none
       integer lun, len
       double precision x(len)
       
       read(lun) x
       return
       end




       subroutine ga_rowprint( text, g_a )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
       character*(*) text
       integer g_a
       integer gtype,vlen,rlen
       integer l_tmp, k_tmp
       integer i
       
       if (ga_nodeid().eq.0) then
         call ga_inquire(g_a,gtype,vlen,rlen)
         if (.not.ma_push_get(MT_DBL,vlen,'tmp',l_tmp,k_tmp))
     $     call errquit('ga_rowprint: cannot allocate',0)
         call ga_get(g_a,1,vlen,1,1,dbl_mb(k_tmp),1)
         write(6,909) text
 909     format(/,a)
         write(6,910) (dbl_mb(k_tmp+i-1),i=1,vlen)
 910     format(8F16.12)
         if (.not.ma_pop_stack(l_tmp))
     $     call errquit('ga_rowprint: cannot pop stack',0)
       endif
       call ga_sync()
       return
       end





      subroutine mcscf_twopdm_print(n,dm2)
      implicit none
      integer n
      double precision dm2(n,n,n,n)
      double precision TOL
      integer i,j,k,l,ltop
      data TOL/1.d-8/
      
      do i=1,n
        do j=1,i
          do k=1,i
            ltop = k
            if (k.gt.j) ltop = j
            do l=1,ltop
              if (abs(dm2(i,j,k,l)).gt.TOL) then
                write(6,771) i,j,k,l,dm2(i,j,k,l)
 771            format(4i5,f18.10)
              endif
            enddo
          enddo
        enddo
      enddo
      return
      end
