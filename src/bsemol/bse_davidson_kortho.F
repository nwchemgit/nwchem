      subroutine bse_davidson_kortho(pars,v,kv,mkv,vkv,work,npoles,n,m,
     $                               first)

      implicit none
#include "bse.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      type(bse_params_t) :: pars
      logical first
      integer npoles,n,m

      integer l_work,k_work,info,i,j,k

      double precision rdum
      double precision v(npoles,n),kv(npoles,n),mkv(npoles,n)
      double precision vkv(n,n),eig(n),work(npoles,n)
      double precision sqrtvkv(n,n)


      ! V**T  K  V
      call dgemm('t','n',n,n,npoles,1d0,v,npoles,kv,npoles,0d0,vkv,n)
      call ga_dgop((/14/),vkv,n*n,'+')

      ! Allocate work field
      call dsyev('v','l',n,vkv,n,eig,rdum,-1,info)
      if(.not.ma_push_get(mt_dbl,int(rdum),'work',l_work,k_work))
     $  call errquit('bse_kortho: failed to allocate work',199,MA_ERR)

      ! Diagonalize matrix
      call dsyev('v','l',n,vkv,n,eig,dbl_mb(k_work),int(rdum),info)
      if (info.ne.0)
     $  call errquit("bse_davidson_kortho: diago failed",0,0)

      ! Deallocate field
      if (.not.ma_chop_stack(l_work))
     $  call errquit('bse_kortho: failed to deallocate work',199,MA_ERR)

      ! Remove small eigenpairs
      m = 0
      do i=n,1,-1
        if(eig(i).lt.1d-10) exit
        m = m + 1
      enddo

      if(m.ne.n) then
        j = 1
        do i=n-m+1,n
          eig(j) = eig(i)
          vkv(:,j) = vkv(:,i)
          j = j + 1
        enddo
      endif  
        
c      m = n
c      k = 1
c  100 continue
c      do i=k,m
c        if (eig(i).lt.1d-9) then
c          do j=i+1,n
c            eig(j-1) = eig(j)
c            vkv(:,j-1) = vkv(:,j)
c          enddo
c          m = m - 1
c          k = i
c          goto 100
c        endif
c      enddo

      ! Obtain VKV^(-1/2)
      do i=1,m
        vkv(:,i) = vkv(:,i)/dsqrt(eig(i))
      enddo

      if (m.eq.0) return

      ! Obtain orthogonalized vectors
      call dgemm('n','n',npoles,m,n,1d0,v,npoles,vkv,n,0d0,work,npoles)
      call dcopy(npoles*m,work,1,v,1)

      if (first) return

      call dgemm('n','n',npoles,m,n,1d0,kv,npoles,vkv,n,0d0,work,npoles)
      call dcopy(npoles*m,work,1,kv,1)

      call dgemm('n','n',npoles,m,n,1d0,mkv,npoles,vkv,n,0d0,work,
     $            npoles)
      call dcopy(npoles*m,work,1,mkv,1)

      end subroutine
