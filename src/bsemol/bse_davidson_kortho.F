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
      m = n
      k = 1
  100 continue
      do i=k,m
        if (eig(i).lt.1d-9) then
          do j=i+1,n
            eig(j-1) = eig(j)
            vkv(:,j-1) = vkv(:,j)
          enddo
          m = m - 1
          k = i
          goto 100
        endif
      enddo

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
