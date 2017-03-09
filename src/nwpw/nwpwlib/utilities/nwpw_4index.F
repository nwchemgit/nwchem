
*     ******************************************
*     *                                        *
*     *          nwpw_4index_init              *
*     *                                        *
*     ******************************************
      subroutine nwpw_4index_init(ispin0,nocc0,nvirt0)
      implicit none
      integer ispin0,nocc0(*),nvirt0(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_4indx common block ****
      integer ispin,nocc(2),nvirt(2),maxrs(2),A4indx(2)
      common /nwpw_4indx/ nocc,nvirt,maxrs,A4indx,ispin

*     **** local variables ****
      integer ms,np,taskid,mysize

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      ispin  = ispin0
      mysize = 0
      do ms=1,ispin
         nocc(ms)  = nocc0(ms)
         nvirt(ms) = nvirt0(ms)
         maxrs(ms) = Nvirt(ms)*Nvirt(ms)/np
         if (mod((Nvirt(ms)*Nvirt(ms)),np).ge.taskid) 
     >      maxrs(ms) = maxrs(ms) + 1
         mysize = mysize + maxrs(ms)*nocc(ms)**2 
      end do

      mysize = 2*mysize
      if (.not.BA_alloc_get(mt_dbl,mysize,'A4indx',A4indx(2),A4indx(1)))
         call errquit('nwpw_4indx_init:get heap',0,MA_ERR)

      return
      end 


*     ******************************************
*     *                                        *
*     *          nwpw_4index_end               *
*     *                                        *
*     ******************************************
      subroutine nwpw_4index_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_4indx common block ****
      integer ispin,nocc(2),nvirt(2),maxrs(2),A4indx(2)
      common /nwpw_4indx/ nocc,nvirt,maxrs,A4indx,ispin

      if (.not.BA_free_heap(A4indx(2)))
         call errquit('nwpw_4indx_end:get heap',0,MA_ERR)

      return
      end

*     ******************************************
*     *                                        *
*     *          nwpw_4index_set               *
*     *                                        *
*     ******************************************
      subroutine nwpw_4index_set(ms,a,b,r,s,v,denom)
      implicit none
      integer ms,a,b,r,s
      real*8 v,denom

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_4indx common block ****
      integer ispin,nocc(2),nvirt(2),maxrs(2),A4indx(2)
      common /nwpw_4indx/ nocc,nvirt,maxrs,A4indx,ispin

*     **** local variables ****
      integer taskid,np

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      abindx = (a-1)+(b-1)*nocc(ms)
      rsindx = (r-1)+(s-1)*nvirt(ms)

      if (mod(rsindx,np).eq.taskid) then
         rs = (rsindx/np)
         indx = 2*abindx + 2*rs*nocc(ms)*nocc(ms)
         dbl_mb(A4indx(1)+indx)    = v
         dbl_mb(A4indx(1)+indx+1)  = v/denom
      end if

      return
      end

*     ******************************************
*     *                                        *
*     *          nwpw_4index_emp2              *
*     *                                        *
*     ******************************************
      real*8 function nwpw_4index_emp2()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_4indx common block ****
      integer ispin,nocc(2),nvirt(2),maxrs(2),A4indx(2)
      common /nwpw_4indx/ nocc,nvirt,maxrs,A4indx,ispin

*     **** local variables ****
      integer a,b,rs,abindx,abindx1,indx,indx1
      real*8  emp2

*     **** external functions ****
      real*8   nwpw_4index_emp2_sub
      external nwpw_4index_emp2_sub

      emp2 = 0.0d0
      if (ispin.eq.1) then
         emp2 = nwpw_4indx_emp2_sub(nocc(1),maxrs(1),dbl_mb(A4indx(1)))
      end if 
      call Parallel_SumAll(emp2)

      nwpw_4indx_emp2 = emp2
      return
      end

      real*8 function nwpw_4index_emp2_sub(nocc,maxrs,A4indx)
      implicit none
      integer nocc,maxrs
      real*8 A4indx(2,nocc,nocc,*)

*     **** local variables ****
      integer a,b,rs
      real*8 emp2

      emp2 = 0.0
      do rs = 1,maxrs
         do b=1,nocc
         do a=1,nocc
            emp2 = emp2 
     >           + A4indx(2,a,b,rs)*(2*A4inx(1,a,b,rs)-A4indx(1,b,a,rs))
         end do
         end do
      end do

      nwpw_4index_emp2_sub = emp2
      return
      end 









