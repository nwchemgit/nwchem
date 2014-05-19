      block data lazymrcc_data
#include "tce_mrcc.fh"
      data l_tmp_on /.false./
      end

      SUBROUTINE tce_mrcc_c2_offset(k_t2_offsetm,d_t2m,
     1 iref,irefnew,d_c,ijblock)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "util.fh"
#include "msgids.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_mrcc.fh"

      integer nprocs
      logical nodezero
      integer isize
      integer size,sizenew
      integer iref
      integer k_t2_offsetm(maxref)
      integer l_c2,k_c2
      integer p1b,p2b,h3b,h4b
      integer mems
      integer i,j,k
      integer m,n
      integer orbindex(4)
      integer orbindexnew(4)
      integer orbspin(4)
      integer irefnew,h3new,p1new
      integer p2new,h4new
      integer h3off,p1off
      integer p2off,h4off
      integer totaloff
      integer ihash,ilength
      integer inoab,invab
      integer inoabn,invabn
      integer ioff
      integer l_t2,k_t2
      integer l_tmp,k_tmp
      integer d_t2m(maxref)
      integer counter
      integer k_a,l_a
      integer d_c
      integer t
      double precision sign
      INTEGER NXTASK
      EXTERNAL NXTASK
      INTEGER NXTASKsub
      EXTERNAL NXTASKsub
      integer next,count
      integer hs,ps
c      integer nref
      integer ijblock(nref,noab+nvab)
      integer ihashold
      integer newp1,newp2,newh3,newh4
      integer iswap
      integer isw1,isw2
c
c     keeping compilers quiet
c
      p1off = 0
      p2off = 0
      h3off = 0
      h4off = 0
c
      if(lusesub) then

      call ga_pgroup_sync(mypgid)
      nprocs = GA_pgroup_NNODES(mypgid)
      count = 0
      next = NXTASKsub(nprocs, 1,mypgid)

      else

      call ga_sync()
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)

      endif

      nodezero=(ga_nodeid().eq.0)

      inoab = nblcks(1,iref)+nblcks(2,iref)
      invab = nblcks(3,iref)+nblcks(4,iref)

      DO p1b = inoab+1,inoab+invab
      DO p2b = p1b,inoab+invab
      DO h3b = 1,inoab
      DO h4b = h3b,inoab

      IF (next.eq.count) THEN

      IF (int_mb(k_spinm(iref)+p1b-1)+int_mb(k_spinm(iref)+p2b-1)
     1.eq.int_mb(k_spinm(iref)+h3b-1)+int_mb(k_spinm(iref)+h4b-1))THEN
      IF (ieor(int_mb(k_symm(iref)+p1b-1),ieor(int_mb(k_symm(iref)+p2b
     1-1),ieor(int_mb(k_symm(iref)+h3b-1),int_mb(k_symm(iref)+h4b-1))))
     1  .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spinm(iref)+p1b-1)+
     1 int_mb(k_spinm(iref)+p2b-1)+int_mb(k_spinm(iref)+h3b-1)
     1 +int_mb(k_spinm(iref)+h4b-1).ne.8)) THEN

         size = int_mb(k_rangem(iref)+p1b-1) *
     1 int_mb(k_rangem(iref)+p2b-1)*int_mb(k_rangem(iref)+h3b-1)*
     2 int_mb(k_rangem(iref)+h4b-1)


              if (.not.ma_push_get(mt_dbl,size,'c2',l_a,k_a))
     1          call errquit('tce_c2_offs: MA problem',11,MA_ERR)
c              write(6,"('C_a',i4)")ga_nodeid()

c         write(6,"('Working on ',I4,I4,I4,I4,I10)")
c     1 p1b,p2b,h3b,h4b,size 
c         call util_flush(6)

c         if(log_mb(k_isactive(iref)+p1b-1).or.
c     1       log_mb(k_isactive(iref)+p2b-1).or.
c     1       log_mb(k_isactive(iref)+h3b-1).or.
c     1       log_mb(k_isactive(iref)+h4b-1)) then
            if(1.eq.1) then

         counter = 0
         ihash = -1
         ihashold = -1
         CALL DFILL(size,0.0d0,dbl_mb(k_a),1)

         do i=1,int_mb(k_rangem(iref)+p1b-1)
           do j=1,int_mb(k_rangem(iref)+p2b-1)
             do m=1,int_mb(k_rangem(iref)+h3b-1)
              do n=1,int_mb(k_rangem(iref)+h4b-1)
               orbspin(4) = int_mb(k_spinm(iref)+h4b-1)-1
               orbspin(3) = int_mb(k_spinm(iref)+h3b-1)-1
               orbspin(2) = int_mb(k_spinm(iref)+p2b-1)-1
               orbspin(1) = int_mb(k_spinm(iref)+p1b-1)-1

c             if(nodezero)
c       write(6,"('SWAPS mu: ',I5)")
        isw1 = (int_mb(k_offsetm(iref)+h3b-1)+m
     2 +int_mb(k_offsetm(iref)+h4b-1)+n)

             dbl_mb(k_a+counter) = 0.0d0
             counter = counter + 1
             
             orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+p1b-1)+i-1))/2
             orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+p2b-1)+j-1))/2
             orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+h3b-1)+m-1))/2
             orbindex(4) = (1 - orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+h4b-1)+n-1))/2

        orbindexnew(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
        orbindexnew(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
        orbindexnew(3) = moindexes(orbindex(3),orbspin(3)+1,iref)
        orbindexnew(4) = moindexes(orbindex(4),orbspin(4)+1,iref)
       

        orbindexnew(1) = moindexes(orbindexnew(1),orbspin(1)+1,irefnew)
        orbindexnew(2) = moindexes(orbindexnew(2),orbspin(2)+1,irefnew)
        orbindexnew(3) = moindexes(orbindexnew(3),orbspin(3)+1,irefnew)
        orbindexnew(4) = moindexes(orbindexnew(4),orbspin(4)+1,irefnew)

cjb ?????????????????????????
c         orbindexnew(1) = orbindex(1)
c         orbindexnew(2) = orbindex(2)
c         orbindexnew(3) = orbindex(3)
c         orbindexnew(4) = orbindex(4)

c        write(6,"(I4,I4,I4,I4,'--->',I4,I4,I4,I4)")
c     1 orbindex(1),orbindex(2),orbindex(3),orbindex(4),
c     2 orbindexnew(1),orbindexnew(2),orbindexnew(3),
c     3 orbindexnew(4)

        inoabn = nblcks(1,irefnew)+nblcks(2,irefnew)
        invabn = nblcks(3,irefnew)+nblcks(4,irefnew)

        p1new = orbinblck(orbindexnew(1),orbspin(1)+1,irefnew)        
        p2new = orbinblck(orbindexnew(2),orbspin(2)+1,irefnew)
        h3new = orbinblck(orbindexnew(3),orbspin(3)+1,irefnew)
        h4new = orbinblck(orbindexnew(4),orbspin(4)+1,irefnew)
c       write(6,"('OLD,..:',I4,I4,I4,I4)")p1b,p2b,h3b,h4b
c       write(6,"('NEW,..:',I4,I4,I4,I4)")p1new,p2new,h3new,h4new

        totaloff = -1
        hs = 0

        if(p1new.le.inoabn)hs=hs+1
        if(p2new.le.inoabn)hs=hs+1
        if(h3new.le.inoabn)hs=hs+1
        if(h4new.le.inoabn)hs=hs+1

        sign = 1.0d0

        if(hs.ne.2)goto 111  

        iswap = 0

        if(p1new.gt.p2new) then
         t = p1new
         p1new = p2new
         p2new = t
         t = orbindexnew(1)
         orbindexnew(1) = orbindexnew(2)
         orbindexnew(2) = t
         t = orbspin(1)
         orbspin(1) = orbspin(2)
         orbspin(2) = t
         iswap = iswap + 1
         sign = sign * (-1.0d0)
        endif

        if(h3new.gt.h4new) then
         t = h3new
         h3new = h4new
         h4new = t
         t = orbindexnew(3)
         orbindexnew(3) = orbindexnew(4)
         orbindexnew(4) = t
         t = orbspin(3)
         orbspin(3) = orbspin(4)
         orbspin(4) = t
         iswap = iswap + 1
         sign = sign * (-1.0d0)
        endif

        p1off = offsetinblck(orbindexnew(1),orbspin(1)+1,irefnew)
        p2off = offsetinblck(orbindexnew(2),orbspin(2)+1,irefnew)
        h3off = offsetinblck(orbindexnew(3),orbspin(3)+1,irefnew)
        h4off = offsetinblck(orbindexnew(4),orbspin(4)+1,irefnew)

        ihash = h4new - 1 + inoabn * (h3new - 1 + inoabn * (p2new-
     &inoabn - 1 + invabn * (p1new - inoabn - 1))) 

        ilength = int_mb(k_t2_offsetm(irefnew))
        ioff = 0
        totaloff = -1

        do k = 1, ilength
          if(int_mb(k_t2_offsetm(irefnew)+k).eq.ihash) then
c            totaloff = int_mb(k_t2_offsetm(irefnew)+k+ilength)
             totaloff = 1
c             if((sign.lt.0.0d0).and.nodezero)write(6,*)'SIGN'
            goto 111
          endif
c          ioff = int_mb(k_t2_offsetm(irefnew)+k+ilength)
        enddo       

 111    continue

        if((p1new.le.inoabn).or.
     1     (p2new.le.inoabn).or.
     2     (h3new.gt.inoabn).or.
     3     (h4new.gt.inoabn)) then
c           write(6,"('BLOCK:',I4,I4,I4,I4)")p1new,p2new,
c     1 h3new,h4new
           totaloff=-1
        endif

        if((totaloff.ne.-1)) then
         ioff = p1off*int_mb(k_rangem(irefnew)+h4new-1)*
     1 int_mb(k_rangem(irefnew)+h3new-1)*
     2 int_mb(k_rangem(irefnew)+p2new-1)+
     2 p2off*int_mb(k_rangem(irefnew)+h4new-1)*
     3 int_mb(k_rangem(irefnew)+h3new-1)+
     4 h3off*int_mb(k_rangem(irefnew)+h4new-1)+h4off

       isw2 = (
     1 int_mb(k_offsetm(irefnew)+h3new-1)+h3off
     2 +int_mb(k_offsetm(irefnew)+h4new-1)+h4off)

c       write(6,"('SWAPS : ',I5,I5,I2)")isw1,isw2,
c     1 mod((mod(isw1,2)+mod(isw2,2)),2)

c          write(6,*)ioff,totaloff

         sizenew = int_mb(k_rangem(irefnew)+p1new-1) *
     1 int_mb(k_rangem(irefnew)+p2new-1)*
     1 int_mb(k_rangem(irefnew)+h3new-1)*
     1 int_mb(k_rangem(irefnew)+h4new-1)

          if((l_tmp_on).and.(ihash.ne.ihashold)) then
             if (.not.ma_pop_stack(l_tmp))
     1            call errquit('tce_c2_offs: MA problem',1,MA_ERR)
             l_tmp_on=.false.
          endif

          if(.not.l_tmp_on) then

          if (.not.ma_push_get(mt_dbl,sizenew,'tmp',l_tmp,k_tmp))
     1       call errquit('tce_c2_offs: MA problem',0,MA_ERR)
c              write(6,"('C_t',i4)")ga_nodeid()
          l_tmp_on=.true.

          call get_hash_block(d_t2m(irefnew),dbl_mb(k_tmp),sizenew,
     1    int_mb(k_t2_offsetm(irefnew)),ihash)
 
          endif

c          if(mod((mod(isw1,2)+mod(isw2,2)),2).ne.0)sign=-sign
c          write(6,*)sign
              dbl_mb(k_a+counter-1)=dbl_mb(k_tmp+ioff)*sign

          ihashold = ihash

        endif
 
c         endif !nonzero C1

c           endif
c           enddo !irefnew
         enddo
         enddo
         enddo
         enddo

         if(l_tmp_on) then
            if (.not.ma_pop_stack(l_tmp))
     1           call errquit('tce_c2_offs: 1MA problem',1,MA_ERR)
            l_tmp_on = .false.
         endif

         else

         inoabn = nblcks(1,irefnew)+nblcks(2,irefnew)
         invabn = nblcks(3,irefnew)+nblcks(4,irefnew)

        newp1 = ijblock(irefnew,p1b)
        newp2 = ijblock(irefnew,p2b)
        newh3 = ijblock(irefnew,h3b)
        newh4 = ijblock(irefnew,h4b)

c        if(newp1.gt.newp2) then
c          t = newp1
c          newp1 = newp2
c          newp2 = t
c        endif
                
c        if(newh3.gt.newh4) then
c          t = newh3
c          newh3 = newh4
c          newh4 = t
c        endif

        ihash = newh4 - 1 + inoabn *
     1  (newh3 - 1 + inoabn *
     1  (newp2 - inoabn - 1 + invabn *
     1  (newp1 - inoabn - 1)))

         CALL DFILL(size,0.0d0,dbl_mb(k_a),1) 

         call get_hash_block(d_t2m(irefnew),dbl_mb(k_a),size,
     1   int_mb(k_t2_offsetm(irefnew)),ihash)

         endif

         call put_hash_block(d_c,dbl_mb(k_a),size,
     1   int_mb(k_t2_offsetm(iref)),(h4b - 1 +
     2  inoab * (h3b - 1 + inoab * (p2b-
     &inoab - 1 + invab * (p1b - inoab - 1)))))

c             call ma_print(dbl_mb(k_a),size,1,'C1')
c              write(6,"('D_a',i4)")ga_nodeid()
             if (.not.ma_pop_stack(l_a))
     Q        call errquit('tce_c2_offs: MA problem',2,MA_ERR)

c         write(6,"('Finished: ',I4,I4,I4,I4)")
c     1 p1b,p2b,h3b,h4b
c         call util_flush(6)

c             if (.not.ma_pop_stack(l_t1))
c     1          call errquit('tce_c1_offs: MA problem',1,MA_ERR)

      endif
      endif
      endif
      if(lusesub) then
       next = NXTASKsub(nprocs,1,mypgid)
      else
       next = NXTASK(nprocs, 1)
      endif
      END IF
      count = count + 1

          enddo ! h6b
          enddo ! p5b
          enddo ! p5b
          enddo ! p5b

c computing offsets, offset files size equal size of amplitude files

c         do
c           do
              
c           enddo
c         enddo 

c purge memory

c         if (.not.ma_pop_stack(l_c1))
c     1       call errquit('tce_mrcc_c1: MA problem',1,MA_ERR)

c delete file
c        do iref=1,nref
c           call deletefile(d_off1m(iref))
c        enddo

c       if(nodezero) then
c          write(6,*)'End of procedure tce_mrcc_c1_offset'
c       endif
c
c --------------------
c  deallocate arrays
c --------------------
c

c deleted

      if(lusesub) then
       next = NXTASKsub(-nprocs,1,mypgid)
       call GA_pgroup_SYNC(mypgid)
      else
       next = NXTASK(-nprocs, 1)
       call GA_SYNC()
      endif
      l_tmp_on=.false.
      RETURN
      END
c $Id$
