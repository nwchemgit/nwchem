c      SUBROUTINE tce_mrcc_c1_offset(size_t1m,nref,k_t1_offsetm,d_t1m)
      SUBROUTINE tce_mrcc_c1_offset(k_t1_offsetm,d_t1m,
     1 iref,irefnew,d_c)
ctce_mrcc_c1_offset(k_t1_offsetm,d_t1m,iref,irefnew,d_t1c)
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
      double precision ga_dble
      double precision ma_dble
      integer d_off1m(maxref)
      integer nodezero
      integer size,sizenew
      character*255 modelname
      character*255 filename
      character*3 namechunk
      integer iref
      integer size_t1m(maxref)
      integer k_t1_offsetm(maxref)
      integer l_c1,k_c1
      integer p5b,h6b
      integer mems
      integer i,j,k
      integer orbindex(2)
      integer orbindexnew(2)
      integer orbspin(2)
      integer irefnew,hnew,pnew
      integer hoff,poff
      integer totaloff
      integer ihash,ilength
      integer ioff
      integer l_t1,k_t1
      integer l_tmp,k_tmp
      integer d_t1m(maxref)
      integer counter
      integer k_a,l_a
      integer d_c
      integer inoabn
      integer x,y,z
      character*1  s,r
      INTEGER NXTASK
      EXTERNAL NXTASK
      INTEGER NXTASKsub
      EXTERNAL NXTASKsub
      integer next,count
      integer t
      double precision sign
      integer isw1,isw2

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
c
c ------------------
c  allocate arrays
c ------------------
c
c       if(nodezero) then
c          write(6,*)'Start of tce_mrcc_c1_offset'
c       endif
c print input arrays
c       do i=1,nref*nref
c         write(6,*)dbl_mb(k_sqc+i-1),dbl_mb(k_heff+i-1)
c       enddo
c create file

c        do iref=1,nref
c           write(namechunk,"(I3.3)")iref
c           call tce_filename('off1'//namechunk,filename)
c           call createfile(filename,d_off1m(iref),size_t1m(iref))
c           call reconcilefile(d_off1m(iref),size_t1m(iref))
c        enddo
c working arrays

c        iref = 1
c        mems = 0

          do p5b = nblcks(1,iref)+nblcks(2,iref)+1,nblcks(1,iref)+
     1 nblcks(2,iref)+nblcks(3,iref)+nblcks(4,iref)
          do h6b = 1,nblcks(1,iref)+nblcks(2,iref)

      IF (next.eq.count) THEN
         
      if (int_mb(k_spinm(iref)+p5b-1) .eq.
     1 int_mb(k_spinm(iref)+h6b-1)) then
      if (ieor(int_mb(k_symm(iref)+p5b-1),int_mb(k_symm(iref)+h6b-1))
     1 .eq. irrep_t) then
      if ((.not.restricted).or.(int_mb(k_spinm(iref)+p5b-1)+
     1 int_mb(k_spinm(iref)+h6b-1).ne.4)) THEN
 
c        write(6,"('Block assigned',I4,I4,' TO ',I3)")
c     1p5b,h6b,ga_nodeid()

         size = int_mb(k_rangem(iref)+p5b-1) * 
     1 int_mb(k_rangem(iref)+h6b-1)

c              if (.not.ma_push_get(mt_dbl,size,'t1',l_t1,k_t1))
c     1          call errquit('tce_c1_offs: MA problem',0,MA_ERR)

              if (.not.ma_push_get(mt_dbl,size,'c1',l_a,k_a))
     1          call errquit('tce_c1_offs: MA problem',1,MA_ERR)

c              call get_hash_block(d_t1m(iref),dbl_mb(k_t1),size,
c     1          int_mb(k_t1_offsetm(iref)),((p5b-noab-1)*noab+h6b-1))

         counter = 0         

         do i=1,int_mb(k_rangem(iref)+p5b-1)
           do j=1,int_mb(k_rangem(iref)+h6b-1)
           orbspin(2) = int_mb(k_spinm(iref)+h6b-1)-1
           orbspin(1) = int_mb(k_spinm(iref)+p5b-1)-1

             dbl_mb(k_a+counter) = 0.0d0
             counter = counter + 1
             
c             do irefnew=1,nref
c             if(irefnew.ne.iref) then

        isw1 = int_mb(k_offsetm(iref)+h6b-1)+j

             orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p5b-1)+i-1))/2
             orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h6b-1)+j-1))/2
c

        orbindexnew(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
        orbindexnew(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
 
        orbindexnew(1) = moindexes(orbindexnew(1),orbspin(1)+1,irefnew)
        orbindexnew(2) = moindexes(orbindexnew(2),orbspin(2)+1,irefnew)
        
c         orbindexnew(1) = orbindex(1)
c         orbindexnew(2) = orbindex(2)

c        if((orbindex(1).ne.orbindexnew(1)).or.
c     1     (orbindex(2).ne.orbindexnew(2))) then
c        write(6,"(I2,'/',I2,'(',I5,',',I5,')','->','(',I5,',',I5,')')")
c     1iref,irefnew,orbindex(1),orbindex(2),orbindexnew(1),orbindexnew(2)

        inoabn = nblcks(1,irefnew)+nblcks(2,irefnew)

        hnew = orbinblck(orbindexnew(2),orbspin(2)+1,irefnew)        
        pnew = orbinblck(orbindexnew(1),orbspin(1)+1,irefnew)

          totaloff=-1

        hoff = offsetinblck(orbindexnew(2),orbspin(2)+1,irefnew)
        poff = offsetinblck(orbindexnew(1),orbspin(1)+1,irefnew)

        ihash = hnew - 1 + inoabn * (pnew - inoabn - 1)
        ilength = int_mb(k_t1_offsetm(irefnew))
        ioff = 0
        totaloff = -1

        do k = 1, ilength
          if(int_mb(k_t1_offsetm(irefnew)+k).eq.ihash) then
            totaloff = ioff
            goto 111
          endif
c          ioff = int_mb(k_t1_offsetm(irefnew)+k+ilength)
        enddo       

 111    continue

        if((pnew.le.inoabn).or.
     1  (hnew.gt.inoabn)) then
          totaloff=-1
        endif

        if(totaloff.ne.-1) then
         ioff = offsetinblck(orbindexnew(1),orbspin(1)+1,irefnew)*
     1 int_mb(k_rangem(irefnew)+hnew-1)+
     2 offsetinblck(orbindexnew(2),orbspin(2)+1,irefnew)

       isw2 = int_mb(k_offsetm(irefnew)+hnew-1)+hoff

c          if(orbspin(1).eq.0) then
c            s='a'
c          else
c            s='b'
c          endif
c          if(orbspin(2).eq.0) then
c            r='a'
c          else
c            r='b'
c          endif

c          write(6,"(I4,I4)")iref,irefnew
c          write(6,"('[',I4,I4,']','(',I4,A1,I4,A1,')-->',
c     1'(',I4,I4,')')")
c     1 p5b,h6b,
c     1 orbindex(1),s,orbindex(2),r,orbindexnew(1),orbindexnew(2)
c      call util_flush(6)

         sizenew = int_mb(k_rangem(irefnew)+pnew-1) *
     1 int_mb(k_rangem(irefnew)+hnew-1)

              if (.not.ma_push_get(mt_dbl,sizenew,'tmp',l_tmp,k_tmp))
     1          call errquit('tce_c1_offs: MA problem',0,MA_ERR)
              call get_hash_block(d_t1m(irefnew),dbl_mb(k_tmp),sizenew,
     1          int_mb(k_t1_offsetm(irefnew)),
     1 ((pnew-inoabn-1)*inoabn+hnew-1))
c             write(6,*)(-dbl_mb(k_t1+counter-1)+dbl_mb(k_tmp+ioff))
c k_heff is global, c(iref) missing! 
          sign=1.0d0
c          if(mod((mod(isw1,2)+mod(isw2,2)),2).ne.0)sign=-1.0d0

              dbl_mb(k_a+counter-1)=dbl_mb(k_tmp+ioff)*sign
c              write(6,"('Counter ',I4,' pnew/hnew',I4,I4,2F16.12)")
c     1 counter,pnew,hnew,dbl_mb(k_tmp+ioff)
c      call util_flush(6)

c*dbl_mb(k_heff+irefnew+(irefnew-1)*nref)*
cdbl_mb(k_sqc+irefnew+(irefnew-1)*nref)

              if (.not.ma_pop_stack(l_tmp))
     1          call errquit('tce_c1_offs: MA problem',1,MA_ERR)
        endif
 
c         endif !nonzero C1

c           endif
c           enddo !irefnew
         enddo
         enddo

         call put_hash_block(d_c,dbl_mb(k_a),size,
     1   int_mb(k_t1_offsetm(iref)),((p5b-nblcks(1,iref)-nblcks(2,iref)
     2 -1)*(nblcks(1,iref)+nblcks(2,iref))+h6b-1))

c             call ma_print(dbl_mb(k_a),size,1,'C1')
             if (.not.ma_pop_stack(l_a))
     1          call errquit('tce_c1_offs: MA problem',2,MA_ERR)

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

      if(lusesub) then
       next = NXTASKsub(-nprocs,1,mypgid)
       call GA_pgroup_SYNC(mypgid)
      else
       next = NXTASK(-nprocs, 1)
       call GA_SYNC()
      endif

c        write(6,"('Doubles:',I8,' Bytes:',I8)")mems,mems*8

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

      RETURN
      END
c $Id$
