c
c ==============================================
c     Create effective Hamiltonian
c ==============================================
c
       subroutine tce_heff(d_em,k_e_offsetm,k_r1_offsetm,
     1 k_r2_offsetm,k_r3_offsetm,k_r4_offsetm,d_r1m,d_r2m,d_r3m,d_r4m,
     2 needt1,needt2,needr3act,needr4act,rtdb)
        implicit none
#include "tce.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "sym.fh"
#include "tce_mrcc.fh"
#include "global.fh"
#include "tce_main.fh"

       integer d_em(maxref)
       integer k_e_offsetm(maxref)
       integer iref,i
       double precision corr
       logical nodezero
       integer k_r1_offsetm(maxref)
       integer k_r2_offsetm(maxref)
       integer k_r3_offsetm(maxref)
       integer k_r4_offsetm(maxref)
       integer d_r1m(maxref),d_r2m(maxref)
       integer d_r3m(maxref),d_r4m(maxref)
       logical needt1,needt2,needr3act
       logical needr4act
       integer rtdb

       nodezero = (ga_nodeid().eq.0)

       if(lusesub)call ga_zero(g_heff)

       do i=1,nref*nref
         dbl_mb(k_heff+i-1) = 0.0d0
       enddo

       do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +ga_nnodes()+ga_nodeid())).or.(.not.lusesub)) then

          call get_block(d_em(iref),corr,1,0)
          dbl_mb(k_heff+iref-1+(iref-1)*nref) = corr+duhfens(iref)

        if(lusesub) then

c          write(6,*)ga_nodeid(),corr+duhfens(iref),iref
          call ga_put(g_heff,nref*(iref-1)+iref,nref*(iref-1)+iref,1,1,
     1    corr+duhfens(iref),1)
        endif

        endif

       enddo

       call tce_heff_offdiagonal(k_r1_offsetm,k_r2_offsetm,
     1 k_r3_offsetm,k_r4_offsetm,d_r1m,d_r2m,d_r3m,d_r4m,needt1,
     2 needt2,needr3act,needr4act,rtdb) 

c       if(nodezero) then
c         call ma_print(dbl_mb(k_heff),nref,nref,'Heff')
c       endif

       return
       end
c
c ==============================================
c     Add offdiagonal elements of Heff
c ==============================================
c
       subroutine tce_heff_offdiagonal(k_r1_offsetm,
     1 k_r2_offsetm,k_r3_offsetm,k_r4_offsetm,d_r1m,d_r2m,d_r3m,
     2 d_r4m,needt1,needt2,needr3act,needr4act,rtdb)
        implicit none
#include "tce.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "sym.fh"
#include "tce_mrcc.fh"
#include "global.fh"
#include "tce_main.fh"

       integer rtdb
       logical nodezero
       integer k_r1_offsetm(maxref)
       integer k_r2_offsetm(maxref)
       integer k_r3_offsetm(maxref)
       integer k_r4_offsetm(maxref)
       integer d_r1m(maxref),d_r2m(maxref)
       integer d_r3m(maxref)
       integer d_r4m(maxref)
       integer iref
       integer i,j,p5b,h6b,k
       integer size,l,m,n,o
       integer l_r1,k_r1,l_r2,k_r2
       integer l_r3,k_r3
       integer l_r4,k_r4
       integer p1b,p2b,h3b,h4b
       integer p3b,p7b,p8b
       INTEGER p4b
       INTEGER p6b
       INTEGER h1b
       INTEGER h2b
       integer h1,h2,h3,p4,p5,p6
       integer h4,p7,p8
       integer orbindex(8)
       integer orbspin(8)
       integer ioccnew(maxorb,2)
       integer iu
       logical isfound
       logical needt1,needt2,needr3act
       logical needr4act
c       logical lusescffv
       integer is,iocc0(maxorb,2)
       integer i1,i2,dist
       double precision dsmult
c       logical limprovet
       EXTERNAL NXTASKsub
       EXTERNAL NXTASK
       INTEGER NXTASKsub
       INTEGER NXTASK
       INTEGER nxt
       INTEGER nprocs
       INTEGER count

       nodezero = (ga_nodeid().eq.0)
       isfound = .false.

c       if (.not.rtdb_get(rtdb,'mrcc:usescffermiv',mt_log,1,lusescffv))
c     1 lusescffv = .false.
c       if (.not.rtdb_get(rtdb,'mrcc:improvetiling',mt_log,1,limprovet))
c     1 limprovet = .false.

      do is=1,2
        do i=1,nmo(is)
          if(i.le.nocc(is)) then
            iocc0(i,is) = 1
          else
            iocc0(i,is) = 0
          end if
        end do
      end do

       do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +ga_nnodes()+ga_nodeid())).or.(.not.lusesub)) then

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
         k_active = k_active_tmpm(iref)

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb
c
c ---------------
c    R1 active
c ---------------
c
         do p5b = noab+1,noab+nvab
         do h6b = 1,noab

      k = 0

      if (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) then
      if (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)).eq.irrep_t)then
      if ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) then
      if(log_mb(k_isactive(iref)+p5b-1).and.
     &log_mb(k_isactive(iref)+h6b-1)) then

        size = int_mb(k_range+h6b-1) * int_mb(k_range+p5b-1)

        if (.not.ma_push_get(mt_dbl,size,'r1mi',l_r1,k_r1))
     1   call errquit('tce_mrcc_iface_r1: MA problem',0,MA_ERR)

        call get_hash_block(d_r1m(iref),dbl_mb(k_r1),size,
     1   int_mb(k_r1_offsetm(iref)),h6b-1+noab*(p5b-noab-1))

        k=0
        do i=1,int_mb(k_range+p5b-1)
        do j=1,int_mb(k_range+h6b-1)
          k = k + 1

        orbspin(1) = int_mb(k_spin+p5b-1)-1
        orbspin(2) = int_mb(k_spin+h6b-1)-1

        orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p5b-1)+i-1))/2
        orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h6b-1)+j-1))/2

        orbindex(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
        orbindex(2) = moindexes(orbindex(2),orbspin(2)+1,iref)

cjb ===========================================================================

        if(isactive(orbindex(1),orbspin(1)+1).and.
     1 isactive(orbindex(2),orbspin(2)+1).or.(.not.limprovet)) then

ccc =====
c       if(nodezero)write(6,"('ACTIVITY: ',2L2)")
c     1 isactive(orbindex(1),orbspin(1)+1),
c     1 isactive(orbindex(2),orbspin(2)+1)
c       if(nodezero)write(6,"('DEBUG: ',5I4)")
c     1 orbindex(1),orbspin(1),
c     1 orbindex(2),orbspin(2),iref
ccc =====

       do iu=1,2
        do n=1,nmo(iu)
          ioccnew(n,iu) = iocc(n,iref,iu)
        enddo
       enddo

      if(iocc(orbindex(1),iref,orbspin(1)+1).eq.
     1 iocc(orbindex(2),iref,orbspin(2)+1)) then
          goto 200
       endif

         ioccnew(orbindex(1),orbspin(1)+1) = iocc(orbindex(2),iref,
     1 orbspin(2)+1)
         ioccnew(orbindex(2),orbspin(2)+1) = iocc(orbindex(1),iref,
     1 orbspin(1)+1)

       do n=1,nref
       isfound = .true.
        do iu=1,2
         do o=1,nmo(iu)
          isfound = isfound.and.(iocc(o,n,iu).eq.ioccnew(o,iu))
         enddo
        enddo
       if(isfound) then
c          write(LuOut,"('Internal amplitude',I4,'->',I4)")iref,n
c        write(LuOut,"('1Internal amplitude',I4,'->',I4,2F16.12)")iref,n,
c     1 dbl_mb(k_r1+m-1)
            if(iref.ne.n) then

             dist = 0
       do iu=1,2
        do i1=1,nmo(iu)
          ioccnew(i1,iu) = iocc(i1,iref,iu)
        enddo
       enddo

             i2 = 0
             do i1=min(orbindex(1),orbindex(2)),
     1 max(orbindex(1),orbindex(2))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(1)-orbindex(2))) then
                  if(iocc(i1+1,iref,orbspin(1)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

         dsmult = 1.0d0

         if(mod(dist,2).eq.1) dsmult = -dsmult

         dbl_mb(k_heff+n-1+(iref-1)*nref) = dbl_mb(k_r1+k-1)*dsmult
c     1 dbl_mb(k_heff+n-1+(iref-1)*nref)

           if(lusesub) then
         call ga_put(g_heff,nref*(iref-1)+n,nref*(iref-1)+n,1,1,
     1   dbl_mb(k_r1+k-1)*dsmult,1)
           endif

            endif
          goto 200
       endif
       enddo

 200   continue

        if((.not.isfound).and.(abs(dbl_mb(k_r1+k-1)).gt.1.0d-5)) then
         if(nodezero) then
           write(LuOut,"('DEBUG: ',4F16.12)")dbl_mb(k_r1+k-1)
           write(LuOut,"('YOU ARE USING INCOMPLETE MODEL SPACE!')")
         endif
c         call errquit('YOU ARE USING INCOMPLETE MODEL SPACE!',1,MA_ERR)
        endif

        endif ! active

        enddo
        enddo

        if (.not.ma_pop_stack(l_r1))
     1   call errquit('tce_mrcc_iface_r1: MA problem',1,MA_ERR)

      endif
      endif
      endif
      endif
      enddo
      enddo
c
c ---------------
c    R2 active
c ---------------
c
      nxt = 0
      count = 0
      if(limprovet) then
       if(lusesub) then
         nprocs=GA_pgroup_NNODES(mypgid)
         nxt=NXTASKsub(nprocs,1,mypgid)
       else
         nprocs=GA_NNODES()
         nxt=NXTASK(nprocs,1)
       endif
      count = 0
      endif

      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO h3b = 1,noab
      DO h4b = h3b,noab

      IF ((nxt.eq.count).or.(.not.limprovet)) THEN

      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN

      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep_t) THEN

      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN

      if(log_mb(k_isactive(iref)+p1b-1).and.
     1 log_mb(k_isactive(iref)+p2b-1).and.
     2 log_mb(k_isactive(iref)+h3b-1).and.
     3 log_mb(k_isactive(iref)+h4b-1)) then

      size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)

        if (.not.ma_push_get(mt_dbl,size,'r2mi',l_r2,k_r2))
     1   call errquit('tce_mrcc_iface_r2: MA problem',0,MA_ERR)

        call get_hash_block(d_r2m(iref),dbl_mb(k_r2),size,
     1   int_mb(k_r2_offsetm(iref)),h4b-1+noab*(h3b-1+noab*(p2b-
     &noab-1+nvab*(p1b - noab - 1))))

c         write(LuOut,"(I4,L3,L3,L3,L3)")
c     1 iref,log_mb(k_isactive(iref)+p1b-1),
c     1 log_mb(k_isactive(iref)+p2b-1),log_mb(k_isactive(iref)+h3b-1),
c     1 log_mb(k_isactive(iref)+h4b-1)
       m = 0

        do i=1,int_mb(k_range+p1b-1)
        do j=1,int_mb(k_range+p2b-1)
        do k=1,int_mb(k_range+h3b-1)
        do l=1,int_mb(k_range+h4b-1)
         m = m + 1
c         write(LuOut,"(I4,'(',I4,I4,I4,I4,'):',2F16.12)")
c     1 iref,i,j,k,l,dbl_mb(k_r2+m-1)
c         write(LuOut,*)int_mb(k_spin+p1b-1)

        orbspin(1) = int_mb(k_spin+p1b-1)-1
        orbspin(2) = int_mb(k_spin+p2b-1)-1
        orbspin(3) = int_mb(k_spin+h3b-1)-1
        orbspin(4) = int_mb(k_spin+h4b-1)-1

        orbindex(1) = (1 - orbspin(1)+ 
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p1b-1)+i-1))/2
        orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p2b-1)+j-1))/2
        orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h3b-1)+k-1))/2
        orbindex(4) = (1 - orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h4b-1)+l-1))/2

        orbindex(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
        orbindex(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
        orbindex(3) = moindexes(orbindex(3),orbspin(3)+1,iref)
        orbindex(4) = moindexes(orbindex(4),orbspin(4)+1,iref)

cjb ===========================================================================

        if(isactive(orbindex(1),orbspin(1)+1).and.
     1 isactive(orbindex(2),orbspin(2)+1).and.
     2 isactive(orbindex(3),orbspin(3)+1).and.
     3 isactive(orbindex(4),orbspin(4)+1).or.(.not.limprovet)) then

c        write(LuOut,"('Real indexes: [',I4,I4,I4,I4,']',
c     1 '[',I2,I2,I2,I2,']')")
c     1 orbindex(1),orbindex(2),orbindex(3),orbindex(4),
c     1 orbspin(1),orbspin(2),orbspin(3),orbspin(4)

       do iu=1,2
        do n=1,nmo(iu)
          ioccnew(n,iu) = iocc(n,iref,iu) 
        enddo
       enddo 
        
       if(((orbindex(1).eq.orbindex(2)).and.(orbspin(1).eq.orbspin(2)))
     1 .or.((orbindex(3).eq.orbindex(4)).and.(orbspin(3).eq.
     2 orbspin(4)))) then
         goto 100
       endif

         ioccnew(orbindex(1),orbspin(1)+1) = iocc(orbindex(3),iref,
     1 orbspin(3)+1)
         ioccnew(orbindex(2),orbspin(2)+1) = iocc(orbindex(4),iref,
     1 orbspin(4)+1)
         ioccnew(orbindex(3),orbspin(3)+1) = iocc(orbindex(1),iref,
     1 orbspin(1)+1)
         ioccnew(orbindex(4),orbspin(4)+1) = iocc(orbindex(2),iref,
     1 orbspin(2)+1)

       do n=1,nref
       isfound = .true.
        do iu=1,2
         do o=1,nmo(iu)
          isfound = isfound.and.(iocc(o,n,iu).eq.ioccnew(o,iu))
         enddo
        enddo
       if(isfound) then
c        write(LuOut,"('2Internal amplitude',I4,'->',I4,2F16.12)")iref,n,
c     1 dbl_mb(k_r2+m-1)
            if(iref.ne.n) then

             dist = 0
       do iu=1,2
        do i1=1,nmo(iu)
          ioccnew(i1,iu) = iocc(i1,iref,iu)
        enddo
       enddo

             i2 = 0
             do i1=min(orbindex(1),orbindex(3)),
     1 max(orbindex(1),orbindex(3))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(1)-orbindex(3))) then
                  if(iocc(i1+1,iref,orbspin(1)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

         ioccnew(orbindex(1),orbspin(1)+1) = iocc(orbindex(3),iref,
     1 orbspin(3)+1)
         ioccnew(orbindex(3),orbspin(3)+1) = iocc(orbindex(1),iref,
     1 orbspin(1)+1)

             i2 = 0
             do i1=min(orbindex(2),orbindex(4)),
     1 max(orbindex(2),orbindex(4))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(2)-orbindex(4))) then
                  if(ioccnew(i1+1,orbspin(2)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

            dsmult = 1.0d0
         if(mod(dist,2).eq.1) dsmult = -dsmult

 
         dbl_mb(k_heff+n-1+(iref-1)*nref) = dbl_mb(k_r2+m-1)*dsmult
c     1 dbl_mb(k_heff+n-1+(iref-1)*nref)

           if(lusesub) then
         call ga_put(g_heff,nref*(iref-1)+n,nref*(iref-1)+n,1,1,
     1   dbl_mb(k_r2+m-1)*dsmult,1)
           endif


            endif
          goto 100
       endif
       enddo

 100    continue

        if((.not.isfound).and.(abs(dbl_mb(k_r2+m-1)).gt.1.0d-5)) then
         if(nodezero) then
           write(LuOut,"('DEBUG: ',4F16.12)")dbl_mb(k_r2+m-1)
           write(LuOut,"('YOU ARE USING INCOMPLETE MODEL SPACE!')")
         endif
c         call errquit('YOU ARE USING INCOMPLETE MODEL SPACE!',2,MA_ERR)
        endif

        endif ! act

        enddo
        enddo
        enddo
        enddo

        if (.not.ma_pop_stack(l_r2))
     1   call errquit('tce_mrcc_iface_r2: MA problem',1,MA_ERR)
      END IF
      END IF
      END IF
      endif
      if(limprovet) then
       if(lusesub) then
        nxt=NXTASKsub(nprocs,1,mypgid)
       else
        nxt=NXTASK(nprocs,1)
       endif
      endif
      endif
       if(limprovet)count = count + 1
      END DO
      END DO
      END DO
      END DO

      if(limprovet) then
      if(lusesub) then
      nxt=NXTASKsub(-nprocs,1,mypgid)
      call GA_Pgroup_SYNC(mypgid)
      else
      nxt=NXTASK(-nprocs,1)
      call GA_SYNC()
      endif
      endif
c
c ---------------
c    R3 active
c ---------------
c
      if(needr3act) then
      DO p4b = noab+1,noab+nvab
      DO p5b = p4b,noab+nvab
      DO p6b = p5b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h3b = h2b,noab
      IF ((.not.restricted).or.(int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1
     &)+int_mb(k_spin+p6b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1).ne.12)) THEN
      IF (int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p4b-1),ieor(int_mb(k_sym+p5b-1),ieor(int_mb(
     &k_sym+p6b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),int
     &_mb(k_sym+h3b-1)))))) .eq. ieor(irrep_v,irrep_t)) THEN
      IF ((log_mb(k_active+p4b-1).eqv..true.).and.(log_mb(k_active+p5b-1
     &).eqv..true.).and.(log_mb(k_active+p6b-1).eqv..true.).and.(log_mb(
     &k_active+h1b-1).eqv..true.).and.(log_mb(k_active+h2b-1).eqv..true.
     &).and.(log_mb(k_active+h3b-1).eqv..true.)) THEN

      size = int_mb(k_range+p4b-1) * int_mb(k_range+p5b-1) * int_
     &mb(k_range+p6b-1) * int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1)
     &* int_mb(k_range+h3b-1)

      if (.not.ma_push_get(mt_dbl,size,'r3mi',l_r3,k_r3))
     1   call errquit('tce_mrcc_iface_r3: MA problem',0,MA_ERR)

      call get_hash_block(d_r3m(iref),dbl_mb(k_r3),size,
     1   int_mb(k_r3_offsetm(iref)),h3b - 1 + noab * (h2b - 1 + noab *
     1 (h1b - 1 + noab * (p6b - noab - 1 + nvab * (p5b - noab - 1 + 
     1 nvab * (p4b - noab - 1))))))

       m = 0

        orbspin(1) = int_mb(k_spin+p4b-1)-1
        orbspin(2) = int_mb(k_spin+p5b-1)-1
        orbspin(3) = int_mb(k_spin+p6b-1)-1
        orbspin(4) = int_mb(k_spin+h1b-1)-1
        orbspin(5) = int_mb(k_spin+h2b-1)-1
        orbspin(6) = int_mb(k_spin+h3b-1)-1

        do p4=1,int_mb(k_range+p4b-1)
        do p5=1,int_mb(k_range+p5b-1)
        do p6=1,int_mb(k_range+p6b-1)
        do h1=1,int_mb(k_range+h1b-1)
        do h2=1,int_mb(k_range+h2b-1)
        do h3=1,int_mb(k_range+h3b-1)

        m = m + 1

        orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p4b-1)+p4-1))/2
        orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p5b-1)+p5-1))/2
        orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p6b-1)+p6-1))/2
        orbindex(4) = (1 - orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h1b-1)+h1-1))/2
        orbindex(5) = (1 - orbspin(5)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h2b-1)+h2-1))/2
        orbindex(6) = (1 - orbspin(6)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h3b-1)+h3-1))/2

        orbindex(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
        orbindex(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
        orbindex(3) = moindexes(orbindex(3),orbspin(3)+1,iref)
        orbindex(4) = moindexes(orbindex(4),orbspin(4)+1,iref)
        orbindex(5) = moindexes(orbindex(5),orbspin(5)+1,iref)
        orbindex(6) = moindexes(orbindex(6),orbspin(6)+1,iref)

c        write(LuOut,"('Real indexes: [',I4,I4,I4,I4,I4,I4,']')")
c     1 orbindex(1),orbindex(2),orbindex(3),orbindex(4),orbindex(5),
c     1 orbindex(6)
c        write(LuOut,"('Spin indexes : [',I4,I4,I4,I4,I4,I4,']')")
c     1 orbspin(1),orbspin(2),orbspin(3),orbspin(4),orbspin(5),orbspin(6)

       do iu=1,2
        do n=1,nbf
          ioccnew(n,iu) = iocc(n,iref,iu)
        enddo
       enddo

      if((iocc(orbindex(1),iref,orbspin(1)+1).eq.
     1 iocc(orbindex(4),iref,orbspin(4)+1)).or.
     2 (iocc(orbindex(2),iref,orbspin(2)+1).eq.
     3 iocc(orbindex(5),iref,orbspin(5)+1)).or.
     4 (iocc(orbindex(3),iref,orbspin(3)+1).eq.
     1 iocc(orbindex(6),iref,orbspin(6)+1))) then
          goto 300
       endif

      if((orbspin(1).ne.orbspin(4)).or.
     1   (orbspin(2).ne.orbspin(5)).or.
     2   (orbspin(3).ne.orbspin(6))) then
          goto 300
       endif

       if(
     1 ((orbindex(1).eq.orbindex(2)).and.(orbspin(1).eq.orbspin(2))).or.
     1 ((orbindex(1).eq.orbindex(3)).and.(orbspin(1).eq.orbspin(3))).or.
     1 ((orbindex(2).eq.orbindex(3)).and.(orbspin(2).eq.orbspin(3))).or.
     1 ((orbindex(4).eq.orbindex(5)).and.(orbspin(4).eq.orbspin(5))).or.
     1 ((orbindex(4).eq.orbindex(6)).and.(orbspin(4).eq.orbspin(6))).or.
     1 ((orbindex(5).eq.orbindex(6)).and.(orbspin(5).eq.orbspin(6)))
     1 ) then
          goto 300
       endif

         ioccnew(orbindex(1),orbspin(1)+1) = iocc(orbindex(4),iref,
     1 orbspin(4)+1)
         ioccnew(orbindex(2),orbspin(2)+1) = iocc(orbindex(5),iref,
     1 orbspin(5)+1)
         ioccnew(orbindex(3),orbspin(3)+1) = iocc(orbindex(6),iref,
     1 orbspin(6)+1)
         ioccnew(orbindex(4),orbspin(4)+1) = iocc(orbindex(1),iref,
     1 orbspin(1)+1)
         ioccnew(orbindex(5),orbspin(5)+1) = iocc(orbindex(2),iref,
     1 orbspin(2)+1)
         ioccnew(orbindex(6),orbspin(6)+1) = iocc(orbindex(3),iref,
     1 orbspin(3)+1)


       do n=1,nref
       isfound = .true.
        do iu=1,2
         do o=1,nbf
          isfound = isfound.and.(iocc(o,n,iu).eq.ioccnew(o,iu))
         enddo
        enddo
       if(isfound) then
c         write(LuOut,"('Internal amplitude',I4,'->',I4,2F16.12)")iref,n,
c     1 dbl_mb(k_r3+m-1)
            if(iref.ne.n) then
cJB START
             dist = 0
       do iu=1,2
        do i1=1,nmo(iu)
          ioccnew(i1,iu) = iocc(i1,iref,iu)
        enddo
       enddo

             i2 = 0
             do i1=min(orbindex(1),orbindex(4)),
     1 max(orbindex(1),orbindex(4))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(1)-orbindex(4))) then
                  if(iocc(i1+1,iref,orbspin(1)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

         ioccnew(orbindex(1),orbspin(1)+1) = iocc(orbindex(4),iref,
     1 orbspin(4)+1)
         ioccnew(orbindex(4),orbspin(4)+1) = iocc(orbindex(1),iref,
     1 orbspin(1)+1)

             i2 = 0
             do i1=min(orbindex(2),orbindex(5)),
     1 max(orbindex(2),orbindex(5))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(2)-orbindex(5))) then
                  if(ioccnew(i1+1,orbspin(5)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

         ioccnew(orbindex(2),orbspin(2)+1) = iocc(orbindex(5),iref,
     1 orbspin(5)+1)
         ioccnew(orbindex(5),orbspin(5)+1) = iocc(orbindex(2),iref,
     1 orbspin(2)+1)

             i2 = 0
             do i1=min(orbindex(3),orbindex(6)),
     1 max(orbindex(3),orbindex(6))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(3)-orbindex(6))) then
                  if(ioccnew(i1+1,orbspin(6)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

            dsmult = 1.0d0
         if(mod(dist,2).eq.1) dsmult = -dsmult

cJB END
              dbl_mb(k_heff+n-1+(iref-1)*nref) = 
     1 dbl_mb(k_r3+m-1)*dsmult

            endif
          goto 300
       endif
       enddo

 300   continue

        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

        if (.not.ma_pop_stack(l_r3))
     1   call errquit('tce_mrcc_iface_r3: MA problem',1,MA_ERR)

      END IF
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      endif
c
c ---------------
c    R4 active
c ---------------
c
      if(needr4act) then
c      DO p5b = noab+1,noab+nvab
c      DO p6b = p5b,noab+nvab
c      DO p7b = noab+1,noab+nvab
c      DO p8b = p7b,noab+nvab
c      DO h1b = 1,noab
c      DO h2b = 1,noab
c      DO h3b = h2b,noab
c      DO h4b = h3b,noab

      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO p8b = p7b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h3b = h2b,noab
      DO h4b = h3b,noab

      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)+int_mb(k_spin+h1b-1)+i
     &nt_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.1
     &6)) THEN
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     &+int_mb(k_spin+p8b-1) .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-
     &1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+p8b-1),ieor(int_mb(k_sym+h1b-1),ieo
     &r(int_mb(k_sym+h2b-1),ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1)
     &))))))) .eq. ieor(irrep_v,ieor(irrep_t,irrep_t))) THEN
      IF ((log_mb(k_active+p5b-1).eqv..true.).and.(log_mb(k_active+p6b-1
     &).eqv..true.).and.(log_mb(k_active+p7b-1).eqv..true.).and.(log_mb(
     &k_active+p8b-1).eqv..true.).and.(log_mb(k_active+h1b-1).eqv..true.
     &).and.(log_mb(k_active+h2b-1).eqv..true.).and.(log_mb(k_active+h3b
     &-1).eqv..true.).and.(log_mb(k_active+h4b-1).eqv..true.)) THEN

      size = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_
     &mb(k_range+p7b-1) * int_mb(k_range+p8b-1) * int_mb(k_range+h1b-1)
     &* int_mb(k_range+h2b-1) * int_mb(k_range+h3b-1) * int_mb(k_range+h
     &4b-1)

      if (.not.ma_push_get(mt_dbl,size,'r4mi',l_r4,k_r4))
     1   call errquit('tce_mrcc_iface_r4: MA problem',0,MA_ERR)

      call get_hash_block(d_r4m(iref),dbl_mb(k_r4),size,
     1   int_mb(k_r4_offsetm(iref)),(h4b - 1 + noab * (h3b - 1 + noab *
     1(h2b - 1 + noab * (h1b - 1 + noab * (p8b - noab - 1 + nvab * (p7b
     1 - noab - 1 + nvab * (p6b - noab - 1 + nvab * (p5b - noab - 1))))
     1)))))

       m = 0

        orbspin(1) = int_mb(k_spin+p5b-1)-1
        orbspin(2) = int_mb(k_spin+p6b-1)-1
        orbspin(3) = int_mb(k_spin+p7b-1)-1
        orbspin(4) = int_mb(k_spin+p8b-1)-1
        orbspin(5) = int_mb(k_spin+h1b-1)-1
        orbspin(6) = int_mb(k_spin+h2b-1)-1
        orbspin(7) = int_mb(k_spin+h3b-1)-1
        orbspin(8) = int_mb(k_spin+h4b-1)-1

        do p5=1,int_mb(k_range+p5b-1)
        do p6=1,int_mb(k_range+p6b-1)
        do p7=1,int_mb(k_range+p7b-1)
        do p8=1,int_mb(k_range+p8b-1)
        do h1=1,int_mb(k_range+h1b-1)
        do h2=1,int_mb(k_range+h2b-1)
        do h3=1,int_mb(k_range+h3b-1)
        do h4=1,int_mb(k_range+h4b-1)

        m = m + 1

        orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p5b-1)+p5-1))/2
        orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p6b-1)+p6-1))/2
        orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p7b-1)+p7-1))/2
        orbindex(4) = (1 - orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p8b-1)+p8-1))/2
        orbindex(5) = (1 - orbspin(5)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h1b-1)+h1-1))/2
        orbindex(6) = (1 - orbspin(6)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h2b-1)+h2-1))/2
        orbindex(7) = (1 - orbspin(7)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h3b-1)+h3-1))/2
        orbindex(8) = (1 - orbspin(8)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h4b-1)+h4-1))/2

        orbindex(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
        orbindex(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
        orbindex(3) = moindexes(orbindex(3),orbspin(3)+1,iref)
        orbindex(4) = moindexes(orbindex(4),orbspin(4)+1,iref)
        orbindex(5) = moindexes(orbindex(5),orbspin(5)+1,iref)
        orbindex(6) = moindexes(orbindex(6),orbspin(6)+1,iref)
        orbindex(7) = moindexes(orbindex(7),orbspin(7)+1,iref)
        orbindex(8) = moindexes(orbindex(8),orbspin(8)+1,iref)

       do iu=1,2
        do n=1,nbf
          ioccnew(n,iu) = iocc(n,iref,iu)
        enddo
       enddo

      if((iocc(orbindex(1),iref,orbspin(1)+1).eq.
     1 iocc(orbindex(5),iref,orbspin(5)+1)).or.
     2 (iocc(orbindex(2),iref,orbspin(2)+1).eq.
     3 iocc(orbindex(6),iref,orbspin(6)+1)).or.
     4 (iocc(orbindex(3),iref,orbspin(3)+1).eq.
     1 iocc(orbindex(7),iref,orbspin(7)+1)).or.
     2 (iocc(orbindex(4),iref,orbspin(4)+1).eq.
     3 iocc(orbindex(8),iref,orbspin(8)+1))) then
          goto 400
       endif

      if((orbspin(1).ne.orbspin(5)).or.
     1   (orbspin(2).ne.orbspin(6)).or.
     2   (orbspin(3).ne.orbspin(7)).or.
     3   (orbspin(4).ne.orbspin(8))) then
          goto 400
       endif

      if(
     1 ((orbindex(1).eq.orbindex(2)).and.(orbspin(1).eq.orbspin(2))).or.
     1 ((orbindex(1).eq.orbindex(3)).and.(orbspin(1).eq.orbspin(3))).or.
     1 ((orbindex(1).eq.orbindex(4)).and.(orbspin(1).eq.orbspin(4))).or.
     1 ((orbindex(2).eq.orbindex(3)).and.(orbspin(2).eq.orbspin(3))).or.
     1 ((orbindex(2).eq.orbindex(4)).and.(orbspin(2).eq.orbspin(4))).or.
     1 ((orbindex(3).eq.orbindex(4)).and.(orbspin(3).eq.orbspin(4))).or.
     1 ((orbindex(5).eq.orbindex(6)).and.(orbspin(5).eq.orbspin(6))).or.
     1 ((orbindex(5).eq.orbindex(7)).and.(orbspin(5).eq.orbspin(7))).or.
     1 ((orbindex(5).eq.orbindex(8)).and.(orbspin(5).eq.orbspin(8))).or.
     1 ((orbindex(6).eq.orbindex(7)).and.(orbspin(6).eq.orbspin(7))).or.
     1 ((orbindex(6).eq.orbindex(8)).and.(orbspin(6).eq.orbspin(8))).or.
     1 ((orbindex(7).eq.orbindex(8)).and.(orbspin(7).eq.orbspin(8)))
     1 ) then
          goto 400
       endif

         ioccnew(orbindex(1),orbspin(1)+1) = iocc(orbindex(5),iref,
     1 orbspin(5)+1)
         ioccnew(orbindex(2),orbspin(2)+1) = iocc(orbindex(6),iref,
     1 orbspin(6)+1)
         ioccnew(orbindex(3),orbspin(3)+1) = iocc(orbindex(7),iref,
     1 orbspin(7)+1)
         ioccnew(orbindex(4),orbspin(4)+1) = iocc(orbindex(8),iref,
     1 orbspin(8)+1)
         ioccnew(orbindex(5),orbspin(5)+1) = iocc(orbindex(1),iref,
     1 orbspin(1)+1)
         ioccnew(orbindex(6),orbspin(6)+1) = iocc(orbindex(2),iref,
     1 orbspin(2)+1)
         ioccnew(orbindex(7),orbspin(7)+1) = iocc(orbindex(3),iref,
     1 orbspin(3)+1)
         ioccnew(orbindex(8),orbspin(8)+1) = iocc(orbindex(4),iref,
     1 orbspin(4)+1)

       do n=1,nref
       isfound = .true.
        do iu=1,2
         do o=1,nbf
          isfound = isfound.and.(iocc(o,n,iu).eq.ioccnew(o,iu))
         enddo
        enddo
       if(isfound) then
c         write(LuOut,"('Internal amplitude',I4,'->',I4,2F16.12)")iref,n,
c     1 dbl_mb(k_r4+m-1)
            if(iref.ne.n) then
c ===============
             dist = 0
       do iu=1,2
        do i1=1,nmo(iu)
          ioccnew(i1,iu) = iocc(i1,iref,iu)
        enddo
       enddo

             i2 = 0
             do i1=min(orbindex(1),orbindex(5)),
     1 max(orbindex(1),orbindex(5))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(1)-orbindex(5))) then
                  if(iocc(i1+1,iref,orbspin(1)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

         ioccnew(orbindex(1),orbspin(1)+1) = iocc(orbindex(5),iref,
     1 orbspin(5)+1)
         ioccnew(orbindex(5),orbspin(5)+1) = iocc(orbindex(1),iref,
     1 orbspin(1)+1)

             i2 = 0
             do i1=min(orbindex(2),orbindex(6)),
     1 max(orbindex(2),orbindex(6))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(2)-orbindex(6))) then
                  if(iocc(i1+1,iref,orbspin(2)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

         ioccnew(orbindex(2),orbspin(2)+1) = iocc(orbindex(6),iref,
     1 orbspin(6)+1)
         ioccnew(orbindex(6),orbspin(6)+1) = iocc(orbindex(2),iref,
     1 orbspin(2)+1)

             i2 = 0
             do i1=min(orbindex(3),orbindex(7)),
     1 max(orbindex(3),orbindex(7))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(3)-orbindex(7))) then
                  if(iocc(i1+1,iref,orbspin(3)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

         ioccnew(orbindex(3),orbspin(3)+1) = iocc(orbindex(7),iref,
     1 orbspin(7)+1)
         ioccnew(orbindex(7),orbspin(7)+1) = iocc(orbindex(3),iref,
     1 orbspin(3)+1)

             i2 = 0
             do i1=min(orbindex(4),orbindex(8)),
     1 max(orbindex(4),orbindex(8))
               i2 = i2 + 1
               if(i2.lt.abs(orbindex(4)-orbindex(8))) then
                  if(iocc(i1+1,iref,orbspin(4)+1).eq.1) then
                     dist=dist+1
                  endif
               endif
             enddo

            dsmult = 1.0d0
         if(mod(dist,2).eq.1) dsmult = -dsmult

c         if(nodezero)then
c          write(6,"('T4 iref/n:',2I4,f4.1)")iref,n,dsmult
c         endif
c ===============
              dbl_mb(k_heff+n-1+(iref-1)*nref) =
     1 dbl_mb(k_r4+m-1)*dsmult
            endif
          goto 400
       endif
       enddo

 400  continue

        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

        if (.not.ma_pop_stack(l_r4))
     1   call errquit('tce_mrcc_iface_r4: MA problem',1,MA_ERR)

      END IF
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

      endif

       endif !sub

       enddo !iref

       return
       end
c
c ==============================================
c     Diagonalize effective Hamiltonian
c ==============================================
c
        subroutine tce_diagonalize_heff(rtdb,iter)
        implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_main.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"

c        integer nref
        double precision heff(nref,nref)
        double precision vl(nref,nref)
        double precision vr(nref,nref)
        double precision er(nref)
        double precision ei(nref)
        double precision work(4*nref)
        integer info
        integer i,j,l
        double precision x
        integer k
        logical nodezero
        integer rtdb,itarget
c        integer nrootmuc
        double precision dvalue,dsum
        double precision dfin
        character*4 ds
        integer l_buff,k_buff
c        integer iignore
        double precision isum
        integer ddblsize,inntsize

ckbn rootfromoverlap
        double precision rootoverlap,rootoverlapmax
        integer stateselectedbyoverlap
        integer iter

        nodezero = (ga_nodeid().eq.0)  
        ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
        inntsize=MA_sizeof(MT_INT,1,MT_BYTE)

        call ga_sync()

        if(lusesub) then
        do i=1,nref*nref
          dbl_mb(k_heff+i-1) = 0.0d0
        enddo
         call ga_sync()
         call ga_get(g_heff,1,nref*nref,1,1,
     1   dbl_mb(k_heff),1)
c         call ga_print(g_heff)
        endif

        do i=1,nref
         do j=1,nref
          x = dbl_mb(k_heff+j-1+(i-1)*nref)
          heff(j,i) = x
         end do
        end do

        do i=1,nref
          er(i) = 0.0d0
          ei(i) = 0.0d0
          do j=1,nref
            vl(i,j)=0.0d0
            vr(i,j)=0.0d0
          enddo
        enddo

       if(nodezero.and.(nref.lt.21)) then
c           call ma_print(dbl_mb(k_heff),nref,nref,'Heff')

        write(LuOut,"(/,'Heff',/,
     1 '=============================================')")
        do i=1,nref
         write(LuOut,"(i5,i5,100F14.8)")ga_nodeid(),i,
     1 (dbl_mb(k_heff+(j-1)*nref+i-1),j=1,nref)
        enddo

       endif
c       call ga_sync()
c       call util_flush(LuOut)
c       if((ga_nodeid().eq.5).and.(nref.lt.21)) then
c           call ma_print(dbl_mb(k_heff),nref,nref,'Heff')
c
c        write(LuOut,"(/,'Heff 5',/,
c     1 '=============================================')")
c        do i=1,nref
c         write(LuOut,"(i5,i5,100F14.8)")ga_nodeid(),i,
c     1 (dbl_mb(k_heff+(j-1)*nref+i-1),j=1,nref)
c        enddo

c       endif
c       call util_flush(LuOut)
c       call ga_sync()

c         call util_flush(LuOut)
c        write(6,*)'BEFORE',ga_nodeid()
c         call util_flush(LuOut)
        call ga_sync()
c       if(nodezero)write(6,*)'TEST 3'
        if(nodezero) then
c        call DGEEV('V','V',nref,heff,nref,er,ei,vl,nref,vr,
        call util_dgeev('V','V',nref,heff,nref,er,ei,vl,nref,vr,
     $                  nref,work,4*nref,info)
        if(info .ne. 0) call errquit('Heff diagonalization',0,CALC_ERR)
        call amp_stabilization(vl,vr,nref)
        endif
c      if(nodezero)write(6,*)'TEST 4'
        call ga_sync()

c         call util_flush(LuOut)
c        write(6,*)'AFTER',ga_nodeid()
c         call util_flush(LuOut)
         
c       if(nodezero.and..not.lconverged) then
c           call ma_print(dbl_mb(k_heff),nref,nref,'Heff')
c       endif

c       if(lconverged.and.nodezero) then
       if(nodezero) then
        if(nref.lt.21) then

        write(LuOut,"(/,'Eigenvalues (real and imaginary)',/,
     1 '=============================================')")
        do i=1,nref
           write(LuOut,"(F18.12,100F14.8)")er(i),ei(i)
        enddo

        write(LuOut,"(/,'Left eigenvectors',/,
     1 '=============================================')")
        do i=1,nref
           write(LuOut,"(i5,100F14.8)")i,(vl(i,j),j=1,nref)
        enddo

c       write(LuOut,"(/,'Left eigenvectors - squares',/,
c     1 '=============================================')")
c        do i=1,nref
c           write(LuOut,"(i5,35f18.12)")i,((vl(i,j)*vl(i,j)),j=1,nref)
c        enddo
        
c        endif
c        endif      

        write(LuOut,"(/,'Right eigenvectors',/,
     1 '=============================================')")
        do i=1,nref
           write(LuOut,"('VR',i5,100F14.8)")i,(vr(i,j),j=1,nref)
        enddo

        endif
        endif

        do i=1,nref
         do j=1,nref
          dbl_mb(k_sqc+(i-1)*nref+j-1)=vr(i,j)
          dbl_mb(k_sqcl+(i-1)*nref+j-1)=vl(i,j)
         enddo
        enddo

       call ga_brdcst(Msg_Vec_EVal+21,dbl_mb(k_sqc),
     1 ddblsize*nref*nref, 0)
       call ga_sync()
       call ga_brdcst(Msg_Vec_EVal+20,dbl_mb(k_sqcl),
     1 ddblsize*nref*nref, 0)
       call ga_sync()

       if (.not.rtdb_get(rtdb,'mrcc:zignore',mt_int,1,iignore))
     1 iignore = 0

        if(.not.nodezero) then
        do i=1,nref
         do j=1,nref
          vr(i,j) = dbl_mb(k_sqc+(i-1)*nref+j-1)
          vl(i,j) = dbl_mb(k_sqcl+(i-1)*nref+j-1)
         enddo
        enddo
        endif

c--------

        epsilon = 0.0d0
        do i=1,nref
           isum = 0.0d0
           do j=1,nref
             isum = isum + vr(j,i)
           enddo
           if((iignore.eq.0).or.(abs(isum).gt.1.0d-5)) then
             if(epsilon.gt.min(epsilon,er(i)))mkroot=i
             epsilon = min(epsilon,er(i))
           endif
        enddo

c--------

       if (.not.rtdb_get(rtdb,'mrcc:rootmuc',mt_int,1,nrootmuc))
     1 nrootmuc = 0
c ------
      if(nrootmuc.gt.0) then ! 1
        dfin = 0.0d0

        do j=1,nref
        dsum = 0.0d0
        isum = 0.0d0
        do i=1,nref
           isum = isum + vr(i,j)
        enddo
        if((iignore.eq.0).or.(abs(isum).gt.1.0d-5)) then
        do i=1,nref
          write(ds,"(I3.3)")i
       if (.not.rtdb_get(rtdb,'mrcc:rootmuc'//ds,mt_dbl,1,dvalue))
     1 dvalue = 0.0d0
          if(dvalue.lt.0.0d0) then
            if(abs(vr(i,j)).lt.abs(dvalue))
     1     dsum = dsum + (dvalue)*(vr(i,j))
          else
           dsum = dsum + (dvalue)*(vr(i,j))
          endif
        enddo
c            write(6,"('SUM:',I4,4F16.12)")j,
c     1 abs(dsum)
          if(dfin.lt.abs(dsum)) then
c            write(6,"('I am watching reference #',I4,4F16.12)")j
            mkroot=j
            epsilon = er(j)
            dfin = abs(dsum)
          endif
        endif
        enddo

      else ! 1

      if (rtdb_get(rtdb,'bwcc:targetroot',mt_int,1,itarget)) then
c       mkroot = itarget
       do i=1,nref

           isum = 0.0d0
           do j=1,nref
c             if(abs(vr(j,i)).lt.1.0d-8) 
              isum = isum + vr(j,i)
           enddo
           if((iignore.eq.0).or.(abs(isum).gt.1.0d-5)) then

        dvalue = er(i)
        k = 0     
        do j=1,nref
        if(j.ne.i) then

           isum = 0.0d0
           do l=1,nref
c             if(abs(vr(l,j)).lt.1.0d-8) 
             isum = isum + vr(l,j)
           enddo
           if((iignore.eq.0).or.(abs(isum).gt.1.0d-5)) then


         if(er(j).lt.er(i))k=k+1

           endif

        endif
        enddo
        if(k.eq.(itarget-1)) then
         itarget = i
         goto 63454
        endif

            endif

       enddo

63454  continue

         mkroot = itarget
         epsilon = er(itarget)

      endif
      endif


ckbn  find the requested root according to overlap instead
      call ga_sync()
      if(iter .gt. iroottooverlapiter) then
c      if(.not.lconverged) then
        if(lrootfromoverlap) then
         stateselectedbyoverlap=0
         rootoverlapmax=-100.0d0
         do i=1,nref
          rootoverlap=0.0d0
          isum = 0.0d0
          do j=1,nref
c          if(nodezero)write(*,'(A,2I5,A,F17.10)')"vr(",j,i,")",vr(j,i)
           rootoverlap=rootoverlap+bwcoefwanted(j)*vr(j,i)
           isum = isum + vr(j,i)
          enddo
          rootoverlap=abs(rootoverlap)
c          write(*,'(A,F17.10)') "Overlap ", rootovelap
c          write(*,*) "i1 ",i," rootoverlapmax ",rootoverlapmax,stateselectedbyoverlap,i
          if(nodezero) write(*,'(A,I5,F17.10,A,F17.10)') "Overlap ",i,
     +     rootoverlap, " isum ",isum
          if((iignore.eq.0).or.(abs(isum).gt.1.0d-5)) then
           if(rootoverlap .gt. rootoverlapmax) then
            rootoverlapmax=rootoverlap
            stateselectedbyoverlap=i
c           write(*,*) "i2 ",i," rootovelapmax ",rootovelapmax,stateselectedbyoverlap,i
           endif
          endif
         enddo
         mkroot= stateselectedbyoverlap
         epsilon = er(stateselectedbyoverlap)
         if(nodezero) write(LuOut,'(A,I5)') 
     +    'MRCC root selected by overlap is ',stateselectedbyoverlap
        endif ! lrootfromoverlap
c      else ! lconverged
c       if(nodezero) write(*,'(A,L3,A,I5)') "Converged ",lconverged,
c    +   " so keeping old root ", mkrootold
c      endif !lconverged
      endif 
      call ga_sync()






        if (nodezero.and.(nref.lt.21)) then
           write(6,"('Target root: ',I4)")mkroot
        end if

ckbn	Introduce some checks before proceeding further.

ckbn    check of mkroot, crash if mkroot gt nref
        if(nodezero) then 
         if(mkroot.gt.nref) call errquit
     +    ('root followed is greater than total number of references',0,
     +    CALC_ERR)
        endif


ckbn check whether there is imaginary eigen values
        do i=1,nref
c         write(LuOut,'(A,F17.10)')'Imaginary eigenvalue',ei(i)
c         if (nodezero) call util_flush(LuOut)
         if(abs(ei(i)).ge.toleiimag) then
          write(LuOut,'(A,F15.10)')
     +            'Warning: complex Heff eigenvalue detected',ei(i)
          if(i.eq.mkroot) then
           write(LuOut,*) "ignorecomplex1 ", ignorecomplex
c          if (rtdb_get(rtdb,'mrcc:ignorecomplex',mt_log,1,
c     +     ignorecomplex)) ignorecomplex = .true.
c           call errquit('Complex root found and no ignorecomplex option',
c     +      0,RTDB_ERR)
c          else
c           if(nodezero) write(*,*) "ignorecomplex1. ", ignorecomplex
c           ignorecomplex = .true.
c          endif
c           if(nodezero) write(*,*) "ignorecomplex2 ", ignorecomplex
           if(ignorecomplex) then
            write(LuOut,'(A,F15.10)')
     +       'Warning: Proceeding with complex Heff eigenvalue ',ei(i)
           else
            call errquit('Warning:complex Heff eigenvalue detected',0,
     +                    UNKNOWN_ERR)
           endif
          endif
         endif
        enddo




        if(nref.gt.20) then
         if(nodezero)write(LuOut,"(/,'Right eigenvector for target',I7,
     1 /)")mkroot
c      if (.not.ma_push_get(mt_dbl,nref,'buff',l_buff,k_buff))
c     1   call errquit('tce_mrcc_iface_buff: MA problem',0,MA_ERR)
         do i=1,nref
          if(abs(vr(i,mkroot)).gt.0.05d0) then
            if(nodezero)write(LuOut,"(I7,'  ',F11.8)")i,vr(i,mkroot)
          endif
         enddo
c        if (.not.ma_pop_stack(l_buff))
c     1   call errquit('tce_mrcc_iface_buff: MA problem',1,MA_ERR)

        endif

        if (nodezero) call util_flush(LuOut)

cjb broadcasts

        call ga_brdcst(Msg_Vec_EVal+MSGINT+30,mkroot,inntsize, 0)
        call ga_brdcst(Msg_Vec_EVal+MSGDBL+31,epsilon,ddblsize, 0)
        call ga_sync()

c        write(LuOut,"(/,'Epsilon: ',2F16.12,/)") epsilon

        return
        end

c
c ==============================================
c     Clean internal amplitudes
c ==============================================
c
        subroutine tce_internal_t_zero(d_t1m,d_t2m,k_t1_offsetm,
     1 k_t2_offsetm,lneedt3,d_t3m,k_t3_offsetm,rtdb)
c     1 k_t2_offsetm,nref,lneedt3,d_t3m,k_t3_offsetm,rtdb)
        implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_main.fh"
#include "rtdb.fh"

        integer d_t1m(maxref),d_t2m(maxref)
        integer d_t3m(maxref)
        integer k_t1_offsetm(maxref),k_t2_offsetm(maxref)
        integer k_t3_offsetm(maxref)
        integer size,p5b,h6b
c        integer nref
        integer l_t1,k_t1
        integer l_t2,k_t2
        integer l_t3,k_t3
        integer i,j,k,l,m
        integer iref,rtdb
        integer p1b,p2b,h3b,h4b
        integer p3b,p4b,h5b,p6b
        integer h1b,h2b
c        logical lneedt3,limprovet
        logical lneedt3
        integer orbindex(6),orbspin(6)
        EXTERNAL NXTASKsub
        EXTERNAL NXTASK
        INTEGER NXTASKsub
        INTEGER NXTASK
        INTEGER nxt
        INTEGER nprocs
        INTEGER count

c       if (.not.rtdb_get(rtdb,'mrcc:improvetiling',mt_log,1,limprovet))
c     1 limprovet = .false.

      do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +ga_nnodes()+ga_nodeid())).or.(.not.lusesub)) then

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

         do p5b = noab+1,noab+nvab
         do h6b = 1,noab

      if (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) then
      if (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)).eq.irrep_t)then
      if ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) then
      if(log_mb(k_isactive(iref)+p5b-1).and.
     &log_mb(k_isactive(iref)+h6b-1)) then

        size = int_mb(k_range+h6b-1) * int_mb(k_range+p5b-1)
        if (.not.ma_push_get(mt_dbl,size,'t1',l_t1,k_t1))
     1   call errquit('tce_mrcc_iface_t1: MA problem',0,MA_ERR)

c         call dfill(size,0.0d0,dbl_mb(k_t1),1)
         do i=1,size
          dbl_mb(k_t1+i-1)=0.0d0
         enddo

cjb ============================

         if(limprovet) then

        call get_hash_block(d_t1m(iref),dbl_mb(k_t1),size,
     1   int_mb(k_t1_offsetm(iref)),h6b-1+noab*(p5b-noab-1))

        k=0
        do i=1,int_mb(k_range+p5b-1)
        do j=1,int_mb(k_range+h6b-1)
          k = k + 1

        orbspin(1) = int_mb(k_spin+p5b-1)-1
        orbspin(2) = int_mb(k_spin+h6b-1)-1

        orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p5b-1)+i-1))/2
        orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h6b-1)+j-1))/2

        orbindex(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
        orbindex(2) = moindexes(orbindex(2),orbspin(2)+1,iref)

        if(isactive(orbindex(1),orbspin(1)+1).and.
     1 isactive(orbindex(2),orbspin(2)+1).or.(.not.limprovet)) then

          dbl_mb(k_t1+k-1) = 0.0d0

          endif
          enddo
          enddo

          endif ! limprovet

c =============================

c         do i=1,size
c           dbl_mb(k_t1+i-1) = 0.0d0
c         enddo

         call put_hash_block(d_t1m(iref),dbl_mb(k_t1),size,
     1   int_mb(k_t1_offsetm(iref)),((p5b-noab-1)*noab+h6b-1))

        if (.not.ma_pop_stack(l_t1))
     1   call errquit('tce_mrcc_iface_t1: MA problem',1,MA_ERR)

      endif
      endif
      endif
      endif
      enddo
      enddo

cjb new in parallel

      nxt = 0
      if(limprovet) then
       if(lusesub) then
         nprocs=GA_pgroup_NNODES(mypgid)
         nxt=NXTASKsub(nprocs,1,mypgid)
       else
         nprocs=GA_NNODES()
         nxt=NXTASK(nprocs,1)
       endif
      count = 0
      endif

      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO h3b = 1,noab
      DO h4b = h3b,noab

      IF ((nxt.eq.count).or.(.not.limprovet)) THEN

      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      if(log_mb(k_isactive(iref)+p1b-1).and.
     1 log_mb(k_isactive(iref)+p2b-1).and.
     2 log_mb(k_isactive(iref)+h3b-1).and.
     3 log_mb(k_isactive(iref)+h4b-1)) then

      size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)

        if (.not.ma_push_get(mt_dbl,size,'t2',l_t2,k_t2))
     1   call errquit('tce_mrcc_iface_t2: MA problem',0,MA_ERR)

c         call dfill(size,0.0d0,dbl_mb(k_t2),1)
         do i=1,size
          dbl_mb(k_t2+i-1)=0.0d0
         enddo

c ===============================================================

         if(limprovet) then

        call get_hash_block(d_t2m(iref),dbl_mb(k_t2),size,
     1   int_mb(k_t2_offsetm(iref)),h4b-1+noab*(h3b-1+noab*(p2b-
     &noab-1+nvab*(p1b - noab - 1))))

       m = 0
        do i=1,int_mb(k_range+p1b-1)
        do j=1,int_mb(k_range+p2b-1)
        do k=1,int_mb(k_range+h3b-1)
        do l=1,int_mb(k_range+h4b-1)
         m = m + 1

        orbspin(1) = int_mb(k_spin+p1b-1)-1
        orbspin(2) = int_mb(k_spin+p2b-1)-1
        orbspin(3) = int_mb(k_spin+h3b-1)-1
        orbspin(4) = int_mb(k_spin+h4b-1)-1

        orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p1b-1)+i-1))/2
        orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+p2b-1)+j-1))/2
        orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h3b-1)+k-1))/2
        orbindex(4) = (1 - orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offset+h4b-1)+l-1))/2

        orbindex(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
        orbindex(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
        orbindex(3) = moindexes(orbindex(3),orbspin(3)+1,iref)
        orbindex(4) = moindexes(orbindex(4),orbspin(4)+1,iref)

        if(isactive(orbindex(1),orbspin(1)+1).and.
     1 isactive(orbindex(2),orbspin(2)+1).and.
     2 isactive(orbindex(3),orbspin(3)+1).and.
     3 isactive(orbindex(4),orbspin(4)+1).or.(.not.limprovet)) then

         dbl_mb(k_t2+m-1)=0.0d0

          endif
         enddo
         enddo
         enddo
         enddo

         endif

c ===============================================================

c        do i=1,size
c           write(LuOut,*)dbl_mb(k_t2+i-1),'->','0.00000000'
c           dbl_mb(k_t2+i-1) = 0.0d0
c        enddo

        call put_hash_block(d_t2m(iref),dbl_mb(k_t2),size,
     1   int_mb(k_t2_offsetm(iref)),((((p1b-noab-1)*nvab+p2b-noab-1)
     2   *noab+h3b-1)*noab+h4b-1))

        if (.not.ma_pop_stack(l_t2))
     1   call errquit('tce_mrcc_iface_t2: MA problem',1,MA_ERR)

      END IF
      END IF
      END IF
      endif

      if(limprovet) then
       if(lusesub) then
        nxt=NXTASKsub(nprocs,1,mypgid)
       else
        nxt=NXTASK(nprocs,1)
       endif
      endif

      endif

       if(limprovet)count = count + 1

      END DO
      END DO
      END DO
      END DO

      if(limprovet) then
      if(lusesub) then
      nxt=NXTASKsub(-nprocs,1,mypgid)
      call GA_Pgroup_SYNC(mypgid)
      else
      nxt=NXTASK(-nprocs,1)
      call GA_SYNC()
      endif
      endif

      if(lneedt3) then

      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
      IF (int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),ieor(int_mb(
     &k_sym+p4b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h5b-1),int
     &_mb(k_sym+h6b-1)))))) .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1
     &)+int_mb(k_spin+p4b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+i
     &nt_mb(k_spin+h6b-1).ne.12)) THEN
      IF ((log_mb(k_isactive(iref)+p2b-1).eqv..true.).and.
     1 (log_mb(k_isactive(iref)+p3b-1).eqv..true.).and.
     2 (log_mb(k_isactive(iref)+p4b-1).eqv..true.).and.
     3 (log_mb(k_isactive(iref)+h1b-1).eqv..true.).and.
     4 (log_mb(k_isactive(iref)+h5b-1).eqv..true.).and.
     5 (log_mb(k_isactive(iref)+h6b-1).eqv..true.)) THEN

      size = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) * int_
     &mb(k_range+p4b-1) * int_mb(k_range+h1b-1) * int_mb(k_range+h5b-1)
     &* int_mb(k_range+h6b-1)

        if (.not.ma_push_get(mt_dbl,size,'t3',l_t3,k_t3))
     1   call errquit('tce_mrcc_iface_t3: MA problem',0,MA_ERR)

        do i=1,size
           dbl_mb(k_t3+i-1) = 0.0d0
        enddo
        call put_hash_block(d_t3m(iref),dbl_mb(k_t3),size,
     1   int_mb(k_t3_offsetm(iref)),(h6b - 1 + noab * 
     2 (h5b - 1 + noab * (h1b
     &- 1 + noab * (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab * (p2
     &b - noab - 1)))))))

        if (.not.ma_pop_stack(l_t3))
     1   call errquit('tce_mrcc_iface_t3: MA problem',1,MA_ERR)

      END IF
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

      endif !needt3

c      write(6,"('CPU BEFORE',I4,F16.12)")ga_nodeid(),util_cpusec()
c      call ga_sync()
      if(lusesub) call ga_pgroup_sync(
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
c      write(6,"('CPU AFTER',I4,F16.12)")ga_nodeid(),util_cpusec()

      endif !sub

      enddo !iref

        return
        end


c $Id$
