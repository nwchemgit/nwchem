ckbn         subroutine tce_mrcc_print_summary(iter,dsummary)
ckbn         implicit none
ckbn#include "tce.fh"
ckbn#include "mafdecls.fh"
ckbn#include "stdio.fh"
ckbn#include "rtdb.fh"
ckbn#include "errquit.fh"
ckbn#include "sym.fh"
ckbn#include "tce_mrcc.fh"
ckbn#include "global.fh"
ckbn#include "tce_main.fh"
ckbn
ckbn      integer iter
ckbn      double precision dsummary(5000,2)
ckbn      logical nodezero
ckbn      integer i,j
ckbn
ckbn      nodezero = (ga_nodeid().eq.0)
ckbn 
ckbn      dsummary(1,2) = 0.0d0
ckbn 
ckbn      if(nodezero) then 
ckbn        write(LuOut,"(/)")
ckbn        call util_print_centered
ckbn     1  (LuOut,'======================================================',
ckbn     2  40,.false.)
ckbn       call util_print_centered
ckbn     1  (LuOut,'Summary output from iterations',
ckbn     2  40,.false.)
ckbn        call util_print_centered
ckbn     1  (LuOut,'======================================================',
ckbn     2  40,.false.)
ckbn        write(LuOut,*)        
ckbn
ckbn        write(LuOut,9100)
ckbn        call util_print_centered
ckbn     1  (LuOut,'------------------------------------------------------',
ckbn     2  40,.false.)
ckbn        do i=1,iter
ckbn         write(LuOut,9200)i,(dsummary(i,j),j=1,2)
ckbn        enddo
ckbn        call util_print_centered
ckbn     1  (LuOut,'======================================================',
ckbn     2  40,.false.)
ckbn      endif
ckbn
ckbn      if (nodezero) call util_flush(LuOut)
ckbn
ckbn 9100 format(15x,'  Iteration  ',3x,'Energy (a.u.)',9x,'Corr. energy')
ckbn 9200 format(17x,1x,i4,40f22.12,7x)
ckbn
ckbn      return
ckbn      end


         subroutine tce_mrcc_print_t1type(d_t1,k_t1_offset,iref,irefnew)
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

      integer p5b,h6b
      logical nodezero
      integer d_t1,k_t1_offset
      integer l_t1,k_t1
      integer size
      integer iref,inoab,counter
      integer orbspin(2),orbindex(2)
      integer i,j
      character*2 s,r
      integer irefnew,lastparam

      nodezero = (ga_nodeid().eq.0)
      inoab = nblcks(1,iref)+nblcks(2,iref)
c DEBUG
c      nodezero=.true.
c
      lastparam=nblcks(1,iref)+nblcks(2,iref)+nblcks(3,iref)+
     1 nblcks(4,iref)

      if (nodezero) write(LuOut,'(A,I4,A)')
     + "Printing T1/R1 array when T1 >0.1 for ",iref," if any"

      DO p5b = nblcks(1,iref)+nblcks(2,iref)+1,lastparam
      DO h6b = 1,nblcks(1,iref)+nblcks(2,iref)
      IF (int_mb(k_spinm(iref)+p5b-1) .eq. int_mb(k_spinm(iref)+
     1h6b-1)) THEN
      IF (ieor(int_mb(k_symm(iref)+p5b-1),int_mb(k_symm(iref)+
     1h6b-1)) .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spinm(iref)+p5b-1)+
     1int_mb(k_spinm(iref)+h6b-1).ne.4)) THEN


         size = int_mb(k_rangem(iref)+p5b-1) *
     1 int_mb(k_rangem(iref)+h6b-1)

         if (.not.ma_push_get(mt_dbl,size,'t1',l_t1,k_t1))
     1   call errquit('tce_c1_offs: MA problem',0,MA_ERR)
         call get_hash_block(d_t1,dbl_mb(k_t1),size,
     1   int_mb(k_t1_offset),((p5b-inoab-1)*inoab+h6b-1))
 
         counter = 0

         do i=1,int_mb(k_rangem(iref)+p5b-1)
           orbspin(1) = int_mb(k_spinm(iref)+p5b-1)-1
           do j=1,int_mb(k_rangem(iref)+h6b-1)
            orbspin(2) = int_mb(k_spinm(iref)+h6b-1)-1

             counter = counter + 1

             orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+p5b-1)+i-1))/2
             orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+h6b-1)+j-1))/2

c       if(nodezero)write(6,*)orbindex(1),orbindex(2)


c       if(nodezero)write(6,*)orbindex(1),orbindex(2)

              if(orbspin(1).eq.0) then 
               s='Pa'
              else 
               s='Pb'
              endif
                
              if(orbspin(2).eq.0) then
               r='Ha' 
              else 
               r='Hb'
              endif

      if(nodezero .and. (abs(dbl_mb(k_t1+counter-1)) .gt. 0.1d0)) then
               write(LuOut,"('(',I5,a2,I5,a2,')=',2F16.12)")

     +moindexes(orbindex(2),orbspin(2)+1,iref),r,
     +moindexes(orbindex(1),orbspin(1)+1,iref),s,
     +dbl_mb(k_t1+counter-1)
      endif

           enddo
          enddo

         if (.not.ma_pop_stack(l_t1))
     1   call errquit('tce_c1_offs: MA problem',2,MA_ERR)

      endif
      endif
      endif

      enddo
      enddo

      call util_flush(LuOut)

       return
       end

         subroutine tce_mrcc_print_t2type(d_t2,k_t2_offset,iref,irefnew)
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

      integer p1b,p2b,h3b,h4b
      logical nodezero
      integer d_t2,k_t2_offset
      integer l_a,k_a
      integer size
      integer iref,inoab,counter
      integer orbspin(4),orbindex(4)
      integer i,j,m,n
      character*2 s,r,t,u
      integer irefnew
      integer ihash,invab

      nodezero = (ga_nodeid().eq.0)
c DEBUG
c      nodezero=.true.
c

      if (nodezero) write(LuOut,'(A,I4,A)')
     + "Printing T2/R2 array when T2 >0.1 for ",iref," if any"

      inoab = nblcks(1,iref)+nblcks(2,iref)
      invab = nblcks(3,iref)+nblcks(4,iref)

      DO p1b = inoab+1,inoab+invab
      DO p2b = p1b,inoab+invab
      DO h3b = 1,inoab
      DO h4b = h3b,inoab

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

         counter = 0

        ihash = (h4b - 1 +
     2  inoab * (h3b - 1 + inoab * (p2b-
     &inoab - 1 + invab * (p1b - inoab - 1))))


        call get_hash_block(d_t2,dbl_mb(k_a),size,
     1  int_mb(k_t2_offset),ihash)

         do i=1,int_mb(k_rangem(iref)+p1b-1)
          orbspin(1) = int_mb(k_spinm(iref)+p1b-1)-1
           do j=1,int_mb(k_rangem(iref)+p2b-1)
            orbspin(2) = int_mb(k_spinm(iref)+p2b-1)-1
             do m=1,int_mb(k_rangem(iref)+h3b-1)
             orbspin(3) = int_mb(k_spinm(iref)+h3b-1)-1
              do n=1,int_mb(k_rangem(iref)+h4b-1)
              orbspin(4) = int_mb(k_spinm(iref)+h4b-1)-1

             counter = counter + 1

             orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+p1b-1)+i-1))/2
             orbindex(2) = (1 - orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+p2b-1)+j-1))/2
             orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+h3b-1)+m-1))/2
             orbindex(4) = (1 - orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+h4b-1)+n-1))/2
              if(orbspin(1).eq.0) then
               s='Pa'
              else
               s='Pb'
              endif

              if(orbspin(2).eq.0) then
               r='Pa'
              else
               r='Pb'
              endif

              if(orbspin(3).eq.0) then
               t='Ha'
              else
               t='Hb'
              endif

              if(orbspin(4).eq.0) then
               u='Ha'
              else
               u='Hb'
              endif

      if(nodezero .and. (abs(dbl_mb(k_a+counter-1)) .gt. 0.1d0)) then
               write(LuOut,"('(',I5,a2,I5,a2,I5,a2,I5,a2,')=',2F16.12)")
     +moindexes(orbindex(3),orbspin(3)+1,iref),t,
     +moindexes(orbindex(4),orbspin(4)+1,iref),u,
     +moindexes(orbindex(1),orbspin(1)+1,iref),s,
     +moindexes(orbindex(2),orbspin(2)+1,iref),r,
     +dbl_mb(k_a+counter-1)
      endif


           enddo
           enddo
           enddo
          enddo

         if (.not.ma_pop_stack(l_a))
     1   call errquit('tce_c1_offs: MA problem',2,MA_ERR)

      endif
      endif
      endif

      enddo
      enddo
           enddo
           enddo

      if (nodezero) call util_flush(LuOut)

       return
       end

         subroutine tce_mrcc_debug_pfile(d_a,size_a,fname,iter,iref)
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

      integer d_a,size_a,iter
      character*4 fname 
      character*3 sname,siter,sipg,sref
      double precision dbuff
      integer i,ipg,iref

      write(sname,"(I3.3)")ga_nodeid()
      write(siter,"(I3.3)")iter
      write(sref,"(I3.3)")iref

      if(lusesub) then
      ipg =int_mb(k_innodes+ga_nnodes()+ga_nodeid())
      else
      ipg = 1
      endif

      write(sipg,"(I3.3)")ipg

      open(unit=20+ga_nodeid(),file='/mscf/home/brab894/JOBS/TESTS/H2O/'
     1 //fname//sname//siter//sipg//sref//'.file',status='unknown')

      do i=1,size_a
        call ga_get(d_a,i,i,1,1,dbuff,1)
        write(20+ga_nodeid(),"(F16.12)")dbuff
      enddo

      close(20+ga_nodeid())

      return
      end

c         subroutine tce_mrcc_denomstats(iter,iref)
c         implicit none
c#include "tce.fh"
c#include "mafdecls.fh"
c#include "stdio.fh"
c#include "rtdb.fh"
c#include "errquit.fh"
c#include "sym.fh"
c#include "tce_mrcc.fh"
c#include "global.fh"
c#include "tce_main.fh"
c 
c         integer iter,iref,nodezero
c
c         nodezero = (ga_nodeid().eq.0)
c
c         return
c         end
c $Id$
