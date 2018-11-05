      SUBROUTINE tce_read_t(d_a,k_a_offset,d_b,k_b_offset)
C     $Id$
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
c     (JB) Interface for Tailored CCSD method. Read active amplitudes in T12DUMP format from an external source
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce_main.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"

      INTEGER d_a, d_b
      INTEGER k_a_offset, k_b_offset
      INTEGER k_a, l_a, i, j, a, b
      logical nodezero
      integer size

      integer h1b,h2b,p3b,p4b

      character*4 irrep_p3,irrep_p4,irrep_h1,irrep_h2
      character*5 spinnames(2),spin_p3,spin_p4,spin_h1,spin_h2
      data spinnames/'alpha','beta '/

      double precision val

      integer r1,r2,r3,r4
      integer r5,r6,r7,r8

      integer spincase
      integer amp_case

      integer ihash, l_r2, k_r2
      integer h1s,h2s,p3s,p4s

      integer off, ddblsize
      integer ihash_old

      nodezero=ga_nodeid().eq.0

      spincase = 0 ! 0 = aa, 1 = bb 2 = ab 3 = a 4 = b

       call ga_zero(d_a)
       call ga_zero(d_b)

       call ga_sync()

       if(nodezero) then

       OPEN(67,FILE='T12DUMP',STATUS='OLD',action='READ')

       read(67,*)r1,r2,r3,r4,r5,r6,r7,r8

       write(6,"('TCC space for alpha: ',I4,' - ',I4,' virt.: ',I4,
     1 ' - ',I4)")r1,r2,r5,r6
       write(6,"('TCC space for beta: ',I4,' - ',I4,' virt.: ',I4,
     1 ' - ',I4)")r3,r4,r7,r8

       tcc_spaces(1) = r1
       tcc_spaces(2) = r2
       tcc_spaces(3) = r3
       tcc_spaces(4) = r4
       tcc_spaces(5) = r5
       tcc_spaces(6) = r6
       tcc_spaces(7) = r7
       tcc_spaces(8) = r8

       ihash_old = -1

       do while(1.eq.1) 
                read(67,FMT=*,END=200)val,j,i,b,a

              if((i.eq.0).and.(j.eq.0).and.(a.eq.0).and.(b.eq.0)) then
                        spincase = spincase + 1
              else

                if(spincase.eq.0) then
                        call tce_get_block_num(i,1,h1b,h1s)
                        call tce_get_block_num(j,1,h2b,h2s)
                        call tce_get_block_num(a,1,p3b,p3s)
                        call tce_get_block_num(b,1,p4b,p4s)
                elseif(spincase.eq.1) then
                        call tce_get_block_num(i,2,h1b,h1s)
                        call tce_get_block_num(j,2,h2b,h2s)
                        call tce_get_block_num(a,2,p3b,p3s)
                        call tce_get_block_num(b,2,p4b,p4s)
                elseif(spincase.eq.2) then
                        call tce_get_block_num(i,1,h1b,h1s)
                        call tce_get_block_num(j,2,h2b,h2s)
                        call tce_get_block_num(a,1,p3b,p3s)
                        call tce_get_block_num(b,2,p4b,p4s)
                elseif(spincase.eq.3) then
                        call tce_get_block_num(j,1,h1b,h1s)
                        call tce_get_block_num(i,1,h2b,h2s)
                elseif(spincase.eq.4) then
                        call tce_get_block_num(j,2,h1b,h1s)
                        call tce_get_block_num(i,2,h2b,h2s)
                endif

c                write(6,"('Reading: ',F18.12,4I4,' Blocks: ',4I4)")
c     1 val,i,j,a,b,h1b,h2b,p3b,p4b

                if(spincase.lt.3) then

                if (int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)
     1            .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)+
     2            int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p3b-1),ieor(int_mb(k_sym+p4b-1),
     1            ieor(int_mb(k_sym+h1b-1),int_mb(k_sym+h2b-1))))
     2            .eq. 0) then

                  size = int_mb(k_range+h1b-1)*int_mb(k_range+h2b-1)*
     1 int_mb(k_range+p3b-1)*int_mb(k_range+p4b-1)

                  ihash = h2b - 1 + noab * (h1b - 1 + noab * (p4b-
     & noab - 1 + nvab * (p3b - noab - 1)))

                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_t2: MA problem',0,MA_ERR)

                  call dfill(size, 0.0d0, dbl_mb(k_r2), 1)

c                  call get_hash_block(d_a,dbl_mb(k_r2),size,
c     1              int_mb(k_a_offset),ihash)

                  off = h2s + int_mb(k_range+h2b-1)*
     1 (h1s+int_mb(k_range+h1b-1)*(p4s+int_mb(k_range+p4b-1)*p3s))
                  dbl_mb(k_r2+off) = val

                  if(h1b.eq.h2b) then
                    off = h1s + int_mb(k_range+h2b-1)*
     1 (h2s+int_mb(k_range+h1b-1)*(p4s+int_mb(k_range+p4b-1)*p3s))
                    dbl_mb(k_r2+off) = -val
                  endif

                  if(p3b.eq.p4b) then
                    off = h2s + int_mb(k_range+h2b-1)*
     1 (h1s+int_mb(k_range+h1b-1)*(p3s+int_mb(k_range+p4b-1)*p4s))
                    dbl_mb(k_r2+off) = -val
                  endif
                  
                  if((h1b.eq.h2b).and.(p3b.eq.p4b)) then
                    off = h1s + int_mb(k_range+h2b-1)*
     1 (h2s+int_mb(k_range+h1b-1)*(p3s+int_mb(k_range+p4b-1)*p4s))
                    dbl_mb(k_r2+off) = val
                  endif

                  call add_hash_block(d_a,dbl_mb(k_r2),size,
     1              int_mb(k_a_offset),ihash)

                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_t2: MA problem',1,MA_ERR)

                  endif
                  endif
                  endif

                elseif(spincase.ge.3) then

            if (int_mb(k_spin+h2b-1) .eq. int_mb(k_spin+h1b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+h2b-1)
     1        +int_mb(k_spin+h1b-1).ne.4)) then
            if (ieor(int_mb(k_sym+h2b-1),int_mb(k_sym+h1b-1)).eq.0) then

                  size = int_mb(k_range+h1b-1)*int_mb(k_range+h2b-1)
                  ihash = ((h2b-noab-1)*noab+h1b-1)
                  off = h1s + int_mb(k_range+h1b-1)*h2s 

                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_t2: MA problem',0,MA_ERR)
                  call get_hash_block(d_b,dbl_mb(k_r2),size,
     1              int_mb(k_b_offset),ihash)

                    dbl_mb(k_r2+off) = val 

                  if((h1b.eq.h2b)) then
                    write(6,"('T1 symmetric')")
                    off = h2s + int_mb(k_range+h2b-1)*h1s
                    dbl_mb(k_r2+off) = val
                  endif


                  call put_hash_block(d_b,dbl_mb(k_r2),size,
     1              int_mb(k_b_offset),ihash)

                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_t2: MA problem',1,MA_ERR)

                endif
                endif
                endif
                endif

              endif
       enddo

 200   continue

       CLOSE(67,STATUS='KEEP')

       endif

       ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
       call ga_brdcst(Msg_Vec_EVal+403,tcc_spaces,
     1 ddblsize*8, 0)


       call ga_sync()

c       call ga_print(d_a)
c       call ga_print(d_b)

      RETURN
      END



      SUBROUTINE tce_get_block_num(index,spin,blck,shift)
C     $Id$
C     Copyright (c) Battelle & Pacific Northwest National Laboratory
C     (2002)
c     (JB) Auxiliary routine to determine the block number
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce_main.fh"


        integer index, spin, blck
        integer i,j

        integer offset
        integer p1b, idx1

        integer shift, p1

        blck = -1

        DO p1b = 1,noab+nvab
                if(int_mb(k_spin+p1b-1).eq.spin) then

        do p1 = 1, int_mb(k_range+p1b-1)

       idx1 = (int_mb(k_mo_index+int_mb(k_offset+p1b-1)+p1-1)+
     1             2-int_mb(k_spin+p1b-1))/2

cjb         write(6,*)'IDX: ',idx1

        if(index.eq.idx1) then
            blck = p1b
            shift = p1 - 1
            goto 300
        endif 

        enddo

                endif
        enddo

 300    continue

        if((blck.lt.0).and.(index.gt.0))
     1 call errquit('Block was not found',index,CALC_ERR)

        return
        end
