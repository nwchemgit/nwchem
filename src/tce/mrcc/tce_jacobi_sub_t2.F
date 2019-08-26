      subroutine tce_jacobi_sub_t2(d_r2,d_t2,k_t2_offset,iter,iref)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_mrcc.fh"
      integer d_r2
      integer d_t2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_t2_offset
      integer size
      integer l_r2,k_r2
      integer l_t2,k_t2
      integer i
      integer nprocs
      integer count
      integer nex
cc      integer nxtask
c
      integer iter
      integer iref
      double precision denom 
c
cc      external nxtask
      INTEGER NXTASKsub
      EXTERNAL NXTASKsub
      logical nodezero
      logical noloadbalance
c *** shift ***
      double precision shift
c *************
      if(iter.le.100) then
       shift=-(2.0d0)*zlshift
      else
       shift=0.0d0
      end if
c
c     ================
c     Loop over blocks
c     ================
c
      nodezero = (ga_nodeid().eq.0)
c      noloadbalance = ((ioalg.eq.4).or. 
c     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs=GA_pgroup_NNODES(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      count = 0
      nex=NXTASKsub(nprocs,1,int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
c      call GA_Pgroup_SYNC(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (nex.eq.count) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)

                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_t2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
c                if(lsubterm) then
                  if (.not.ma_push_get(mt_dbl,size,'t2',l_t2,k_t2))
     1              call errquit('tce_jacobi_t2: MA problem',0,MA_ERR)
                  call get_hash_block(d_t2,dbl_mb(k_t2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
c                endif

                  i = 0
c *** shift added ***
ccc               shift=0.50d0
c *******************
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1

              denom = (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1)
     5              +mrccshift)

       if((abs(denom).lt.0.01d0).and.
     1 (abs(dbl_mb(k_r2+i-1)/denom).gt.1.0d0))
     1  then
       write(6,"('2DENOM CLOSE TO ZERO: ',F16.10,F16.10,F16.10)")
     1 denom,dbl_mb(k_t2+i-1),dbl_mb(k_r2+i-1)
c       if(denom.lt.0.0d0) then
c        shift=-1.0d10 
c       else 
c        shift=1.0d10
c       endif
       endif

              if(.not.lsubterm) then

         if(forcedegen) then

                          dbl_mb(k_r2+i-1) = (dbl_mb(k_r2+i-1))
     1           / (orbdegenenergy(iref,h3b)
     1             +orbdegenenergy(iref,h4b)
     1             -orbdegenenergy(iref,p1b)
     1             -orbdegenenergy(iref,p2b)+mrccshift)
         else
                          dbl_mb(k_r2+i-1) = (dbl_mb(k_r2+i-1))
     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1)
     5              +shift+mrccshift)
         endif
              else
  
         if(forcedegen) then

                          dbl_mb(k_r2+i-1) = (dbl_mb(k_r2+i-1)
     1             -mrccshift*dbl_mb(k_t2+i-1))
     1           / (orbdegenenergy(iref,h3b)
     1             +orbdegenenergy(iref,h4b)
     1             -orbdegenenergy(iref,p1b)
     1             -orbdegenenergy(iref,p2b)+mrccshift)

        else

                          dbl_mb(k_r2+i-1) = (dbl_mb(k_r2+i-1)
     1             -mrccshift*dbl_mb(k_t2+i-1))
     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1)
     5              +shift+mrccshift)

c                          dbl_mb(k_r2+i-1) = ((dbl_mb(k_r2+i-1)
c     1             -mrccshift*dbl_mb(k_t2+i-1))*denom)
c     1           / (denom*denom+shift*shift)
          endif
              endif
c       shift=0.0d0

c       if(abs(dbl_mb(k_r2+i-1)).gt.20.0d0) then
c           if(dbl_mb(k_r2+i-1).lt.0.0d0) then
c               dbl_mb(k_r2+i-1)=-20.0d0
c           else
c               dbl_mb(k_r2+i-1)=20.0d0
c           endif
c           write(6,"('2RESIDUE HAS BEEN CUTED')")
c       endif
      if(iter .lt. 4) then
       if(.not. lreadt) then
        if(abs(dbl_mb(k_r2+i-1)).gt.0.1d0) then
         if(dbl_mb(k_r2+i-1).lt.0.0d0) then
               dbl_mb(k_r2+i-1)=-0.01d0
         else
               dbl_mb(k_r2+i-1)=0.01d0
         endif
c        if(nodezero) 
c          write(luout,*)"t2 residue has been modified to 0.01"
c        if(nodezero) call util_flush(LuOut)
        endif
       endif
      endif




                        enddo
                      enddo
                    enddo
                  enddo
                  call add_hash_block(d_t2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
c update of the res.-double vector to the form of increment used in DIIS proc.
c                  call put_hash_block(d_r2,dbl_mb(k_r2),size,
c     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
c     2              *noab+h3b-1)*noab+h4b-1))
c ----------------------------------------------------------------------------
                  if (nodezero.and.util_print('t2',print_debug)) 
     1              then
                    call get_hash_block(d_t2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                    call ma_print_compact
     2              (dbl_mb(k_r2),size,1,'t2')
                  endif
c                 if(lsubterm) then
                  if (.not.ma_pop_stack(l_t2))
     1              call errquit('tce_jacobi_t2: MA problem',1,MA_ERR)
c                 endif
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_t2: MA problem',1,MA_ERR)
                endif
                endif
                endif
      nex=NXTASKsub(nprocs,1,int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      nex=NXTASKsub(-nprocs,1,int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      call GA_Pgroup_SYNC(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      call util_flush(LuOut)
      return
      end
c
c
c
c
c $Id$
