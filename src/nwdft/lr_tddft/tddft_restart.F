      subroutine tddft_restart(ipol,nroots,ntrials,nov,algorithm,
     1  g_trials,g_buffer1,g_buffer2,d_trials,fn_trials,tda)
c
c Read trial vectors from disk
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "util.fh"
#include "stdio.fh"
#include "dra.fh"
c
      integer ipol             ! = 1 (RDFT); =2 (UDFT)
      integer nroots           ! Number of roots sought
      integer ntrials          ! Number trial vectors
      integer nov(2)           ! Number of occupied virtual pairs
      integer algorithm        ! =1 multiple fock build; =2: conventional
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_buffer1(2)     ! GA handle for buffer vectors
      integer g_buffer2(2)     ! GA handle for buffer vectors
      integer d_trials(2)      ! DRA handle for trial vectors
      character(*) fn_trials   ! File name for trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
c
      integer i,j,k,l,m,n
      integer l_trsl, k_trsl
      integer nlen(2)
      integer ios
      integer request
      integer dummy
      integer inntsize
      integer shift
      integer ntrials_old
      integer ipol_verify
      integer nov_verify(2)
c
      double precision r
      double precision tga_ddot
      character*32 pname
c
c     IO UNIT NUMBER
      integer ioun
      parameter (ioun=3258)
c
c     ==================================================================
c
      pname='tddft_restart: '
c     
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
c
      nlen(1)=max(1,nov(1))
      nlen(2)=max(1,nov(2))
c
      if (ga_nodeid().eq.0) then
c
c     Open output file.
        open(ioun,status='old',form='unformatted',file=fn_trials,
     1    iostat=ios)
      end if
c
      call ga_brdcst(Msg_Vec_NBF+MSGINT, ios, inntsize, 0)
      if (ios.ne.0) then 
        if (ga_nodeid().eq.0)
     1    write(LuOut,*) 'Error opening trials, using default guess'
        return
      end if
c
      if (ga_nodeid().eq.0) then
c
c     Read how many trial vectors we have (TDA: nroots; RPA: 2*nroots)
        read(ioun) ntrials_old
        read(ioun) ipol_verify
        read(ioun) nov_verify(1), nov_verify(2)
      end if
      call ga_brdcst(Msg_Vec_stat+MSGINT, ntrials_old, inntsize, 0)
      call ga_brdcst(Msg_Vec_stat+MSGINT, ipol_verify, inntsize, 0)
      call ga_brdcst(Msg_Vec_stat+MSGINT, nov_verify, 2*inntsize, 0)
      if ((ipol_verify.ne.ipol).or.(nov_verify(1).ne.nov(1)).or.
     1  (nov_verify(2).ne.nov(2))) then
        if (ga_nodeid().eq.0)
     1    write(LuOut,*) 'ipol/nov inconsistency, using default guess'
        return
      end if
c
c     Check if we are increasing number of roots from previous
c      calculation. If so, overwrite only the unit vector guesses
c      for which we already have solutions and save the rest.
c      For RPA this means we have to shift nroots-(ntrials_old/2)
c      vectors to the "right" in g_trials.
c      For TDA we do not have to shift anything.
      shift=0
      if (.not.tda) then
        ntrials=ntrials_old
c       ntrial_old has to be even because we store 2*nroots vectors for
c        RPA calculations
        if (nroots.gt.ntrials_old/2) then
          shift=ntrials_old/2
          do i=1,ipol
            do j=nroots,(ntrials_old/2)+1,-1
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                call ga_copy_patch('n',g_trials(i),1,nlen(i),j,j,
     1            g_trials(i),1,nlen(i),j+shift,j+shift)
              else
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nlen(i),1,1,d_trials(i),1,nlen(i),j,j,
     2            request)
                dummy=dra_wait(request)
                dummy=dra_write_section(.false.,g_buffer1(i),1,nlen(i),
     1            1,1,d_trials(i),1,nlen(i),j+shift,j+shift,request)
                dummy=dra_wait(request)
              end if
            end do ! nroots, -1
          end do ! ipol
c         New number of trial vectors
          ntrials=nroots+shift
        end if
      end if
c
      do i=1,ipol
c
c     Allocate memory
        if (.not.ma_push_get(mt_dbl,nlen(i),"mat slice",l_trsl,k_trsl))
     1      call errquit(trim(pname)//"failed to alloc mat slice",0,0)
c
c     Read trials vectors from file and store
       do j=1,ntrials_old
         if (ga_nodeid().eq.0) call sread(ioun,dbl_mb(k_trsl),nlen(i))
         if ((algorithm.eq.1).or.(algorithm.eq.2)) then
           if (ga_nodeid().eq.0)
     1       call ga_put(g_trials(i),1,nlen(i),j,j,dbl_mb(k_trsl),
     2         nlen(i))
           call ga_sync()
         else
           if (ga_nodeid().eq.0)
     1       call ga_put(g_buffer1(i),1,nlen(i),1,1,dbl_mb(k_trsl),
     2         nlen(i))
           call ga_sync()
           dummy=dra_write_section(.false.,g_buffer1(i),1,nlen(i),1,1,
     1       d_trials(i),1,nlen(i),j,j,request)
           dummy=dra_wait(request)
         end if
       end do ! ntrials_old
c
c     Deallocate memory
       if (.not.ma_pop_stack(l_trsl))
     1   call errquit(trim(pname)//"failed to pop stack",0,0)
      end do ! ipol
c
c     Close output file
      if (ga_nodeid().eq.0) close(ioun)
c
c     Gram-Schmidt Process
      do n=1,ntrials
c
c     Normalize
        r=0.d0
        do i=1,ipol
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            r=r+ga_ddot_patch(g_trials(i),'n',1,nlen(i),n,n,
     1        g_trials(i),'n',1,nlen(i),n,n)
          else
            dummy=dra_read_section(.false.,g_buffer1(i),1,nlen(i),1,1,
     1        d_trials(i),1,nlen(i),n,n,request)
            dummy=dra_wait(request)
            r=r+ga_ddot_patch(g_buffer1(i),'n',1,nlen(i),1,1,
     1        g_buffer1(i),'n',1,nlen(i),1,1)
          end if
        end do ! ipol
        if (r.lt.1.d-15) 
     1    call errquit(trim(pname)//'unable to normalize',0,0)
        r=1.d0/dsqrt(r)
        do i=1,ipol
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            call ga_scale_patch(g_trials(i),1,nlen(i),n,n,r)
          else
            call ga_scale_patch(g_buffer1(i),1,nlen(i),1,1,r)
            dummy=dra_write_section(.false.,g_buffer1(i),1,nlen(i),1,1,
     1        d_trials(i),1,nlen(i),n,n,request)
            dummy=dra_wait(request)
          end if
        end do ! ipol
c
c     Orthogonalize
        do m=n+1,ntrials
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            r=0.d0
            do i=1,ipol
              r=r+ga_ddot_patch(g_trials(i),'n',1,nlen(i),n,n,
     1          g_trials(i),'N',1,nlen(i),m,m)
            end do
            do i=1,ipol
              call ga_add_patch(1.d0,g_trials(i),1,nlen(i),m,m,-r,
     1          g_trials(i),1,nlen(i),n,n,g_trials(i),1,nlen(i),m,m)
            end do
          else
            r=0.d0
            do i=1,ipol
              dummy=dra_read_section(.false.,g_buffer1(i),1,nlen(i),1,1,
     1          d_trials(i),1,nlen(i),n,n,request)
              dummy=dra_wait(request)
              dummy=dra_read_section(.false.,g_buffer2(i),1,nlen(i),1,1,
     1          d_trials(i),1,nlen(i),m,m,request)
              dummy=dra_wait(request)
              r=r+ga_ddot_patch(g_buffer1(i),'n',1,nlen(i),1,1,
     1          g_buffer2(i),'n',1,nlen(i),1,1)
            end do
            do i=1,ipol
              call ga_add_patch(1.d0,g_buffer2(i),1,nlen(i),1,1,-r,
     1          g_buffer1(i),1,nlen(i),1,1,g_buffer2(i),1,nlen(i),1,1)
              dummy=dra_write_section(.false.,g_buffer2(i),1,nlen(i),1,
     1          1,d_trials(i),1,nlen(i),m,m,request)
              dummy=dra_wait(request)
            end do
          end if
        end do ! ntrials (m)
      end do ! ntrials (n)
c
      return
      end
