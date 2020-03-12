      subroutine tddft_checkpoint(rtdb,ipol,ispintag,nroots,ntrials,nov,
     1  algorithm,apbvec,ambvec,g_trials,g_buffer,d_trials,tda)
c
c Output current trial vectors so calculation can be restarted
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
      integer rtdb             ! Runtime database handle
      integer ipol             ! = 1 (RDFT); =2 (UDFT)
      integer ispintag         ! 1 = Restricted Singlet, 2 = Restricted Triplet, 3=Unrestricted
      integer nroots           ! Number of roots sought
      integer ntrials          ! Number of trial vectors
      integer nov(2)           ! Number of occupied virtual pairs
      integer algorithm        ! =1 multiple fock build; =2: conventional
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_buffer(2)      ! GA handle for buffer vectors
      integer d_trials(2)      ! DRA handle for trial vectors
      double precision apbvec(ntrials,ntrials) ! (X+Y)
      double precision ambvec(ntrials,ntrials) ! (X-Y)
      logical tda              ! True if Tamm-Dancoff approximation
c
      integer g_x, g_y
      integer i,j,k
      integer m,n
      integer l_trsl, k_trsl
      integer nlen(2)
      integer ios
      integer request
      integer dummy
      integer inntsize
c
      character*32 pname
      character*(nw_max_path_len) fname
c
c     IO UNIT NUMBER
      integer ioun
      parameter (ioun=3258)
c
c     ==================================================================
c
      pname='tddft_checkpoint: '
c     
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
c
      nlen(1)=max(1,nov(1))
      nlen(2)=max(1,nov(2))
c
      if (.not.rtdb_cget(rtdb,'tddft:trials',1,fname))
     1  call util_file_name('trials',.false.,.false.,fname)
c
      if (ga_nodeid().eq.0) then
c
c     Open output file.
c         
        open(ioun,status='replace',form='unformatted',file=fname,
     1    iostat=ios)
      end if
c
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ios, inntsize, 0)
      if (ios.ne.0) then
        if (ga_nodeid().eq.0)
     1    write(luout, *) trim(pname)//"failed to open: "//trim(fname)
        return
      end if
c
      if (.not.rtdb_cput(rtdb,'tddft:trials',1,fname))
     1  call errquit(trim(pname)//'failed writing to rtdb',0,RTDB_ERR)
c
      if (ga_nodeid().eq.0) then
c
        if (tda) then
          write(ioun) nroots
        else
          write(ioun) 2*nroots
        end if
        write(ioun) ipol
        write(ioun) ispintag
        write(ioun) nov(1), nov(2)
      end if
c
      do i=1,ipol
c
c     Allocate memory
        if (.not.ma_push_get(mt_dbl,nlen(i),"mat slice",l_trsl,k_trsl))
     1    call errquit(trim(pname)//"failed to alloc mat slice",0,0)
        if (.not.ga_create(mt_dbl,nlen(i),nroots,'x',-1,nroots,g_x))
     1    call errquit(trim(pname)//'failed to create g_x',0,GA_ERR)
        if (.not.tda) then
          if (.not.ga_create(mt_dbl,nlen(i),nroots,'y',-1,nroots,g_y))
     1      call errquit(trim(pname)//'failed to create g_y',0,GA_ERR)
        end if
c
c     Zero arrays
        call ga_zero(g_x)
        if (.not.tda) then
          call ga_zero(g_y)
        end if
c
c     Construct solution vectors and output to file
        do n=1,nroots
          do m=1,ntrials
            if ((algorithm.eq.1).or.(algorithm.eq.2)) then
              call tga_add(apbvec(m,n),g_trials(i),g_x,m,n,1)
              if (.not.tda) 
     1          call tga_add(ambvec(m,n),g_trials(i),g_y,m,n,1)
            else
              dummy=dra_read_section(.false.,
     1          g_buffer(i),1,nlen(i),1,1,d_trials(i),1,nlen(i),m,m,
     2          request)
              dummy=dra_wait(request)
              call tga_add(apbvec(m,n),g_buffer(i),g_x,1,n,1)
              if (.not.tda)
     1          call tga_add(ambvec(m,n),g_buffer(i),g_y,1,n,1)
            end if 
          end do ! ntrials
          call ga_get(g_x,1,nlen(i),n,n,dbl_mb(k_trsl),nlen(i))
          if (ga_nodeid().eq.0)
     1      call swrite(ioun,dbl_mb(k_trsl),nlen(i))
          if (.not.tda) then
            call ga_get(g_y,1,nlen(i),n,n,dbl_mb(k_trsl),nlen(i))
            if (ga_nodeid().eq.0)
     1        call swrite(ioun,dbl_mb(k_trsl),nlen(i))
          end if
        end do ! nroots
c
c     Deallocate memory
        if (.not.ma_pop_stack(l_trsl))
     $    call errquit(pname//"failed to pop stack",0,0)
        if (.not.ga_destroy(g_x))
     $    call errquit(pname//"failed to destroy g_x",0,0)
        if (.not.tda) then
          if (.not.ga_destroy(g_y))
     $      call errquit(pname//"failed to destroy g_y",0,0)
        end if
      enddo  ! ipol
c
c     Close output file
c
      if (ga_nodeid().eq.0) close(ioun)
c
      return
      end
