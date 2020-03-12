c
c     Load vectors (X+Y) and (X-Y) resulting from the TDDFT
c     calculations. Load only the vectors corresponding to states
c     for which the gradients are to be evaluated.
c     HvD 5/2006, NG 11/2012
c
      subroutine tddft_grad_load_xy(tda,ipol,naoc,nav,nov,
     +           nroots,iroots,fn_civecs,
     +           nsroots,isroots,fn_civecs_singlet,
     +           ntroots,itroots,fn_civecs_triplet,
     +           g_xpy,g_xmy,Eroot)
c
c     The workings of this routine need to be aligned with:
c     - tddft_analysis: which stores the vectors
c     - tddft_iguess:  which restores the vectors in case of restarts
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "tddft_grad_util.fh"
#include "stdio.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "mafdecls.fh"
c
c     Input:
c
      logical tda              ! .true. if the Tamm-Dancoff 
                               ! approximation applies
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer naoc(2)          ! Number of active occupied orbitals
      integer nav(2)           ! Number of active virtual orbitals
      integer nov(2)           ! Number of occupied-virtual pairs
      integer nroots           ! The number of unrestricted roots
      integer iroots(nroots)   ! The unrestricted root numbers
      integer nsroots          ! The number of singlet roots
      integer isroots(nsroots) ! The singlet root numbers
      integer ntroots          ! The number of triplet roots
      integer itroots(ntroots) ! The triplet root numbers
c
      character*(*) fn_civecs  ! The unrestricted vectors file
      character*(*) fn_civecs_singlet ! Singlet vectors file
      character*(*) fn_civecs_triplet ! Triplet vectors file
c
c     Output:
c
      integer g_xpy(2)         ! The (X+Y) vectors
      integer g_xmy(2)         ! The (X-Y) vectors
c
      double precision Eroot(nroots+nsroots+ntroots) ! The root energies
c
c     Local:
c
      logical vector_exists ! .true. if the file exists
      logical tda_verify    ! to check consistency of tda
      integer ipol_verify   ! to check consistency of ipol
      integer nroots_verify
      integer nocc_verify(2)
      integer nmo_verify(2)
      integer nfc_verify(2)
      integer nfv_verify(2)
      integer nov_verify(2)
      integer nrt           ! the number of roots in the current file
      integer i, j, k, l, m ! counters
      double precision r    ! single word buffer for data transfer
      integer ilo(3),ihi(3) ! to specify the patch to set
      integer ld(3)         ! to specify leading dimensions
      integer nrts          ! total number of roots
      logical oroot
c
      integer itrial,ii
      integer l_trials, k_trials
      character*32 pname
c
      integer inntsize,ddblsize,ok
c
c Daniel (2-19-13): Variables for checking the <S**2> expectation value.
      integer ihdl_s2saved, iptr_s2saved
      integer tot_nroots
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
c
      pname="tddft_grad_load_xy: "    
c
c     Code:
c
      ld(1) = 1
      ld(2) = 1
      ld(3) = 0
      oroot = ga_nodeid().eq.0
c
c Daniel (2-19-13): Code for testing <S**2> expectation values.
      tot_nroots = nroots + nsroots + ntroots
      if (.not.ma_push_get(mt_dbl,tot_nroots,'s2saved',
     1                     ihdl_s2saved,iptr_s2saved))
     2  call errquit(pname//'failed to allocate s2saved', 0,
     3    MA_ERR)
c
      call ga_sync()
c
      if (oroot) then
        if (ipol.eq.1) then
c
c         Singlet states
c
          if (nsroots.gt.0) then
c
c           Restore singlet state ci vectors from file
c
            inquire(file=fn_civecs_singlet,exist=vector_exists)
            if (vector_exists) then
              open(unit=69,file=fn_civecs_singlet,form='unformatted',
     +             status='unknown',err=90000)
              rewind(69)
              read(69,err=90000) tda_verify
              read(69,err=90000) ipol_verify
              read(69,err=90000) nroots_verify
              read(69,err=90000) nocc_verify(1),nocc_verify(2)
              read(69,err=90000) nmo_verify(1),nmo_verify(2)
              read(69,err=90000) nfc_verify(1),nfc_verify(2)
              read(69,err=90000) nfv_verify(1),nfv_verify(2)
              read(69,err=90000) nov_verify(1),nov_verify(2)
              read(69,err=90000)
c
              nrt = nroots_verify
c
              if ((tda.eqv.tda_verify).and.(ipol.eq.ipol_verify).and.
     +         (nov(1).eq.nov_verify(1)).and.(nov(2).eq.nov_verify(2))) 
     +             then
c Daniel (2-19-13): This part searches for the largest singlet root a
c gradient is requested for.  This prevents us from doing a lot of work
c by searching through all roots if, for instance, the excitations 
c calculation uses 10 roots, but the gradient is requested for root 2.  
c If no roots are requested, the calculation will die here.
                j = 0
                do i = 1, nsroots
                 j = max(j,isroots(i))
                 if (isroots(i).lt.1.or.isroots(i).gt.nrt) 
     +            call errquit(pname//'singlet root out of range',0,
     +                  INPUT_ERR)
                enddo
                nrt = j
c
c               Now for the actual reading of the roots
c
                do i = 1, nrt

c                 Allocate memory
                  if (.not.ma_push_get(mt_dbl,naoc(1)*nav(1),"slice",
     &                l_trials, 
     &                k_trials))call errquit(trim(pname)//"slice",0,0)


                  k = 0
c Daniel (2-19-13): Here we need to account for reordering of the 
c excited states.  The <S**2> value is irrelevant for singlet roots,
c since it is always 0.
                  do j = 1, nsroots
                    if (i.eq.isroots(j)) k = j
                  enddo
                  read(69,err=90000) r  ! energy line
                  if (k.gt.0) Eroot(k) = r
                  read(69,err=90000) r  ! s2 line
c
c                 Read in a chunk of the trials
                  call sread(69,dbl_mb(k_trials),naoc(1)*nav(1))
                  itrial = 0
                  if (k.gt.0) then
                  do l = 1, naoc(1) 
                    do m = 1, nav(1)
                       itrial = itrial+1 
                       r = dbl_mb(k_trials+itrial-1)
                        ilo(1) = k
                        ilo(2) = l
                        ilo(3) = m
                        ihi(1) = k
                        ihi(2) = l
                        ihi(3) = m
                        call nga_put(g_xpy(1),ilo,ihi,r,ld)
                    enddo
                  enddo
               endif
c
                  if (.not.tda) then
c
c                   Read in a chunk of the trials
                    call sread(69,dbl_mb(k_trials),naoc(1)*nav(1))
                    itrial = 0
                    if (k.gt.0) then
                    do l = 1, naoc(1) 
                      do m = 1, nav(1)
                        itrial = itrial+1 
                        r = dbl_mb(k_trials+itrial-1)
                          ilo(1) = k
                          ilo(2) = l
                          ilo(3) = m
                          ihi(1) = k
                          ihi(2) = l
                          ihi(3) = m
                          call nga_put(g_xmy(1),ilo,ihi,r,ld)
                      enddo
                    enddo
                 endif
                  endif   ! .not.tda
c
c                 Deallocate memory
                  if (.not.ma_pop_stack(l_trials))
     &             call errquit(trim(pname)//"failed to pop stack",0,0)
c
                enddo   ! nrt
              else
                call errquit(pname//'vectors file corrupt',0,DISK_ERR)
              endif
              close(69)
            else
              call errquit(pname//'singlet vectors absent',0,DISK_ERR)
            endif
          endif
c
c         Triplet states
c
          if (ntroots.gt.0) then
c
c           Restore triplet state ci vectors from file
c
            inquire(file=fn_civecs_triplet,exist=vector_exists)
            if (vector_exists) then
              open(unit=69,file=fn_civecs_triplet,form='unformatted',
     +             status='unknown',err=90010)
              rewind(69)
              read(69,err=90000) tda_verify
              read(69,err=90000) ipol_verify
              read(69,err=90000) nroots_verify
              read(69,err=90000) nocc_verify(1),nocc_verify(2)
              read(69,err=90000) nmo_verify(1),nmo_verify(2)
              read(69,err=90000) nfc_verify(1),nfc_verify(2)
              read(69,err=90000) nfv_verify(1),nfv_verify(2)
              read(69,err=90000) nov_verify(1),nov_verify(2)
              read(69,err=90000)
c
              nrt = nroots_verify
c
              if ((tda.eqv.tda_verify).and.
     +            (ipol.eq.ipol_verify).and.
     +            (nov(1).eq.nov_verify(1)).and.
     +            (nov(2).eq.nov_verify(2))) then
                j = 0
                do i = 1, ntroots
                  j = max(j,itroots(i))
                  if (itroots(i).lt.1.or.itroots(i).gt.nrt) 
     +             call errquit(pname//'triplet root out of range',0,
     +                      INPUT_ERR)
                enddo
                nrt = j
c
c               Now for the actual reading of the roots
c
                do i = 1, nrt
c
c                 Allocate memory
                  if (.not.ma_push_get(mt_dbl,naoc(1)*nav(1),"slice",
     &                l_trials, 
     &                k_trials))call errquit(trim(pname)//"slice",0,0)
c
                  k = 0
c Daniel (2-19-13): Here we need to account for reordering of the 
c excited states.  The <S**2> value is irrelevant for triplet roots,
c since it is always 2.
                  do j = 1, ntroots
                    if (i.eq.itroots(j)) k = j+nsroots
                  enddo
c
                  read(69,err=90000) r  ! energy line
                  if (k.gt.0) Eroot(k) = r
                  read(69,err=90000) r   ! s2 line
c
c                 Read in a chunk of the trials
                  call sread(69,dbl_mb(k_trials),naoc(1)*nav(1))
                  itrial = 0
                  if (k.gt.0) then
                  do l = 1, naoc(1) 
                    do m = 1, nav(1) 
                      itrial = itrial+1 
                      r = dbl_mb(k_trials+itrial-1)
                        ilo(1) = k
                        ilo(2) = l
                        ilo(3) = m
                        ihi(1) = k
                        ihi(2) = l
                        ihi(3) = m
                        call nga_put(g_xpy(1),ilo,ihi,r,ld)
                    enddo
                  enddo
                  endif
c
                  if (.not.tda) then
c                   Read in a chunk of the trials
                    call sread(69,dbl_mb(k_trials),naoc(1)*nav(1))
                    itrial = 0
                    if (k.gt.0) then
                    do l = 1, naoc(1) 
                      do m = 1, nav(1) 
                        itrial = itrial+1 
                        r = dbl_mb(k_trials+itrial-1)
                          ilo(1) = k
                          ilo(2) = l
                          ilo(3) = m
                          ihi(1) = k
                          ihi(2) = l
                          ihi(3) = m
                          call nga_put(g_xmy(1),ilo,ihi,r,ld)
                      enddo
                    enddo
                    endif
                  endif    ! .not. tda
c
c                 Deallocate memory
                  if (.not.ma_pop_stack(l_trials))
     &             call errquit(trim(pname)//"failed to pop stack",0,0)

                enddo  ! nrt

              else
                call errquit(pname//'vectors file corrupt',0,DISK_ERR)
              endif
              close(69)
            else
              call errquit(pname//'triplet vectors absent',0,DISK_ERR)
            endif
          endif
c
        else if (ipol.eq.2) then
c
c         Unrestricted
c
          if (nroots.gt.0) then
c
c           Restore unrestricted state ci vectors from file
c
            inquire(file=fn_civecs,exist=vector_exists)
            if (vector_exists) then
              open(unit=69,file=fn_civecs,form='unformatted',
     &           status='unknown',err=90020)
              rewind(69)
              read(69,err=90000) tda_verify
              read(69,err=90000) ipol_verify
              read(69,err=90000) nroots_verify
              read(69,err=90000) nocc_verify(1),nocc_verify(2)
              read(69,err=90000) nmo_verify(1),nmo_verify(2)
              read(69,err=90000) nfc_verify(1),nfc_verify(2)
              read(69,err=90000) nfv_verify(1),nfv_verify(2)
              read(69,err=90000) nov_verify(1),nov_verify(2)
              read(69,err=90000)
c
              nrt = nroots_verify
c
              if ((tda.eqv.tda_verify).and.
     +            (ipol.eq.ipol_verify).and.
     +            (nov(1).eq.nov_verify(1)).and.
     +            (nov(2).eq.nov_verify(2))) then
                j = 0
                do i = 1, nroots
                  j = max(j,iroots(i))
                  if (iroots(i).lt.1.or.iroots(i).gt.nrt) call errquit
     +            (pname//'unrestricted root out of range',0,INPUT_ERR)
                enddo
                nrt = j
c
c               Now for the actual reading of the roots
c
                do i = 1, nrt
c
                  k = 0
c Daniel (2-19-13): Here we need to account for reordering of the 
c excited states and also the <S**2> value.
                  do j = 1, nroots
                    if (i.eq.iroots(j)) k = j
                  enddo
c
                  read(69,err=90000) r   ! energy line
                  if (k.gt.0) Eroot(k) = r 
                  read(69,err=90000) r   ! s2 line
c
                  do j = 1, ipol
c
c                   Allocate memory
                    if (.not.ma_push_get(mt_dbl,naoc(j)*nav(j),"slice",
     &                l_trials, 
     &                k_trials))call errquit(trim(pname)//"slice",0,0)
c
c                   Read in a chunk of the trials
                    call sread(69,dbl_mb(k_trials),naoc(j)*nav(j))
                    itrial = 0
                    if (k.gt.0) then
                    do l = 1, naoc(j) 
                      do m = 1, nav(j) 
                        itrial = itrial + 1
                        r = dbl_mb(k_trials+itrial-1)
                          ilo(1) = k
                          ilo(2) = l
                          ilo(3) = m
                          ihi(1) = k
                          ihi(2) = l
                          ihi(3) = m
                          call nga_put(g_xpy(j),ilo,ihi,r,ld)
                      enddo
                    enddo 
                    endif
c
                    if (.not.tda) then
c                     Read in a chunk of the trials
                      call sread(69,dbl_mb(k_trials),naoc(j)*nav(j))
                      itrial = 0
                      if (k.gt.0) then
                      do l = 1, naoc(j)
                        do m = 1, nav(j)
                          itrial = itrial + 1
                          r = dbl_mb(k_trials+itrial-1)
                            ilo(1) = k
                            ilo(2) = l
                            ilo(3) = m
                            ihi(1) = k
                            ihi(2) = l
                            ihi(3) = m
                            call nga_put(g_xmy(j),ilo,ihi,r,ld)
                        enddo
                      enddo
                      endif
                     endif ! .not.tda
c
c                   Deallocate memory
                    if (.not.ma_pop_stack(l_trials))
     &              call errquit(trim(pname)//"failed to pop stack",0,0)
c
                  enddo  ! j = 1,ipol
                enddo  ! nrt
c
              else
                call errquit(pname//'vector file corrupt',0,DISK_ERR)
              endif
              close(69)
            else
              call errquit(pname//'no unrestricted vectors present!?',
     +        0,DISK_ERR)
            endif
          endif
        else
          call errquit(pname//'illegal value for ipol',0,UNKNOWN_ERR)
        endif
      endif ! oroot
c
c     Broadcast status
      ok =1
      call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)
c
c Daniel (2-19-13): Deallocate local arrays.
      if (.not.ma_pop_stack(ihdl_s2saved)) 
     1  call errquit(pname//'failed to deallocate s2saved', 0, MA_ERR)
      return
c
90000 continue
      call errquit(pname//'error reading singlet ci vectors',0,DISK_ERR)
      return
c
90010 continue
      call errquit(pname//'error reading triplet ci vectors',0,DISK_ERR)
      return
c
90020 continue
      call errquit(pname//'error reading unrestricted ci vectors',0,
     &      DISK_ERR)
c
      return
      end
c $Id$
