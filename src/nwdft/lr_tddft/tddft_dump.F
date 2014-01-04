      subroutine tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     1  g_prod_p,g_prod_m,g_buffer1,d_trials,d_prod_p,d_prod_m,
     2  ivector,ntrials,nproducts,tda,algorithm)
c
c $Id$
c
c Print all current trial and product vectors for debug
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
c
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer nov(2)           ! Number of occupied virtual pairs
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for A+B product vectors
      integer g_prod_m(2)      ! GA handle for A-B product vectors
      integer g_buffer1(2)     ! GA handle for buffer vectors
      integer d_trials(2)      ! DRA handle for trial vectors
      integer d_prod_p(2)      ! DRA handle for A+B product vectors
      integer d_prod_m(2)      ! DRA handle for A-B product vectors
      integer ivector          ! Vectors to print
      integer ntrials          ! Current number of trial vectors
      integer nproducts        ! Current number of product vectors
      logical tda              ! True if Tamm-Dancoff approximation
      integer algorithm        ! Algorithm
c
      double precision t,p,q   ! Work variables
      integer i,j,k,l
      logical nodezero
      integer request
      integer dummy
c
      nodezero=(ga_nodeid().eq.0)
c
c -------------
c Print vectors
c -------------
c
c Determine the length of the trial vector: Occ*Virt
      call tddft_lentrialvec(ipol,nocc,nmo,nfc,nfv,nov)
c
      if (ivector.le.ntrials) then
        if ((algorithm.eq.1).or.(algorithm.eq.2)) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9000) ivector
            if (ivector.le.nproducts) then
              if (.not.tda) then
                write(LuOut,9030)
              else
                write(LuOut,9020)
              endif
            else
              write(LuOut,9010)
            endif
            do i=1,ipol
              l=0
              do j=nfc(i)+1,nocc(i)
                do k=nocc(i)+1,nmo(i)-nfv(i)
                  l=l+1
                  call ga_get(g_trials(i),l,l,ivector,ivector,t,1)
                  if (ivector.le.nproducts) then
                    call ga_get(g_prod_p(i),l,l,ivector,ivector,p,1)
                    if (.not.tda) then
                      call ga_get(g_prod_m(i),l,l,ivector,ivector,q,1)
                      write(LuOut,9040) i,j,k,t,p,q
                    else
                      write(LuOut,9040) i,j,k,t,p
                    endif
                  else
                    write(LuOut,9040) i,j,k,t
                  endif
                enddo
              enddo
            enddo
          endif
        else
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9050) ivector
          endif
          do i=1,ipol
            dummy=dra_read_section(.false.,g_buffer1(i),1,nov(i),1,1,
     1        d_trials(i),1,nov(i),ivector,ivector,request)
            dummy=dra_wait(request)
            if (nodezero) then
              l=0
              do j=nfc(i)+1,nocc(i)
                do k=nocc(i)+1,nmo(i)-nfv(i)
                  l=l+1
                  call ga_get(g_buffer1(i),l,l,1,1,t,1)
                  write(LuOut,9040) i,j,k,t
                enddo
              enddo
            endif
          enddo
          if (ivector.le.nproducts) then
            if (nodezero) write(LuOut,9060) ivector
            do i=1,ipol
              dummy=dra_read_section(.false.,g_buffer1(i),1,nov(i),1,1,
     1          d_prod_p(i),1,nov(i),ivector,ivector,request)
              dummy=dra_wait(request)
              if (nodezero) then
                l=0
                do j=nfc(i)+1,nocc(i)
                  do k=nocc(i)+1,nmo(i)-nfv(i)
                    l=l+1
                    call ga_get(g_buffer1(i),l,l,1,1,p,1)
                    write(LuOut,9040) i,j,k,p
                  enddo
                enddo
              endif
            enddo
            if (tda) then
              if (nodezero) write(LuOut,9070) ivector
              do i=1,ipol
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_prod_m(i),1,nov(i),ivector,ivector,request)
                dummy=dra_wait(request)
                if (nodezero) then
                  l=0
                  do j=nfc(i)+1,nocc(i)
                    do k=nocc(i)+1,nmo(i)-nfv(i)
                      l=l+1
                      call ga_get(g_buffer1(i),l,l,1,1,q,1)
                      write(LuOut,9040) i,j,k,q
                    enddo
                  enddo
                endif
              enddo
            endif
          endif
        endif
      endif
 9000 format('Vector',i5)
 9010 format(' Spin  Occ  Vir     X+Y')
 9020 format(' Spin  Occ  Vir     X+Y       A+B')
 9030 format(' Spin  Occ  Vir     X+Y       A+B       A-B')
 9040 format(3i5,3f10.5)
 9050 format('X+Y vector',i5)
 9060 format('(A+B)(X+Y) vector',i5)
 9070 format('(A-B)(X+Y) vector',i5)
c           
c ------
c Return
c ------
c
      return
      end
