      subroutine tddft_dump(ipol,nocc,nmo,g_trials,
     1  g_prod_p,g_prod_m,ivector,ntrials,nproducts,tda)
c
c $Id: tddft_dump.F,v 1.1 2002-03-22 02:10:36 sohirata Exp $
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
c
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for A+B product vectors
      integer g_prod_m(2)      ! GA handle for A-B product vectors
      integer ivector          ! Vectors to print
      integer ntrials          ! Current number of trial vectors
      integer nproducts        ! Current number of product vectors
      logical tda              ! True if Tamm-Dancoff approximation
c
      double precision t,p,q   ! Work variables
      integer i,j,k,l
c
c -------------
c Print vectors
c -------------
c
      if (ivector.le.ntrials) then
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
          do j=1,nocc(i)
            do k=nocc(i)+1,nmo(i)
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
 9000 format('Vector',i5)
 9010 format(' Spin  Occ  Vir     X+Y')
 9020 format(' Spin  Occ  Vir     X+Y       A+B')
 9030 format(' Spin  Occ  Vir     X+Y       A+B       A-B')
 9040 format(3i5,3f10.5)
c           
c ------
c Return
c ------
c
      return
      end
