* $Id$
*


*     ***********************************
*     *					*
*     *	        Dnex_init		*	
*     *					*
*     ***********************************

      subroutine Dnex_init(ispin_in,ne_in,map_in)
      implicit none
      integer ispin_in
      integer ne_in(2)
      integer map_in

#include "Dnex.fh"

      ispin   = ispin_in
      ne(1)   = ne_in(1)
      ne(2)   = ne_in(2)

      neq(1) = 0
      neq(2) = 0
      call D1dB_init(3,ne(1),map_in)
      call D1dB_nq(3,neq(1))
      if (ispin.eq.2) then
        call D1dB_init(4,ne(2),map_in)
        call D1dB_nq(4,neq(2))
      end if

      return
      end


*     ***********************************
*     *					*
*     *	          Dnex_end   		*	
*     *					*
*     ***********************************

      subroutine Dnex_end()
      implicit none

#include "Dnex.fh"

      call D1dB_end(3)
      if (ispin.eq.2) call D1dB_end(4)

      return
      end
         

*     ***********************************
*     *					*
*     *	         Dnexall_ntoqp 		*	
*     *					*
*     ***********************************

      subroutine Dnexall_ntoqp(n,q,p)      
      implicit none
      integer n,q,p

#include "Dnex.fh"

      if (n.le.ne(1)) then
        call D1dB_ktoqp(3,n,q,p)
      else
        call D1dB_ktoqp(4,n,q,p)
        q = q + neq(1)
      end if

      return
      end



*     ***********************************
*     *					*
*     *	        Dnexall_neq		*	
*     *					*
*     ***********************************

      subroutine Dnexall_neq(nqtmp)
      implicit none
      integer nqtmp(2)

#include "Dnex.fh"

      nqtmp(1) = neq(1)
      nqtmp(2) = neq(2)
      return 
      end


c     ****************************************
c     *                                      *
c     *        Dnexall_m_diagonalize         *
c     *                                      *
c     ****************************************

*   This routine diagonalizes hml = v*eig*v'.  On exit hml is replaced by v.
*   ---- Does not use Dnexall common block! 

      subroutine Dnexall_m_diagonalize(mb,ispin,ne,hml,eig,assending)
      implicit none
      integer    mb,ispin,ne(2)
      real*8     hml(*),eig(*)
      logical    assending

#include "bafdecls.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer MASTER
      parameter (MASTER=0)
      logical value
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,ierr
      integer tmp1(2)


      call nwpw_timing_start(17)

         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
            call dcopy((ne(1)+ne(2)),0.0d0,0,eig,1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            call dcopy(ne(mb),0.0d0,0,eig,1)
         end if

         value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),
     >                       'tmp1',tmp1(2),tmp1(1))
         if (.not. value)
     >      call errquit('Dnexall_m_diagonalize:out of stack',
     >                    1,MA_ERR)

*        ***** diagonalize the matrix *****
         do ms=ms1,ms2
            shift1 = 1+(ms-1)*ishift1
            shift2 = 1+(ms-1)*ishift2
            if (ne(ms).le.0) go to 30

!$OMP MASTER
            ierr = 0
            call DSYEV('V','U',ne(ms),
     >                 hml(shift2),ne(ms),
     >                 eig(shift1),
     >                 dbl_mb(tmp1(1)),2*ne(1)*ne(1),
     >                 ierr)

            if (.not.assending)
     >         call EIGSRT(eig(shift1),
     >                   hml(shift2),
     >                   ne(ms),ne(ms))
!$OMP END MASTER
!$OMP BARRIER

  30        continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,ne(1)*ne(1)+ne(2)*ne(2),
     >                                  hml)
         else
            call Parallel_Brdcst_values(MASTER,ne(mb)*ne(mb),hml)
         end if

         value = BA_pop_stack(tmp1(2))
         if (.not. value)
     >    call errquit('error popping stack in Dnexall_m_diagonalize',
     >                 0,MA_ERR)


      call nwpw_timing_end(17)
      return
      end



c     ****************************************
c     *                                      *
c     *        Dnexall_fmf_Multiply          *
c     *                                      *
c     ****************************************

*  uses rotation algorithm

      subroutine Dnexall_fmf_Multiply(mb,ispin,ne,Ain,npack1,
     >                                   hml,alpha,
     >                                   Aout,beta)
      implicit none
      integer    mb,ispin,ne(2)
      complex*16 Ain(*)
      integer    npack1
      real*8     hml(*)
      real*8     alpha
      complex*16 Aout(*)
      real*8     beta

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,shift3,ishift2,ishift3

      call nwpw_timing_start(16)

         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift  = 1 + (ms-1)*ne(1)*npack1
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM_OMP('N','N',2*npack1,n,n,
     >                (alpha),
     >                Ain(shift),  2*npack1,
     >                hml(shift2),    n,
     >                (beta),
     >                Aout(shift),2*npack1)
   30       continue
         end do

      call nwpw_timing_end(16)
      return
      end







