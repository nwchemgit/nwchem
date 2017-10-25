      subroutine dimqm_f2d(rtdb, fld, muind, muold, xyz, icmplx, dir,
     $                     pm, doefield)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c     Input variables
      integer rtdb
      integer icmplx
      double precision fld(3, nDIM, icmplx)
      double precision muind(3, nDIM, icmplx)
      double precision muold(3, nDIM, icmplx)
      double precision xyz(3, nDIM)
      character dir, pm
      logical doefield
c
c     Common variables used
c      logical lcpim, ldebug
c      integer alg, nDIM, nOrder
c      double precision totcharge
c
c     Local variables
      integer i, i3, i0, id, iAtom, iOrder
      double precision tempArr(nOrder, icmplx)
      logical stat
      logical iscmplx
c
c     Initallize
c      write(luout,*) "i3: ", 3*nDIM
c      write(luout,*) "nORDER: ", nOrder
      id = ga_nodeid()
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "Begin dimqm_f2d"
      i3 = 3 * nDIM
      iscmplx = .false.
      if(icmplx == 2) iscmplx = .true.
c
c     Assemble b array
      tempArr = ZERO
c     Turn 3XNXicmplx matrix into 3Nxicmplx vector
      do i = 1, icmplx
        do iAtom = 1, nDIM
          i0 = 3 * (iAtom - 1)
          tempArr(i0+1,i) = fld(1,iAtom,i)
          tempArr(i0+2,i) = fld(2,iAtom,i)
          tempArr(i0+3,i) = fld(3,iAtom,i)
        end do
      end do

c
c     If local fields, include the external field.  Real only.
c jbecca: START. We are missing contributions to the dipoles from the 
c        external applied electric field. This needs to be done carefully
c        as to not accidentally pull in and add induced dipoles from the
c        external field repeatedly. adding an additional logical variable
c        to check for this.
c     NOTE: Since this has (.not. iscmplx), it will only work for static
c           cases. FD with efield will not and should not work. 
c      if(lefield .and. doefield .and. lclFld .and. .not. iscmplx) then
c         write(luout,*)'finitefieldX',finitefieldX
c         write(luout,*)'finitefieldY',finitefieldY
c         write(luout,*)'finitefieldZ',finitefieldZ
c         do i = 1, icmplx
c            do iAtom = 1, nDIM
c               i0 = 3 * (iAtom - 1)
cc          tempArr(i0+1:i0+3,1) = tempArr(i0+1:i0+3,1) ! + external field
c                  tempArr(i0+1,i) = tempArr(i0+1,i) + finiteFieldX
c                  tempArr(i0+2,i) = tempArr(i0+2,i) + finiteFieldY
c                  tempArr(i0+3,i) = tempArr(i0+3,i) + finiteFieldZ
c            enddo
c         end do
c      end if
c
c  ---------------------------------------------------------------
c  If the equations are complex, turn tempArr into a complex array
c  so that we can pass the data to the complex solvers.
c  ---------------------------------------------------------------
c
c      if(iscmplx) then
c        if(.not.ma_push_get(mt_dcpl, nOrder,'dimqm:cmplxArr',
c     $                      l_cmplxArr, k_cmplxArr))
c     $  call errquit('dimqm_f2d malloc k_cmplxArr failed', 1, MA_ERR)
c        do iOrder = 1, nOrder
c          dcpl_mb(k_cmplxArr+iOrder-1) = DCMPLX(tempArr(iOrder,1),
c     $                                    tempArr(iOrder,2))
c        end do
c      end if
c
c  ---------------------------------------------------------------------------
c  Solve equation a*x=b. Before calling, tempArr is b. After, tempArr is x.
c  Use either a UDU decomposed matrix to solve directly or solve iteratively.
c  ---------------------------------------------------------------------------
c
      if(alg == 0) then
c       Exact Solver.  Might be done another time.
      else
c       Iterative solver!
        call dimqm_iterative(rtdb, tempArr(:,:), muold, xyz, icmplx,
     $                       iscmplx, dir, pm)
      end if
c      write(luout,*) "Dipoles: "
c      write(luout,*) tempArr
c
c  ------------------------------------------------------
c  Collect dipoles from x vector.
c  Neglect a contribution if explicitly not including it.
c  ------------------------------------------------------
c
      call dfill(i3*icmplx, ZERO, muind, 1)
c
c     Turn from a [3N x icmplx] vector back to a [3 x N xicmplx] matrix
       do i = 1, icmplx
         do iAtom = 1, nDIM
c            write(luout,*) iAtom
            i0 = 3 * ( iAtom - 1 )
c            write(luout,*) i0
            muind(1,iAtom, i) = tempArr(i0+1,i)
c            muind(1,iAtom) = ONE
c            write(luout,*) muind(1,iAtom)
c            stat = MA_verify_allocator_stuff()
            muind(2,iAtom,i) = tempArr(i0+2,i)
c            write(luout,*) muind(2,iAtom)
c            stat = MA_verify_allocator_stuff()
            muind(3,iAtom,i) = tempArr(i0+3,i)
c            write(luout,*) muind(3,iAtom)
c            stat = MA_verify_allocator_stuff()
         end do
       end do
c      end if

      

      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "End dimqm_f2d"
      end subroutine dimqm_f2d
