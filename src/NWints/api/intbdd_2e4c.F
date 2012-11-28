c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Computes the 2nd derivative integrals of a block of shell quartets
C> with respect to the nuclear positions
C>
C> Compute the 2nd derivative integrals of a block of shell quartets. Only
C> the non-zero integrals and the corresponding labels are returned. This
C> routine cannot be called without a preceding call to `intbdd_init4c` with
C> the same block information.
C>
C> the integrals evaluated are defined as
C> \f{eqnarray*}{
C> ({\mu}{\rho}|{\nu}{\lambda}) &=& \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2} \\\\
C> eri(1,*) &=& \frac{\partial^2 ({\mu}{\rho}|{\nu}{\lambda})}{\partial (X_\mu^{x})^2} \\\\
C> &\ldots& \\\\
C> eri(78,*) &=& \frac{\partial^2 ({\mu}{\rho}|{\nu}{\lambda})}{\partial (X_\lambda^{z})^2} \\\\
C> \f}
C>
C> \return Returns .true. if not all integrals fitted in the buffer
C> provided and an additional call to `intbdd_2e4c` is required, and
C> returns .false. if all integrals were retrieved.
C>
      logical function intbdd_2e4c(brain, icl, jcl, ketin, kcl, lcl,
     &     num_q, q4, use_q4, zerotol, canonicalize,
     &     ilab, jlab, klab, llab, eri,
     &     l_erilab, nint, lscr, scr)
      implicit none
#include "numb_qP.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
c::functions
      logical intbdd_nw_2e4c
      external intbdd_nw_2e4c
c::passed
c:tex-\begin{verbatim}
      integer brain              !< [Input]  basis set handle for bra 
      integer ketin              !< [Input]  basis set handle for ket 
      integer num_q              !< [Input]  number of quartets input
      integer icl(num_q)         !< [Input]  i-contraction quartet labels
      integer jcl(num_q)         !< [Input]  j-contraction quartet labels
      integer kcl(num_q)         !< [Input]  k-contraction quartet labels
      integer lcl(num_q)         !< [Input]  l-contraction quartet labels
      logical use_q4             !< [Input]  logical for use of q4
      double precision q4(num_q) !< [Input] symmetry prefactors for ints
*     
      integer l_erilab           !< [Input]  length of eri and label arrays
      integer ilab(l_erilab)     !< [Output] integral labels for ``i'' 
      integer jlab(l_erilab)     !< [Output] integral labels for ``j'' 
      integer klab(l_erilab)     !< [Output] integral labels for ``k'' 
      integer llab(l_erilab)     !< [Output] integral labels for ``l'' 
      integer nint               !< [Output] number of integrals returned
      double precision eri(78,l_erilab) !< [Output] integral values
      integer lscr               !< [Input] length of scratch array
      double precision scr(lscr) !< [Input] scratch array
      double precision zerotol   !< [Input] zero integral threshold
      logical canonicalize       !< [Input] Canonicalize integral labels?
c:tex-\end{verbatim}
      
c::local
      double precision roff(3)
c
      intbdd_2e4c = .false.
      if (num_q.eq.0) then
        nint = 0
        return
      endif
      if (txs_b2e4c) then
        call texas_hf2_m(
     &      brain, icl(num_q_start), jcl(num_q_start),
     &      ketin, kcl(num_q_start), lcl(num_q_start), 
     &      num_q-num_q_start+1,
     &      q4(num_q_start), use_q4,
     &      roff,roff,roff,roff,.false.,
     &      eri,l_erilab, ilab, jlab, klab, llab, nint, .true.,
     &      more_ints_txs, scr, lscr, zerotol,'der2_int')
        intbdd_2e4c = more_ints_txs
      else
c
c currently, only texas 2 e- 2nd derivative integrals are available.  The rest
c of these mechanics need to be added in if another integral package becomes
c available (i.e. HONDO)
c
c       intbdd_2e4c = intbdd_nw_2e4c(
c    &      brain, icl, jcl, ketin, kcl, lcl,
c    &      num_q, q4, use_q4, zerotol, canonicalize,
c    &      ilab, jlab, klab, llab, eri,
c    &      l_erilab, nint, lscr, scr)
        write(luout,*)'intbdd_2e4c: could not use the texas '
        write(luout,*)'            integral second derivatives'
        write(luout,*)' Please notify nwchem-support@emsl.pnl.gov with'
        write(luout,*)'        the input deck and output available'
        call errquit('intbdd_2e4c: fatal error',911, INT_ERR)
      endif
c
      return
      end
C>
C> \brief Compute the 2nd derivative integrals of a block of shell quartets
C> using the NWChem McMurchie-Davidson code.
C>
C> Compute and return the non-zero 2nd derivative integrals and corresponding
C> labels based on the input of a group of shell quartets (a block). This 
C> routine explicitly uses the McMurchie-Davidson code, it should \b NEVER
C> be called directly by an application module! Please call `intbdd_2e4c`
C> instead.
C>
      logical function intbdd_nw_2e4c(
     &    brain, icl, jcl, ketin, kcl, lcl,
     &    num_q, q4, use_q4, zerotol, canonicalize,
     &     ilab, jlab, klab, llab, eri,
     &     l_erilab, nint, lscr, scr)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "numb_qP.fh"
c::functions
      integer int_nint_cart
      external int_nint_cart
c::   passed
c:tex-\begin{verbatim}
      integer brain              !< [Input]  basis set handle for bra 
      integer ketin              !< [Input]  basis set handle for ket 
      integer num_q              !< [Input]  number of quartets input
      integer icl(num_q)         !< [Input]  i-contraction quartet labels
      integer jcl(num_q)         !< [Input]  j-contraction quartet labels
      integer kcl(num_q)         !< [Input]  k-contraction quartet labels
      integer lcl(num_q)         !< [Input]  l-contraction quartet labels
      logical use_q4             !< [Input]  logical for use of q4
      double precision q4(num_q) !< [Input] symmetry prefactors for ints
*     
      integer l_erilab           !< [Input]  length of eri and label arrays
      integer ilab(l_erilab)     !< [Output] integral labels for ``i'' 
      integer jlab(l_erilab)     !< [Output] integral labels for ``j'' 
      integer klab(l_erilab)     !< [Output] integral labels for ``k'' 
      integer llab(l_erilab)     !< [Output] integral labels for ``l'' 
      integer nint               !< [Output] number of integrals returned
      double precision eri(78,l_erilab) !< [Output] integral values
      integer lscr               !< [Input] length of scratch array
      double precision scr(lscr) !< [Input] scratch array
      double precision zerotol   !< [Input] zero integral threshold
      logical canonicalize       !< [Input] Canonicalize integral labels?
c:tex-\end{verbatim}
*
c:local
      integer i_q 
      integer i_scr_ptr
      integer newlscr
      integer idatom(4)
      integer nint_s
      integer numint
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi
*
c
      intbdd_nw_2e4c = .true.
*
      if (num_q_start.eq.num_q_last) then
         write(luout,*)' l_erilab: too small ',l_erilab
         call errquit(' intbdd_nw_2e4c: infinite blocking loop ',911,
     &             INT_ERR)
      endif
c
*... only do one quartet of derivative integrals at a time!!
      i_q = num_q_start
      num_q_last = num_q_start
*......
      nint_s = 1
*
      if (.not.bas_cn2bfr(brain,icl(i_q),ilo,ihi))
     &    call errquit(' intbdd_nw_2e4c: cn2bfr failed for ish',911,
     &           INT_ERR)
      nint_s = nint_s*(ihi-ilo+1)       ! i partial block size
*
      if (.not.bas_cn2bfr(brain,jcl(i_q),jlo,jhi))
     &    call errquit(' intbdd_nw_2e4c: cn2bfr failed for jsh',911,
     &              INT_ERR)
      nint_s = nint_s*(jhi-jlo+1)       ! j partial block size
*
      if (.not.bas_cn2bfr(brain,kcl(i_q),klo,khi))
     &    call errquit(' intbdd_nw_2e4c: cn2bfr failed for ksh',911,
     &          INT_ERR)
      nint_s = nint_s*(khi-klo+1)       ! k partial block size
*
      if (.not.bas_cn2bfr(brain,lcl(i_q),llo,lhi))
     &    call errquit(' intbdd_nw_2e4c: cn2bfr failed for lsh',911,
     &           INT_ERR)
      nint_s = nint_s*(lhi-llo+1)       ! l partial block size
*......

      if (nint_s.gt.l_erilab) then
        write(luout,*)' nint_s    = ',nint_s
        write(luout,*)' l_erilab  = ',l_erilab
        call errquit(' intbdd_nw_2e4c: buffer size too small ',911,
     &            INT_ERR)
      endif
      numint = int_nint_cart(                    ! get size of cartesian block
     &    brain,icl(i_q),brain,jcl(i_q),
     &    ketin,kcl(i_q),ketin,lcl(i_q))
      numint = numint*78
      i_scr_ptr = numint + 1
      newlscr   = lscr - numint
      if (newlscr .lt. 0) then
        write(luout,*)' newlscr too small      :',newlscr
        write(luout,*)' lscr                   :',lscr
        write(luout,*)' numint                 :',(numint/78)
        write(luout,*)' numint spherical       :',nint_s
        write(luout,*)' numint hess            :',numint
        write(luout,*)' numint hess spherical  :',(nint_s*78)
        write(luout,*)' i bf size is ',(ihi-ilo+1)
        write(luout,*)' j bf size is ',(jhi-jlo+1)
        write(luout,*)' k bf size is ',(khi-klo+1)
        write(luout,*)' l bf size is ',(lhi-llo+1)
        call errquit('intbdd_nw_2e4c: fatal error',911, INT_ERR)
      endif
*
* use memory copy instead of first part of scratch.  Otherwise scratch is 
* way to big to do derivativeswith sphericals
*
      call intdd_2e4c(
     &    brain,icl(i_q),jcl(i_q),
     &    ketin,kcl(i_q),lcl(i_q),
     &    lscr,scr,
     &    numint,eri,idatom)
      call dcopy (numint,eri,1,scr,1)
      if (nint_s.gt.numint) then
        write(luout,*)' spherical/real size larger than cartesean size'
        write(luout,*)' nint_s = ',nint_s
        write(luout,*)' numint = ',numint
        call errquit('intbdd_nw_2e4c: fatal error ',911, INT_ERR)
      endif
* scr has integral derivatives in nwchem order (nint_s,3,4)      
      nint = 0
      call intbdd_nw_label_copy(scr,eri,nint_s,nint,idatom,
     &    ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &    ilab,jlab,klab,llab,
     &    zerotol,canonicalize)
      if (use_q4.and.nint.gt.0) then
        call dscal((nint*78),q4(i_q),eri,1)
      endif
      if (i_q.eq.num_q) then
        num_q_start = 1
        intbdd_nw_2e4c = .false.
      else
        num_q_start = i_q + 1
      endif
c
      end
      subroutine intbdd_nw_label_copy(eri_nworder,eri,nint_s,nint,
     &    idatom,
     &    ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &    ilab,jlab,klab,llab,
     &    zerotol,canonicalize)
      implicit none
c::passed
      integer nint_s
      integer nint
      integer idatom(4)
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi
      integer ilab(nint_s)
      integer jlab(nint_s)
      integer klab(nint_s)
      integer llab(nint_s)
      double precision zerotol
      logical canonicalize
      double precision eri_nworder(llo:lhi,klo:khi,jlo:jhi,ilo:ihi,3,4)
      double precision eri(3,4,nint_s)
c::local
      integer i, xyz, ii, jj, kk, ll
      double precision val_max
      double precision vals(3,4)
      integer numgt1
*
      call dcopy((78*nint_s),0.0d00,0,eri,1)
      nint = 0
      numgt1 = 0
      do i = 1,4
        if (idatom(i).ge.1) numgt1 = numgt1+1
      enddo
      if (numgt1.eq.0) return
*
c right now there doesn't seem to be any logic in intdd_2e4c to set idatom
c correctly.  May have to zero the buffer before making the call to that
c routine.
c     do i = 1,4
c       if (idatom(i).lt.1) then   ! ensure that unused buffers from idatom are zero
c         call dcopy((3*nint_s),0.0d00,0,
c    &        eri_nworder(llo,klo,jlo,ilo,1,i),1)
c       endif
c     enddo
c
c  The do loop structure needs to be fixed to handle getting 78 matrices
c  The current setup is for the gradient code.
c
      do ii = ilo,ihi
        do jj = jlo,jhi
          do kk = klo,khi
            do ll = llo,lhi
              val_max = 0.0d00
              do i = 1,4
                do xyz = 1,3
*assume: memory copy is cheaper than 2 memory access with stride nint_s ?
*check this out
                  vals(xyz,i) = eri_nworder(ll,kk,jj,ii,xyz,i)
                  val_max = max(val_max,abs(vals(xyz,i)))
                enddo
              enddo
              if (val_max.gt.zerotol) then
                nint = nint + 1
                if (canonicalize) then
                  call int_canon(
     &                ii,        jj,        kk,        ll,
     &                ilab(nint),jlab(nint),klab(nint),llab(nint))
                else
                  ilab(nint) = ii
                  jlab(nint) = jj
                  klab(nint) = kk
                  llab(nint) = ll
                endif
                do i = 1,4
                  do xyz = 1,3
                    eri(xyz,i,nint) = vals(xyz,i)
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
*
      return
      end
C> @}
