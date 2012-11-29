c $Id$  
*
C> \ingroup nwint
C> @{
C>
C> \brief Initialize the blocking integral API based on the input
C> of a group of shell quartets (a block)
C>
C> Initialize the blocking integral API based on the input
C> of a group of shell quartets (a block). This routine must be called
C> before attempting to calculate the integrals of the block.
C>
C> The shell quartet information may be reordered
C> since this routine pulls out blocks that consist of only \f$s\f$, \f$p\f$, 
C> and \f$sp\f$ functions.  These blocks are computed using the \f$sp\f$ rotated
C> axis code since it is faster than even the Texas integral code.
C>
C> \return Returns .true. if the initialization was successfull and the
C> integrals are ready to be calculated by `intb_2e4c`, and .false.
C> otherwise (i.e. a problem was detected).
C>
C> The following code excerpt describes the proper use of the blocking
C> API routines:
C> 
C> * begin atom/shell loops
C> 
C> \code
C>       call collect_group_of_shells()
C>       okay = intb_init4c($\cdots$) ! with group of shells info     
C>       if (.not.okay) stop ' error setting up blocking interface '
C> 00001 continue
C>       more = intb_2e4c($\cdots$) ! with group of shells info
C>       call use_integrals_in_some_way()
C>       if (more) goto 00001
C> \endcode
C> 
C> * end atom/shell loops
C> 
c:tex-\subsection{intb\_init4c}
c:tex-This logical function sets up the blocking integral API based on the
c:tex-input of a group of shell quartets (a block).  
c:tex-
c:tex-{\it Return Values:}
c:tex-\begin{tabular}{|c|p{60mm}|} \hline
c:tex- .true.  & blocking API initialization okay\\ \hline
c:tex- .false. & blocking API detected a problem \\ \hline
c:tex-\end{tabular}
c:tex-
c:tex-{\it Side Effects:}  The shell quartet information may be reordered
c:tex-since this routine pulls out blocks that consist of only $s$, $p$, 
c:tex-and $sp$ functions.  Thes blocks are computed using the $sp$ rotated
c:tex-axis code since it is faster than even the Texas integral code.
c:tex-
c:tex-The following code excerpt describes the proper use of the blocking
c:tex-API routines.\\
c:tex-{\it Pseudo Code:}
c:tex-\begin{verbatim}
c:tex-
c:tex-* begin atom/shell loops
c:tex-
c:tex-      call collect_group_of_shells()
c:tex-      okay = intb_init4c($\cdots$) ! with group of shells info     
c:tex-      if (.not.okay) stop ' error setting up blocking interface '
c:tex-00001 continue
c:tex-      more = intb_2e4c($\cdots$) ! with group of shells info
c:tex-      call use_integrals_in_some_way()
c:tex-      if (more) goto 00001
c:tex-
c:tex-* end atom/shell loops
c:tex-
c:tex-\end{verbatim}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function intb_init4c(brain, icl, jcl, ketin, kcl, lcl,
     &    num_q, q4, use_q4, lscr, scr, l_erilab, block_eff)
c:tex-\end{verbatim}
      implicit none
*
c
c blocking initialization routine for 2 electron 4 center integrals
c
*
c:tex-\begin{verbatim}
      integer brain      !< [Input] basis set handle for bra basis
      integer ketin      !< [Input] basis set handle for ket basis
      integer num_q      !< [Input] number of quartets
      integer icl(num_q) !< [Input] i-contraction labels for quartets
      integer jcl(num_q) !< [Input] j-contraction labels for quartets
      integer kcl(num_q) !< [Input] k-contraction labels for quartets
      integer lcl(num_q) !< [Input] l-contraction labels for quartets
      double precision q4(num_q) !< [Input] scaling factors
      logical use_q4            !< [Input] true if scaling
      integer l_erilab   !< [Input] size of eri and label arrays that 
                         !< will be used in `intb_2e4c`
      integer lscr       !< [Input] length of scratch array
      double precision scr(lscr)  !< [Input] scratch array
      double precision block_eff  !< [Output] blocking efficiency 
c:tex-\end{verbatim}
#include "mafdecls.fh"
#include "candoP.fh"
#include "numb_qP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c
      logical status
c
      logical cando_txs
      external cando_txs
c
      integer iu, ju, ku, lu
      integer itype, jtype, ktype, ltype
      integer inprim, jnprim, knprim, lnprim
      integer ingen, jngen, kngen, lngen
      integer nprimtest
      integer iq, bra, ket
      double precision tmp
#include "bas_ibs_sfn.fh"
c
      num_q_start = 1           
      num_q_last  = -1
c
      more_ints_txs = .false. 
      intb_init4c = .true.
      block_eff = -1.0d00
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
c     The SP integrals are faster than texas for S and SP shells with
c     contraction of 3 or greater.  Select quartets which have only
c     S or SP shells of which all have contraction greater than or
c     equal to 3.  These shells are put onto the end of the list and
c     then evaluated first using the existing NW blocking interface.
c
#if (defined(SUN) || defined(SOLARIS))
c     scalar machines ... older SUNS use 1 ... ultra-2 wants 3
*      nprimtest = 1
      nprimtest = 3
#else 
c     pipelined/vector machines
      nprimtest = 3
#endif
      nsp = 0
      if (.not. (user_cando_sp.and.(.not.def_cando_sp))) then
         do 10 iq = 1, num_q
            iu = (sf_ibs_cn2ucn(icl(iq),bra))
            itype = infbs_cont(CONT_TYPE,  iu, bra)
            if (itype .gt. 0) goto 10
            inprim= infbs_cont(CONT_NPRIM, iu, bra)
            if (inprim .lt. nprimtest) goto 10
            ingen = infbs_cont(CONT_NGEN,  iu, bra)
            if (itype.eq.0 .and. ingen.gt.1) goto 10
c     
            ju = (sf_ibs_cn2ucn(jcl(iq),bra))
            jtype = infbs_cont(CONT_TYPE,  ju, bra)
            if (jtype .gt. 0) goto 10
            jnprim= infbs_cont(CONT_NPRIM, ju, bra)
            if (jnprim .lt. nprimtest) goto 10
            jngen = infbs_cont(CONT_NGEN,  ju, bra)
            if (jtype.eq.0 .and. jngen.gt.1) goto 10
c     
            ku = (sf_ibs_cn2ucn(kcl(iq),bra))
            ktype = infbs_cont(CONT_TYPE,  ku, ket)
            if (ktype .gt. 0) goto 10
            knprim= infbs_cont(CONT_NPRIM, ku, ket)
            if (knprim .lt. nprimtest) goto 10
            kngen = infbs_cont(CONT_NGEN,  ku, ket)
            if (ktype.eq.0 .and. kngen.gt.1) goto 10
c     
            lu = (sf_ibs_cn2ucn(lcl(iq),bra))
            ltype = infbs_cont(CONT_TYPE,  lu, ket)
            if (ltype .gt. 0) goto 10
            lnprim= infbs_cont(CONT_NPRIM, lu, ket)
            if (lnprim .lt. nprimtest) goto 10
            lngen = infbs_cont(CONT_NGEN,  lu, ket)
            if (ltype.eq.0 .and. lngen.gt.1) goto 10
c     
c     Now have a shell quartet that fits the criteria.  Swap to beginning
c     of the list
c     
            nsp = nsp + 1
            if (iq .ne. nsp) then
               iu       = icl(iq)
               ju       = jcl(iq)
               ku       = kcl(iq)
               lu       = lcl(iq)
               icl(iq)  = icl(nsp)
               jcl(iq)  = jcl(nsp)
               kcl(iq)  = kcl(nsp)
               lcl(iq)  = lcl(nsp)
               icl(nsp) = iu
               jcl(nsp) = ju
               kcl(nsp) = ku
               lcl(nsp) = lu
               if (use_q4) then
                  tmp     = q4(iq)
                  q4(iq)  = q4(nsp)
                  q4(nsp) = tmp
               endif
            endif
 10      enddo
#if defined(DEBUG_TXS)
         write(luout,*) ' num_q nsp ', num_q, nsp
#endif
      endif
c

*.... check shells for txs only after intb_init is called and established        
      txs_b2e4c = .true.
      status = .true.
*.... check first 20 quartets but every other one
      do iq = 1,min(20,num_q),2
        status = status.and. (cando_txs(brain,icl(iq),jcl(iq)))
        status = status.and. (cando_txs(ketin,kcl(iq),lcl(iq)))
        if (.not.status) then
          txs_b2e4c = .false.
          return
        endif
      enddo
      end
C> @}
