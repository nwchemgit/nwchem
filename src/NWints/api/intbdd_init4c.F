c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Initialize the blocking 2nd derivative integral API based on the input
C> of a group of shell quartets (a block)
C>
C> Initialize the blocking 2nd derivative integral API based on the input
C> of a group of shell quartets (a block). This routine must be called
C> before attempting to calculate the integrals of the block.
C>
C> The shell quartet information may be reordered
C> since this routine pulls out blocks that consist of only \f$s\f$, \f$p\f$, 
C> and \f$sp\f$ functions.  These blocks are computed using the \f$sp\f$ rotated
C> axis code since it is faster than even the Texas integral code.
C>
C> \return Returns .true. if the initialization was successfull and the 2nd
C> derivative integrals are ready to be calculated by `intbdd_2e4c`, and .false.
C> otherwise (i.e. a problem was detected).
C>
c:tex-\subsection{intbdd\_init4c}
c:tex-This logical function sets up the blocking integral second derivative 
c:tex-API based on the input of a group of shell quartets (a block).  
c:tex-
c:tex-{\it Return Values:}
c:tex-\begin{tabular}{|c|p{60mm}|} \hline
c:tex- .true.  & blocking API initialization okay\\ \hline
c:tex- .false. & blocking API detected a problem \\ \hline
c:tex-\end{tabular}
c:tex-
c:tex-
c:tex-The following code excerpt describes the proper use of the blocking
c:tex-API routines.\\
c:tex-{\it Pseudo Code:}
c:tex-\begin{verbatim}
c:tex-
c:tex-* begin atom/shell loops
c:tex-
c:tex-      call collect_group_of_shells()
c:tex-      okay = intbdd_init4c($\cdots$) ! with group of shells info     
c:tex-      if (.not.okay) stop ' error setting up blocking interface '
c:tex-00001 continue
c:tex-      more = intbdd_2e4c($\cdots$) ! with group of shells info
c:tex-      call use_integrals_in_some_way()
c:tex-      if (more) goto 00001
c:tex-
c:tex-* end atom/shell loops
c:tex-
c:tex-\end{verbatim}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function intbdd_init4c(brain, icl, jcl, ketin, kcl, lcl,
     &    num_q, q4, use_q4, lscr, scr, l_erilab, block_eff)
c:tex-\end{verbatim}
      implicit none
*
c
c blocking initialization routine for 2 electron 4 center integrals
c
*
c:tex-\begin{verbatim}
      integer brain      ! [input] basis set handle for bra basis
      integer ketin      ! [input] basis set handle for ket basis
      integer num_q      ! [input] number of quartets
      integer icl(num_q) ! [input] i-contraction labels for quartets
      integer jcl(num_q) ! [input] j-contraction labels for quartets
      integer kcl(num_q) ! [input] k-contraction labels for quartets
      integer lcl(num_q) ! [input] l-contraction labels for quartets
      double precision q4(num_q) ! [input] scaling factors
      logical use_q4            ! [input] true if scaling
      integer l_erilab   ! [input] size of eri and label arrays that 
*................................. will be used in intbdd_2e4c.F
      integer lscr       ! [input] length of scratch array
      double precision scr(lscr)  ! [input] scratch array
      double precision block_eff  ! [output] blocking efficiency 
c:tex-\end{verbatim}
#include "candoP.fh"
#include "numb_qP.fh"
c
      logical status
c
      logical cando_txs
      external cando_txs
c
      integer iq
c
      num_q_start = 1           
      num_q_last  = -1
c
      more_ints_txs = .false. 
      intbdd_init4c = .true.
      if (num_q.eq.0) return
      block_eff = -1.0d00
c
*.... check shells for txs only after intbdd_init is called and established        
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
