c $Id: int_hf1sp.F,v 1.19 2003-10-17 22:53:42 carlfahl Exp $
*
c:tex-% part of the internal API routines.
c:tex-\subsection{int\_hf1sp}
c:tex-This is a layer routine that calls the 
c:tex-McMurchie-Davidson one electron routine.  This layer
c:tex-routine handles all permutations to compute sp integrals.
c:tex-This routine should never be called by an application module.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_hf1sp(
     &       xyzi,expi,coefi, i_nprim, i_ngen, Li, ictri,
     &       xyzj,expj,coefj, j_nprim, j_ngen, Lj, ictrj,
     &       xyz,zan,exinv,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &       dryrun,scr,lscr,msg)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
#include "util.fh"
c::passed
c:tex-For an integral $<i|Operator|j>$
c:tex-\begin{verbatim}
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i  
      integer ictri    ! [input] lexical atom index for function i 
      integer j_nprim  ! [input] num. prims on function j     
      integer j_ngen   ! [input] num general conts on func. j 
      integer Lj       ! [input] angular momentum of func. j  
      integer ictrj    ! [input] lexical atom index for function j 
      integer nat      ! [input] number of atoms
      integer lscr     ! [input] size of scratch array
      integer lstv     ! [input] size of any integral buffer
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision xyz(3,nat)  ! [input] all atom positions
      double precision zan(nat)  ! [input] charges on all atoms
      double precision exinv(nat)   ! [input] inverse nuclear exponents
      double precision scr(lscr)  ! [scratch] scratch buffers
      double precision S(lstv)  ! [output] overlap integrals
      double precision T(lstv)  ! [output] kinetic energy integrals
      double precision V(lstv)  ! [output] potential integrals
      logical doS     ! [input] compute overlap (True/False)
      logical doT     ! [input] compute kinetic (True/False)
      logical doV     ! [input] compute potential (True/False)
      logical canAB   ! [input] compute only canonical ints (false only)
      logical dryrun  ! [input] true means only compute required memory
      character*(*) msg ![input] calling func. identification message
c:tex-\end{verbatim}
c::local
      integer offset
      integer nintx, nint_xsp
c
*debug_print      integer iii
c
      if (Li.gt.-1.and.Lj.gt.-1) then
        call int_hf1sp_ecp(
     &         xyzi,expi,coefi,i_nprim,i_ngen,Li,ictri,
     &         xyzj,expj,coefj,j_nprim,j_ngen,Lj,ictrj,
     &         xyz,zan,exinv,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
c
      elseif (Li.eq.-1.and.Lj.eq.-1) then
c.............................................................. (SP|SP)
c________________________________ (S|S)
        call int_hf1sp_ecp(
     &         xyzi,expi,coefi,i_nprim,1,0,ictri,
     &         xyzj,expj,coefj,j_nprim,1,0,ictrj,
     &         xyz,zan,exinv,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
c________________________________ (S|P)
        offset = 1
        call int_hf1sp_ecp(
     &         xyzi,
     &         expi,coefi(1,1),
     &         i_nprim,1,0,ictri,
     &         xyzj,
     &         expj,coefj(1,2),
     &         j_nprim,1,1,ictrj,
     &         xyz,zan,exinv,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
c________________________________ (P|S)
        offset = 1 + 3 
        call int_hf1sp_ecp(
     &         xyzi,
     &         expi,coefi(1,2),
     &         i_nprim,1,1,ictri,
     &         xyzj,
     &         expj,coefj(1,1),
     &         j_nprim,1,0,ictrj,
     &         xyz,zan,exinv,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
c________________________________ (P|P)
        offset = 1 + 3 + 3
        call int_hf1sp_ecp(
     &         xyzi,
     &         expi,coefi(1,2),
     &         i_nprim,1,1,ictri,
     &         xyzj,
     &         expj,coefj(1,2),
     &         j_nprim,1,1,ictrj,
     &         xyz,zan,exinv,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
*
        if (doS) call int_1spsp(S,1)
        if (doT) call int_1spsp(T,1)
        if (dov) call int_1spsp(V,1)
      elseif (Li.eq.-1) then
c.............................................................. (SP|X)
c________________________________ (S|X)
        call int_hf1sp_ecp(
     &         xyzi,expi,coefi,i_nprim,1,0,ictri,
     &         xyzj,expj,coefj,j_nprim,j_ngen,Lj,ictrj,
     &         xyz,zan,exinv,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
c________________________________ (P|X)
        offset = 1*(Lj+1)*(Lj+2)/2*j_ngen
        call int_hf1sp_ecp(
     &         xyzi,
     &         expi,coefi(1,2),
     &         i_nprim,1,1,ictri,
     &         xyzj,
     &         expj,coefj,j_nprim,j_ngen,Lj,ictrj,
     &         xyz,zan,exinv,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
      elseif (Lj.eq.-1) then
c.............................................................. (X|SP)
c________________________________ (X|S)
        call int_hf1sp_ecp(
     &         xyzi,expi,coefi,i_nprim,i_ngen,Li,ictri,
     &         xyzj,expj,coefj,j_nprim,1,0,ictrj,
     &         xyz,zan,exinv,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
c________________________________ (X|P)
        offset = 1*(Li+1)*(Li+2)/2*i_ngen
        call int_hf1sp_ecp(
     &         xyzi,
     &         expi,coefi,i_nprim,i_ngen,Li,ictri,
     &         xyzj,
     &         expj,coefj(1,2),
     &         j_nprim,1,1,ictrj,
     &         xyz,zan,exinv,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr,msg)
c
        if (li.eq.0.and.i_ngen.eq.1) then
c * for i_ngen = 1
c * nothing needed since (s|s), (s|p) is same order as (s|sp)
c * (s|s), (s|x), (s|y), (s|z)
          continue
        else if (li.eq.1.and.i_ngen.eq.1) then
          if (doS) call int_1psp(S,1)
          if (doT) call int_1psp(T,1)
          if (doV) call int_1psp(V,1)
        else if (li.eq.2.and.i_ngen.eq.1) then
          if (doS) call int_1dsp(S,1)
          if (doT) call int_1dsp(T,1)
          if (doV) call int_1dsp(V,1)
        else
          nintx   = (Li+1)*(Li+2)/2*i_ngen  ! size of X
          nint_xsp = nintx * 4              ! size of X * size of sp
          if (doS) then
            call dcopy(nint_xsp,S,1,scr,1)
            call int_sp1b(S,scr,scr(nintx+1),nintx,1)
          endif
          if (doT) then
            call dcopy(nint_xsp,T,1,scr,1)
            call int_sp1b(T,scr,scr(nintx+1),nintx,1)
          endif
          if (doV) then
            call dcopy(nint_xsp,V,1,scr,1)
            call int_sp1b(V,scr,scr(nintx+1),nintx,1)
          endif
        endif
      else
        write(luout,*)'int_hf1sp called by ',msg
        call errquit('int_hf1sp: unknown case finished',911, INT_ERR)
      endif
c
c
      end
c:tex-% part of the internal API routines.
c:tex-\subsection{int\_hf1sp\_ecp}
c:tex-This is a layer routine that calls the 
c:tex-McMurchie-Davidson one electron routine.  This layer
c:tex-routine handles all options for computing ecp integrals.
c:tex-This routine should never be called by an application module.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_hf1sp_ecp(
     &       xyzi,expi,coefi, i_nprim, i_ngen, Li, ictri,
     &       xyzj,expj,coefj, j_nprim, j_ngen, Lj, ictrj,
     &       xyz,zan,exinv,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &       dryrun,scr,lscr,msg)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
#include "util.fh"
c::passed
c:tex-For an integral $<i|Operator|j>$
c:tex-\begin{verbatim}
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i  
      integer ictri    ! [input] lexical atom index for function i 
      integer j_nprim  ! [input] num. prims on function j     
      integer j_ngen   ! [input] num general conts on func. j 
      integer Lj       ! [input] angular momentum of func. j  
      integer ictrj    ! [input] lexical atom index for function j 
      integer nat      ! [input] number of atoms
      integer lscr     ! [input] size of scratch array
      integer lstv     ! [input] size of any integral buffer
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision xyz(3,nat)  ! [input] all atom positions
      double precision zan(nat)  ! [input] charges on all atoms
      double precision exinv(nat)   ! [input] inverse nuclear exponents
      double precision scr(lscr)  ! [scratch] scratch buffers
      double precision S(lstv)  ! [output] overlap integrals
      double precision T(lstv)  ! [output] kinetic energy integrals
      double precision V(lstv)  ! [output] potential integrals
      logical doS     ! [input] compute overlap (True/False)
      logical doT     ! [input] compute kinetic (True/False)
      logical doV     ! [input] compute potential (True/False)
      logical canAB   ! [input] compute only canonical ints (false only)
      logical dryrun  ! [input] true means only compute required memory
      character*(*) msg ![input] calling func. identification message
c:tex-\end{verbatim}
c::local
      integer nintV
      integer offset
*
      call hf1(
     &      xyzi,expi,coefi,i_nprim,i_ngen,Li,
     &      xyzj,expj,coefj,j_nprim,j_ngen,Lj,
     &      xyz,zan,exinv,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &      dryrun,scr,lscr)
*
      if (any_ecp.and.doV.and.msg.ne.'int_1eefc') then
        nintV = int_nbf_x(Li)*i_ngen
        nintV = nintV*int_nbf_x(Lj)*j_ngen
        offset = nintV + 1
*       write(luout,*)' lscr to ecp_hf1:',(lscr-nintV)
        if (ictri.lt.0.or.ictrj.lt.0) then
          write(luout,*)' int_hf1sp_ecp: ictri = ',ictri
          write(luout,*)' int_hf1sp_ecp: ictrj = ',ictrj
          write(luout,*)'int_hf1sp_ecp called by ',msg
          call errquit(
     &        'int_hf1sp: unknown center for ECPs',
     &        911, INT_ERR)
        endif
        call int_ecp_hf1(
     &      xyzi,expi,coefi,i_nprim,i_ngen,Li,ictri,
     &      xyzj,expj,coefj,j_nprim,j_ngen,Lj,ictrj,
     &      scr,nintV,scr(offset),(lscr-nintV),
     &      dryrun)
*... sum ecp into V block
        call daxpy(nintV, 1.0d00, scr,1,V,1)
      endif
*
      end
c:tex-% part of the internal API routines
c:tex-\subsection{int\_1psp}
c:tex-This routine transforms integrals from the way they
c:tex-were computed $(p|s)$, $(p|p)$ to $(p|sp)$.
c:tex-The transformation is done in place as follows:
c:tex-\begin{tabular}{rccc}
c:tex- & computed & transformed & \\
c:tex- &  order   &   order     & \\
c:tex-~1 & (x$|$s)  &  (x$|$s)  & \\
c:tex-~2 & (y$|$s)  &  (x$|$x)  & 2 $\rightarrow$ 5\\
c:tex-~3 & (z$|$s)  &  (x$|$y)  & 3 $\rightarrow$ 9\\
c:tex-~4 & (x$|$x)  &  (x$|$z)  & 4 $\rightarrow$ 2\\
c:tex-~5 & (x$|$y)  &  (y$|$s)  & 5 $\rightarrow$ 3\\
c:tex-~6 & (x$|$z)  &  (y$|$x)  & 6 $\rightarrow$ 4\\
c:tex-~7 & (y$|$x)  &  (y$|$y)  & 7 $\rightarrow$ 6\\
c:tex-~8 & (y$|$y)  &  (y$|$z)  & 8 $\rightarrow$ 7\\
c:tex-~9 & (y$|$z)  &  (z$|$s)  & 9 $\rightarrow$ 8\\
c:tex-10 & (z$|$x)  &  (z$|$x)  &  \\
c:tex-11 & (z$|$y)  &  (z$|$y)  &  \\
c:tex-12 & (z$|$z)  &  (z$|$z)  &  \\
c:tex-\end{tabular}
c:tex- 
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1psp(block,num_blocks)
c:tex-\end{verbatim}
c
c transforms a (p|sp) block to correct order in place
c integrals in block were calculated (p|s),(p|p).
c      computed transformed 
c       order     order
c  1.  (x|s)     (x|s)     
c  2.  (y|s)     (x|x)     2 -> 5
c  3.  (z|s)     (x|y)     3 -> 9
c  4.  (x|x)     (x|z)     4 -> 2
c  5.  (x|y)     (y|s)     5 -> 3
c  6.  (x|z)     (y|x)     6 -> 4
c  7.  (y|x)     (y|y)     7 -> 6
c  8.  (y|y)     (y|z)     8 -> 7
c  9.  (y|z)     (z|s)     9 -> 8
c 10.  (z|x)     (z|x)     
c 11.  (z|y)     (z|y)     
c 12.  (z|z)     (z|z)     
c
      implicit none
c:tex-\begin{verbatim}
      integer num_blocks  ! [input] num. blocks to transform
      double precision block(12,num_blocks) ! [input/output] 
c. . . . . . . . . . . . . . . . . . . . . .! integral block
c:tex-\end{verbatim}
c
      integer ib
      double precision temp1, temp2
      do 00100 ib = 1,num_blocks
        temp1       = block(2,ib)  ! 2 -> temp1
        temp2       = block(3,ib)  ! 3 -> temp2
        block(2,ib) = block(4,ib)  ! 4 -> 2
        block(3,ib) = block(5,ib)  ! 5 -> 3
        block(4,ib) = block(6,ib)  ! 6 -> 4
        block(6,ib) = block(7,ib)  ! 7 -> 6
        block(7,ib) = block(8,ib)  ! 8 -> 7
        block(8,ib) = block(9,ib)  ! 9 -> 8
        block(9,ib) = temp2        ! temp2 -> 9 or 3 -> 9
        block(5,ib) = temp1        ! temp1 -> 5 or 2 -> 5
00100 continue
      end
c:tex-% part of the internal API routines
c:tex-\subsection{int\_1dsp}
c:tex-This routine transforms integrals from the way they
c:tex-were computed $(d|s)$, $(d|p)$ to $(d|sp)$.
c:tex-The transformation is done in place as follows:
c:tex-\begin{tabular}{rccc}      
c:tex-  &  computed  & transformed & \\
c:tex-  &  ~order    &  order & \\
c:tex-~1 & (xx$|$s)  &  (xx$|$s)  & \\
c:tex-~2 & (xy$|$s)  &  (xx$|$x)  &  ~2 $\rightarrow$ 5~~ \\
c:tex-~3 & (xz$|$s)  &  (xx$|$y)  &  ~3 $\rightarrow$ 9~~ \\
c:tex-~4 & (yy$|$s)  &  (xx$|$z)  &  ~4 $\rightarrow$ 13~ \\
c:tex-~5 & (yz$|$s)  &  (xy$|$s)  &  ~5 $\rightarrow$ 17~ \\
c:tex-~6 & (zz$|$s)  &  (xy$|$x)  &  ~6 $\rightarrow$ 21~ \\
c:tex-~7 & (xx$|$x)  &  (xy$|$y)  &  ~7 $\rightarrow$ 2~~ \\
c:tex-~8 & (xx$|$y)  &  (xy$|$z)  &  ~8 $\rightarrow$ 3~~ \\
c:tex-~9 & (xx$|$z)  &  (xz$|$s)  &  ~9 $\rightarrow$ 4~~ \\
c:tex-10 & (xy$|$x)  &  (xz$|$x)  &  10 $\rightarrow$ 6~~ \\
c:tex-11 & (xy$|$y)  &  (xz$|$y)  &  11 $\rightarrow$ 7~~ \\
c:tex-12 & (xy$|$z)  &  (xz$|$z)  &  12 $\rightarrow$ 8~~ \\
c:tex-13 & (xz$|$x)  &  (yy$|$s)  &  13 $\rightarrow$ 10~ \\
c:tex-14 & (xz$|$y)  &  (yy$|$x)  &  14 $\rightarrow$ 11~ \\
c:tex-15 & (xz$|$z)  &  (yy$|$y)  &  15 $\rightarrow$ 12~ \\
c:tex-16 & (yy$|$x)  &  (yy$|$z)  &  16 $\rightarrow$ 14~ \\
c:tex-17 & (yy$|$y)  &  (yz$|$s)  &  17 $\rightarrow$ 15~ \\
c:tex-18 & (yy$|$z)  &  (yz$|$x)  &  18 $\rightarrow$ 16~ \\
c:tex-19 & (yz$|$x)  &  (yz$|$y)  &  19 $\rightarrow$ 18~ \\
c:tex-20 & (yz$|$y)  &  (yz$|$z)  &  20 $\rightarrow$ 19~ \\
c:tex-21 & (yz$|$z)  &  (zz$|$s)  &  21 $\rightarrow$ 20~ \\
c:tex-22 & (zz$|$x)  &  (zz$|$x)  & \\
c:tex-23 & (zz$|$y)  &  (zz$|$y)  & \\
c:tex-24 & (zz$|$z)  &  (zz$|$z)  & \\
c:tex-\end{tabular}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1dsp(block,num_blocks)
c:tex-\end{verbatim}
c
c transforms a (d|sp) block to correct order in place
c integrals in block were calculated (d|s),(d|p).
c      computed transformed 
c       order     order
c  1.  (xx|s)     (xx|s)
c  2.  (xy|s)     (xx|x)  2 -> 5   *
c  3.  (xz|s)     (xx|y)  3 -> 9   *
c  4.  (yy|s)     (xx|z)  4 -> 13  *
c  5.  (yz|s)     (xy|s)  5 -> 17  *
c  6.  (zz|s)     (xy|x)  6 -> 21  *
c  7.  (xx|x)     (xy|y)  7 -> 2   *
c  8.  (xx|y)     (xy|z)  8 -> 3   *
c  9.  (xx|z)     (xz|s)  9 -> 4   *
c 10.  (xy|x)     (xz|x) 10 -> 6   *
c 11.  (xy|y)     (xz|y) 11 -> 7   *
c 12.  (xy|z)     (xz|z) 12 -> 8   *
c 13.  (xz|x)     (yy|s) 13 -> 10  *
c 14.  (xz|y)     (yy|x) 14 -> 11  *
c 15.  (xz|z)     (yy|y) 15 -> 12  *
c 16.  (yy|x)     (yy|z) 16 -> 14  *
c 17.  (yy|y)     (yz|s) 17 -> 15  *
c 18.  (yy|z)     (yz|x) 18 -> 16  *
c 19.  (yz|x)     (yz|y) 19 -> 18  *
c 20.  (yz|y)     (yz|z) 20 -> 19  *
c 21.  (yz|z)     (zz|s) 21 -> 20  *
c 22.  (zz|x)     (zz|x)
c 23.  (zz|y)     (zz|y)
c 24.  (zz|z)     (zz|z)
c
      implicit none
c:tex-\begin{verbatim}
      integer num_blocks  ! [input] num. blocks to transform
      double precision block(24,num_blocks) ! [input/output] 
c . . . . . . . . . . . . . . . . . . . . . ! integral block
c:tex-\end{verbatim}
c
      double precision temp1, temp2, temp3
      integer ib
c
      do 00100 ib = 1, num_blocks
        temp1        = block(2,ib)    ! 2 -> temp1
        block(2,ib)  = block(7,ib)    ! 7 -> 2
        temp2        = block(3,ib)    ! 3 -> temp2
        block(3,ib)  = block(8,ib)    ! 8 -> 3
        temp3        = block(4,ib)    ! 4 -> temp3
        block(4,ib)  = block(9,ib)    ! 9 -> 4
        block(9,ib)  = temp2          ! temp2 -> 9 or 3 -> 9
        temp2        = block(6,ib)    ! 6 -> temp2
        block(6,ib)  = block(10,ib)   ! 10 -> 6
        block(7,ib)  = block(11,ib)   ! 11 -> 7
        block(8,ib)  = block(12,ib)   ! 12 -> 8
        block(10,ib) = block(13,ib)   ! 13 -> 10
        block(11,ib) = block(14,ib)   ! 14 -> 11
        block(12,ib) = block(15,ib)   ! 15 -> 12
        block(14,ib) = block(16,ib)   ! 16 -> 14
        block(15,ib) = block(17,ib)   ! 17 -> 15
        block(17,ib) = block(5,ib)    ! 5 -> 17
        block(5,ib)  = temp1          ! temp1 -> 5 or 2 -> 5
        block(16,ib) = block(18,ib)   ! 18 -> 16
        block(18,ib) = block(19,ib)   ! 19 -> 18
        block(19,ib) = block(20,ib)   ! 20 -> 19
        block(20,ib) = block(21,ib)   ! 21 -> 20
        block(13,ib) = temp3          ! temp3 -> 13 or 4 -> 13
        block(21,ib) = temp2          ! temp2 -> 21 or 6 -> 21
00100 continue
      end
c:tex-% part of the internal API routines
c:tex-\subsection{int\_1spsp}
c:tex-This routine transforms integrals from the way they
c:tex-were computed $(s|s)$, $(s|p)$, $(p|s)$, $(p|p)$  
c:tex-to $(sp|sp)$.
c:tex-The transformation is done in place as follows:
c:tex-\begin{tabular}{rccc}      
c:tex-  &  computed & transformed   &   \\
c:tex-  &   order   &  order  &   \\
c:tex-~1 & (s$|$s)  &  (s$|$s)  &   \\
c:tex-~2 & (s$|$x)  &  (s$|$x)  &   \\
c:tex-~3 & (s$|$y)  &  (s$|$y)  &   \\
c:tex-~4 & (s$|$z)  &  (s$|$z)  &   \\
c:tex-~5 & (x$|$s)  &  (x$|$s)  &   \\
c:tex-~6 & (y$|$s)  &  (x$|$x)  &  ~6 $\rightarrow$ 9\\
c:tex-~7 & (z$|$s)  &  (x$|$y)  &  ~7 $\rightarrow$ 13\\
c:tex-~8 & (x$|$x)  &  (x$|$z)  &  ~8 $\rightarrow$ 6\\
c:tex-~9 & (x$|$y)  &  (y$|$s)  &  ~9 $\rightarrow$ 7\\
c:tex-10 & (x$|$z)  &  (y$|$x)  &  10 $\rightarrow$ 8 \\ 
c:tex-11 & (y$|$x)  &  (y$|$y)  &  11 $\rightarrow$ 10\\
c:tex-12 & (y$|$y)  &  (y$|$z)  &  12 $\rightarrow$ 11\\
c:tex-13 & (y$|$z)  &  (z$|$s)  &  13 $\rightarrow$ 12\\
c:tex-14 & (z$|$x)  &  (z$|$x)  &   \\
c:tex-15 & (z$|$y)  &  (z$|$y)  &   \\
c:tex-16 & (z$|$z)  &  (z$|$z)  &   \\
c:tex-\end{tabular}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1spsp(block,num_blocks)
c:tex-\end{verbatim}
c
c transforms a 1 electron integral block holding the 16
c (sp|sp) integrals calculated from 4 calls:
c           (s|s),(s|p),(p|s),(p|p) to the (sp|sp) order
c
c      computed transformed 
c       order      order
c  1.   (s|s)      (s|s)
c  2.   (s|x)      (s|x)  
c  3.   (s|y)      (s|y)
c  4.   (s|z)      (s|z)
c  5.   (x|s)      (x|s)
c  6.   (y|s)      (x|x)   6 -> 9
c  7.   (z|s)      (x|y)   7 -> 13
c  8.   (x|x)      (x|z)   8 -> 6
c  9.   (x|y)      (y|s)   9 -> 7
c 10.   (x|z)      (y|x)  10 -> 8 
c 11.   (y|x)      (y|y)  11 -> 10
c 12.   (y|y)      (y|z)  12 -> 11
c 13.   (y|z)      (z|s)  13 -> 12
c 14.   (z|x)      (z|x)
c 15.   (z|y)      (z|y)
c 16.   (z|z)      (z|z)
      implicit none
c:tex-\begin{verbatim}
      integer num_blocks        ! [input] num. blocks to transform
      double precision block(16,num_blocks) ! [input/output] 
c . . . . . . . . . . . . . . . . . . . . . ! integral block
c:tex-\end{verbatim}
c
      double precision temp
      integer ib
c
      do 00100 ib = 1,num_blocks
        temp         = block(6, ib)
        block(6, ib) = block(8, ib)
        block(8, ib) = block(10,ib)
        block(10,ib) = block(11,ib)
        block(11,ib) = block(12,ib)
        block(12,ib) = block(13,ib)
        block(13,ib) = block(7, ib)
        block(7, ib) = block(9, ib)
        block(9, ib) = temp
00100 continue
      end
c:tex-% part of the internal API routines
c:tex-\subsection{int\_1spa}
c:tex-This routine transforms integrals from the way they
c:tex-were computed $(s|X)$, $(p|X)$, to $(sp|X)$.
c:tex-The transformation is {\bf NOT} done in place:
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_sp1a(sp_block,s_block,p_block,sizeb,num_blocks)
c:tex-\end{verbatim}
c
*..  sp 1 electron transformation routine 
*.. (ash|bsh) block of integrals where a is an sp shell
c
      implicit none
c::passed
c:tex-\begin{verbatim}
      integer sizeb  ! [input] size of non sp block
      integer num_blocks ! [input] num of blocks to transform
*
* . . . . . . . . . . . . . . .! [output] (sp|X) transformed integral block
      double precision sp_block(sizeb,1:4,num_blocks) 
*
* . . . . . . . . . . . . . . . . . . . . . .! [input] computed (s|X) block
      double precision s_block(sizeb,num_blocks)
*
* . . . . . . . . . . . . . . . . . . . . . .! [input] computed (p|X) block
      double precision p_block(sizeb,2:4,num_blocks)
c:tex-\end{verbatim}
c::local
      integer block
      integer ib
c
      do 00100 block=1,num_blocks
        do 00200 ib = 1,sizeb
          sp_block(ib,1,block) = s_block(ib,block)
          sp_block(ib,2,block) = p_block(ib,2,block)
          sp_block(ib,3,block) = p_block(ib,3,block)
          sp_block(ib,4,block) = p_block(ib,4,block)
00200   continue
00100 continue
      end
c:tex-% part of the internal API routines
c:tex-\subsection{int\_sp1b}
c:tex-This routine transforms integrals from the way they
c:tex-were computed $(X|s)$, $(X|p)$, to $(X|sp)$.
c:tex-The transformation is {\bf NOT} done in place:
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_sp1b(sp_block,s_block,p_block,sizea,num_blocks)
c:tex-\end{verbatim}
c
*..  sp 1 electron transformation routine 
*.. (ash|bsh) block of integrals where b is an sp shell
c
      implicit none
c::passed
c:tex-\begin{verbatim}
      integer sizea ! [input] size of non sp block
      integer num_blocks ! [input] num of blocks to transform
*
* . . . . . . . . . . . . . . .! [output] (X|sp) transformed integral block
      double precision sp_block(1:4,sizea,num_blocks)
*
* . . . . . . . . . . . . . . . . . . . . . .! [input] computed (X|s) block
      double precision s_block(sizea,num_blocks)
*
* . . . . . . . . . . . . . . . . . . . . . .! [input] computed (X|p) block
      double precision p_block(2:4,sizea,num_blocks)
c:tex-\end{verbatim}
c::local
      integer block
      integer ia
c
      do 00100 block = 1,num_blocks
        do 00200 ia = 1,sizea
          
          sp_block(1,ia,block) = s_block(ia,block)
          sp_block(2,ia,block) = p_block(2,ia,block)
          sp_block(3,ia,block) = p_block(3,ia,block)
          sp_block(4,ia,block) = p_block(4,ia,block)
00200   continue
00100 continue
      end
