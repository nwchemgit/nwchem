c $Id: int_mpole.F,v 1.9 1997-04-25 00:40:53 d3e129 Exp $
*
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsubsection{int\_mpole}
c:tex-This routine returns multipole integrals up to the level $lmax$
c:tex-
c:tex- The general form is $< shell | pole | shell >$. 
c:tex- Integrals are returned in shell blocks of $<ish|L|jsh>$ $L$=0 to $lmax$
c:tex- one block for each $L$ value.  
c:tex- for example $ish$ = p and $L$ = 1 and $jsh$ = p you would get 
c:tex- (3*1*3)+(3*3*3)=36 integrals.
c:tex- The order would be:\\
c:tex-  $<x|x>, <x|y>, <x|z>, <y|x> \ldots <z|z>$ (first nine) \\
c:tex-  $<x|x|x>, <x|x|y>, <x|x|z>, <x|y|x>, <x|y|y> \dots <x|z|z>$ (second nine)\\
c:tex-  $<y|x|x>, <y|x|y>, <y|x|z>, <y|y|x>, <y|y|y> ... <y|z|z>$ (third nine)\\
c:tex-  $<z|x|x>, <z|x|y>, <z|x|z>, <z|y|x>, <z|y|y> ... <z|z|z>$ (fourth nine).
c:tex-
c:tex-
c:tex-
c:tex- See the int_order code inside $\ldots$/NWints/int for specific order of a 
c:tex- set of shells and dipole order. 
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mpole(i_basis, ish, j_basis, jsh, lmax, centerl,
     &    lscr, scr, lmpint, MP)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
c::functions
      logical int_chk_init
      integer int_nint_cart
      external int_chk_init
      external int_nint_cart
c::passed
c:tex-\begin{verbatim}
      integer i_basis             ! [input] basis set handle for ish
      integer ish                 ! [input] i shell/contraction
      integer j_basis             ! [input] basis set handle for jsh
      integer jsh                 ! [input] j shell/contraction
      integer lmax                ! [input] maximum lvalue for multipole integrals in this batch
      double precision centerl(3) ! [input] coordinates of multipole
      integer lscr                ! [input] length of scratch array
      double precision scr(lscr)  ! [input] scratch array
      integer lmpint              ! [input] length of multipole integrals array
      double precision MP(lmpint) ! [output] multipole integrals
c:tex-\end{verbatim}
c::local      
      logical shells_ok
      integer ibas, Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom
      integer jbas, Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_geom
      integer ucont, lval
      integer l_int, ij_int, num_int
      integer int_point, int_size
      logical inline_chk_sh      
      integer lpole
      integer lmpsiz            ! Local copy of lmpint
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c... statement function for int_chk_sh
      inline_chk_sh(ibas,ish) =
     $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(ibas)))
c
c check initialization
c      
      if (.not.int_chk_init('int_mpole'))
     &       call errquit('int_mpole: int_init was not called' ,0)
c
c  check if spherical
c
      call int_nospherical_check(i_basis,'int_mpole:i_basis')
      call int_nospherical_check(j_basis,'int_mpole:j_basis')
      call int_nogencont_check(i_basis,'int_mpole:i_basis')
      call int_nogencont_check(j_basis,'int_mpole:j_basis')
      call int_nospshell_check(i_basis,'int_mpole:i_basis')
      call int_nospshell_check(j_basis,'int_mpole:j_basis')
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
c
      lmpsiz = lmpint
c
      shells_ok = inline_chk_sh(ibas,ish)
      shells_ok = shells_ok .and. inline_chk_sh(jbas,jsh)
      if (.not. shells_ok)
     &       call errquit('int_mpole: invalid contraction/shell',0)
c
***   set defNxyz such that it can handle the maximum multi-pole
c
      lpole = lmax/4 + 1
      call defNxyz(lpole)
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(6,*)'int_mpole: WARNING: possible geometry inconsistency'
        write(6,*)'i_basis geometry handle:',i_geom
        write(6,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
      if (i_gen.gt.1 .or. j_gen.gt.1) then
        write(6,*)' hf3ois does not handle general contractions yet'
        call errquit('int_mpole: general contraction error ',911)
      endif
c
      num_int = 0
      do 00100 lval = 0,lmax
        l_int  = (lval+1)*(lval+2)/2
        ij_int = int_nint_cart(i_basis, ish, j_basis, jsh, 0,0, 0,0)
        int_point = num_int + 1
        int_size  = lmpsiz - num_int
        num_int = num_int + l_int*ij_int
        if (num_int.gt.lmpsiz) then
          write(6,*)' int_mpole: lmpsiz  = ',lmpsiz
          write(6,*)' int_mpole: num_int = ',num_int
          call errquit('int_mpole: lmpint too small ',911)
        endif
        call hf3ois(
     &       coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &       dbl_mb(mb_exndcf(i_icfp,ibas)),i_prim, Li,
     &       centerl,DCexp,
     &       DCcoeff, 1,  lval,
     &       coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &       dbl_mb(mb_exndcf(j_icfp,jbas)),j_prim, Lj,
c...................................TriDiag..DryRun
     &       MP(int_point),int_size,.false.,.false.,scr,lscr)
00100 continue
*     lmpsiz = num_int
      end
