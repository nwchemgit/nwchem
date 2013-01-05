c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute projected multipole integrals up to the specified level
C>
C> This routine computes projected multipole integrals up to level 
C> \f$l_{max}\f$ (\f$0 \rightarrow l_{max}\f$).
C> The general form is \f$\langle pole | shell \rangle\f$.
C>  
C> Integrals are returned in \f$\langle pole | shell \rangle\f$ blocks
C> one block for each \f$L\f$ value \f$0 \rightarrow l_{max}\f$. 
C> 
C> For example, a multipole, L = 1 and a D shell 
C> would yield (1+3)*6 = 24 integrals.
C> 
C> The order would be:
C> <table>
C> <tr><td> <0|xx> </td><td> <0|xy> </td><td> <0|xz> </td><td> <0|yy> </td><td> \f$\ldots\f$ </td><td> <0|zz> </td><td> first six </td></tr>
C> <tr><td> <x|xx> </td><td> <x|xy> </td><td> <x|xz> </td><td> <x|yy> </td><td> \f$\ldots\f$ </td><td> <x|zz> </td><td> second six </td></tr>
C> <tr><td> <y|xx> </td><td> <y|xy> </td><td> <y|xz> </td><td> <y|yy> </td><td> \f$\ldots\f$ </td><td> <y|zz> </td><td> third six </td></tr>
C> <tr><td> <z|xx> </td><td> <z|xy> </td><td> <z|xz> </td><td> <z|yy> </td><td> \f$\ldots\f$ </td><td> <z|zz> </td><td> third six </td></tr>
C> </table>
C>
c:tex-% This is part of the standard API routines
c:tex-\subsection{int\_projpole}
c:tex-This routine computes projected multipole integrals up to level 
c:tex-lmax (0$\rightarrow$lmax):\\
c:tex-The general form is $<$pole$|$shell$>$
c:tex- 
c:tex-Integrals are returned in $<$pole$|$shell$>$ blocks one block for each 
c:tex-L value 0$\rightarrow$lmax. 
c:tex-
c:tex-For example, a multipole, L = 1 and a d shell 
c:tex-would yield (1+3)*6 = 24 integrals.
c:tex-
c:tex-The order would be:\\
c:tex-\begin{tabular}{ccccccc}
c:tex-$<0|xx>$ & $<0|xy>$ & $<0|xz>$ & $<0|yy>$ & $\ldots$ & $<0|zz>$ & first six\\
c:tex-$<x|xx>$ & $<x|xy>$ & $<x|xz>$ & $<x|yy>$ & $\ldots$ & $<x|zz>$ & second six\\
c:tex-$<y|xx>$ & $<y|xy>$ & $<y|xz>$ & $<y|yy>$ & $\ldots$ & $<y|zz>$ & third six \\
c:tex-$<z|xx>$ & $<z|xy>$ & $<z|xz>$ & $<z|yy>$ & $\ldots$ & $<z|zz>$ & fourth six
c:tex-\end{tabular}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_projpole(i_basis, ish, centerl, lmax, 
     &    lscr, scr, lmpint, MP)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
#include "int_nbf.fh"
c
c::functions
      logical int_chk_init
      logical cando_nw_1e, cando_nw
      integer int_nint_cart, int_nint
      external int_chk_init
      external cando_nw_1e, cando_nw
      external int_nint_cart, int_nint
c::passed
c:tex-\begin{verbatim}
      integer i_basis             !< [Input] basis set handle for ish
      integer ish                 !< [Input] i shell/contraction
      integer lmax                !< [Input] maximum lvalue for 
                                  !< multipole ints in this batch
      double precision centerl(3) !< [Input] coordinates of multipole
      integer lscr                !< [Input] length of scratch array
      double precision scr(lscr)  !< [Input] scratch array
      integer lmpint              !< [Input/Output] length of 
                                  !< multipole integrals array
      double precision MP(lmpint) !< [Output] multipole integrals
c:tex-\end{verbatim}
c::local      
      logical shells_ok
      integer ibas, Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom
      integer i_gen_save
      integer ucont, lval
      integer l_int, i_int, i_int_x,num_int
      integer int_point, int_size
      integer lpole
      logical inline_chk_sh      
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c... statement function for int_chk_sh
      inline_chk_sh(ibas,ish) =
     $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(ibas)))
c
c check initialization
c      
      if (.not.int_chk_init('int_projpole'))
     &       call errquit('int_projpole: int_init was not called' ,0, 
     &          INT_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
c
      shells_ok = inline_chk_sh(ibas,ish)
      if (.not. shells_ok)
     &       call errquit('int_projpole: invalid contraction/shell',0, 
     &            BASIS_ERR)
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
      num_int = 0
      do 00100 lval = 0,lmax
        l_int   = (lval+1)*(lval+2)/2
        i_int_x = int_nint_cart(i_basis, ish, 0, 0, 0,0, 0,0)
        i_int   = int_nint(i_basis, ish, 0, 0, 0,0, 0,0)
        int_point = num_int + 1
        int_size  = lmpint - num_int
        if ((num_int+l_int*i_int_x).gt.lmpint) then
          write(luout,*)
     &        ' int_projpole: lmpint  = ',lmpint
          write(luout,*)
     &        ' int_projpole: (num_int+l_int*i_int_x) = ',
     &        (num_int+l_int*i_int_x)
          call errquit('int_projpole: lmpint too small ',911, INT_ERR)
        endif
        num_int = num_int + l_int*i_int
        if (cando_nw_1e(i_basis,ish,0)) then
          call int_hf1sp(
     &        centerl,DCexp,
*                                   ! no ECPs allowed when projpole?
     &        DCcoeff, 1, 1, lval, -1, 
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li, i_cent,
     &        coords(1,1,i_geom),charge(1,i_geom),
     &        geom_invnucexp(1,i_geom),ncenter(i_geom),
     &        MP(int_point),scr,scr,int_size,
c............. doS    doT     doV     canonical
     &        .true.,.false.,.false.,.false.,
c.............. dryrun
     &        .false.,scr,lscr,'int_projpole')
          if (bas_spherical(ibas)) then    ! only one basis set
            if (Li.eq.-1) then
              i_gen_save = i_gen
              i_gen = 1
            endif
            call spcart_2ctran(MP(int_point),scr,lscr,
     &          int_nbf_x(lval),int_nbf_x(lval),lval,1,.false.,
     &          int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,
     &          bas_spherical(ibas),.false.)
            if (Li.eq.-1) i_gen = i_gen_save
          endif
        else
          call errquit
     &        ('int_projpole: could not do sp or nw integrals',0, 
     &           INT_ERR)
        endif
c
00100 continue
      end
C> @}
