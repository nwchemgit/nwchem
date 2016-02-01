c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute the multipole integrals up to order `lmax`
C>
C> The general form is \f$\langle shell | pole | shell \rangle\f$. 
C> Integrals are returned in shell blocks of \f$<ish|L|jsh>\f$ where
C> \f$L\f$ runs from 0 to \f$l_{max}\f$
C> one block for each \f$L\f$ value.  
C>
C> The integral for each \f$L\f$ value computed is:
C> \f{eqnarray*}{
C> ({\mu}|{\hat{L}}|{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X,r_{1})\hat{L}g_{\lambda}(X,r_{1})dr_{1}
C> \f}
C>
C> For example, if `ish` = \f$p\f$ and \f$L\f$ = 1 and `jsh` = \f$p\f$ you
C> would get 
C> (3*1*3)+(3*3*3)=36 integrals.
C> The order would be:
C> <table>
C> <tr><td> \f$<x|x>\f$ </td><td> \f$<x|y>\f$ </td><td> \f$<x|z>\f$ </td><td> \f$<y|x>\f$ </td><td> \f$<y|y>\f$ </td><td> \f$\ldots\f$ </td><td> \f$<z|z>\f$ </td><td> (first nine) </td></tr>
C> <tr><td> \f$<x|x|x>\f$ </td><td> \f$<x|x|y>\f$ </td><td> \f$<x|x|z>\f$ </td><td> \f$<y|x|x>\f$ </td><td> \f$<y|x|y>\f$ </td><td> \f$\ldots\f$ </td><td> \f$<z|x|z>\f$ </td><td> (second nine) </td></tr>
C> <tr><td> \f$<x|y|x>\f$ </td><td> \f$<x|y|y>\f$ </td><td> \f$<x|y|z>\f$ </td><td> \f$<y|y|x>\f$ </td><td> \f$<y|y|y>\f$ </td><td> \f$\ldots\f$ </td><td> \f$<z|y|z>\f$ </td><td> (third nine) </td></tr>
C> <tr><td> \f$<x|z|x>\f$ </td><td> \f$<x|z|y>\f$ </td><td> \f$<x|z|z>\f$ </td><td> \f$<y|z|x>\f$ </td><td> \f$<y|z|y>\f$ </td><td> \f$\ldots\f$ </td><td> \f$<z|z|z>\f$ </td><td> (fourth nine) </td></tr>
C> </table>
C>


c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mpole}
c:tex-This routine returns multipole integrals up to the level $lmax$
c:tex-
c:tex- The general form is $< shell | pole | shell >$. 
c:tex- Integrals are returned in shell blocks of $<ish|L|jsh>$ $L$=0 to $lmax$
c:tex- one block for each $L$ value.  
c:tex-
c:tex- For example $ish$ = p and $L$ = 1 and $jsh$ = p you would get 
c:tex- (3*1*3)+(3*3*3)=36 integrals.
c:tex- The order would be:\\
c:tex-\begin{tabular}{cccccccc}
c:tex-  $<x|x>$   & $<x|y>$   & $<x|z>$   & $<y|x>$   & $<y|y>$  & $\ldots$ & $<z|z>$   & (first nine) \\
c:tex-  $<x|x|x>$ & $<x|x|y>$ & $<x|x|z>$ & $<x|y|x>$ & $<x|y|y>$ & $\ldots$ & $<x|z|z>$ & (second nine)\\
c:tex-  $<y|x|x>$ & $<y|x|y>$ & $<y|x|z>$ & $<y|y|x>$ & $<y|y|y>$ & $\ldots$ & $<y|z|z>$ & (third nine)\\
c:tex-  $<z|x|x>$ & $<z|x|y>$ & $<z|x|z>$ & $<z|y|x>$ & $<z|y|y>$ & $\ldots$ & $<z|z|z>$ & (fourth nine)
c:tex-\end{tabular}
c:tex-
c:tex-The integral for each L value computed is:
c:tex-\begin{eqnarray*}
c:tex-({\mu}{\hat{L}}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X,r_{1})\hat{L}g_{\lambda}(X,r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-
c:tex-
c:tex- See the int\_order code inside $\ldots$/NWints/int for specific order of a 
c:tex- set of shells and dipole order. 
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mpole(i_basis, ish, j_basis, jsh, lmax, centerl,
     &    lscr, scr, lmpint, MP)
c:tex-\end{verbatim}
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "apiP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "rel_nwc.fh"
c
c::functions
      logical int_chk_init
      integer int_nint_cart, int_nint
      external int_chk_init
      external int_nint_cart, int_nint
c::passed
c:tex-\begin{verbatim}
      integer i_basis             !< [Input] basis set handle for ish
      integer ish                 !< [Input] i shell/contraction
      integer j_basis             !< [Input] basis set handle for jsh
      integer jsh                 !< [Input] j shell/contraction
      integer lmax                !< [Input] maximum lvalue for 
                                  !< multipole integrals in this batch
      double precision centerl(3) !< [Input] coordinates of multipole
      integer lscr                !< [Input] length of scratch array
      double precision scr(lscr)  !< [Input] scratch array
      integer lmpint              !< [Input] length of multipole ints
      double precision MP(lmpint) !< [Output] multipole integrals
c:tex-\end{verbatim}
c::local      
      logical any_spherical
      logical shells_ok, orel, oirel, ojrel, oNR
      integer ibas, isbas, jbas, jsbas, lbas, sbas, abas
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom, i_icfpS
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_geom, j_icfpS
      integer i_gen_save, j_gen_save
      integer ucont, uconts, lval
      integer l_int, ij_int, ij_int_x, num_int
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
     &       call errquit('int_mpole: int_init was not called' ,0, 
     &             INT_ERR)
c
c  check if spherical
c
      call int_nogencont_check(i_basis,'int_mpole:i_basis')
      call int_nogencont_check(j_basis,'int_mpole:j_basis')
      call int_nospshell_check(i_basis,'int_mpole:i_basis')
      call int_nospshell_check(j_basis,'int_mpole:j_basis')
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
c
      oNR = .true.
      oirel = .false.
      ojrel = .false.
      orel = .false.
c
      if (dyall_mod_dir) then
c
c     get basis set handles; relativistic integral option only valid
c     if both ibas and jbas are the ao basis.
c
        lbas = lc_bsh + BASIS_HANDLE_OFFSET
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        orel = ibas .eq. abas .and. jbas .eq. abas
      end if
c
      lmpsiz = lmpint
c
      shells_ok = inline_chk_sh(ibas,ish)
      shells_ok = shells_ok .and. inline_chk_sh(jbas,jsh)
      if (.not. shells_ok)
     &       call errquit('int_mpole: invalid contraction/shell',0, 
     &        BASIS_ERR)
c
***   set defNxyz such that it can handle the maximum multi-pole
c
      lpole = lmax/4 + 1
      if (dyall_mod_dir) lpole = lpole+1
      call defNxyz(lpole)
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,ibas) .ne. 0)) then
        oirel = .true.
        isbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_mpole: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          ibas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        isbas = ibas
      end if
c
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
      i_icfpS = infbs_cont(CONT_ICFP ,uconts,isbas)
c
c   j shell
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,jbas) .ne. 0)) then
        ojrel = .true.
        jsbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1estv: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          jbas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        jsbas = jbas
      end if
c
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
      j_icfpS = infbs_cont(CONT_ICFP ,uconts,jsbas)
c
      oNR = .not.(oirel.and.ojrel)
      orel = oirel.or.ojrel
c
      if (i_geom.ne.j_geom) then
#if defined(USE_OPENMP)
!$omp critical(serialize_warnp)
#endif
        if (WarnP.eq.0) then
          write(luout,*)
     &      'int_mpole: WARNING: possible geometry inconsistency'
          write(luout,*)'i_basis geometry handle:',i_geom
          write(luout,*)'j_basis geometry handle:',j_geom
          WarnP = 1
        endif
#if defined(USE_OPENMP)
!$omp end critical(serialize_warnp)
#endif
      endif
c
      if (i_gen.gt.1 .or. j_gen.gt.1) then
        write(luout,*)
     &      ' hf3ois does not handle general contractions yet'
        call errquit('int_mpole: general contraction error ',911,
     &            INT_ERR)
      endif
c
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
c
      num_int = 0
      do 00100 lval = 0,lmax
        l_int  = int_nbf_x(lval) ! (lval+1)*(lval+2)/2
        ij_int_x = int_nint_cart(i_basis, ish, j_basis, jsh, 0,0, 0,0)
        ij_int   = int_nint(i_basis, ish, j_basis, jsh, 0,0, 0,0)
        int_point = num_int + 1
        int_size  = lmpsiz - num_int
        if ((num_int+(l_int*ij_int_x)).gt.lmpsiz) then
          write(luout,*)' int_mpole: lmpsiz                     = ',
     &        lmpsiz
          write(luout,*)' int_mpole: l_int                      = ',
     &        l_int
          write(luout,*)' int_mpole: ij_int_x                   = ',
     &        ij_int_x
          write(luout,*)' int_mpole: ij_int                     = ',
     &        ij_int
          write(luout,*)' int_mpole: lval                       = ',
     &        lval
          write(luout,*)' int_mpole: (num_int+(l_int*ij_int_x)) = ',
     &        (num_int+(l_int*ij_int_x))
          call errquit('int_mpole: lmpint too small ',911, INT_ERR)
        endif
        num_int = num_int + l_int*ij_int
        if (orel) then
          call rel_mpole (
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfpS,isbas)),i_prim, Li,
     &        centerl,DCexp,
     &        DCcoeff, 1,  lval,
     &        coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfps,jsbas)),j_prim, Lj,
c....................................TriDiag......DryRun
     &        MP(int_point),int_size,.false.,oNR,.false.,
     &        scr,lscr,rel_dbg,rel_typ)
        else
          call hf3ois(
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),i_prim, Li,
     &        centerl,DCexp,
     &        DCcoeff, 1,  lval,
     &        coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfp,jbas)),j_prim, Lj,
c...................................TriDiag..DryRun
     &        MP(int_point),int_size,.false.,.false.,scr,lscr)
        end if
        if (any_spherical) then
          if (Li.eq.-1) then
            i_gen_save = i_gen
            i_gen = 1
          endif
          if (Lj.eq.-1) then
            j_gen_save = i_gen
            j_gen = 1
          endif
          call spcart_3ctran(MP(int_point),scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x(lval),int_nbf_x(lval),lval,1,.false.,
     &        int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &        .false.)
          if (Li.eq.-1) i_gen = i_gen_save
          if (Lj.eq.-1) j_gen = j_gen_save
        endif
00100 continue
      end
C> @}
