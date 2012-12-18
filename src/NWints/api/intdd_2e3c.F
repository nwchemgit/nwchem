C> \ingroup nwint
C> @{ 
C>
C> \brief Compute 3-center 2-electron integral 2nd derivatives
C>       
C> Compute the 3-center 2-electron integral 2nd derivatives as given by
C> \f{eqnarray*}{
C> \frac{\partial^2({\mu}|{\nu}{\lambda})}{\partial X_x\partial X_y} = \int_{-\infty}^{\infty} \frac{\partial^2 g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}
C> g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})}{\partial X_x\partial X_y}dr_{1}dr_{2}
C> \f}
C> The integral derivatives are stored in an order that is consistent with
C> the declaration `ERI(natcrd*(natcrd+1)/2,nint)`, where `nint` is the number
C> of integrals in the shell triplet, `natcrd` is the number of Cartesian
C> coordinates in the shell triplet, i.e. 9 ( = 3 atoms * 3 coordinates).
C> The actual lexical indeces of the atoms on which the shells `ish`,
C> `jsh` and `ksh` are centered are returned in `idatom`.
C> 
      subroutine intdd_2e3c(brain, ish, ketin, jsh, ksh,
     &       lscr, scr, leri, eri, idatom)
c $Id$
      implicit none
c
* basic api routine to generate 4 center two electron 
* integral second derivatives 
c
* The buffer comes out as:
*      (upper triangle block derivative, basis indecies)
c------------------------------------------------------------
c construct all 10 blocks of sec.der. (output) from 6 blocks:
c
c          AA AB AC AD                AA AB AC 
c             BB BC BD      from         BB BC
c                CC CD                      CC
c                   DD  
c      1-6, 7-15,16-24,25-33         1-6, 7-15,16-24
c          34-39,40-48,49-57             25-30,31-39
c                58-63,64-72                   40-45
c                      73-78
c------------------------------------------------------------
* block AA: (1-6)  | block BB: (34-39)| block CC: (58-63)   |
c       axax=1     |       bxbx=34    |       cxcx=58       |
c       axay=2     |       bxby=35    |       cxcy=59       |
c       axaz=3     |       bxbz=36    |       cxcz=60       |
c       ayay=4     |       byby=37    |       cycy=61       |
c       ayaz=5     |       bybz=38    |       cycz=62       |
c       azaz=6     |       bzbz=39    |       czcz=63       |
* block AB: (7-15) | block BC: (40-48)| block CD: (64-72)   |
c       axbx=7     |       bxcx=40    |       cxdx=64       |
c       axby=8     |       bxcy=41    |       cxdy=65       |
c       axbz=9     |       bxcz=42    |       cxdz=66       |
c       aybx=10    |       bycx=43    |       cydx=67       |
c       ayby=11    |       bycy=44    |       cydy=68       |
c       aybz=12    |       bycz=45    |       cydz=69       |
c       azbx=13    |       bzcx=46    |       czdx=70       |
c       azby=14    |       bzcy=47    |       czdy=71       |
c       azbz=15    |       bzcz=48    |       czdz=72       |
* block AC: (16-24)| block BD: (49-57)| block DD: (73-78)   |
c       axcx=16    |       bxdx=49    |       dxdx=73       |
c       axcy=17    |       bxdy=50    |       dxdy=74       |
c       axcz=18    |       bxdz=51    |       dxdz=75       |
c       aycx=19    |       bydx=52    |       dydy=76       |
c       aycy=20    |       bydy=53    |       dydz=77       |
c       aycz=21    |       bydz=54    |       dzdz=78       |
c       azcx=22    |       bzdx=55    |                     |
c       azcy=23    |       bzdy=56    |                     |
c       azcz=24    |       bzdz=57    |                     |
* block AD: (25-33)|                  |                     |
c       axdx=25    |                  |                     |
c       axdy=26    |                  |                     |
c       axdz=27    |                  |                     |
c       aydx=28    |                  |                     |
c       aydy=29    |                  |                     |
c       aydz=30    |                  |                     |
c       azdx=31    |                  |                     |
c       azdy=32    |                  |                     |
c       azdz=33    |                  |                     |
c------------------------------------------------------------
#include "stdio.fh"
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
c::external subroutines used
c errquit
c::functions
      integer int_nint_cart
      external int_nint_cart
ckw
      integer int_nint
      external int_nint
c
      logical cando_nw
      logical cando_txs
      external cando_nw
      external cando_txs
ckw
c::passed
      integer brain !< [Input] basis set handle for bra basis
      integer ish   !< [Input] lexical contraction index 
      integer jsh   !< [Input] lexical contraction index 
      integer ketin !< [Input] basis set handle for ket basis
      integer ksh   !< [Input] lexical contraction index 
      integer lscr  !< [Input] length of scratch array
      integer leri  !< [Input] length of eri array
      double precision scr(lscr) !< [Scratch] scratch array for integral code.
      double precision eri(leri) !< [Output]  array for two electron integral derivatives.
c NOTE: length of idatom is always 4 because there can be at most 4 centers involved
      integer idatom(4)          !< [Output]  array identifying centers for derivatives
c
c::local
      integer nint, ucont
      integer bra, ket, ab_geom, cd_geom
      integer iatom
      integer jatom
      integer katom
c
      double precision roff(3)
      integer txs_i, txs_j, txs_k, txs_l
      logical status_nw, status_txs
      logical dum_log
      integer nintzero, num_quart, dummy_lab
      double precision q4
ckw
c
      logical used_nw
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      used_nw = .false.
      nint = int_nint_cart(brain,ish,0, 0,ketin,jsh,ketin,ksh)
*
      if (nint*78.gt.leri) then
        write(luout,*) 'nint*78 = ',nint*78
        write(luout,*) 'leri     = ',leri
        call errquit('intdd_2e3c: nint>leri error',911, INT_ERR)
      endif
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
      if (ab_geom.ne.cd_geom) then
        write(luout,*)
     &      'intdd_2e3c.F: two different geometries for',
     &         ' derivatives?'
        call errquit('intdd_2e3c: geom error ',911, GEOM_ERR)
      endif
c
      ucont = (sf_ibs_cn2ucn(ish,bra))
      iatom = (sf_ibs_cn2ce(ish,bra))
      idatom(1)=iatom
      idatom(2)=0
c
      ucont = (sf_ibs_cn2ucn(jsh,ket))
      jatom = (sf_ibs_cn2ce(jsh,ket))
      idatom(3)=jatom
c
      ucont = (sf_ibs_cn2ucn(ksh,ket))
      katom = (sf_ibs_cn2ce(ksh,ket))
      idatom(4)=katom
c
      if (iatom.eq.jatom.and.jatom.eq.katom) then
        call dcopy((nint*78),0.0d00,0,eri,1)
        call ifill(4,-1,idatom,1)
        return
      endif
c
      status_nw  = .false.
      status_txs =.true.
      status_txs = status_txs .and.
     &    cando_txs(brain,ish,0).and.cando_txs(ketin,jsh,ksh)
c
      if (status_txs) then
        call dcopy(3,0.0d00,0,roff,1)
        q4 = 1.0d00
        txs_i = ish
        txs_j = 0
        txs_k = jsh
        txs_l = ksh
        num_quart=1
        dum_log=.false.
c
        call texas_hf2_m(
     &      brain,txs_i,txs_j,
     &      ketin,txs_k,txs_l,num_quart,
     &      q4,.false.,
c...............................use roff set false
     &      roff,roff,roff,roff,.false.,
     &      eri, leri, dummy_lab, dummy_lab, dummy_lab, dummy_lab,
c...............gen labs .. more_integrals
     &      nint, .false., dum_log, scr, lscr, 0.0d0,'der2_int')
        write(6,*)' intdd:2e3c:txshfm:nint = ',nint
        if (nint .eq. 0) then
          nintzero = int_nint(brain,ish,0,0,ketin,jsh,ketin,ksh)
          nintzero = nintzero*78
          call dcopy(nintzero, 0.0d0, 0, eri, 1)
          call ifill(4,-9,idatom,1)
        endif
c
      else
        write(luout,*)'intdd_2e3c: could not use the texas '
        write(luout,*)'           integral second derivatives'
        write(luout,*)' Please notify nwchem-support@emsl.pnl.gov with'
        write(luout,*)'        the input deck and output available'
        call errquit('intdd_2e3c: fatal error',911, INT_ERR)
      endif
c
      end
C> @}
