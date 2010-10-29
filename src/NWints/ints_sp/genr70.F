      subroutine genr70(bra, ish, xi, jsh, xj, ket, ksh, xk, lsh, xl,
     $     eri)
c $Id$
c
c     Generic interface to the gaussian two electron integral package.
c     Evaluates those integrals involving only s, p and sp functions.
c
c     Internally has to order shells to map to routines
c     sp0000  sp0001  sp0011  sp0101  sp0111  sp1111
c
c     sinfo obtains geometrical information about the four centers
c
      implicit none
#include "errquit.fh"
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "shllfo.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
      integer bra, ket          ! [input] bra/ket basis handles
      integer ish, jsh, ksh, lsh ![input] indices of contractions
      double precision xi(3), xj(3), xk(3), xl(3) ! [input] coords of centers
      double precision eri(*)   ! [output] returns 4-d array of integrals
c
      integer itype, jtype      ! used to id 0000, 0001, ...
      integer ngout             ! no. integrals produced
      integer ib(4)             ! the map from input order to internal
      integer ibinv(4)          ! the inverse map to ib
      integer bas(4)            ! used to permute basis handles
      integer typ(4)            ! used to permute shell types
      integer ush(4)            ! used to permute unique shell indices
      integer i, iush, jush, kush, lush, ibas, jbas, kbas, lbas
      integer ngi, ngj, ngk, ngl, iptc, jptc, kptc, lptc, ngij, ngkl
      integer ipte, jpte, kpte, lpte, isht, jsht, ksht, lsht
      double precision x(3,4)   ! used to permute coordinates
      double precision gout(256) ! Internal integral buffer
c
#ifdef BASIS_DEBUG
      logical bas_check_handle
      external bas_check_handle
#endif
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c
c     check basis handles, prepare to get data from basis structures
c
#ifdef BASIS_DEBUG
      if (.not. (bas_check_handle(bra,'genr70 bra') .and.
     $     bas_check_handle(ket,'genr70 ket'))) call errquit
     $     ('genr70: invalid bra/ket basis ', 0, BASIS_ERR)
#endif
      bas(1) = bra + BASIS_HANDLE_OFFSET
      bas(2) = bra + BASIS_HANDLE_OFFSET
      bas(3) = ket + BASIS_HANDLE_OFFSET
      bas(4) = ket + BASIS_HANDLE_OFFSET
c
c     Map from user shell/contractions to internal unique contractions
c     and get the shell types
c
      ush(1) = (sf_ibs_cn2ucn(ish,bas(1)))
      ush(2) = (sf_ibs_cn2ucn(jsh,bas(2)))
      ush(3) = (sf_ibs_cn2ucn(ksh,bas(3)))
      ush(4) = (sf_ibs_cn2ucn(lsh,bas(4)))
      typ(1) = infbs_cont(CONT_TYPE,ush(1),bas(1))
      typ(2) = infbs_cont(CONT_TYPE,ush(2),bas(2))
      typ(3) = infbs_cont(CONT_TYPE,ush(3),bas(3))
      typ(4) = infbs_cont(CONT_TYPE,ush(4),bas(4))
      ib(1)  = 1
      ib(2)  = 2
      ib(3)  = 3
      ib(4)  = 4
      do i = 1, 3
         x(i,1) = xi(i)
         x(i,2) = xj(i)
         x(i,3) = xk(i)
         x(i,4) = xl(i)
      enddo
c
c     Put data into arrays so that 
c
c     Figure out the necessary permutation
c
      itype = 8*abs(typ(1)) + 4*abs(typ(2)) + 2*abs(typ(3)) +
     $     abs(typ(4)) + 1
c
      go to (120,120,160,120,180,120,160,120,220,140,200,140,180,180,
     +     240,120),itype
      call errquit('genr70: invalid shell types', itype, BASIS_ERR)
c
c     !! Note that some of these permutations are DIFFERENT from
c     !! from those in GAMESS/GAUSSIAN due to ishell etc. not being used
c     !! and the use of ibinv below
c
c     types 0000,0001,0101,0011,0111,1111 are unaltered
  120 go to 260
c     types 1001,1011 have ij switched
 140  ib(1) = 2
      ib(2) = 1
      go to 260
c     types 0010,0110 have kl switched
 160  ib(3) = 4
      ib(4) = 3
      go to 260
c     types 0100,1100,1101 have pairs ij and kl switched
 180  ib(1) = 3
      ib(2) = 4
      ib(3) = 1
      ib(4) = 2
      go to 260
c     type 1010 has ij switched and kl switched
 200  ib(1) = 2
      ib(2) = 1
      ib(3) = 4
      ib(4) = 3
      go to 260
c     type 1000  has pairs ij and kl switched followed by kl switch
 220  ib(1) = 3
      ib(2) = 4
      ib(3) = 2
      ib(4) = 1
      go to 260
c     type 1110 has pairs ij and kl switched followed by ij switch
 240  ib(1) = 4
      ib(2) = 3
      ib(3) = 1
      ib(4) = 2
  260 continue
c
c     Want ij contraction level to be the greatest
c
      ngi  = infbs_cont(CONT_NPRIM,ush(ib(1)),bas(ib(1)))
      ngj  = infbs_cont(CONT_NPRIM,ush(ib(2)),bas(ib(2)))
      ngk  = infbs_cont(CONT_NPRIM,ush(ib(3)),bas(ib(3)))
      ngl  = infbs_cont(CONT_NPRIM,ush(ib(4)),bas(ib(4)))
      ngij = ngi*ngj
      ngkl = ngk*ngl
c
c     only 6 standard types remain. 0000,0001,0011,0101,0111,1111
c     specify these by jtype
c
      go to (280,300,300,320,300,340,340,360,300,340,340,360,320,360,
     +     360,380),itype
  280 jtype = 1
      ngout = 1
      if (ngij .lt. ngkl) then
         ibinv(1) = ib(1)
         ibinv(2) = ib(2)
         ib(1) = ib(3)
         ib(2) = ib(4)
         ib(3) = ibinv(1)
         ib(4) = ibinv(2)
      endif
      go to 400
  300 jtype = 2
      ngout = 4
      go to 400
  320 jtype = 3
      ngout = 16
      go to 400
  340 jtype = 4
      ngout = 64
      go to 400
  360 jtype = 5
      ngout = 64
      go to 400
  380 jtype = 6
      ngout = 256
      if (ngij .lt. ngkl) then
         ibinv(1) = ib(1)
         ibinv(2) = ib(2)
         ib(1) = ib(3)
         ib(2) = ib(4)
         ib(3) = ibinv(1)
         ib(4) = ibinv(2)
      endif
  400 continue
c
c     Zero the internal integral buffer
c
      call dfill(ngout, 0.0d0, gout, 1)
c
c     Permute indices and get data from basis data structure
c
      iush = ush(ib(1))         ! Unique shell labels
      jush = ush(ib(2))
      kush = ush(ib(3))
      lush = ush(ib(4))
      ibas = bas(ib(1))         ! Basis set indexes
      jbas = bas(ib(2))
      kbas = bas(ib(3))
      lbas = bas(ib(4))
      isht = typ(ib(1))         ! Unique shell types
      jsht = typ(ib(2))
      ksht = typ(ib(3))
      lsht = typ(ib(4))
      ngi  = infbs_cont(CONT_NPRIM,iush,ibas) ! No. of primtives
      ngj  = infbs_cont(CONT_NPRIM,jush,jbas)
      ngk  = infbs_cont(CONT_NPRIM,kush,kbas)
      ngl  = infbs_cont(CONT_NPRIM,lush,lbas)
      iptc = infbs_cont(CONT_ICFP ,iush,ibas) ! Coeff. ptr.
      jptc = infbs_cont(CONT_ICFP ,jush,jbas)
      kptc = infbs_cont(CONT_ICFP ,kush,kbas)
      lptc = infbs_cont(CONT_ICFP ,lush,lbas)
      ipte = infbs_cont(CONT_IEXP ,iush,ibas) ! Exp. ptr
      jpte = infbs_cont(CONT_IEXP ,jush,jbas)
      kpte = infbs_cont(CONT_IEXP ,kush,kbas)
      lpte = infbs_cont(CONT_IEXP ,lush,lbas)
c
c     Sinfo computes info about axis rotation etc.
c
      call sinfo(
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     $      iush, isht, dbl_mb(mb_exndcf(iptc,ibas)),
     &      dbl_mb(mb_exndcf(ipte,ibas)),
     $      ngi, x(1,ib(1)),
     $      jush, jsht, dbl_mb(mb_exndcf(jptc,jbas)),
     &      dbl_mb(mb_exndcf(jpte,jbas)),
     $      ngj, x(1,ib(2)),
     $      kush, ksht, dbl_mb(mb_exndcf(kptc,kbas)),
     &      dbl_mb(mb_exndcf(kpte,kbas)),
     $      ngk, x(1,ib(3)),
     $      lush, lsht, dbl_mb(mb_exndcf(lptc,lbas)),
     &      dbl_mb(mb_exndcf(lpte,lbas)),
     $      ngl, x(1,ib(4)), jtype)
c
      go to(440,480,500,520,540,560),jtype
  440 call sp0000(eri)          ! Fast return for (ss|ss)
      return
 480  call sp0001(gout)
      go to 1080
 500  call sp0011(gout)
      go to 1080
 520  call sp0101(gout)
      go to 1080
 540  call sp0111(gout)
      go to 1080
 560  call sp1111(gout)
      goto 1080
c
1080  continue
      isht = typ(1)
      jsht = typ(2)
      ksht = typ(3)
      lsht = typ(4)
      ibinv(ib(1)) = 1
      ibinv(ib(2)) = 2
      ibinv(ib(3)) = 3
      ibinv(ib(4)) = 4
c
      call qout70(gout, eri, isht, jsht, ksht, lsht,
     $     ibinv(1), ibinv(2), ibinv(3), ibinv(4))
c
      return
      end
