* $Id$
c===================================================================
      subroutine transfor(bl,nbls,ixint,
     *                    ityps,jtyps,ktyps,ltyps,
     *                    ilenx,jlenx,klenx,llenx,
     *                    ilens,jlens,klens,llens)
c--------------------------------------------------------------------
c This routine is called only if PNL basis set was specified as 
c                       SPHERICAL
c
c Note: Despite the comment above it is still necessary to check 
c       whether lens.ne.lenx. The reason is that in resolution-of-the-
c       density calculations (e.g. RI-SCF, RI-MP2 and other similar
c       methods) we have 2 basis sets, the AO basis and the fitting
c       basis. These may have different specifications for the 
c       angular momenta, e.g. spherical harmonic AO basis combined
c       with a cartesian fitting basis. In these cases we must transform
c       only the bra or the ket of the integrals but not both. Failing
c       to appreciate this causes rogue data to be pulled in leading to
c       invalid results.
c
c ityps-ltyps PNL type of spherical harmonic shells
c ilenx-llenx are cartesian sizes (see below)
c ilens-llens are spher.harmon. sizes (see below)
c--------------------------------------------------------------------
c Shell  s  p  l   d5  d6  f7  f10  g9 g15  h11 h21  i13 i28
c--------------------------------------------------------------------
c size = 1  3  4   5   6   7   10   9  15   11  21   13  28

c--------------------------------------------------------------------
c itype= 1  2  3   4   5   6   7        8        9       10 ! texas
c size = 1  3  4   5   6   7   10   9  15   11  21   13  28
c--------------------------------------------------------------------
c ityp1= 1  2  3   4   4   5   5    6   6    7   7    8   8 ! texas
c size1= 1  3  4   6   6   10  10   15 15   21  21   28  28

c--------------------------------------------------------------------
c PNL    0  1  -1  2   2   3   3    4   4    5   5    6   6 ! pnl 
c--------------------------------------------------------------------
c siz_c  1  3  4   6   6   10  10   15  15  21   21  28  28
c siz_s  1  3  4   5   5    7   7    9   9  11   11  13  13
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
      dimension bl(*)
c--------------------------------------------------------------------
c bl(ixint) is :
c xint(nbls,ilenx*jlenx*klenx*llenx) - input (cart) integrals
c xint(nbls,           *           ) - output(sphe) integrals
c
c--------------------------------------------------------------------
c check if transformation is nedded for this block of integrals:
c
      if( max(ityps,jtyps,ktyps,ltyps).le.1 ) RETURN
c
c     write(6,*)'PNL types:',ityps,jtyps,ktyps,ltyps
c     write(6,*)'PNL sizeS:',ilens,jlens,klens,llens
c     write(6,*)'PNL sizeC:',ilenx,jlenx,klenx,llenx
c
c--------------------------------------------------------------------
      ndim_bra=ilenx*jlenx
      ndim_ket=klenx*llenx
      n_braket=ndim_bra*ndim_ket
c--------------------------------------------------------------------
c allocate memory for scratch :
c
      call getmem(n_braket, iscrt)
c--------------------------------------------------------------------
c
c transpose integral array :  X(nbls,n_braket) --> Y(n_braket,nbls)
c
      call getmem(nbls*n_braket, iyint)
c
      lda=nbls
      ldb=n_braket
      call trspmo(bl(ixint),lda,  bl(iyint),ldb)
c--------------------------------------------------------------------
c do transformation for BRA <ij| :
c
      max_bra=max(ityps,jtyps)
      if( max_bra.gt.1.and.(jlenx+ilenx).ne.(jlens+ilens) ) then
         iaddress=iyint
         do 100 nn=1,nbls
            call spcart_bra2etran(bl(iaddress),bl(iscrt),
     *                            jlenx,ilenx,
     *                            jlens,ilens,
     *                            jtyps,ityps,
     *                            1    ,1    ,
     *                            ndim_ket,   
     *                            .false.)  
ccc  *                            .true. )  ! print 
            iaddress=iaddress+n_braket
  100    continue
         ndim_bra=ilens*jlens
      endif
c--------------------------------------------------------------------
c do transformation for KET |kl> :
c
      max_ket=max(ktyps,ltyps)
      if( max_ket.gt.1.and.(llenx+klenx).ne.(llens+klens) ) then
         iaddress=iyint
         do 200 nn=1,nbls
            call spcart_ket2etran(bl(iaddress),bl(iscrt),
     *                            llenx,klenx,
     *                            llens,klens,
     *                            ltyps,ktyps,
     *                            1    ,1    ,
     *                            ndim_bra, 
     *                            .false.)  
ccc  *                            .true. )  ! print 
            iaddress=iaddress+n_braket
  200    continue
      endif
c--------------------------------------------------------------------
c transpose back : Y(n_braket,nbls) --> X(nbls,n_braket)
c
c     lda=nbls
c     ldb=n_braket
c
      call trspmo(bl(iyint),ldb,  bl(ixint),lda)
c--------------------------------------------------------------------
c release allocated memory :
c
      call retmem(2)
c--------------------------------------------------------------------
      end
c===================================================================
      subroutine get_spher_pnl_type(itypx,ityps)
c--------------------------------------------------------------------
c This routine returns PNL-type of SPHERICAL HARMONIC shell's (ityps)
c using TXS CARTESIAN types itypx
c
c--------------------------------------------------------------------
c Shell  s  p  l   d    f    g    h    i  
c--------------------------------------------------------------------
c PNL_S  0  1  -1  2    3    4    5    6  ! bas_spherical=.true.
c PNL_C  0  1  -1  2    3    4    5    6  ! bas_spherical=.false.
c--------------------------------------------------------------------
c
c Shell  s  p  l   d5  d6  f7   f10  g15  h21  i28
c--------------------------------------------------------------------
c TXS-   1  2  3   4   5   6    7    8    9    10   ! type
c TXS-1  1  2  3   4   4   5    5    6    7     8   ! type1
c--------------------------------------------------------------------
c siz_c  1  3  4   6   6   10   10   15   21   28 
c siz_s  1  3  4   5   5    7   7    9    11   13 
c--------------------------------------------------------------------
      dimension ish_type(8)
c
c TXS shells         s  p  l   d    f    g    h    i
c TXS_types          1  2  3   4    5    6    7    8 
c
      data ish_type/ 0, 1, -1, 2,   3,   4,   5,    6  /
c--------------------------------------------------------------------
      if(itypx.gt.8) then
         call errquit('texas tried to transform orbitals higher than I',
     $        0, INT_ERR)
      endif
c
      ityps=ish_type(itypx)
c--------------------------------------------------------------------
      end
c===================================================================
