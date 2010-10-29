        subroutine mxtrn(mxin,mxout,idim,jdim)
C     $Id$
        implicit none
        integer idim,jdim
        double precision mxin(idim,jdim),mxout(jdim,idim)
        integer i,j
c
        do i=1,idim
         do j=1,jdim
          mxout(j,i)=mxin(i,j)
         enddo
        enddo
c
        return
        end


    
      subroutine prtfmat(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,nvir)
      implicit none
      integer nvir
      double precision f1n(nvir,nvir),f1t(nvir,nvir),f2n(nvir,nvir),
     &                 f2t(nvir,nvir),f3n(nvir,nvir),f3t(nvir,nvir),
     &                 f4n(nvir,nvir),f4t(nvir,nvir)
c
      integer c,b
c
      print *,' printing out f matrices'
      do c=1,nvir
      do b=1,nvir
      write(6,100)b,c,f1n(b,c),f1t(b,c),f2n(b,c),f2t(b,c),
     &                f3n(b,c),f3t(b,c),f4n(b,c),f4t(b,c)
  100 format(2i5,4e14.5,/,10x,4e14.5)
      enddo
      enddo
c
      return
      end      



      integer function cart_2e4c(basis,ish,jsh,ksh,lsh)
      implicit none
c
#include "bas.fh"
#include "nwc_const.fh"
#include "apiP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
c
      integer basis ! [input] basis set handle
      integer ish   ! [input] shell/contraction index
      integer jsh   ! [input] shell/contraction index
      integer ksh   ! [input] shell/contraction index
      integer lsh   ! [input] shell/contraction index
c
      integer La, a_nbf, Lb, b_nbf, Lc, c_nbf, Ld, d_nbf
      integer usebas, ucont
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      usebas  = basis + BASIS_HANDLE_OFFSET
      ucont   = (sf_ibs_cn2ucn(ish,usebas))
      La      = infbs_cont(CONT_TYPE ,ucont,usebas)
      ucont   = (sf_ibs_cn2ucn(jsh,usebas))
      Lb      = infbs_cont(CONT_TYPE ,ucont,usebas)
      ucont   = (sf_ibs_cn2ucn(ksh,usebas))
      Lc      = infbs_cont(CONT_TYPE ,ucont,usebas)
      ucont   = (sf_ibs_cn2ucn(lsh,usebas))
      Ld      = infbs_cont(CONT_TYPE ,ucont,usebas)
c
      a_nbf = int_nbf_x(La)
      b_nbf = int_nbf_x(Lb)
      c_nbf = int_nbf_x(Lc)
      d_nbf = int_nbf_x(Ld)
c
      cart_2e4c=a_nbf*b_nbf*c_nbf*d_nbf
c
      end
