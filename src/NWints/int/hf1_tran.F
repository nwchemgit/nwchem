      subroutine hf1_tran_shift(Gab,Gabcp,ngen,la,lb,nca,ncb)
* $Id$
      implicit none
      integer la, lb
      integer nca, ncb
      integer ngen
      double precision
     &      Gabcp(1:((lb+1)*(lb+2)/2),1:((la+1)*(la+2)/2),ngen)
      double precision
     &      Gab(1:ncb*((lb+1)*(lb+2)/2),1:nca*((la+1)*(la+2)/2))
c
      integer kgen, igen, la2, lb2, ninta, nintb
      integer ica, icb, ia, ib
c
      la2 = (la+1)*(la+2)/2
      lb2 = (lb+1)*(lb+2)/2
*      write(6,*)' nca     = ',nca
*      write(6,*)' ncb     = ',ncb
*      write(6,*)' la/la2  = ',la,la2
*      write(6,*)' lb/lb2  = ',lb,lb2
      igen = 0
      do ica = 1,nca
        do icb = 1,ncb
          kgen = igen + 1
          do ia = 1,la2
            ninta=(ica-1)*la2+ia
            do ib = 1,lb2
              nintb=(icb-1)*lb2+ib
*              write(6,*)'nintb,ninta,ib,ia,kgen,igen',
*     &              nintb,ninta,ib,ia,kgen,igen,ica,icb
              Gab(nintb,ninta) = Gabcp(ib,ia,kgen)
            enddo
          enddo
          igen = igen + 1
        enddo
      enddo
c
*      igen = 0
*      do ica = 1,nca
*        do icb = 1,ncb
*          igen = igen + 1
*          write(6,*)' ica ',ica,' icb ',icb,' igen ',igen
*        enddo
*      enddo
*      do igen = 1,ngen
*        write(6,*)' gabcp matrix for general contraction ',igen
*        call output(gabcp(1,1,igen),1,lb2,1,la2,lb2,la2,1)
*      enddo
*      igen = 0
*      kgen = 0
*      do ica = 1,nca
*        do icb = 1,ncb
*          igen = igen + 1
*          do ia = 1,la2
*            do ib = 1,lb2
*              if (abs(gabcp(ib,ia,igen)).gt.1.0d-08) then
*                kgen = kgen + 1
*              endif
*            enddo
*          enddo
*        enddo
*      enddo
*      write(6,*)' number of non-zero gabcp elements ',kgen
*      write(6,*)' gab transformed '
*      call output(gab,1,(lb2*ncb),1,(la2*nca),(lb2*ncb),(la2*nca),1)
*      kgen = 0
*      do ia = 1,(nca*la2)
*        do ib = 1,(ncb*lb2)
*          if (abs(gab(ib,ia)).gt.1.0d-08) then
*            kgen = kgen + 1
*          endif
*        enddo
*      enddo
*      write(6,*)' number of non-zero gab elements ',kgen
c
      end
      subroutine hf1_tran_gen(gctrana,gctranb,PAIRp,NPP,
     &      Acoefs,NPA,NCA,Bcoefs,NPB,NCB)
      implicit none
c
c      form the general contraction transformation matrix for 1e (2 index quantities)
c
#include "stdio.fh"
c::passed
      integer NPP  ! [input] number of primitive pairs used on center P (=A,B)
      integer NPA  ! [input] number of primitives on center A
      integer NPB  ! [input] number of primitives on center B
      integer NCA  ! [input] number of contractions on center A
      integer NCB  ! [input] number of contractions on center B
      integer PAIRp (2,NPP) ! [input] kept indices of contraction coeffs
      double precision Acoefs(NPA,NCA) ! [input] conctraction coeffs on center A
      double precision Bcoefs(NPB,NCB) ! [input] conctraction coeffs on center B
      double precision gctrana(NPP,NCA) ! [output] general contraction coefs for A multiply
      double precision gctranb(NCB,NPP) ! [output] general contraction coefs for B multiply
c::local
      integer ica,icb,ipp
c
      do 00100 ica = 1,NCA
        do 00200 ipp = 1,NPP
*          write(6,*)' ica,ipp,pair = ',ica, ipp,PAIRP(1,ipp)
          gctrana(ipp,ica) = Acoefs(PAIRP(1,ipp),ica)
00200   continue
00100 continue
      do 00300 icb = 1,NCB
        do 00400 ipp = 1,NPP
*          write(6,*)' icb,ipp,pair = ',icb, ipp,PAIRP(2,ipp)
          gctranb(icb,ipp) = Bcoefs(PAIRP(2,ipp),icb)
00400   continue
00300 continue
      end
