* $Id$
      subroutine conbuf2_der_ne(firstc,nbls,nbls1,xt1,lt1,lt2,buf2,indx,
     *     aax,bbx,ccx)
c     
c     this is called only for where.eq.'forc'
c     
      implicit real*8 (a-h,o-z)
      logical firstc
c     
      common/obarai/
     *     lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     *     nqi,nqj,nqk,nql,nsij,nskl,
     *     nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c     
#include "texas_lpar.fh"
c     
      dimension indx(*)
      dimension xt1(nbls1*lt1*lt2)
c     
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
c     
      dimension buf2(4*nbls*lt1*lt2)
c
c-----------------------------------------------------------
c     buf2(1,nbls,lt1,lt2) - ordinary contraction
c     buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c     buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c     buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
c     
c     Temporarily adjust array to speed adressing
c
      do i=1,nbls1
         indx(i)=(indx(i)-1)*4
      enddo
c
      if (firstc) then
         do kl=nfu(nqkl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               ijkl1   = (ij-1 + lt1*(kl-1))
               ijklbuf = 4*nbls*ijkl1
               ijklxt  =  nbls1*ijkl1
               do i=1,nbls1
                  ijkl=indx(i)+ijklbuf
                  xxt1 = xt1(i+ijklxt)
                  buf2(1+ijkl)=xxt1
                  buf2(2+ijkl)=xxt1*aax(i)
                  buf2(3+ijkl)=xxt1*bbx(i)
                  buf2(4+ijkl)=xxt1*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nsij)+1,nfu(nsij+1)
               ijkl1   = (ij-1 + lt1*(kl-1))
               ijklbuf = 4*nbls*ijkl1
               ijklxt  =  nbls1*ijkl1
               do i=1,nbls1
                  ijkl=indx(i)+ijklbuf
                  xxt1 = xt1(i+ijklxt)
                  buf2(1+ijkl)=xxt1
                  buf2(2+ijkl)=xxt1*aax(i)
                  buf2(3+ijkl)=xxt1*bbx(i)
                  buf2(4+ijkl)=xxt1*ccx(i)
               enddo
            enddo
         enddo
         firstc=.false.
      else
         do kl=nfu(nqkl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               ijkl1   = (ij-1 + lt1*(kl-1))
               ijklbuf = 4*nbls*ijkl1
               ijklxt  =  nbls1*ijkl1
               do i=1,nbls1
                  ijkl=indx(i)+ijklbuf
                  xxt1 = xt1(i+ijklxt)
                  buf2(1+ijkl)=buf2(1+ijkl)+xxt1
                  buf2(2+ijkl)=buf2(2+ijkl)+xxt1*aax(i)
                  buf2(3+ijkl)=buf2(3+ijkl)+xxt1*bbx(i)
                  buf2(4+ijkl)=buf2(4+ijkl)+xxt1*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nsij)+1,nfu(nsij+1)
               ijkl1   = (ij-1 + lt1*(kl-1))
               ijklbuf = 4*nbls*ijkl1
               ijklxt  =  nbls1*ijkl1
               do i=1,nbls1
                  ijkl=indx(i)+ijklbuf
                  xxt1 = xt1(i+ijklxt)
                  buf2(1+ijkl)=buf2(1+ijkl)+xxt1
                  buf2(2+ijkl)=buf2(2+ijkl)+xxt1*aax(i)
                  buf2(3+ijkl)=buf2(3+ijkl)+xxt1*bbx(i)
                  buf2(4+ijkl)=buf2(4+ijkl)+xxt1*ccx(i)
               enddo
            enddo
         enddo
      endif
c     
c     Undo temporary shift
c
      do i=1,nbls1
         indx(i) = indx(i)/4+1
      enddo
c     
      end
      subroutine conbuf2_der_eq(firstc,nbls,nbls1,xt1,lt1,lt2,buf2,indx,
     *     aax,bbx,ccx)
c     
cccc  this is called only for where.eq.'forc' and nbls=nbls1
c     
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
      logical firstc
c     
      common/obarai/
     *     lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     *     nqi,nqj,nqk,nql,nsij,nskl,
     *     nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c     
#include "texas_lpar.fh"
c     
      dimension indx(*)
      dimension xt1(nbls*lt1*lt2) ! Was nbls1
c     
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
c     
      dimension buf2(4,nbls*lt1*lt2)
c-----------------------------------------------------------
c     buf2(1,nbls,lt1,lt2) - ordinary contraction
c     buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c     buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c     buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
c     
c
      if (nbls1 .ne. nbls) call errquit('conbuf2_der_eq: ne!',nbls,
     &       INT_ERR)
c
      if (firstc) then
         do kl=nfu(nqkl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               ijkl = nbls*(ij-1 + lt1*(kl-1))
               do i=1,nbls
                  buf2(1,i+ijkl)=xt1(i+ijkl)
                  buf2(2,i+ijkl)=xt1(i+ijkl)*aax(i)
                  buf2(3,i+ijkl)=xt1(i+ijkl)*bbx(i)
                  buf2(4,i+ijkl)=xt1(i+ijkl)*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nsij)+1,nfu(nsij+1)
               ijkl = nbls*(ij-1 + lt1*(kl-1))
               do i=1,nbls
                  buf2(1,i+ijkl)=xt1(i+ijkl)
                  buf2(2,i+ijkl)=xt1(i+ijkl)*aax(i)
                  buf2(3,i+ijkl)=xt1(i+ijkl)*bbx(i)
                  buf2(4,i+ijkl)=xt1(i+ijkl)*ccx(i)
               enddo
            enddo
         enddo
         firstc=.false.
      else
         do kl=nfu(nqkl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               ijkl = nbls*(ij-1 + lt1*(kl-1))
               do i=1,nbls
                  buf2(1,i+ijkl)=buf2(1,i+ijkl)+xt1(i+ijkl)
                  buf2(2,i+ijkl)=buf2(2,i+ijkl)+xt1(i+ijkl)*aax(i)
                  buf2(3,i+ijkl)=buf2(3,i+ijkl)+xt1(i+ijkl)*bbx(i)
                  buf2(4,i+ijkl)=buf2(4,i+ijkl)+xt1(i+ijkl)*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nsij)+1,nfu(nsij+1)
               ijkl = nbls*(ij-1 + lt1*(kl-1))
               do i=1,nbls
                  buf2(1,i+ijkl)=buf2(1,i+ijkl)+xt1(i+ijkl)
                  buf2(2,i+ijkl)=buf2(2,i+ijkl)+xt1(i+ijkl)*aax(i)
                  buf2(3,i+ijkl)=buf2(3,i+ijkl)+xt1(i+ijkl)*bbx(i)
                  buf2(4,i+ijkl)=buf2(4,i+ijkl)+xt1(i+ijkl)*ccx(i)
               enddo
            enddo
         enddo
      endif
c     
      end

