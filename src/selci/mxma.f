C> \ingroup selci
C> @{
      subroutine selci_mxma(a,iac,iar, b,ibc,ibr, c,icc,icr,
     $     nar,nac,nbc)
*
* $Id$
*
      implicit integer (a-z)
      real*8 a(*),b(*),c(*)
      real*8     zero,     bkj
      parameter (zero=0d0)
c
c  initialize the output matrix c().
c
c      if (iac.le.0 .or. iar.le.0 .or. ibc.le.0 .or. ibr.le.0 .or.
c     $    icc.le.0 .or. icr.le.0 .or. nar.le.0 .or. nac.le.0 .or.
c     $    nbc.le.0) then
c         write(6,*) iac,iar,ibc,ibr,icc,icr,nar,nac,nbc
c         call errquit('mxma: debug test failed',iac)
c      endif
c
      i1j=1
      do 200 j=1,nbc
         ij=i1j
         do 100 i=1,nar
            c(ij)=zero
            ij=ij+icc
100      continue
         i1j=i1j+icr
200   continue
c
      i1j=1
      k1j=1
      do 500 j=1,nbc
         i1k=1
         kj=k1j
         do 400 k=1,nac
            bkj=b(kj)
            kj=kj+ibc
            if(bkj.ne.zero)then
               ij=i1j
               ik=i1k
               do 300 i=1,nar
                  c(ij)=c(ij)+a(ik)*bkj
                  ij=ij+icc
                  ik=ik+iac
300            continue
            endif
            i1k=i1k+iar
400      continue
         i1j=i1j+icr
         k1j=k1j+ibr
500   continue
c
      return
      end
C> @}
