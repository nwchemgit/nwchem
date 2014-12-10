      function i2addr(iorb,jorb,korb,lorb,ijklof,noccsym,no12sym)
*
* obtain address of integral (iorb jorb ! korb lorb) in molcas order 
* iorb jorb korb lorb corresponds to symmetry ordered indeces !!
* integrals assumed in core 
*
* ITRA_ROUTE switch added May 2011
*
*
      implicit real*8(a-h,o-z)
*
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'cintfo.inc'
      include 'multd2h.inc'
      INCLUDE 'crun.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      
*
      dimension ijklof(nsmob,nsmob,nsmob)
      logical isymj,ksyml,isymk,jsyml,ijsymkl,iksymjl
*. 
      ntest = 00
      noccsym_l = noccsym
      if (i12s.eq.0.and.i34s.eq.0) noccsym_l = 1
*
      iabs = iorb
      ism = ismfto(ireost(iorb))
      ioff = ibso(ism)
*
      jabs = jorb
      jsm = ismfto(ireost(jorb))
      joff = ibso(jsm)
*
      kabs = korb
      ksm = ismfto(ireost(korb))
      koff = ibso(ksm)
*
      labs = lorb
      lsm = ismfto(ireost(lorb))
      loff = ibso(lsm)
*
      if( ntest.ge. 100) then
        write(6,*) ' gmijkl at your service '
        write(6,*) ' iorb iabs ism ioff ',iorb,iabs,ism,ioff
        write(6,*) ' jorb jabs jsm joff ',jorb,jabs,jsm,joff
        write(6,*) ' korb kabs ksm koff ',korb,kabs,ksm,koff
        write(6,*) ' lorb labs lsm loff ',lorb,labs,lsm,loff
      end if
*
c test
      ijsm = multd2h(ism,jsm)
      klsm = multd2h(ksm,lsm)
      ijklsm = multd2h(ijsm,klsm)
      if (ijklsm.ne.1) then
        i2addr = -1
        return
      end if
c test end
*
      if (noccsym_l.eq.0.and.
     &     (jsm.gt.ism .or. ( ism.eq.jsm .and. jabs.gt.iabs))) then
        isym=jsm
        jsym=ism
        i = jabs - joff + 1
        j = iabs - ioff + 1
      else
        isym=ism
        jsym=jsm
        i = iabs - ioff + 1
        j = jabs - joff + 1
      end if
      if (noccsym_l.eq.0) then
        ijblk=jsym+isym*(isym-1)/2
      else
        ijblk = (isym-1)*nsmob + jsym
      end if
      if ( noccsym_l.eq.0 .and.
     &     (lsm.gt.ksm  .or. ( ksm.eq.lsm .and. labs.gt.kabs)) ) then
        ksym=lsm
        lsym=ksm
        k = labs -loff + 1
        l = kabs - koff + 1
      else
        ksym=ksm
        lsym=lsm
        k = kabs - koff + 1
        l = labs -loff + 1
      end if
      if ( noccsym_l.eq.0) then
        klblk=lsym+ksym*(ksym-1)/2
      else
        klblk = (ksym-1)*nsmob + lsym
      end if
*
      if ( klblk.gt.ijblk .and. no12sym.eq.0 ) then
        itemp=isym
        isym=ksym
        ksym=itemp
        itemp=jsym
        jsym=lsym
        lsym=itemp
        itemp=ijblk
        ijblk=klblk
        klblk=itemp
*
        itemp = i
        i = k
        k = itemp
        itemp = j
        j = l
        l = itemp
      end if
      if(ntest .ge. 100 ) then
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' isym,jsym,ksym,lsym',isym,jsym,ksym,lsym
      end if
*
*  define offset for given symmetry block
      ibloff = ijklof(isym,jsym,ksym)
      if(ntest .ge. 100 )
     &write(6,*) ' ibloff isym jsym ksym ', ibloff,isym,jsym,ksym
      if (noccsym_l.eq.0) then
        isymj=isym.eq.jsym
        ksyml=ksym.eq.lsym
        isymk=(no12sym.eq.0).and.isym.eq.ksym
        jsyml=(no12sym.eq.0).and.jsym.eq.lsym
        ijsymkl=isymj.and.ksyml
        iksymjl=isymk.and.jsyml
      else
        isymj=.false.
        ksyml=.false.
        isymk=.false.
        jsyml=.false.
        ijsymkl=.false.
        iksymjl=(no12sym.eq.0).and.(isym.eq.ksym).and.(jsym.eq.lsym)
      end if
*
      itorb=ntoobs(isym)
      jtorb=ntoobs(jsym)
      ktorb=ntoobs(ksym)
      ltorb=ntoobs(lsym)
c?    print *,' itorb,jtorb,ktorb,ltorb',itorb,jtorb,ktorb,ltorb
      if ( isymj ) then
        ijpairs=itorb*(itorb+1)/2
        ij=j+i*(i-1)/2
      else
        ijpairs=itorb*jtorb
        IF(ITRA_ROUTE.EQ.1) THEN
          ij=j + (i-1)*jtorb
        ELSE
          ij=i + (j-1)*itorb
        END IF
      end if 
*
      if(ksyml ) then
        klpairs=ktorb*(ktorb+1)/2
        kl=l+k*(k-1)/2
      else
        klpairs=ktorb*ltorb
        IF(ITRA_ROUTE.EQ.1) THEN
          kl=l+(k-1)*ltorb
        ELSE
          KL = K + (L-1)*KTORB
        END IF
      end if
c?    print *,' ijpairs,klpairs',ijpairs,klpairs
*
      if ( iksymjl ) then
        if ( ij.gt.kl ) then
          kloff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijpairs-kloff
        else
          ijoff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klpairs-ijoff
        end if
      else
        ijkl=ij+(kl-1)*ijpairs
      end if
      if( ntest .ge. 100 )
     & write(6,*) ' ijkl ', ijkl
*
      i2addr = ibloff-1+ijkl
      if ( ntest .ge. 100 ) then
      write(6,*) 'i j k l ', i,j,k,l
      write(6,*) ' ibloff ijkl ',ibloff,ijkl
        write(6,*) ' i2addr  = ', i2addr
      end if
*
      return
      end 
*
      function i2addr2(iorb,jorb,korb,lorb,ijklof,i12,i34,i1234)
*
* obtain address of 4index quantity (iorb,jorb;korb,lorb)  
* iorb jorb korb lorb corresponds to symmetry ordered indeces !!
*
* i12   (0,-1,1) symmetry wrt permutation of iorb and jorb
* i34   (0,-1,1) symmetry wrt permutation of korb and lorb
* i1234 (0,-1,1) symmetry wrt permutation of pairs iorb,jorb and korb,lorb
*
* ITRA_ROUTE switch added May 2011
*
      implicit real*8(a-h,o-z)
*
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'cintfo.inc'
      include 'multd2h.inc'
      include 'frorbs.inc'
      INCLUDE 'crun.inc'
      
*
      dimension ijklof(nsmob,nsmob,nsmob)
      logical isymj,ksyml,isymk,jsyml,ilsymjk,iksymjl,swapij,swapkl
*. 
      ntest = 00
*
      iabs = iorb
      isym = ismfto(ireost(iorb))
      ioff = ibso(isym)
      i = iabs - ioff + 1 - nfrobs(isym)
*
      jabs = jorb
      jsym = ismfto(ireost(jorb))
      joff = ibso(jsym)
      j = jabs - joff + 1 - nfrobs(jsym)
*
      kabs = korb
      ksym = ismfto(ireost(korb))
      koff = ibso(ksym)
      k = kabs - koff + 1 - nfrobs(ksym)
*
      labs = lorb
      lsym = ismfto(ireost(lorb))
      loff = ibso(lsym)
      l = labs -loff + 1 - nfrobs(lsym)
*
      if (i.le.0.or.j.le.0.or.k.le.0.or.l.le.0) then
        i2addr2 = -1
        return
      end if

*
      if( ntest.ge. 100) then
        write(6,*) ' gmijkl at your service '
        write(6,*) ' iorb iabs isym ioff ',iorb,iabs,isym,ioff
        write(6,*) ' jorb jabs jsym joff ',jorb,jabs,jsym,joff
        write(6,*) ' korb kabs ksym koff ',korb,kabs,ksym,koff
        write(6,*) ' lorb labs lsym loff ',lorb,labs,lsym,loff
      end if
*
c test
      ijsym = multd2h(isym,jsym)
      klsym = multd2h(ksym,lsym)
      ijklsym = multd2h(ijsym,klsym)
      if (ijklsym.ne.1) then
        if (ntest.ge.100)
     &       write(6,*) 'WARNING: symmetry error in i2addr2'
        i2addr2 = -2
        return
      end if
c test end

      if(ntest .ge. 100 ) then
        write(6,*) ' as entered:'
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' isym,jsym,ksym,lsym',isym,jsym,ksym,lsym
      end if
*
*     resort index quadruple to standard sequence
*
      ijblk = (min(isym,jsym)-1)*nsmob + max(isym,jsym)
      klblk = (min(ksym,lsym)-1)*nsmob + max(ksym,lsym)

      ijdx = (min(iabs,jabs)-1)*ntoob + max(iabs,jabs)
      kldx = (min(kabs,labs)-1)*ntoob + max(kabs,labs)
      if ( i1234.ne.0 .and.
     &     (ijblk.gt.klblk .or.
     &      (ijblk.eq.klblk .and. ijdx.gt.kldx) ) ) then
        itemp = isym
        isym = ksym
        ksym = itemp
        itemp = jsym
        jsym = lsym
        lsym = itemp
        itemp = i
        i = k
        k = itemp
        itemp = j
        j = l
        l = itemp
      end if

      swapij = ( i12.ne.0 .and. (isym.gt.jsym .or.
     &           (isym.eq.jsym .and. i.gt.j) )     ) 
      swapkl = ( i34.ne.0 .and. (ksym.gt.lsym .or.
     &           (ksym.eq.lsym .and. k.gt.l) )     ) 
      swapij = swapij.or.( i34.ne.0 .and. i12.eq.0 .and. swapkl )
      swapkl = swapkl.or.( i12.ne.0 .and. i34.eq.0 .and. swapij )

      swapij = swapij.or.( i34.ne.0.and.i12.eq.0 .and.
     &                     ksym.eq.lsym .and. k.eq.l
     &              .and.(isym.gt.jsym .or.(isym.eq.jsym .and. i.gt.j)))

      swapkl = swapkl.or.( i12.ne.0.and.i34.eq.0 .and.
     &                     isym.eq.jsym .and. i.eq.j
     &              .and.(ksym.gt.lsym .or.(ksym.eq.lsym .and. k.gt.l)))

      if (swapij) then
        itemp = isym
        isym = jsym
        jsym = itemp
        itemp = i
        i = j
        j = itemp
      end if

      if (swapkl) then
        itemp = ksym
        ksym = lsym
        lsym = itemp
        itemp = k
        k = l
        l = itemp
      end if
*
      if(ntest .ge. 100 ) then
        write(6,*) ' resorted to:'
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' isym,jsym,ksym,lsym',isym,jsym,ksym,lsym
      end if
*
*  offset for given symmetry block
      ibloff = ijklof(isym,jsym,ksym)

      if(ntest .ge. 100 )
     &     write(6,*) ' ibloff isym jsym ksym ', ibloff,isym,jsym,ksym

      isymj=(i12.ne.0).and.isym.eq.jsym
      ksyml=(i34.ne.0).and.ksym.eq.lsym
      isymk=(i1234.ne.0).and.isym.eq.ksym
      jsyml=(i1234.ne.0).and.jsym.eq.lsym
      ilsymjk=(i12.ne.0.and.i34.eq.0.and.i1234.ne.0).and.
     &         isym.eq.lsym.and.jsym.eq.ksym
      iksymjl=isymk.and.jsyml
*
      itorb=ntaobs(isym)
      jtorb=ntaobs(jsym)
      ktorb=ntaobs(ksym)
      ltorb=ntaobs(lsym)
c?    print *,' itorb,jtorb,ktorb,ltorb',itorb,jtorb,ktorb,ltorb
      if ( isymj.and.i12.eq.1 ) then
        ijpairs=itorb*(itorb+1)/2
        ij=i+j*(j-1)/2
        ijdia = 1
      else if ( isymj.and.i12.eq.-1 ) then
        if (i.eq.j) then
          i2addr2 = -1
          return
        end if
        ijpairs=itorb*(itorb-1)/2
        ij=i+(j-1)*(j-2)/2
        ijdia = -1
      else
        ijpairs=itorb*jtorb
        IF(ITRA_ROUTE.EQ.1) THEN
          ij=i + (j-1)*itorb
        ELSE
          ij=j + (i-1)*jtorb
        END IF
        ijdia = 0
      end if 
*
      if(ksyml.and.i34.eq.1) then
        klpairs=ktorb*(ktorb+1)/2
        kl=k+l*(l-1)/2
        kldia = 1
      else if(ksyml.and.i34.eq.-1) then
        if (k.eq.l) then
          i2addr2 = -1
          return
        end if
        klpairs=ktorb*(ktorb-1)/2
        kl=k+(l-1)*(l-2)/2
        kldia = -1
      else
        klpairs=ktorb*ltorb
        IF(ITRA_ROUTE.EQ.1) THEN
          kl=k+(l-1)*ktorb
        ELSE
          kl=l+(k-1)*ltorb
        END IF
        kldia = 0
      end if
c?    print *,' ijpairs,klpairs',ijpairs,klpairs
*
      if ( iksymjl .and. i1234.eq.1 ) then
        if (ijdia.eq.kldia) then 
          if ( ij.le.kl ) then
            ijkl=(kl-1)*kl/2 + ij
          else
            ijkl=(ij-1)*ij/2 + kl
          end if
        else
          stop 'incomplete i2addr2'
        end if
      else if ( iksymjl .and. i1234.eq.-1 ) then
        if (ijdia.eq.kldia) then
          if ( ij.lt.kl ) then
            ijkl=(kl-2)*(kl-1)/2 + ij
          else if (ij.gt.kl) then
            ijkl=(ij-2)*(ij-1)/2 + kl
          else
            ! not allowed
            i2addr2 = -1
            return
          end if
        else if (ijdia.eq.1.and.kldia.eq.0) then
          ! ij is upper-triangular packed
          ! kl is full matrix
          ! case 1:
          ! i.le.j and k.le.l
          if (k.le.l) then
            kl_tri = k+l*(l-1)/2            
            if (ij.gt.kl_tri) then
              ijkl = (ij-1)*(ij-2)/2 + kl_tri
            else if (ij.lt.kl_tri) then
              ijkl = (kl_tri-1)*(kl_tri-2)/2 + ij
            else
              ! not allowed
              i2addr2 = -1
              return
            end if
          ! case 2:
          ! i.lt.j and k.gt.l
          else if (i.lt.j.and.k.gt.l) then
            ijkl_off = ijpairs*(ijpairs-1)/2
            ij_tri = (j-1)*(j-2)/2+i
            kl_tri = (k-1)*(k-2)/2+l
            if (ij_tri.le.kl_tri) then
              ijkl = ijkl_off+kl_tri*(kl_tri-1)/2+ij_tri
            else
              ijkl = ijkl_off+ij_tri*(ij_tri-1)/2+kl_tri
            end if
          ! forbidden case:
          else
            i2addr2 = -1
            return
          end if
          
        else if (ijdia.eq.0.and.kldia.eq.1) then
          nnmx = itorb
            stop 'adapt'
          ! case 1:
          ! i.le.j and k.le.l
          if (k.le.l) then
            kl_tri = l+k*(k-1)/2            
            if (ij.gt.kl_tri) then
              ijkl = ij*(ij-1)/2 + kl_tri
            else if (ij.lt.kl_tri) then
              ijkl = kl_tri*(kl_tri-1)/2 + ij
            else
              ! not allowed
              i2addr2 = -1
              return
            end if
          ! case 2:
          ! i.lt.j and k.gt.l
          else if (i.lt.j.and.k.gt.l) then
            ijkl_off = ijpairs*(ijpairs-1)/2
            ij_tri = (j-1)*(j-2)/2+i
            kl_tri = (k-1)*(k-2)/2+l
            if (ij_tri.le.kl_tri) then
              ijkl = ijkl_off+kl_tri*(kl_tri-1)/2+ij_tri
            else
              ijkl = ijkl_off+ij_tri*(ij_tri-1)/2+kl_tri
            end if
          ! forbidden case:
          else
            i2addr2 = -1
            return
          end if

        else
          stop 'incomplete i2addr2'
        end if
      else if (ksym.eq.lsym.and.ijdia.eq.1.and.kldia.eq.0) then
        if (i.lt.j) then
          ij_tri = (j-1)*(j-2)/2+i
          ij_tripairs = itorb*(itorb-1)/2
          ijkl = ij_tri + (kl-1)*ij_tripairs
        else if (i.eq.j) then
          ijkl_off = itorb*(itorb-1)/2*klpairs
          kl_tri = l*(l-1)/2+k
          ijkl = ijkl_off + i + (kl_tri-1)*itorb
        else
          stop 'unexpected'
        end if
      else if (isym.eq.jsym.and.ijdia.eq.0.and.kldia.eq.1) then
        if (k.lt.l) then
          kl_tri = (l-1)*(l-2)/2+k
          kl_tripairs = ktorb*(ktorb-1)/2
          ijkl = kl_tri + (ij-1)*kl_tripairs
        else if (k.eq.l) then
          ijkl_off = ktorb*(ktorb-1)/2*ijpairs
          ij_tri = j*(j-1)/2+i
          ijkl = ijkl_off + k + (ij_tri-1)*ktorb
        else
          stop 'unexpected'
        end if
        
      else if (ilsymjk) then
        ! i<j, k>l
        IF(ITRA_ROUTE.EQ.1) THEN
          ij = (j-1)*itorb + i
          kl = (k-1)*ltorb + l
        ELSE
          ij = (i-1)*jtorb + j
          kl = (l-1)*ktorb + k
        END IF
        
        if (ij.le.kl) then
          ijkl = kl*(kl-1)/2+ij
        else
          ijkl = ij*(ij-1)/2+kl
        end if
      else
        ijkl=ij+(kl-1)*ijpairs
      end if
      if( ntest .ge. 100 )
     & write(6,*) ' ijkl ', ijkl
*
      i2addr2 = ibloff-1+ijkl
      if( ntest .ge. 100 ) then
      write(6,*) 'i j k l ', i,j,k,l
      write(6,*) ' ibloff ijkl ',ibloff,ijkl
        write(6,*) ' i2addr2  = ', i2addr2
      end if
*
      return
      end 

      subroutine ijkl2iadr(ijkl,iadr,nadr,ntoob,ireost,
     &                     ijklof,i12,i34,i1234)
c
c     return the orbital quadruples for a list of addresses listed on iadr()
c      
      implicit none

      integer, intent(in) ::
     &     nadr, iadr(nadr), ntoob,
     &     i12, i34, i1234, ijklof(*), ireost(*)

      integer, intent(out) ::
     &     ijkl(4,nadr)

      integer ::
     &     idx, jdx, kdx, ldx, ii, iad

      integer, external ::
     &     i2addr2
      
      ijkl(1:4,1:nadr) = 0

      do idx = 1, ntoob
        do jdx = 1, ntoob
          do kdx = 1, ntoob
            do ldx = 1, ntoob
              iad = i2addr2(idx,jdx,kdx,ldx,ijklof,i12,i34,i1234)
              if (iad.lt.0) cycle
              do ii = 1, nadr
                if (iad.eq.iadr(ii).and.ijkl(1,ii).eq.0) then
                  ijkl(1,ii) = ireost(idx)
                  ijkl(2,ii) = ireost(jdx)
                  ijkl(3,ii) = ireost(kdx)
                  ijkl(4,ii) = ireost(ldx)
                end if
              end do
            end do
          end do
        end do
      end do

      return
      
      end 
      FUNCTION GTIJKL(I,J,K,L)
*
* Obtain integral (I J ! K L )
*
* reads oper and will get similarity transformed integrals if
* I_USE_SIMTRH==1
* if i_unrorb.eq.1 it will obey ispcas (both on oper.inc)
*
* I,J,K L refers to orbitals in  Type ordering
*                         ==============
*
* from GTIJKL_SM_AB, replacing old GTIJKL
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'oper.inc'
*
      XIJKL = 0.0D0
*
      ! set up offset and pointer array
      IF (I_USE_SIMTRH.EQ.0) THEN
        NOCCSYM=0
        IF (I_UNRORB.EQ.0.OR.ISPCAS.EQ.1) THEN
          NO12SYM=0
          K2ADR = KINT2
          KP2ADR= KPINT2
        ELSE IF (I_UNRORB.EQ.1.AND.ISPCAS.EQ.2) THEN
          NO12SYM=0
          K2ADR = KINT2BB
          KP2ADR= KPINT2
        ELSE IF (I_UNRORB.EQ.1.AND.(ISPCAS.EQ.3.OR.ISPCAS.EQ.4)) THEN
          NO12SYM=1
          K2ADR = KINT2AB
          KP2ADR= KPINT2AB
        ELSE
          WRITE(6,*) 'unknown case: i_unrorb, ispcas: ',i_unrorb, ispcas
          STOP 'GTIJKL'
        END IF
      ELSE
        NOCCSYM=1
        IF (I_UNRORB.EQ.0) THEN
          NO12SYM=0
          K2ADR = KINT2_SIMTRH
          KP2ADR= KPINT2_SIMTRH
        ELSE IF (I_UNRORB.EQ.1.AND.ISPCAS.EQ.1) THEN
          NO12SYM=0
          K2ADR = KINT2_SIMTRH_AA
          KP2ADR= KPINT2_SIMTRH
        ELSE IF (I_UNRORB.EQ.1.AND.ISPCAS.EQ.2) THEN
          NO12SYM=0
          K2ADR = KINT2_SIMTRH_BB
          KP2ADR= KPINT2_SIMTRH
        ELSE IF (I_UNRORB.EQ.1.AND.(ISPCAS.EQ.3.OR.ISPCAS.EQ.4)) THEN
          NO12SYM=1
          K2ADR = KINT2_SIMTRH_AB
          KP2ADR= KPINT2_SIMTRH_AB
        END IF
          
      END IF
      IF (.NOT.(I_UNRORB.EQ.1.AND.ISPCAS.EQ.4)) THEN
        IADR = I2ADDR(IREOTS(I),IREOTS(J),
     &              IREOTS(K),IREOTS(L),
     &              int_mb(KP2ADR),NOCCSYM,NO12SYM)
      ELSE
        IADR = I2ADDR(IREOTS(K),IREOTS(L),
     &              IREOTS(I),IREOTS(J),
     &              int_mb(KP2ADR),NOCCSYM,NO12SYM)
      END IF

      IF (IADR.GT.0) THEN
        XIJKL = dbl_mb(K2ADR-1+IADR)
      ELSE
        XIJKL = 0D0
      END IF
*
      GTIJKL = XIJKL
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' 2e integral for I,J,K,L = ', I,J,K,L
        WRITE(6,*) ' is ', XIJKL
      END IF
      RETURN
      END
      SUBROUTINE PTIJKL(I,J,K,L,XINT,XLIST)
*
* Put integral (I J ! K L ) on XINT to its correct place in XLIST
*
*  a quick hack: never tested, never used :-))
*
* reads oper and will behave correspondingly
*
* I,J,K L refers to active orbitals in  Type ordering
*                                      ==============
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'oper.inc'

      DIMENSION XLIST(*)
*
      ! set up offset and pointer array
      IF (I_USE_SIMTRH.EQ.0) THEN
        NOCCSYM=0
        IF (I_UNRORB.EQ.0.OR.ISPCAS.EQ.1) THEN
          NO12SYM=0
          KP2ADR= KPINT2
        ELSE IF (I_UNRORB.EQ.1.AND.ISPCAS.EQ.2) THEN
          NO12SYM=0
          KP2ADR= KPINT2
        ELSE IF (I_UNRORB.EQ.1.AND.(ISPCAS.EQ.3.OR.ISPCAS.EQ.4)) THEN
          NO12SYM=1
          KP2ADR= KPINT2AB
        ELSE
          WRITE(6,*) 'unknown case: i_unrorb, ispcas: ',i_unrorb, ispcas
          STOP 'GTIJKL'
        END IF
      ELSE
        NOCCSYM=1
        IF (I_UNRORB.EQ.0) THEN
          NO12SYM=0
          KP2ADR= KPINT2_SIMTRH
        ELSE IF (I_UNRORB.EQ.1.AND.ISPCAS.EQ.1) THEN
          NO12SYM=0
          KP2ADR= KPINT2_SIMTRH
        ELSE IF (I_UNRORB.EQ.1.AND.ISPCAS.EQ.2) THEN
          NO12SYM=0
          KP2ADR= KPINT2_SIMTRH
        ELSE IF (I_UNRORB.EQ.1.AND.(ISPCAS.EQ.3.OR.ISPCAS.EQ.4)) THEN
          NO12SYM=1
          KP2ADR= KPINT2_SIMTRH_AB
        END IF
          
      END IF
      IF (.NOT.(I_UNRORB.EQ.1.AND.ISPCAS.EQ.4)) THEN
        IADR = I2ADDR(IREOTS(I),IREOTS(J),
     &              IREOTS(K),IREOTS(L),
     &              dbl_mb(KP2ADR),NOCCSYM,NO12SYM)
      ELSE
        IADR = I2ADDR(IREOTS(K),IREOTS(L),
     &              IREOTS(I),IREOTS(J),
     &              dbl_mb(KP2ADR),NOCCSYM,NO12SYM)
      END IF

      IF (IADR.GT.0) THEN
        XLIST(IADR) = XINT
      ELSE
        IF (XINT.GT.1D3*EPSILON(1D0)) THEN
          WRITE(6,*) 'WARNING: INTEGRAL ',I,J,K,L,
     &         ' IS NOT TOTALLY SYMMETRIC, THUS IGNORED!'
          WRITE(6,*) 'VALUE = ',XINT
        END IF
      END IF
*
      RETURN
      END
      FUNCTION GTIJKL_GN(IORB,JORB,KORB,LORB)
*
* Obtain integral (IORB JORB ! KORB LORB) from current
* active array of integrals. It is assumed that the
* current list is a complete integral list
*
* IST allows switching between symmetry(1)- and type(2)-ordered orbitals
* 
* Jeppe Olsen, July 2011
* 
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
      NTEST = 000
      IF(NTEST.GE.200) THEN
       WRITE(6,*)  ' Output from GT_IJKL_GN '
       WRITE(6,'(A,4I4)') ' IORB, JORB, KORB, LORB = ',
     &                      IORB, JORB, KORB, LORB   
      END IF
*
*
      KINT2_LA = KINT2_A(IE2ARRAY_A)
      KPINT2_LA = KPINT2_A(IE2ARRAY_A)
*
      IST = 2
      I2ADDR = I2ADDR_GN(IORB,JORB,KORB,LORB,
     &         int_mb(KPINT2_LA),IST)
*
      IF(I2ADDR.GT.0) THEN
        XINT = dbl_mb(KINT2_LA-1+I2ADDR)
      ELSE
        WRITE(6,*) ' Negative integral address '
        WRITE(6,*) ' IORB, JORB, KORB, LORB, I2ADDR = ',
     &               IORB, JORB, KORB, LORB, I2ADDR
        WRITE(6,*) ' IE2ARRAY_A, KPINT2_LA = ', 
     &               IE2ARRAY_A, KPINT2_LA
        STOP ' Negative integral address '
      END IF
*
      GTIJKL_GN = XINT
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*)  ' Output from GT_IJKL_GN '
       WRITE(6,'(A,4I4)') ' IORB, JORB, KORB, LORB = ',
     &                      IORB, JORB, KORB, LORB   
       WRITE(6,'(A,I8,E15.8)') 
     & 'Address and value of integral ', I2ADDR, GTIJKL_GN
      END IF
*
      RETURN
      END


*
* Obtain 
      FUNCTION I2ADDR_GN(IORB,JORB,KORB,LORB,IJKLOF,IST)
*
* obtain address of integral (iorb jorb ! korb lorb) using 
* new integral order..
* IST = 1 => Symmetry ordered orbital indices
* IST = 2 => Type     ordered orbital indices
*
* Permutational symmetry is defined by I12S_A,I34S_A,I1234S_A coming through
* CINTFO
*
*
      INCLUDE 'implicit.inc'
*
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'cintfo.inc'
      include 'multd2h.inc'
      INCLUDE 'crun.inc'
      
*
      dimension ijklof(nsmob,nsmob,nsmob)
      logical isymj,ksyml,ijsymkl
*. 
      ntest = 00
*. Reform to symmetry ordered indices if required
      IF(IST.EQ.2) THEN
        IABS = IREOTS(IORB)
        JABS = IREOTS(JORB)
        KABS = IREOTS(KORB)
        LABS = IREOTS(LORB)
      ELSE
        IABS = IORB
        JABS = JORB
        KABS = KORB
        LABS = LORB
      END IF
*
      IF(IST.EQ.1) THEN
        ism = ismfto(ireost(iorb))
      ELSE
        ISM = ISMFTO(IORB)
      END IF
      ioff = ibso(ism)
*
      IF(IST.EQ.1) THEN
        jsm = ismfto(ireost(jorb))
      ELSE
        JSM = ISMFTO(JORB)         
      END IF
      joff = ibso(jsm)
*
      IF(IST.EQ.1) THEN
        ksm = ismfto(ireost(korb))
      ELSE
        KSM = ISMFTO(KORB)         
      END IF
      koff = ibso(ksm)
*
      IF(IST.EQ.1) THEN
        lsm = ismfto(ireost(lorb))
      ELSE
        LSM = ISMFTO(LORB)        
      END IF
      loff = ibso(lsm)
*
      if( ntest.ge. 100) then
        write(6,*) ' gmijkl at your service '
        write(6,*) ' iorb iabs ism ioff ',iorb,iabs,ism,ioff
        write(6,*) ' jorb jabs jsm joff ',jorb,jabs,jsm,joff
        write(6,*) ' korb kabs ksm koff ',korb,kabs,ksm,koff
        write(6,*) ' lorb labs lsm loff ',lorb,labs,lsm,loff
*
        WRITE(6,*)  '  I12S_A,I34S_A,I1234S_A = ',
     &                 I12S_A,I34S_A,I1234S_A
      end if
*. Test symmetry
      ijsm = multd2h(ism,jsm)
      klsm = multd2h(ksm,lsm)
      ijklsm = multd2h(ijsm,klsm)
      if (ijklsm.ne.1) then
        i2addr = -1
        return
      end if
* 
      if (I12S_A.EQ.1.and.
     &     (jsm.gt.ism .or. ( ism.eq.jsm .and. jabs.gt.iabs))) then
        isym=jsm
        jsym=ism
        i = jabs - joff + 1
        j = iabs - ioff + 1
      else
        isym=ism
        jsym=jsm
        i = iabs - ioff + 1
        j = jabs - joff + 1
      end if
      if (I12S_A.eq.1) then
        ijblk=jsym+isym*(isym-1)/2
      else
        ijblk = (isym-1)*nsmob + jsym
      end if
      if ( I34S_A.EQ.1 .and.
     &     (lsm.gt.ksm  .or. ( ksm.eq.lsm .and. labs.gt.kabs)) ) then
        ksym=lsm
        lsym=ksm
        k = labs -loff + 1
        l = kabs - koff + 1
      else
        ksym=ksm
        lsym=lsm
        k = kabs - koff + 1
        l = labs -loff + 1
      end if
      if (I34S_A .eq.1) then
        klblk=lsym+ksym*(ksym-1)/2
      else
        klblk = (ksym-1)*nsmob + lsym
      end if
*
      if ( klblk.gt.ijblk .and. I1234S_A.EQ.1 ) then
        itemp=isym
        isym=ksym
        ksym=itemp
        itemp=jsym
        jsym=lsym
        lsym=itemp
        itemp=ijblk
        ijblk=klblk
        klblk=itemp
*
        itemp = i
        i = k
        k = itemp
        itemp = j
        j = l
        l = itemp
      end if
      if(ntest .ge. 100 ) then
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' isym,jsym,ksym,lsym',isym,jsym,ksym,lsym
      end if
*
*  define offset for given symmetry block
      ibloff = ijklof(isym,jsym,ksym)
      if(ntest .ge. 100 )
     &write(6,*) ' ibloff isym jsym ksym ', ibloff,isym,jsym,ksym
      isymj=.false.
      ksyml=.false.
      ijsymkl=.false.
      if (I12S_A.eq.1) isymj=isym.eq.jsym
      IF(I34S_A.EQ.1) KSYML = KSYM.EQ.LSYM
      IF(I1234S_A.EQ.1) IJSYMKL = (ISYM.EQ.KSYM).AND.(JSYM.EQ.LSYM)
*
      itorb=ntoobs(isym)
      jtorb=ntoobs(jsym)
      ktorb=ntoobs(ksym)
      ltorb=ntoobs(lsym)
c?    print *,' itorb,jtorb,ktorb,ltorb',itorb,jtorb,ktorb,ltorb
      if ( isymj ) then
        ijpairs=itorb*(itorb+1)/2
        ij=j+i*(i-1)/2
      else
        ijpairs=itorb*jtorb
        ij=i + (j-1)*itorb
      end if 
*
      if(ksyml ) then
        klpairs=ktorb*(ktorb+1)/2
        kl=l+k*(k-1)/2
      else
        klpairs=ktorb*ltorb
        KL = K + (L-1)*KTORB
      end if
c?    print *,' ijpairs,klpairs',ijpairs,klpairs
*
      if ( IJSYMKL ) then
        if ( ij.gt.kl ) then
          kloff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijpairs-kloff
        else
          ijoff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klpairs-ijoff
        end if
      else
        ijkl=ij+(kl-1)*ijpairs
      end if
      if( ntest .ge. 100 )
     & write(6,*) ' ijkl ', ijkl
*
      I2ADDR_GN = ibloff-1+ijkl
      if( ntest .ge. 100 ) then
      write(6,*) 'i j k l ', i,j,k,l
      write(6,*) ' ibloff ijkl ',ibloff,ijkl
        write(6,*) ' I2ADDR_GN  = ', I2ADDR_GN
      end if
*
      return
      end 
c $Id$
