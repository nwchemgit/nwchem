c
c
c   Generates vertex and arc tables for 
c   a full CI wavefunction
c
c   Vertex table is the standard definition, i.e.
c   number of pathways to the vertex from head. Zero
c   value means vertex cannot be reached for that
c   symmetry.
c
c   Arc weight table is non-negative for any allowed
c   arc between vertices. The arc weight is just the
c   diametrically-opposite vertex value. The strings
c   computed will be relative within each symmetry.
c
c   The absolute string index is obtained by adding
c   on the offset stored in arctab(1,1,isym)
c
c   RJH !! changed to arctab(0,0,isym)
c   
c
      subroutine detci_vatable( norb, nel, nsym, osym,
     $                          vtxtab, arctab )
*
* $Id$
*
      implicit none
#include "detciP.fh"
      integer norb
      integer nel
      integer nsym
      integer osym(norb)
      integer vtxtab(0:norb,0:nel,nsym)
      integer arctab(0:norb,0:nel,nsym)
c
#include "symmdef.fh"
      integer iorb, ello, el, weight, offset, iel
      integer isym, jsym, ksym
#include "bitops.fh"
#include "symmmul.fh"
c
c
c
      call ifill((norb+1)*(nel+1)*nsym,0,vtxtab,1)
      vtxtab(0,0,1) = 1
      do iorb = 1, norb
        ello = iorb + nel - norb
        if (ello.le.0) ello = 0
        do el=ello,nel
          do isym=1,nsym
            weight = vtxtab((iorb-1),el,isym)
C            jsym = ieor((osym(iorb)-1),(isym-1)) + 1
            jsym = MULX(osym(iorb),isym)
            if (el.gt.0) weight = weight + vtxtab((iorb-1),(el-1),jsym)
            vtxtab(iorb,el,isym) = weight
          enddo
        enddo
      enddo

c$$$      WRITE(6,941)
c$$$ 941  FORMAT(//,'VERTEX TABLES',/)
c$$$      DO ISYM=1,NSYM
c$$$        WRITE(6,911) ISYM
c$$$        DO IORB=0,NORB
c$$$          WRITE(6,901) (VTXTAB(IORB,EL,ISYM),EL=0,NEL)
c$$$        ENDDO
c$$$      ENDDO

      offset = 0
      call ifill(((norb+1)*(nel+1)*nsym),-1,arctab,1) 
      do jsym=1,nsym
        do iorb=1,norb
          isym=osym(iorb)
          do iel=1,nel
            if (vtxtab(iorb,iel,jsym).ne.0) then
C              ksym = ieor((jsym-1),(isym-1)) + 1
              ksym = MULX(jsym,isym)
              if (vtxtab(iorb-1,iel-1,ksym).ne.0) then
                arctab(iorb,iel,jsym) = vtxtab((iorb-1),iel,jsym)
              endif
            endif
          enddo
        enddo
c       arctab(1,1,jsym) = offset
        arctab(0,0,jsym) = offset
        offset = offset + vtxtab(norb,nel,jsym)
      enddo

c$$$      WRITE(6,944)
c$$$ 944  FORMAT(//,'ARC WEIGHT TABLE',/)
c$$$      DO ISYM=1,NSYM
c$$$        WRITE(6,911) ISYM
c$$$        DO IORB=0,NORB
c$$$          WRITE(6,901) (ARCTAB(IORB,EL,ISYM),EL=1,NEL)
c$$$        ENDDO
c$$$      ENDDO
c$$$ 911  format('Symmetry:',i5)
c$$$ 901  format(16i5)

      return
      end






c
c
c
c  Compute string index from occupation vector
c  Start at head (vaccuum) and tranverse down graph
c  Each index in elv() gives the traverse, update symmetry and
c  look up weight from arctab()
c
c  Note: symmtery offset is stored in arctab(1,1,isym)
c  
c     RJH now arctab(0,0,isym)
c
c
c
      integer function detci_elv2str( norb, nel, nsym, osym, 
     $                                arctab, elv )
      implicit none
      integer norb
      integer nel
      integer nsym
      integer osym(norb)
      integer arctab(0:norb,0:nel,nsym)
      integer elv(nel)
      integer isym, iel, iorb, wght
#include "symmdef.fh"
#include "bitops.fh"
#include "symmmul.fh"

      detci_elv2str = 1
      isym = 1
      if (nel.gt.0) then
        do iel=1,nel
          iorb = elv(iel)
          isym = MULX(isym,osym(iorb))
c          write(6,*) iorb, iel, isym, arctab(iorb,iel,isym)
          wght = arctab(iorb,iel,isym)

          if (wght.eq.-1) stop 116

          detci_elv2str = detci_elv2str + wght
        enddo
c       detci_elv2str = detci_elv2str + arctab(1,1,isym)
        detci_elv2str = detci_elv2str + arctab(0,0,isym)
      endif
      
      return
      end




c
c  Compute index string from Boolean occupation
c  vector.
c
      integer function detci_occ2str( norb, nel, nsym, osym,
     $                              arctab, occ )
      implicit none
      integer norb
      integer nel
      integer nsym
      integer osym(norb)
      integer arctab(0:norb,0:nel,nsym)
      integer occ(norb)
      integer isym, iel, iorb, wght
#include "symmdef.fh"
#include "bitops.fh"
#include "symmmul.fh"

      detci_occ2str = 1
      isym = 1
      iel = 0
      do iorb=1,norb
        if (occ(iorb).eq.1) then
          iel = iel + 1
C          if (iel.gt.nel) stop 118
          isym = MULX(isym,osym(iorb))
          wght = arctab(iorb,iel,isym)
          if (wght.eq.-1) stop 117
          detci_occ2str = detci_occ2str + wght
        endif
      enddo
c     detci_occ2str = detci_occ2str + arctab(1,1,isym)
      detci_occ2str = detci_occ2str + arctab(0,0,isym)

      return
      end







c
c  Return Boolean occupation from
c  electron occupation
c
      subroutine detci_elv2occ( norb, nel, elv, occ )
      implicit none
      integer norb
      integer nel
      integer elv(nel)
      integer occ(norb)
      integer i, iel

      iel = 1
      do i=1,norb
        if ((iel.le.nel).and.(i.eq.elv(iel))) then
          occ(i) = 1
          iel = iel + 1
        else
          occ(i) = 0
        endif
      enddo
      return
      end


c
c Return electron occupation vector from
c Boolean occupation vector
c

      subroutine detci_occ2elv( norb, nel, occ, elv )
      implicit none
      integer norb
      integer nel
      integer elv(nel)
      integer occ(norb)
      integer i, iel

      iel = 0
      do i=1,norb
        if (occ(i).eq.1) then
          iel = iel + 1
          elv(iel) = i
        endif
      enddo
      return
      end







c
c  Return the hole complement of an
c  electron occupation string
c
      subroutine detci_elv2holv( norb, nel, elv, holv )
      implicit none
      integer norb
      integer nel
      integer elv(nel)
      integer holv(norb-nel)
      integer i, iel, ihole

      ihole = 0
      iel = 1
      do i=1,norb
        if ((iel.le.nel).and.(i.eq.elv(iel))) then
          iel = iel + 1
        else
          ihole = ihole + 1
          holv(ihole) = i
        endif
      enddo
          
      return
      end






      integer function detci_elvsymm( norb, nel, nsym, osym, elv )
      implicit none
      integer norb
      integer nel
      integer nsym
      integer osym(norb)
      integer elv(nel)
      integer i
#include "symmdef.fh"
#include "bitops.fh"
#include "symmmul.fh"

      detci_elvsymm = 1
      do i=1,nel
        detci_elvsymm = MULX(detci_elvsymm,osym(elv(i)))
      enddo

      return
      end



      integer function detci_occsymm( norb, nel, nsym, osym, occ )
      implicit none
      integer norb
      integer nel
      integer nsym
      integer osym(norb)
      integer occ(norb)
      integer i
#include "symmdef.fh"
#include "bitops.fh"
#include "symmmul.fh"

      detci_occsymm = 1
      do i=1,norb
        if (occ(i).eq.1) detci_occsymm = MULX(detci_occsymm,osym(i))
      enddo

      return
      end




c
c                                
c  Compute binomial coefficient
c                              
c      n!
c   ------ ,     m+k <= n         with minimal overflow 
c   m!(n-m)!                   (so theorem says this is always integer div if (m+k) <= n)
c
c
       integer function detci_binomial(n,m)
       implicit none
       integer n,m
       integer m1,k1
       integer i,ii,jj

       m1 = max(max(m,(n-m)),1)
       k1 = max(min(m,(n-m)),1)
       jj = min(2,k1)
       ii = 1
       do i=n,m1+1,-1
         ii = ii*i
 100     continue
         if ((jj.ne.1).and.(mod(ii,jj).eq.0)) then
           ii = ii/jj
           jj = jj + 1
           if (jj.gt.k1) jj = 1
           goto 100
         endif
       enddo
       detci_binomial = ii
       return
       end







c
c  Return next occupation string with respect to input string
c  by "cascading" electrons into successive orbitals.
c  True lexical ordering, *not* compatible with ordering 
c  using arc-table and symmtery.
c  Recursive and non-recursive versions.
c
c
#ifdef FORTRAN_RECURSION
      logical function detci_getnextelv( norb, nel, elv )
      implicit none
      integer norb
      integer nel
      integer elv(nel)
      integer iel
      logical status

      iel = nel
      status = .true.
      call detci_nextelv_( norb, nel, iel, elv, status )
      detci_getnextelv = status

      return
      end



      subroutine detci_nextelv_( norb, nel, iel, elv, status )
      implicit none
      integer norb
      integer nel
      integer iel
      integer elv(nel)
      logical status

      elv(iel) = elv(iel) + 1
      if (elv(iel).le.(norb-nel+iel)) return
      if (iel.eq.1) then
        status = .false.
        return
      endif
      iel = iel - 1
      call detci_nextelv_( norb, nel, iel, elv, status )
      iel = iel + 1
      elv(iel) = elv(iel-1) + 1
      return
      end

#else

      logical function detci_getnextelv( norb, nel, elv )
      implicit none
      integer norb
      integer nel
      integer elv(nel)
      integer iel, i

      detci_getnextelv = .false.
      if (nel.eq.0) return
      iel = nel
 100  continue
      elv(iel) = elv(iel) + 1
      if (elv(iel).gt.(norb-nel+iel)) then
        iel = iel - 1
        if (iel.eq.0) then
          detci_getnextelv = .false.
          return
        endif
        goto 100
      endif
      do i=iel+1,nel
        elv(i) = elv(i-1) + 1
      enddo
      detci_getnextelv = .true.
      return
      end

#endif

        







c
c  Generates excitation operator table
c  Each entry for excitation, E_ji, and string, I, 
c  is encoded as:
c         ex(1,iex,istr) = jstr
c         ex(2,iex,istr) = jsym
c         ex(3,iex,istr) = ij
c         ex(4,iex,istr) = phase
c         ex(5,iex,istr) = j
c         ex(6,iex,istr) = i
c
c
      subroutine detci_excit( norb, nel, nsym, nstr, nex, osym,
     $                        ijmap, atab, ex )
      implicit none
#include "detciP.fh"
#include "detci.fh"
      integer norb
      integer nel
      integer nsym
      integer nstr
      integer nex
      integer osym(norb)
      integer ijmap(norb,norb)
      integer atab(0:norb,0:nel,nsym)
      integer ex(6,nex,nstr)
c
      integer nhole
      integer icreate
      integer elv(detci_maxelec)
      integer occv(detci_maxorb)
      integer holv(detci_maxorb)
      integer istr, jstr, isym, jsym
      integer iel, ihole, iex, iexsym
      integer phase1, phase2
      integer i, j, n1
c
c
c
      if (nel.eq.0) return        ! infinite loop otherwise
      n1 = norb + 1
      nhole = norb - nel
      do i=1,nel
        elv(i) = i
      enddo
c
c  Loop through all strings
c
 100  continue
      istr = detci_elv2str( norb, nel, nsym, osym, atab, elv )
      if (istr.gt.nstr) then 
         write(6,*) ' bad string number? ', istr, '  ', nstr
         write(6,*) ' norb ', norb
         write(6,*) ' nel  ', nel
         write(6,*) ' nsym ', nsym
         write(6,*) ' osym ', (osym(i), i=1, norb)
         write(6,*) ' elv  ', (elv(i), i=1, nel)
         stop 3344
      end if
      isym = detci_elvsymm( norb, nel, nsym, osym, elv )
      call detci_elv2occ( norb, nel, elv, occv )
      call detci_elv2holv( norb, nel, elv, holv )
      iex = 0
c
c
c  Outer symmetry loop to enforce symmtery
c  ordering on excitations
c      
      do iexsym=1,nsym
        do iel=1,nel
          occv(elv(iel)) = 0
          phase1 = (-1)**(1+mod(iel,2))
          do ihole=0,nhole
            if (ihole.eq.0) then
              icreate = elv(iel) ! E_ii 
            else
              icreate = holv(ihole) ! E_ij
            endif
            occv(icreate) = 1
            jsym = detci_occsymm( norb, nel, nsym, osym, occv )
            if (jsym.eq.iexsym) then
              iex = iex + 1
              if (iex.gt.nex) stop 333
              if (istr.gt.nstr) stop 334
              phase2 = phase1
              do j=1,icreate-1
                if (occv(j).eq.1) phase2 = phase2*(-1)
              enddo
              jstr=detci_occ2str(norb, nel, nsym, osym, atab, occv)
              ex(1,iex,istr) = jstr
              ex(2,iex,istr) = jsym
              ex(3,iex,istr) = ijmap(elv(iel),icreate)
              ex(4,iex,istr) = phase2
              ex(5,iex,istr) = icreate
              ex(6,iex,istr) = elv(iel)
            endif
            occv(icreate) = 0
          enddo
          occv(elv(iel)) = 1
        enddo
      enddo
      if (detci_getnextelv( norb, nel, elv )) goto 100
        
      return
      end





c
c Symmetry-blocked triangular sum (ij) mapping
c
      subroutine detci_ijmap( norb, nsym, osym, nij, ijmap )
      implicit none
      integer norb
      integer nsym
      integer osym(norb)
      integer nij
      integer ijmap(norb,norb)
      integer i, j
c
      nij = 0
      do i=1,norb
        do j=1,i
          nij = nij + 1
          ijmap(i,j) = nij
          ijmap(j,i) = nij
        enddo
      enddo

      return
      end








c
c Opaque routine!....
c Given an ordered string pointer array,"ip", for 
c (nel-1) positions, create electron vector "elv" with
c orbital "k" inserted. Phase is also returned.
c
c This is "k" creation operation on (nel-1) string 
c *except* using indirection
c
c
      subroutine detci_ptr2elv(norb,nel,np,nlab,ip,lab,k,elv,iph)
      implicit none
      integer norb, nel
      integer np
      integer nlab
      integer ip(np)
      integer lab(nlab)
      integer k
      integer elv(nel)
      integer iph
      integer iel,i

      iel = 0
      do i=1,np
        if (lab(ip(i)).lt.k) then
          iel = iel + 1
          elv(iel) = lab(ip(i))
          iph = iph*(-1)
        endif
      enddo
      iel = iel + 1
      elv(iel) = k
      do i=1,np
        if (lab(ip(i)).gt.k) then
          iel = iel + 1
          elv(iel) = lab(ip(i))
        endif
      enddo
      return
      end

      



c
c Column gather on CI matrix (with phase)
c
      subroutine detci_cigather(n,m,nv,g_x,iv,ph,y,x,usex)
      implicit none
#include "global.fh"
      integer n,m,nv
      integer iv(nv),ph(nv)
      double precision x(n,m)
      logical usex
      integer g_x
      double precision y(n,*)
      integer i,j,ii
      double precision q

      if (usex) then
        do i=1,nv
          ii = iv(i)
          q = ph(i)
          do j=1,n
            y(j,i) = q*x(j,ii)
          enddo
        enddo
      else
        do i=1,nv
          ii = iv(i)
          q = ph(i)
          call ga_get(g_x,1,n,ii,ii,y(1,i),1)
          do j=1,n
            y(j,i) = q*y(j,i)
          enddo
        enddo
      endif
      return
      end








c
c Column gather on CI matrix no phase
c
      subroutine detci_cigather0(n,m,nv,g_x,iv,y)
      implicit none
#include "global.fh"
      integer n,m,nv
      integer g_x
      integer iv(nv)
      double precision y(n,*)
      integer i,ii

      do i=1,nv
        ii = iv(i)
        call ga_get(g_x,1,n,ii,ii,y(1,i),1)
      enddo
      return
      end








c
c Dot product
c
      subroutine detci_dotab(istr,n,nstr,f,c,s)
      implicit none
      integer istr, n, nstr
      double precision f(nstr)
      double precision c(nstr,n)
      double precision s(nstr,n)
      integer i
      double precision ddot
      external ddot

      do i=1,n
        s(istr,i) = s(istr,i) + ddot(nstr,f,1,c(1,i),1)
      enddo
      return
      end



c
c  Sparse dot product
c
      subroutine detci_dotabx( istr, n, nstr, nex, f, idx, c, s )
      implicit none
      integer istr, n, nex, nstr
      integer idx(nex)
      double precision f(nex)
      double precision c(n,nstr)
      double precision s(n,nstr)
      double precision xx
      integer i, j


      do j=1,nex
        xx = f(j)
	if (xx .ne. 0.0d0) then
          do i=1,n
            s(i,istr) = s(i,istr) + xx*c(i,idx(j))
          enddo
        endif
      enddo
      return
      end





c
c
c CI vector scatter accumulate
c
      subroutine detci_ciscatter(n,m,nv,x,iv,g_y,y,usey)
      implicit none
#include "global.fh"
      integer n,m,nv
      double precision x(n,nv)
      double precision y(n,m)
      logical usey
      integer g_y
      integer iv(nv)
      integer i,ii,j

      do i=1,nv
        ii = iv(i)
	if (usey) then
	  do j = 1, n
	     y(j,ii) = y(j,ii) + x(j,i)
          enddo
	else
          call ga_acc(g_y,1,n,ii,ii,x(1,i),n,1.d0)
	endif
      enddo
      return
      end







c
c
c Return operator indices for
c
c               + + +
c      |lhs> = a a a ... a a a ... |rhs>
c               i j k     r s t
c     
c
      integer function detci_elvop( nel, rhs, lhs, phase, a, c )
      implicit none
      integer nel
      integer rhs(nel)
      integer lhs(nel)
      integer phase
      integer a(nel)
      integer c(nel)
      integer i,j,nop1,nop2

      do i=1,nel
        a(i) = -1
        c(i) = -1
      enddo

      nop1 = 0
      phase = 1
      do i=1,nel
        do j=1,nel
          if (lhs(j).eq.rhs(i)) goto 111
        enddo
        phase = phase * ((-1)**(i-1-nop1))
        nop1 = nop1 + 1
        a(nop1) = rhs(i)
 111    continue
      enddo
      
      nop2 = 0
      do i=1,nel
        do j=1,nel
          if (rhs(j).eq.lhs(i)) goto 112
        enddo
        phase = phase * ((-1)**(i-1))
        nop2 = nop2 + 1
        c(nop2) = lhs(i)
 112    continue
      enddo
      
      if (nop1.ne.nop2) stop 1177
      detci_elvop = nop1
      return
      end
        



c
c  Print out dominant CI vector elements
c  with corresponding configuration info
c
      subroutine detci_civec_print( norb, nsym, nela, nelb,
     $                              nstra, nstrb, osym,
     $                              map, ataba, atabb,
     $                              g_x, thresh )
      implicit none
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
      integer norb
      integer nsym
      integer nela, nelb
      integer nstra, nstrb
      integer osym(norb)
      integer map(norb)
      integer ataba(0:norb,0:nela,nsym), atabb(0:norb,0:nelb,nsym)
      integer g_x
      double precision thresh
c
      integer elva(detci_maxelec), elvb(detci_maxelec)
      integer iaa(detci_maxelec), ibb(detci_maxelec)
      integer i, isym, jsym, istr, jstr, ici
      double precision xx
c
      if (ga_nodeid() .gt. 0) return
c
c
c Loop over pairs of alpha,beta strings = CI index
c
      if (ga_nodeid().eq.0) write(6,909)
 909  format(/3x,'Index',11x,'Coefficient',10x,'Config. Occupation')
      do i=1,nela
        elva(i) = i
      enddo
 101  continue
      isym = detci_elvsymm( norb, nela, nsym, osym, elva )
      istr = detci_elv2str( norb, nela, nsym, osym, ataba, elva )
      do i=1,nelb
        elvb(i) = i
      enddo
 102  continue
      jsym = detci_elvsymm( norb, nelb, nsym, osym, elvb )
      jstr = detci_elv2str( norb, nelb, nsym, osym, atabb, elvb )
      ici = (istr-1)*nstrb + jstr
      call ga_get(g_x,jstr,jstr,istr,istr,xx,1)
      if (abs(xx).ge.thresh) then
        do i=1,nela
          iaa(i) = map(elva(i))
        enddo
        do i=1,nelb
          ibb(i) = map(elvb(i))
        enddo
        call util_isort(nela,iaa)
        call util_isort(nelb,ibb)
        if (ga_nodeid().eq.0) write(6,901) ici,xx,
     $                                     (iaa(i),i=1,nela),
     $                                     (ibb(i),i=1,nelb)
 901    format(i10,5x,f15.8,8x,20i3)
      endif
c
      if (detci_getnextelv( norb, nelb, elvb )) goto 102
      if (detci_getnextelv( norb, nela, elva )) goto 101
      
      return
      end




      double precision function detci_refenergy( norb, nela, nelb, eps )
      implicit none
#include "detciP.fh"
      integer norb
      integer nela
      integer nelb
      double precision eps(norb)
c
      integer i,j,jj,ii(detci_maxorb)
      double precision xx
c
      detci_refenergy = -0.2d0
      call ifill(norb,0,ii,1)
      do i=1,nela
        xx = 100.d0
        jj=-3333
        do j=1,norb
          if ((ii(j).eq.0).and.(eps(j).lt.xx)) then
            xx = eps(j)
            jj = j
          endif
        enddo
        detci_refenergy = detci_refenergy + xx
        if (jj .eq. -3333) stop 543
        ii(jj) = 1
      enddo
      call ifill(norb,0,ii,1)
      do i=1,nelb
        xx = 100.d0
        jj=-3333
        do j=1,norb
          if ((ii(j).eq.0).and.(eps(j).lt.xx)) then
            xx = eps(j)
            jj = j
          endif
        enddo
        detci_refenergy = detci_refenergy + xx
        if (jj .eq. -3333) stop 542
        ii(jj) = 1
      enddo
      return
      end











      subroutine detci_strlist( norb, nsym, nel, nstr, osym, atab )
      implicit none
#include "detciP.fh"
#include "detci.fh"
      integer norb
      integer nsym
      integer nel
      integer nstr
      integer osym(norb)
      integer atab(0:norb,0:nel,nsym)
c
      integer elv(detci_maxelec)
      integer i, isym, istr
c
c
c
      write(6,909)
 909  format(//3x,'Symm. Index',12x,'Occupation')
      do i=1,nel
        elv(i) = i
      enddo
 101  continue
      isym = detci_elvsymm( norb, nel, nsym, osym, elv )
      istr = detci_elv2str( norb, nel, nsym, osym, atab, elv )
      write(6,988) isym,istr,(elv(i),i=1,nel)
 988  format(i3,i6,5x,10i3)
      if (detci_getnextelv( norb, nel, elv )) goto 101
      
      return
      end







c
c  String index range for some irrep
c
      subroutine detci_symstrr( norb, nsym, nel, nstr, atab,
     $                          ksym, slo, shi ,
     .     mxcnt, kci,domxcnt,oscreen)
      implicit none
#include "detciP.fh"
#include "detci.fh"
      integer norb
      integer nsym
      integer nel
      integer nstr
      integer atab(0:norb,0:nel,nsym)
      integer ksym
      integer slo, shi
      logical domxcnt,oscreen
      double precision kci(*)
      double precision mxcnt
      integer jstr
c
c  slo > shi is possible for no orbital of symmetry ksym
c
      slo = 0
      shi = 0
      if (nel.gt.0) then
c     slo = atab(1,1,ksym) + 1
         slo = atab(0,0,ksym) + 1
         if (ksym.eq.nsym) then
            shi = nstr
         else
c     shi = atab(1,1,ksym+1)
            shi = atab(0,0,ksym+1)
         endif
      endif
      
      if(domxcnt) then
         do jstr=slo,shi
            mxcnt = max(mxcnt,abs(kci(jstr)))
            if (oscreen) kci(jstr) = 0.d0
         enddo
      endif
      
      return
      end










c
c  Generate lowest Aufbau configuration
c  for given symmetry
c
      subroutine detci_guesscfg( norb, nsym, nela, nelb,
     $                           osym, eps, symstate, cfg, ncfg )
      implicit none
#include "errquit.fh"
#include "detci.fh"
#include "detciP.fh"
      integer norb
      integer nsym
      integer nela, nelb
      integer osym(norb)
      double precision eps(norb)
      integer symstate
      integer cfg(nela+nelb,*)
      integer ncfg
c
      integer ix(detci_maxorb),kx(detci_maxorb)
      integer jx(2*detci_maxorb)
      integer i,j,ii,jj,nopen,noorb,isym
      double precision elow, xx
#include "symmdef.fh"
#include "bitops.fh"
#include "symmmul.fh"
c
c  Sort orbitals by energy
c
      ncfg = 1
      call ifill(norb,0,ix,1)
      do i=1,norb
        xx = 1000.d0
        ii = 0
        do j=1,norb
          if ((eps(j).lt.xx).and.(ix(j).eq.0)) then
            xx = eps(j)
            ii = j
          endif
        enddo
        ix(ii) = i
        kx(i) = ii
      enddo
c
c  Doubly occupy lowest nelb orbitals 
c
      call ifill((nela+nelb)*2,0,cfg,1)
      jj = 0
      do i=1,nelb
        cfg(i,1) = kx(i)
        cfg(i,2) = kx(i)
        cfg(nela+i,1) = kx(i)
        cfg(nela+i,2) = kx(i)
      enddo
c
c  If closed shell, but the requested symmetry is
c  not symmetryic one, we need to find an
c  open-shell singlet with alpha and beta in
c  different orbitals
c
      if (nela.eq.nelb.and.symstate.gt.1) then
         ncfg = 2
         do i=nela,1,-1
            do j=nela,norb
            isym = osym(kx(i))
            isym = MULX(isym,osym(kx(j)))
            if (isym.eq.symstate) goto 110
            enddo
         enddo
         call errquit('singlet symmetry state not found in CAS',0,
     &       GEOM_ERR)
 110     cfg(nela+i,1) = kx(j)
         cfg(i,2) = kx(j)
      endif
c
c  Loop through all possible permutations of
c  (nela-nelb) electons in (norb-nelb) orbitals
c  Find lowest energy with appropriate symmetry
c
      if (nela.gt.nelb) then
        nopen = nela - nelb
        noorb = norb - nelb
        elow = 10000.d0
        cfg(nelb+1,1) = -1
        do i=1,nopen
          ix(i) = i
        enddo
 111    continue
        isym = 1
        xx = 0.d0
        do i=1,nopen
          jx(i) = kx(nelb + ix(i))
          isym = MULX(isym,osym(jx(i)))
          xx = xx + eps(jx(i))
        enddo
c
c  If correct symmetry & lower energy cfg copy result
c
        if ((isym.eq.symstate).and.(xx.lt.elow)) then
          elow = xx
          do i=1,nopen
            cfg(nelb+i,1) = jx(i)
          enddo
        endif
        if (detci_getnextelv(noorb,nopen,ix)) goto 111
        if (cfg(nelb+1,1).eq.-1)
     $    call errquit('symmetry state not found in CAS',0,
     &       GEOM_ERR)
      endif
c
c  Sort config by ascending indices
c
      call util_isort(nela,cfg(1,1))
      call util_isort(nelb,cfg(nela+1,1))
      if (ncfg.eq.2) then
         call util_isort(nela,cfg(1,2))
         call util_isort(nelb,cfg(nela+1,2))
      endif
c
      return
      end
c
c  Generate some CI guess vector
c  from a list of configurations and 
c  associated coefficients
c
      subroutine detci_ciguess( norb, nsym, nela, nelb, nstra, nstrb,
     $                          osym, ataba, atabb,
     $                          ngs, cfggs, cgs, g_civec )
      implicit none
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
      integer norb                            ! [input] Orbitals
      integer nsym                            ! [input] Irreps
      integer nela                            ! [input] Alpha electrons
      integer nelb                            ! [input] Beta electrons
      integer nstra                           ! [input] Alpha strings
      integer nstrb                           ! [input] Beta strings
      integer osym(norb)                      ! [input] Orbital irreps     
      integer ataba(0:norb,0:nela,nsym)           ! [input] Alpha arc weights
      integer atabb(0:norb,0:nelb,nsym)           ! [input] Beta arc weights
      integer ngs                             ! [input] Number of configurations
      integer cfggs((nela+nelb),ngs)          ! [input] Configurations
      double precision cgs(ngs)               ! [input] Coefficients
      integer g_civec                         ! [output] GA CI-vector
c
      integer astr, bstr, i
      integer elva(detci_maxelec), elvb(detci_maxelec)
      double precision xx
      double precision dnrm2
      external dnrm2
c
c
      call ga_zero(g_civec)
      do i=1,ngs
        call icopy(nela,cfggs(1,i),1,elva,1)
        call icopy(nelb,cfggs(nela+1,i),1,elvb,1)
        astr = detci_elv2str( norb, nela, nsym, osym, ataba, elva )
        bstr = detci_elv2str( norb, nelb, nsym, osym, atabb, elvb )
        call ga_put(g_civec,bstr,bstr,astr,astr,cgs(i),1)
      enddo
      xx = sqrt(ga_ddot(g_civec,g_civec))
      if (xx.gt.0.0001d0) call ga_dscal(g_civec,(1.d0/xx))

      return
      end





      subroutine detci_etrace( norb, ntij, ijmap,
     $                         h, g, onepdm, twopdm )
      implicit none
      integer norb
      integer ntij
      integer ijmap(norb,norb)
      double precision h(ntij)
      double precision g(ntij,ntij)
      double precision onepdm(norb,norb)
      double precision twopdm(norb,norb,norb,norb)
c
      integer i,j,k,l
      integer ij, kl
      double precision xx, x1, x2

      write(6,'(//,a)') 'Test Density matrices'
      xx = 0.d0
      do i=1,norb
        xx = xx + onepdm(i,i)
      enddo

      x1 = 0.d0
      x2 = 0.d0
      do i=1,norb
        do j=1,norb
          ij = ijmap(i,j)
          x1 = x1 + h(ij)*onepdm(i,j)
          do k=1,norb
            do l=1,norb
              kl = ijmap(k,l)
              x2 = x2 + g(ij,kl)*twopdm(i,j,k,l)
            enddo
          enddo
        enddo
      enddo
      write(6,'(10x,a,16x,f15.6)') '1-pdm trace',xx
      write(6,'(10x,a,10x,f20.10)') 'energy trace',(x1+x2)

      return
      end






c
c  Scale the CI vector with sum of orbital energies - ground-state sum
c
c
c                  C
c                   i
c     C  = ------------------
c      i   E    -  sum_(k) e
c           ref             k
c
c
      subroutine detci_diagscale( norb, nsym, nela, nelb, nstra, nstrb,
     $                            osym, ataba, atabb, eref, eps,
     $                            g_civec )
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
#include "mafdecls.fh"
      integer norb                            ! [input] Orbitals
      integer nsym                            ! [input] Irreps
      integer nela                            ! [input] Alpha electrons
      integer nelb                            ! [input] Beta electrons
      integer nstra                           ! [input] Alpha strings
      integer nstrb                           ! [input] Beta strings
      integer osym(norb)                      ! [input] Orbital irreps     
      integer ataba(0:norb,0:nela,nsym)           ! [input] Alpha arc weights
      integer atabb(0:norb,0:nelb,nsym)           ! [input] Beta arc weights
      double precision eref                   ! [input] Reference ground-state
      double precision eps(norb)              ! [input] Orbital energies
      integer g_civec                         ! [input] CI-vector
c$$$      double precision civec(nstrb,nstra)     ! [input] CI-vector
c
      integer elva(detci_maxelec), elvb(detci_maxelec)
      integer i, isym, jsym, istr, jstr
      integer rlo, rhi, clo, chi, myid
      double precision ea, eb, xx
      double precision denominator
      integer l_c, k_c
c
c$$$      DOUBLE PRECISION XRANDOM
c$$$      DOUBLE PRECISION DRAND
c
      myid = ga_nodeid()
      call ga_distribution( g_civec, myid, rlo, rhi, clo, chi )
      if (.not. ma_push_get(mt_dbl, nstrb, 'detci_dscal: beta',
     $     l_c, k_c)) call errquit('detci_diagscale: ma', nstrb, MA_ERR)
c
c Loop over pairs of alpha,beta strings = CI index
c
      do i=1,nela
        elva(i) = i
      enddo
 101  continue
      isym = detci_elvsymm( norb, nela, nsym, osym, elva )
      istr = detci_elv2str( norb, nela, nsym, osym, ataba, elva )
      if ((istr.ge.clo).and.(istr.le.chi)) then
        call ga_get(g_civec,1,nstrb,istr,istr,dbl_mb(k_c),1)
        ea = 0.d0
        do i=1,nela
          ea = ea + eps(elva(i))
        enddo
        do i=1,nelb
          elvb(i) = i
        enddo
 102    continue
        eb = 0.d0
        do i=1,nelb
          eb = eb + eps(elvb(i))
        enddo
        jsym = detci_elvsymm( norb, nelb, nsym, osym, elvb )
        jstr = detci_elv2str( norb, nelb, nsym, osym, atabb, elvb )
*        call ga_get(g_civec,jstr,jstr,istr,istr,xx,1)
        xx = dbl_mb(k_c+jstr-1)
	denominator = (eref - ea - eb )
        if (abs(denominator) .lt. 0.1d0) denominator = 
     $    sign(0.1d0,denominator)
        xx = xx / denominator
*        call ga_put(g_civec,jstr,jstr,istr,istr,xx,1)
        dbl_mb(k_c+jstr-1) = xx
c
C
C ******* TESTING ONLY ********
C  Deliberately introduce spin 
c  contamination to test spin-adaption
C
C      XRANDOM = ((0.5d0-DRAND(0))*0.001)*(EB+EA)
C      CIVEC(JSTR,ISTR) = CIVEC(JSTR,ISTR)/(EREF - EA - EB + XRANDOM)
C
        if (detci_getnextelv( norb, nelb, elvb )) goto 102
        call ga_put(g_civec,1,nstrb,istr,istr,dbl_mb(k_c),1)
      endif
      if (detci_getnextelv( norb, nela, elva )) goto 101
c
      if (.not. ma_pop_stack(l_c)) call errquit
     $   ('detci:dscal ma corrupt',0, MA_ERR)
      call ga_sync
c
      end








c
c  Project out coefficients of given irrep
c
c
      subroutine detci_symmproject( norb, nsym, nela, nelb,
     $                              nstra, nstrb, osym, 
     $                              ataba, atabb, ksym, oscreen,
     $                              mxcnt, g_civec )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
      integer norb                            ! [input] Orbitals
      integer nsym                            ! [input] Irreps
      integer nela                            ! [input] Alpha electrons
      integer nelb                            ! [input] Beta electrons
      integer nstra                           ! [input] Alpha strings
      integer nstrb                           ! [input] Beta strings
      integer osym(norb)                      ! [input] Orbital irreps     
      integer ataba(0:norb,0:nela,nsym)           ! [input] Alpha arc weights
      integer atabb(0:norb,0:nelb,nsym)           ! [input] Beta arc weights
      integer ksym                            ! [input] Requested irrep
      logical oscreen                         ! [input] Toggle projecting
      double precision mxcnt                  ! [output] Maximum symmetry contamination
      integer g_civec                         ! [input/output] CI-vector
c
      integer istr, jstr, asym, bsym
      integer alo, ahi, blo, bhi
      integer myid, rlo, rhi, cilo, cihi
      integer k_ci, ldc
c
#include "symmdef.fh"
#include "bitops.fh"
#include "symmmul.fh"
c
      mxcnt = 0.d0
      myid = ga_nodeid()
      call ga_distribution( g_civec, myid, rlo, rhi, cilo, cihi )
      do asym=1,nsym
        call detci_symstrr( norb, nsym, nela, nstra, ataba, asym,
     $                      alo, ahi ,
     ,        mxcnt, dbl_mb(k_ci), .false.,oscreen)
        do istr=alo,ahi
          if ((istr.ge.cilo).and.(istr.le.cihi)) then
            call ga_access( g_civec, rlo, rhi, istr, istr, k_ci, ldc )
            do bsym=1,nsym
              if (MULX(asym,bsym).ne.ksym) then
                call detci_symstrr( norb, nsym, nelb, nstrb, atabb,
     $                              bsym, blo, bhi ,
     ,        mxcnt, dbl_mb(k_ci), .true., oscreen)
!                do jstr=blo,bhi
!                  mxcnt = max(mxcnt,abs(dbl_mb(k_ci+jstr-1)))
!                  if (oscreen) dbl_mb(k_ci+jstr-1) = 0.d0
!                enddo
              endif
            enddo
            call ga_release( g_civec, rlo, rhi, istr, istr )
          endif
        enddo
      enddo
      call ga_sync()

      return
      end









       subroutine transpose_nw( n, m, x, y )
       implicit none
       integer n, m
       double precision x(n,m)
       double precision y(m,n)
       integer i, j

       do i=1,n
         do j=1,m
           y(j,i) = x(i,j)
         enddo
       enddo
       return
       end
