      logical function  CCT_UHF (rtdb)
c 
c....  development version UHF CCSD(T)
c      ==========================================
c      CHECK
c      c/a symmetry
c      triangular symmetries (ij/ab in summation)
c      are sorted in  cct_uhf_t2  ......
c      can be (but are not)sorted by the factors in the occ's
c      ==========================================
c
c...    SYMMETRY of integrals requires careful CONSIDERATION
c...    SYMMETRY of T's requires careful CONSIDERATION
c...    All needws careful checking
c
c...   dgemmm description'
c
c...   call dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
c...   C:= alpha*A*B + beta*C
c...   transa,transb : 'N' : asis ; 'T' transposed
c...   m : # rows A,C            : nar
c...   n : # columns B,C         : nbc
c...   k : # columns A, rows B   : nac (nlink)
c...   (refers to matrix after possible transpose)
c...   ld : leading dimensions
c
*
* $Id$
*
      implicit none
#include "util.fh"
c
c
c...   symbolic matrix elements are all mapped onto A3A
c     (ij)t(ij,ce)*t(ij,af)  *  (bk) <ba||ek><bc||fk>
c        T(ij)(ceaf)  (Y)           I(bk)(aecf) (X)
c...   ...  integrals and T';s are treated alike   ...
c...   properties :   RESULT  SUM1 SUM2
c                     a e c f  b k  i j 
c     *note*
c     b,k are running indices for first contraction type(*,1)
c     a,e are on same and c,f are on same "integral" for ,1
c     i,j are running indices for second contraction type(*,2)
c     a,f are on same and c,e are on same "integral" for ,2
c     --------------------------------------------------------
c     e,f are complete ranges / a,c get batched
c     e,f (2,4)  are the  unweighted
c     --------------------------------------------------------
c
c      sij=sce=saf / sbk=sae=scf
c      do sbk=1,nirr
c       do sij=1,nirr
c        do sa=1,nirr
c        for "symmetric" batch of c,a ; all e,f in core
c        within summation bk/ij may be batched (saves space)
c
c      storage
c      bk,ea, so for one a all e's and then bk's
c      The result is (ae,cf) / (af,ce)
c      in each column/row the e/f are the fastest running indices
c         
c
      integer rtdb
c
      double precision r_a2a,r_a2b,r_b2b,r_a3a,r_a3b,r_b3b
      double precision r_a4a,r_a4b,r_b4b
      double precision r_a1a,r_a1b,r_b1b,flop
      double precision cct_uhf_T2
c
      integer ncase_a2a,ncase_a2b,ncase_b2b
      integer ncase_a3a,ncase_a3b,ncase_b3b
      integer ncase_a4a,ncase_a4b,ncase_b4b
      integer ncase_11a,ncase_11b,ncase_14a,ncase_14b
      parameter (ncase_a2a=5,ncase_a2b=6,ncase_b2b=5)
      parameter (ncase_a3a=5,ncase_a3b=6,ncase_b3b=5)
      parameter (ncase_a4a=4,ncase_a4b=6,ncase_b4b=4)
      parameter (ncase_11a=4,ncase_11b=4,ncase_14a=6,ncase_14b=6)
c
      integer cases_a2a(9,ncase_a2a),occ_a2a(11)
      character*20 a2a,type_a2a(2,2)
      integer cases_a2b(9,ncase_a2b),occ_a2b(11)
      character*20 a2b,type_a2b(2,2)
      integer cases_b2b(9,ncase_b2b),occ_b2b(11)
      character*20 b2b,type_b2b(2,2)
c
      integer cases_a3a(9,ncase_a3a),occ_a3a(11)
      character*20 a3a,type_a3a(2,2)
      integer cases_a3b(9,ncase_a3b),occ_a3b(11)
      character*20 a3b,type_a3b(2,2)
      integer cases_b3b(9,ncase_b3b),occ_b3b(11)
      character*20 b3b,type_b3b(2,2)
c
      integer cases_a4a(9,ncase_a4a),occ_a4a(11)
      character*20 a4a,type_a4a(2,2)
      integer cases_a4b(9,ncase_a4b),occ_a4b(11)
      character*20 a4b,type_a4b(2,2)
      integer cases_b4b(9,ncase_b4b),occ_b4b(11)
      character*20 b4b,type_b4b(2,2)
c
      integer cases_11a(9,ncase_11a),occ_11a(11)
      character*20 a11,type_11a(2,2)
      integer cases_11b(9,ncase_11b),occ_11b(11)
      character*20 b11,type_11b(2,2)
c
      integer cases_14a(9,ncase_14a),occ_14a(11)
      character*20 a14,type_14a(2,2)
      integer cases_14b(9,ncase_14b),occ_14b(11)
      character*20 b14,type_14b(2,2)
c
c======================================================================
      data a2a/'A2A'/
c                     k e i f  a b  j c  fac
      data cases_a2a/ 1,1,1,1, 1,1, 0,0,  1,
     2                2,2,2,2, 2,2, 0,0,  1,
     3                1,2,1,2, 1,2, 2,1,  2,
     4                2,1,2,1, 1,2, 1,2,  2,
     5                2,1,1,2, 1,2, 0,0,  4/
      data occ_a2a  / 1,2,1,2, 1,2, 1,2, -2,      1,3/
      data type_a2a /'<ab||ek>(abek)','<ab||ek>(abek)',
     1               'T(ij,ab)(aibj)','T(ij,ab)(aibj)'/
c======================================================================
      data a2b/'A2B'/
c                     i n k e  a b  j c
      data cases_a2b/ 1,1,1,1, 1,1, 0,0,  1,
     2                2,2,2,2, 2,2, 0,0,  1,
     3                1,2,2,1, 1,2, 0,0,  2,
     4                2,1,1,2, 1,2, 0,0,  2,
     5                2,1,2,1, 1,2, 1,2,  2,
     6                1,2,1,2, 1,2, 2,1,  2/
      data occ_a2b  / 1,1,1,2, 2,2, 1,2, +2,      1,3/
      data type_a2b /'T(ij,ab)(abij)','<ba||ek>(baek)',
     1               'T(ij,ab)(aibj)','<ij||an>(iajn)'/
c======================================================================
      data b2b/'B2B'/
c                     k m i n  a b  j c
      data cases_b2b/ 1,1,1,1, 1,1, 0,0,  1,
     2                2,2,2,2, 2,2, 0,0,  1,
     3                1,2,1,2, 1,2, 2,1,  2,
     4                2,1,2,1, 1,2, 1,2,  2,
     5                2,1,1,2, 1,2, 0,0,  4/
      data occ_b2b  / 1,1,1,1, 2,2, 1,2, -2,      1,3/
      data type_b2b /'T(ij,ab)(abij)','T(ij,ab)(abij)',
     1               '<ij||an>(iajn)','<ij||an>(iajn)'/
c======================================================================
      data a3a/'A3A'/
c...  cases contains the spin cases and the factor to multiply by
c...  occ contains the occupations (1=occ,2=virt) and the divide factor
c...  and (to be save) the indices to be weighted in the intermediate)
c...  type contains types for intgrals/t's in first and second mult
c                     1 2 3 4  5 6  7 8
c                     a e c f  b k  i j, fac 
      data cases_a3a/ 1,1,1,1, 0,0, 1,1,  1,
     2                2,2,2,2, 0,0, 2,2,  1,
     3                1,2,1,2, 2,1, 1,2,  2,
     4                2,1,2,1, 1,2, 1,2,  2,
     5                1,1,2,2, 0,0, 1,2, -4/
c..    above since we seem to get abab while specifying abba
c    5                1,1,2,2, 0,0, 1,2,  4/
      data occ_a3a  / 2,2,2,2, 2,1, 1,1,  2,   1,3/
      data type_a3a /'<ba||ek>(bkea)','<ba||ek>(bkea)',
     1               'T(ij,ce)(ijec)','T(ij,ce)(ijec)'/
c======================================================================
      data a3b/'A3B'/
c                     c n a e  b k  i j
      data cases_a3b/ 1,1,1,1, 0,0, 1,1,  1,
     2                2,2,2,2, 0,0, 2,2,  1,
     3                1,1,2,2, 0,0, 1,2,  2,
     4                2,2,1,1, 0,0, 1,2,  2,
     5                1,2,1,2, 2,1, 1,2,  2,
     6                2,1,2,1, 1,2, 1,2,  2/
      data occ_a3b  / 2,1,2,2, 2,1, 1,1, -2,   1,3/
      data type_a3b /'T(ij,ab)(aibj)','<ba||ek>(bkea)',
     1               'T(ij,ab)(aibj)','<ij||an>(ijan)'/
c======================================================================
      data b3b/'B3B'/
c                     a m c n  b k  i j
      data cases_b3b/ 1,1,1,1, 0,0, 1,1,  1,
     2                2,2,2,2, 1,1, 2,2,  1,
     3                1,2,1,2, 2,1, 1,2,  2,
     4                2,1,2,1, 1,2, 1,2,  2,
     5                2,2,1,1, 0,0, 1,2,  4/
      data occ_b3b  / 2,1,2,1, 2,1, 1,1, -2,   1,3/
      data type_b3b /'T(ij,ce)(ciej)','T(ij,ce)(ciej)',
     1               '<cm||ij>(ijcm)','<cm||ij>(ijcm)'/
c======================================================================
      data a4a/'A4A'/
c                     1 2 3 4  5 6  7 8
c                     j e b f  a k  c i, fac 
      data cases_a4a/ 1,1,1,1, 0,0, 0,0,  1,
     2                2,2,2,2, 0,0, 0,0,  1,
     3                1,1,2,2, 1,2, 0,0,  2,
     4                2,2,1,1, 2,1, 0,0,  2/
      data occ_a4a  / 1,2,2,2, 2,1, 2,1,  1,   1,3/
      data type_a4a /'T(ij,ce)(ciej)','<ba||ek>(bkea)',
     1               'T(ij,ce)(ciej)','<ba||ek>(bkea)'/
c======================================================================
      data a4b/'A4B'/
c                     b e j n  a k  c i
      data cases_a4b/ 1,1,1,1, 0,0, 0,0,  1,
     2                2,2,2,2, 0,0, 0,0,  1,
     3                2,1,1,2, 1,2, 0,0,  1,
     4                1,2,2,1, 2,1, 0,0,  1,
     5                1,1,2,2, 0,0, 2,1,  1,
     6                2,2,1,1, 0,0, 1,2,  1/
      data occ_a4b  / 2,2,1,1, 2,1, 2,1, -1,   1,3/
      data type_a4b /'<ba||ek>(bkea)','<ba||ek>(bkea)',
     1               'T(ij,ab)(aibj)','T(ij,ab)(aibj)'/
c======================================================================
      data b4b/'B4B'/
c                     b m j n  a k  c i
      data cases_b4b/ 1,1,1,1, 0,0, 0,0,  1,
     2                2,2,2,2, 0,0, 0,0,  1,
     3                2,1,1,2, 1,2, 0,0,  2,
     4                1,2,2,1, 2,1, 0,0,  2/
      data occ_b4b  / 2,1,1,1, 2,1, 2,1,  1,   1,3/
      data type_b4b /'T(ij,ab)(aibj)','<ba||ek>(bkea)',
     1               'T(ij,ab)(aibj)','<ba||ek>(bkea)'/
c======================================================================
c
c     T1 diagrams (cf the T1 comments) "UNCHECKED"
c
c======================================================================
      data a11/'11A'/
c                     1 2 3 4  5 6  7 8
c                     i j f k  a b  c 0, fac 
      data cases_11a/ 1,1,1,1, 1,1, 1,0,  1,
     2                2,2,2,2, 2,2, 2,0,  1,
     3                1,2,2,1, 1,2, 1,0,  4,
     4                2,1,1,2, 1,2, 2,0,  4/
      data occ_11a  / 1,1,2,1, 2,2, 2,0,  4,   0,3/
      data type_11a /'<ij||ab>(abij)','<ij||ab>(abij)',
     1               'T(k,c)(ck)','T(ij,ce)(ceij)'/
c======================================================================
      data b11/'11B'/
c                     1 2 3 4  5 6  7 8
c                     i j n k  a b  c 0, fac 
      data cases_11b/ 1,1,1,1, 1,1, 1,0,  1,
     2                2,2,2,2, 2,2, 2,0,  1,
     3                1,2,2,1, 1,2, 1,0,  4,
     4                2,1,1,2, 1,2, 2,0,  4/
      data occ_11b  / 1,1,1,1, 2,2, 2,0,  4,   0,3/
      data type_11b /'T(ij,ce)(ceij)','<ij||ab>(abij)',
     1               'T(k,c)(ck)','<ij||ab>(abij)'/
c======================================================================
      data a14/'14A'/
c                     1 2 3 4  5 6  7 8
c                     b i f k  a j  c 0, fac 
      data cases_14a/ 1,1,1,1, 0,0, 1,0,  1,
     2                2,2,2,2, 0,0, 2,0,  1,
     3                2,2,1,1, 0,0, 1,0,  1,
     4                1,1,2,2, 0,0, 2,0,  1,
     5                2,1,2,1, 1,2, 1,0,  1,
     6                1,2,1,2, 2,1, 2,0,  1/
      data occ_14a  / 2,1,2,1, 2,1, 2,0, -1,   0,3/
      data type_14a /'T(ij,ce)(ciej)','<ij||ab>(aibj)',
     1               'T(k,c)(ck)','<ij||ab>(aibj)'/
c======================================================================
      data b14/'14B'/
c                     1 2 3 4  5 6  7 8
c                     b i n k  a j  c 0, fac 
      data cases_14b/ 1,1,1,1, 0,0, 1,0,  1,
     2                2,2,2,2, 0,0, 2,0,  1,
     3                2,2,1,1, 0,0, 1,0,  1,
     4                1,1,2,2, 0,0, 2,0,  1,
     5                2,1,2,1, 1,2, 1,0,  1,
     6                1,2,1,2, 2,1, 2,0,  1/
      data occ_14b  / 2,1,1,1, 2,1, 2,0,  1,   0,3/
      data type_14b /'<ij||ab>(aibj)','<an||kj>(ankj)',
     1               'T(k,c)(ck)','T(ij,ce)(ceij)'/
c======================================================================
c
c     r_a2a = cct_uhf_T2(a2a,occ_a2a,cases_a2a,ncase_a2a,type_a2a,
c    1                   r_a1a)
c     r_a2b = cct_uhf_T2(a2b,occ_a2b,cases_a2b,ncase_a2b,type_a2b,
c    1                   r_a1b)
c     r_b2b = cct_uhf_T2(b2b,occ_b2b,cases_b2b,ncase_b2b,type_b2b,
c    1                   r_b1b)
c     print *,' a2a ',r_a2a
c     print *,' a2b ',r_a2b
c     print *,' b2b ',r_b2b
c     print *,' total term2 ',r_a2a+r_a2b+r_b2b
c
      print *,' diagram ',A3A,ncase_a3a
      r_a3a = cct_uhf_T2(a3a,occ_a3a,cases_a3a,ncase_a3a,type_a3a,
     1                   flop)
c     r_a3b = cct_uhf_T2(a3b,occ_a3b,cases_a3b,ncase_a3b,type_a3b,
c    1                   flop)
c     r_b3b = cct_uhf_T2(b3b,occ_b3b,cases_b3b,ncase_b3b,type_b3b,
c    1                   flop)
      print *,' a3a ',r_a3a
c     print *,' a3b ',r_a3b
c     print *,' b3b ',r_b3b
c     print *,' total term3 ',r_a3a+r_a3b+r_b3b
         stop
c
      r_a4a = cct_uhf_T2(a4a,occ_a4a,cases_a4a,ncase_a4a,type_a4a,
     1                   flop)
      r_a4b = cct_uhf_T2(a4b,occ_a4b,cases_a4b,ncase_a4b,type_a4b,
     1                   flop)
      r_b4b = cct_uhf_T2(b4b,occ_b4b,cases_b4b,ncase_b4b,type_b4b,
     1                   flop)
      print *,' a4a ',r_a4a
      print *,' a4b ',r_a4b
      print *,' b4b ',r_b4b
      print *,' total term4 ',r_a4a+r_a4b+r_b4b
c
      CCT_UHF = .true.
c
      return
      end
      subroutine cct_uhf_address(nirr_in,mult_in,iro_in,
     1                           occa,occb,norb)
c
c...  set up cct_uhf adddressing
c     
c...  all alpha's first / beta's second in numbering
c
c...   for t2 vectors, t1 vectors ...
c...  normal ordering is ij ; i.ge.j, so goes for symmetry as well
c...  e.g.
c...  total symmetry 1  orbsym : (1,1) (2,2) (3,3) (4,4)
c...  total symmetry 2  orbsym : (1,2) (3,4)
c...  total symmetry 3  orbsym : (1,3) (2,4)
c...  total symmetry 4  orbsym : (1,4) (2,1)
c...   and (4,1) = -(1,4) transposed
c
c     universal order
c     do j=1,n
c      do i=1,n
c       a(i,j)
c     so if i< j => + (apart from it's own sign)
c
      implicit none
c
      integer nirr_in,mult_in(8,8),iro_in(*),occa,occb,norb
c
c...  cct_table_uhf
c
#include "cct_table_uhf.fh"
c
      integer sij,si,sj,i,spin,k,spinx,spiny
      integer nbk,nij,nab
c
      nirr = nirr_in
c
      do si=1,nirr
       do sj=1,nirr
        mult(si,sj) = mult_in(si,sj)
       end do
      end do
c
      orb = norb
      virt(1) = orb-occa
      virt(2) = orb-occb
      occ(1) = occa
      occ(2) = occb
c
      k = 0
      do spin=1,2
       do i=1,occ(spin)
        k = k + 1
        iocc(k) = 1
        ispin(k) = spin
       end do
       do i=1,virt(spin)
        k = k + 1
        iocc(k) = 2
        ispin(k) = spin
       end do
      end do
c
      do i=1,norb
       iro(i) = iro_in(i)
       iro(i+norb) = iro_in(i+norb)
       do spin=1,2
        k = i+(spin-1)*norb
        if (i.le.occ(spin)) then
          nocc(iro(k),spin) = nocc(iro(k),spin) + 1
        else
          nvirt(iro(k),spin) = nvirt(iro(k),spin) + 1
        end if
       end do
      end do
c
c...  some addressing
c     spin can be deduced from above norb
c
      i = 1
      do spin=1,2
       do si = 1,nirr
        socc(si,spin) = i
        i = i + nocc(si,spin)
       end do
       socc(nirr+1,spin) = i
       do si = 1,nirr
        svirt(si,spin) = i
        i = i + nvirt(si,spin)
       end do
       svirt(nirr+1,spin) = i
      end do
c
c...  len_bk
c...  spinx => virtual, spiny => occupied (if applicable)
c
      do spinx = 1,2
       do spiny = 1,2
        do sij=1,nirr
         nbk = 0
         nij = 0
         nab = 0
         do si=1,nirr
          sj = mult(sij,si)
          nbk = nbk + nocc(si,spiny)*nvirt(sj,spinx)
          if (si.eq.sj.and.spinx.eq.spiny) then
           nij = nij + nocc(si,spinx)*(nocc(si,spinx)+1)/2
           nab = nab + nvirt(si,spinx)*(nvirt(si,spinx)+1)/2
          else if (si.le.sj) then
           nij = nij + nocc(si,spinx)*nocc(sj,spiny)
           nab = nab + nvirt(si,spinx)*nvirt(sj,spiny)
          end if
         end do
         lenbk(sij,spinx,spiny) = nbk
         lenij(sij,spinx,spiny) = nij
         lenab(sij,spinx,spiny) = nab
        end do
       end do
      end do
c
c...  calculate various constants
c
      do i = 1, 2*orb
         fac(i) = 1.0d0
      end do
c
      return
      end
      
      double precision function cct_uhf_T2
     1       (diagram,occupation,cases,ncase,type,
     2        value_X1X)
c
c...  CCT T2 diagram evaluator; modelled for A3A
c
c...  if diagram is A2A,A2B or B2B we can generate the 
c...  corresponding A1A,A1B and B1B simultaneously => value_X1X
c
      implicit none
#include "errquit.fh"
c
c     (ij)t(ij,ce)*t(ij,af)  *  (bk) <ba||ek><bc||fk>
c        T(ij)(ceaf)  (Y)           I(bk)(aecf) (X)
c
c...  scheme 
c     get as many integrals in (Global) core
c     calculate X ; then corresponding T and be done
c...  deciding factor is the aecf (n**4) storage.
c...  storage :
c     sbk : store  <ba||ek> as bk,ae (bk along columns)
c          do sb=1,nirr
c           sk = mult(sbk,sb)
c           if (sk.lt.sb) skip
c          same scheme for ae (along rows)
c         
c
c...   symbolic matrix element (for A3A)
c
c...   output : accumulated value
c...   input : diagram : text string to identify diagram
c...           occuption,cases,ncase,type:diagram definition (see below)
c
c     integer ncase,case
c     parameter (ncase=5)
c     integer cases(9,ncase),occupation(9)
c     integer c,e,a,f,i,j,b,k
c...    (occupations of these)
c     character*20 type(2,2)
c...  cases contains the spin cases and the factor to multiply by
c...  occ contains the occupations (1=occ,2=virt) and the divide factor
c...  type contains types for intgrals/t's in first and second mult
c     (ij)t(ij,ce)*t(ij,af)  *  (bk) <ba||ek><bc||fk>
c                 1 2 3 4  5 6  7 8
c                 a e c f  b k  i j  fac 
c     data cases/ 1,1,1,1, 0,0, 1,1,  1,
c    1            2,2,2,2, 0,0, 2,2,  1,
c    2            1,2,1,2, 2,1, 1,2,  2,
c    3            2,1,2,1, 1,2, 1,2,  2,
c    4            1,1,2,2, 0,0, 1,2,  4/
c     data occupation  / 2,2,2,2, 2,1, 1,1,  2/
c     data type /'<ba||ek>(bkea)','<ba||ek>(bkea)',
c    1           'T(ij,ce)(ijec)','T(ij,ce)(ijec)'/
c
c     *note*
c     b,k are running indices for first contraction type(*,1)
c     a,e are on same and c,f are on same "integral" for ,1
c     i,j are running indices for second contraction type(*,2)
c     a,f are on same and c,e are on same "integral" for ,2
c     --------------------------------------------------------
c     e,f are complete ranges / a,c get batched
c     e,f (2,4)  are the  unweighted
c     --------------------------------------------------------
c
c
      integer occupation(9),ncase,cases(9,ncase)
      double precision value_X1X
      character*(*) diagram,type(2,2)
c
#include "mafdecls.fh"
#include "cct_table_uhf.fh"
c
      logical ga_create,ga_destroy
      integer avail_ga,ga_memory_avail,ga_memory,ga_nnodes
      integer g_RESULT,g_TESULT
c
      integer spinc,spine,spina,spinf,spini,spinj,spinb,spink,spin
      integer a,e,c,f,b,k,i,j
      integer eb,fb,bb,kb,ab,ae,cb,ce,ib,jb
      integer na,nc,ne,nf
      integer nae,ncf,naf,nce
      integer sbk,sij,sa,se,sf,sc
      integer case,nbatch,nbcha,nbchc,size,size_1
      logical symmetry,a2a,a1a
      double precision factor,cct_sum,sum,cct_sum_a1a
      double precision deb1,deb2,deb3,deb4
c
      sum = 0.0d0
      a2a = diagram.eq.'A2A'.or.diagram.eq.'A2B'.or.diagram.eq.'B2B'
c..    (could deduce this from fact that both c and a are occupied)
      if (a2a) value_X1X = 0.0d0
c
c...   get available global memory
c...   be sure to have 2 local ne arrays per processor in final stage
c
       avail_ga = ga_memory_avail()/8
       call ga_igop(13,avail_ga,1,'min')
c...    remember that there is overhead (in space) as well
       ga_memory = ga_nnodes()*avail_ga*0.9
c     
c                 1 2 3 4  5 6  7 8
c                 a e c f  b k  i j  fac 
       print *,' diagram ',diagram
      a = occupation(1)
      e = occupation(2)
      c = occupation(3)
      f = occupation(4)
      b = occupation(5)
      k = occupation(6)
      i = occupation(7)
      j = occupation(8)
c
      if (a2a.and.(a.eq.2.or.c.eq.2)) call errquit('A2A error',0,
     &       UNKNOWN_ERR)
c
      do  case=1,ncase
c
       spina = cases(1,case)
       spine = cases(2,case)
       spinc = cases(3,case)
       spinf = cases(4,case)
       spinb = cases(5,case)
       spink = cases(6,case)
       spini = cases(7,case)
       spinj = cases(8,case)
c
       factor = cases(9,case)/(occupation(9)*1.0D0)
c
c...   handle triangular factors 
c
       if (i.eq.j.and.spini.eq.spinj) factor = factor*2
       if (b.eq.k.and.spinb.eq.spink) factor = factor*2
c
       eb = s_occ_virt(1,spine,e)
       fb = s_occ_virt(1,spinf,f)
c
       do sbk=1,nirr
        do sij=1,nirr
         do sa=1,nirr
          sf = mult(sij,sa)
          se = mult(sbk,sa)
          sc = mult(sbk,sf)
          if (sc.lt.sa.or.sf.lt.se) go to 10
c
          ne = n_occ_virt(se,spine,e)
          nf = n_occ_virt(sf,spinf,f)
c
          a1a = a2a.and.spina.eq.spinc.and.se.eq.sf
          size_1 = 0
          if (a1a) size_1 = 2*ne**2
c
          symmetry = (se.eq.sf).and.spin.eq.spin 
          size = ne*nf*n_occ_virt(sa,spina,a)*n_occ_virt(sc,spinc,c) 
c         if (.not.symmetry) size = size * 2 ?
          nbatch = max(size/(ga_memory-size_1),1)
c...   'squares' for now  CHECK !!!
          nbcha = dsqrt((nbatch-1)*1.0d0) + 1.00001d0 
          nbchc = (nbatch-1)/nbcha + 1 
          na = (n_occ_virt(sa,spina,a)-1)/nbcha+1 
          nc = (n_occ_virt(sc,spinc,c)-1)/nbchc+1 
c
          do ab=s_occ_virt(sa,spina,a),s_occ_virt(sa+1,spina,a)-1,na
           ae = min(ab+na,s_occ_virt(sa+1,spina,a))-1
           na = ae-ab+1
           do cb=s_occ_virt(sc,spinc,c),s_occ_virt(sc+1,spinc,c)-1,nc
            ce = min(cb+nc,s_occ_virt(sc+1,spinc,c))-1
            nc = ce-cb+1
c
            nae = na*ne
            ncf = nc*nf
            naf = na*nf
            nce = nc*ne
c
            if (.not.ga_create(MT_DBL,nae,ncf,'Result',0,0,g_RESULT))
     1          call errquit('no global space for RESULT',0, GA_ERR)
            call ga_zero(g_RESULT)
c
c...    loop over bk spin if sbk = 0 
c
            spin = spinb
            if (spin.eq.0.and.spink.ne.0) 
     1          call errquit('whats this spin',0, UNKNOWN_ERR)
            spinb = max(spinb,1)
            spink = max(spink,1)
5           bb = s_occ_virt(sbk,spinb,b)
            kb = s_occ_virt(sbk,spink,k)
c...   * I(aa+bb)(aaaa)
            call mult_cct(bb,kb,eb,fb,
     1                    type(1,1),ab,ae,
     1                    type(2,1),cb,ce,
     1                    g_RESULT)
            if (spin.eq.0) then
             spin = 1
             spinb = 2
             spink = 2
             go to 5
            end if
c
            call cct_scale_ac(g_RESULT,ab,ne,cb,nf)
c
c...    check if the I and the T are actually identical (A4A)
c
            if (ab.eq.cb.and.ae.eq.ce.and.bb.eq.ib.and.kb.eq.jb.and.
     1         type(1,1).eq.type(1,2).and.type(2,1).eq.type(2,2)) then
             sum = sum + cct_sum(g_RESULT,g_RESULT,na,nc,ne,nf) * factor
             go to  9
            end if
c
c      * T(aa)(aaaa)
            if (.not.ga_create(MT_DBL,naf,nce,'Tesult',0,0,g_TESULT))
     1          call errquit('no global space for TESULt',0, GA_ERR)
            call ga_zero(g_TESULT)
c
            spin = spini
            if (spin.eq.0.and.spinj.ne.0) 
     1          call errquit('whats this spini',0, UNKNOWN_ERR)
            spini = max(spini,1)
            spinj = max(spinj,1)
6           ib = s_occ_virt(sbk,spini,i)
            jb = s_occ_virt(sbk,spinj,j)
            call mult_cct(ib,jb,fb,eb,
     1                    type(1,2),ab,ae,
     1                    type(2,2),cb,ce,
     1                     g_TESULT)
            if (spin.eq.0) then
             spin = 1
             spini = 2
             spinj = 2
             go to 6 
            end if
c
            call cct_scale_ac(g_TESULT,ab,ne,cb,nf)
c           
c           sum = sum + cct_sum(g_RESULT,g_TESULT,na,nc,ne,nf) * factor
            deb1 = cct_sum(g_RESULT,g_TESULT,na,nc,ne,nf) * factor
            print *,' case ',case,' intermediate ',deb1
            sum = sum + deb1
c
            if (a1a.and.ab.eq.cb) then
             value_X1X = value_X1X +  
     1                   cct_sum_a1a(g_RESULT,g_TESULT,na,ne) * factor
            end if
 
c
c           if (cb.eq.ab) then
c             sum = sum + cct_sum(g_RESULT,g_TESULT,na,nc,ne,nf) * factor
c            else
c             call ga_zero(g_TESULT)
c             spin = spini
c             spini = max(spini,1)
c             spinj = max(spinj,1)
c7            ib = s_occ_virt(sbk,spini,i)
c             jb = s_occ_virt(sbk,spinj,j)
c             call mult_cct(ib,jb,eb,fb,
c     1                     type(1,2),cb,ce,
c     1                     type(2,2),ab,ae,
c     1                      g_TESULT)
c             if (spin.eq.0) then
c              spin = 1
c              spini = 2
c              spinj = 2
c              go to 7 
c             end if
c
c            call cct_scale_ac(g_TESULT,ab,ne,cb,nf)
c
c             sum = sum + cct_sum(g_RESULT,g_TESULT,na,nc,ne,nf) * factor
c            end if
c
            if (.not.ga_destroy(g_TESULT)) 
     1          call errquit('unable to destroy TESULT',0, GA_ERR)
            if (.not.ga_destroy(g_RESULT)) 
     1          call errquit('unable to destroy RESULT',0, GA_ERR)
c
9          end do
          end do
10       end do
        end do
       end do
      end do
c
      cct_uhf_T2 = sum
c
      return
      end
      subroutine mult_cct(bb,kb,eb,fb,
     1                    typeA,ab,ae, typeB,cb,ce,
     1                    g_RESULT)
c
      implicit none
#include "errquit.fh"
c
c...   general (:)) routine to do t*t or int*int (asymmetric)
c...   notation taken from A3A Integral*Integral
c...                       (bk)<ba||ek><bc||fk>
c...                        integrals as (bk,ea) 
c
c...   the first bb and bk are used to denote the full range
c...   of that spin/space/occupation "symmetry" => bbb,bbk
c...   the first be and bf are used to denote full range of
c...   that  spin/occupation "symmetry"  => eb,fb
c
c...   IS SYMMETRY FULLY EXPLOITED ?
c
#include "mafdecls.fh"
#include "cct_table_uhf.fh"
c
      integer g_RESULT
c
      integer bb,kb,eb,fb,ab,ae,cb,ce
      character*(*) typeA,typeB
c
      logical ga_create,ga_destroy
      integer ga_memory_avail,ga_nnodes,avail_ga
      integer g_MATA,g_MATB,size
      integer se,sf,sbk,ee,fe,la,lc,lea,lcf,lbk,llbk,len
      integer batch,nbatch,bkb,bke
      integer spina,spine,spinc,spinf,spinb,spink
      integer occe,occf,occb,occk
      logical symmetric
c
c...   get available global memory
c
       avail_ga = ga_memory_avail()/8
       call ga_igop(13,avail_ga,1,'min')
c...    remember that there is overhead (in space) as well
        avail_ga = ga_nnodes()*avail_ga*0.9
c
c...  figure out spins
c
      spina = ispin(ab)
      spinc = ispin(cb)
      spine = ispin(eb)
      spinf = ispin(fb)
      spinb = ispin(bb)
      spink = ispin(kb)
c...  figure occupations (1=occ,2-virt)
      occe = iocc(eb)
      occf = iocc(fb)
      occb = iocc(bb)
      occk = iocc(kb)
c
      symmetric = (ab.eq.cb.and.ae.eq.ce.and.typeA.eq.typeB.and.
     1             spine.eq.spinf)
      size = avail_ga
      if (.not.symmetric) size = size/2
c
c...  can use symmetry to save a get; we have to use sqrt(fac) then
c
      if (iro(ae).ne.iro(ab)) call errquit('sym in mult_cct',0,
     &       UNKNOWN_ERR)
      if (iro(ce).ne.iro(cb)) call errquit('sym in mult_cct',0,
     &       UNKNOWN_ERR)
c
      sbk = mult(iro(bb),iro(kb))
      se = mult(iro(ab),sbk)
      sf = mult(iro(cb),sbk)
      eb = s_occ_virt(se,spine,occe)
      ee = s_occ_virt(se+1,spine,occe)-1
      fb = s_occ_virt(sf,spinf,occf)
      fe = s_occ_virt(sf+1,spinf,occf)-1
      la = ae-ab+1
      lc = ce-cb+1
      lea = (ee-eb+1)*la
      lcf = (fe-fb+1)*lc
      lbk = lenxy(sbk,spinb,spink,occb+occk)
      llbk = lbk
      len = lbk*(la+lc)*(ee-eb+1)*(fe-fb+1)
      nbatch = (len-1)/size + 1
      len = (lbk-1)/nbatch + 1
c
c...  get global integral array's
c
      if (.not.ga_create(MT_DBL,len,lea,'matA(a,e)',0,0,g_MATA))
     1   call errquit('unable to crate global arrray matA',0, GA_ERR)
      if (.not.symmetric) then
       if (.not.ga_create(MT_DBL,len,lcf,'matA(c,f)',0,0,g_MATB))
     1     call errquit('unable to create global arrray matB',0, GA_ERR)
      else
       g_MATB = g_MATA
      end if
c
      bke = 0
      do batch=1,nbatch
       bkb = bke + 1
       bke = min(bke+len,llbk)
       lbk = (bke-bkb+1)
       call get_cct(g_MATA,typeA,bkb,bke,spinb,spink,
     1              eb,ee,spine,ab,ae,spina)
       if (.not.symmetric)
     1 call get_cct(g_MATB,typeB,bkb,bke,spinb,spink,
     2              fb,fe,spinf,cb,ce,spinc)
c
c...  note ;, the integrals are weighted for the summation indices
c   
       call ga_dgemm('T','N',lea,lcf,lbk,1.0d0,g_MATA,g_MATB,
     1               1.0d0,g_RESULT)
c
      end do
c
      if (.not.symmetric) then
       if (.not.ga_destroy(g_MATB)) call errquit('cant destroy MATB',0,
     &       GA_ERR)
      end if
      if (.not.ga_destroy(g_MATA)) call errquit('cant destroy MATA',0,
     &       GA_ERR)
c
      return
      end
      subroutine cct_scale_ac(g_R,ab,ne,cb,nf)
c
c...  subroutine to add scaling factors to the intermediates
c...  in the global array g_result
c...  the a/c are always the slower changing indices
c...  ae and cf are pairs
c...  template taken from ga_screen (from 1995) and RJH
c
       implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
       integer g_R,ab,cb,ne,nf
c
       integer ilo, ihi, jlo, jhi, idim, jdim, ld, ind
c
       call ga_sync()
c
       call ga_distribution(g_R, ga_nodeid(), ilo, ihi, jlo, jhi)
       if (ihi.gt.0 .and. jhi.gt.0) then
          idim = ihi - ilo + 1
          jdim = jhi - jlo + 1
          call ga_access(g_R, ilo, ihi, jlo, jhi, ind, ld)
          call cct_do_scale(dbl_mb(ind), ld, ilo, ihi, jlo, jhi,
     1                      ab,ne,cb,nf)
          call ga_release_update(g_R, ilo, ihi, jlo, jhi)
       endif
c
       call ga_sync()
c
       end
      subroutine cct_do_scale(a, ld, ilo, ihi, jlo, jhi,
     1                        ab,ne,cb,nf)
c
c...  apply scaling factors to intermediate in part of 
c...  global array a; Note the the scaling is sqrt(scale)
c
      implicit none
c
#include "cct_table_uhf.fh"
c
      integer ld, ilo, ihi, jlo, jhi, ab,cb,ne,nf,i,j
      double precision a(ilo:ilo+ld-1,jlo:jhi),facj
c
      do j = jlo, jhi
         facj = fac((j-1)/nf+cb)
         do i = ilo, ihi
            a(i,j) = a(i,j)*fac((i-1)/ne+ab)*facj
         end do
      end do
c
      end
      double precision function cct_sum(g_R,g_T,na,nc,ne,nf)
c
      implicit none
#include "errquit.fh"
c
c...  function to sum the partial contributions in global arrays R,T
c...  sum R(ae,cf)*T(af,ce)
c...  e and f are the fastest moving indices (and are "complete")
c...  this uses a bit of ma_memory (routine above should be aware)
c
      integer g_R,g_T,na,nc,ne,nf
      integer a,c,f
      integer ihi_R,ilo_R,jlh_R,ilh_T,jlo_T,jhi_T
      integer l_R,k_R,l_T,k_T
      integer type_cct_sum
      data type_cct_sum/42/
      double precision sum,ddot
c
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      if (.not. ma_push_get(mt_dbl, ne,'R',l_R, k_R))
     1    call errquit('cannot allocate R in cct_sum',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, ne,'T',l_T, k_T))
     1    call errquit('cannot allocate T in cct_sum',0, MA_ERR)

      sum = 0.0d0
      do a=1,na
       do c=1,nc
        do f=ga_nodeid()+1,nf,ga_nnodes()
c...     get all e's
         ilo_R = (a-1)*ne+1
         ihi_R = ilo_R+ne-1
         jlh_R = (c-1)*nf+f
         jlo_T = (c-1)*ne+1
         jhi_T = jlo_T+ne-1
         ilh_T = (a-1)*nf+f
c
         call ga_get(g_T,ilh_T,ilh_T,jlo_T,jhi_T,dbl_mb(k_T),1)
         call ga_get(g_R,ilo_R,ihi_R,jlh_R,jlh_R,dbl_mb(k_R),ne)
         sum = sum + ddot(ne,dbl_mb(k_R),1,dbl_mb(k_T),1)
c
        end do
       end do
      end do
c
      if (.not. ma_pop_stack(l_T)) call errquit('mapop T in cct_sum',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_R)) call errquit('mapop R in cct_sum',0,
     &       MA_ERR)
c
      call ga_dgop(type_cct_sum,sum,1,'+')
c
      cct_sum = sum
c
      return
      end
      double precision function cct_sum_a1a(g_R,g_T,na,ne)
c
      implicit none
#include "errquit.fh"
c
c...  function to sum the partial contributions in global arrays R,T
c...  for "diagonal" A1A case  (so a=c)
c...  sum R(ae,cf)*T(af,ce)
c...  e and f are the fastest moving indices (and are "complete")
c...  this uses a bit of ga_memory which should be availabele
c
      integer g_R,g_T,na,ne
      integer g_RR,g_TT,a,ihi,ilo
      logical ga_create,ga_destroy
      double precision ga_ddot
c
#include "mafdecls.fh"
c
      if (.not.ga_create(MT_DBL,ne,ne,'A1A-aaRR',0,0,g_RR))
     1   call errquit('cannot create global A1A-RR',0, GA_ERR)
      if (.not.ga_create(MT_DBL,ne,ne,'A1A-aaTT',0,0,g_TT))
     1   call errquit('cannot create global A1A-TT',0, GA_ERR)
      call ga_zero(g_RR)
      call ga_zero(g_TT)
c
      do a=1,na
       ilo = (a-1)*ne+1
       ihi = ilo+ne-1
       call ga_add_patch(1.0d0,g_R,ilo,ihi,ilo,ihi,
     1                   1.0d0,g_RR,1,ne,1,ne,
     2                   1.0d0,g_RR,1,ne,1,ne)
       call ga_add_patch(1.0d0,g_T,ilo,ihi,ilo,ihi,
     1                   1.0d0,g_TT,1,ne,1,ne,
     2                   1.0d0,g_TT,1,ne,1,ne)
      end do
c
      cct_sum_a1a = ga_ddot(g_RR,g_TT)
c
      if (.not.ga_destroy(g_RR)) call errquit('cannot destroy RR',0,
     &       GA_ERR)
      if (.not.ga_destroy(g_TT)) call errquit('cannot destroy TT',0,
     &       GA_ERR)
c
      return
      end
      double precision function cct_uhf_T1
     1       (diagram,occupation,cases,ncase,type,
     2        value_X1X)
c
c...  CCT T1 diagram evaluator; modelled after the T2 version
c...  may well disappear and be replaced by ifs+comments
c
      implicit none
c
c     examples:
c
c     (c)T(k,c)*T<in||cb>  *  (aj)<ab||ij><an||kj>  (14B)
c                        or
c     (kc)T(k,c)T(kj,cf)  *  (ab)<ij||ab><fi||ab>  (12A)
c...  the kc contractions are sort of the diagonal of the c
c     they are ther for each case !!
c
c     the 4-4 contraction is always first
c     so the T1's follow ...
c
c     map to CCT T2 diagram evaluator; modelled for A3A
c     (ij)t(ij,ce)*t(ij,af)  *  (bk) <ba||ek><bc||fk>
c        T(ij)(ceaf)  (Y)           I(bk)(aecf) (X)
c
c     (j<=>k,a<=>b
c     (c)T(j,c)*T<in||ca>  *  (bk)<ba||ik><bn||jk>  (14B)
c     i=>e,n<=>c,j=>f (watch weighting) then n=>i
c     (i)T(f,i)*T<ec||ia>  *  (bk)<ba||ek><bc||fk>  (14B)
c     (and of course occupations not intentional)
c                 1 2 3 4  5 6  7 8
c                 a e c f  b k  i j  fac 
c                 b i n k  a j  c 0
c      rationale :
c        b,k (5,6) : summation indices of first contraction
c        i   (7)   : summation index of 2nd contract ; j=0
c        f   (4)   : other index of T1 (always occ)
c        still a,e / c,f on same
c      **this will allow it all in same routine**
c..     does require a bit of testing on j=0
c
c     *note*
c     b,k are running indices for first contraction type(*,1)
c     a,e are on same and c,f are on same "integral" for ,1
c     n is running index for second contraction type(*,2)
c     the second contraction is quite different (T1)
c     --------------------------------------------------------
c     e,f are complete ranges / a,c get batched
c     weights are applied to all except c (now)
c     --------------------------------------------------------
c     
c     *note* due to the t1 contraction that global  array is weird
c     also it's always a2a but differently
c
      integer occupation(9),ncase,cases(9,ncase)
      double precision value_X1X
      character*(*) diagram,type(2,2)
c
      return
      end
      subroutine mult_cct_t1(bb,kb,ke,eb,fb,
     1                    typeA,ab,ae, typeB,cb,ce,
     1                    g_RESULT)
c
      implicit none
#include "errquit.fh"
c*********************** incomplete develop
c
c...   general (:)) routine to do t1*t or t1*int (asymmetric)
c...   notation taken from 11A T1*T2 but translated a bit 
c...                  (c)T(k,c)*T(ij,cf) =>
c...                  (b)T(k,b)*T(ij,bf)
c...   should be able to do as well
c...                  (bk)T(k,b)*T(kj,bf)
c
c
#include "mafdecls.fh"
#include "cct_table_uhf.fh"
c
      integer g_RESULT
c
      integer bb,kb,ke,eb,fb,ab,ae,cb,ce
      character*(*) typeA,typeB
c
      logical ga_create,ga_destroy
      integer ga_memory_avail,ga_nnodes,avail_ga
      integer g_MATA,g_MATB,size
      integer se,sf,sbk,ee,fe,la,lc,lea,lcf,lbk,llbk,len
      integer batch,nbatch,bkb,bke
      integer spina,spine,spinc,spinf,spinb,spink
      integer occe,occf,occb,occk
      logical symmetric
c
c...   get available global memory
c
       avail_ga = ga_memory_avail()/8
       call ga_igop(13,avail_ga,1,'min')
c...    remember that there is overhead (in space) as well
        avail_ga = ga_nnodes()*avail_ga*0.9
c
c...  figure out spins
c
      spina = ispin(ab)
      spinc = ispin(cb)
      spine = ispin(eb)
      spinf = ispin(fb)
      spinb = ispin(bb)
      spink = ispin(kb)
c...  figure occupations (1=occ,2-virt)
      occe = iocc(eb)
      occf = iocc(fb)
      occb = iocc(bb)
      occk = iocc(kb)
c
      symmetric = (ab.eq.ae.and.cb.eq.ce.and.typeA.eq.typeB.and.
     1             spine.eq.spinf.and.spine.eq.spinf)
      size = avail_ga
      if (.not.symmetric) size = size/2
c
c...  can use symmetry to save a get; we have to use sqrt(fac) then
c
      if (iro(ae).ne.iro(ab)) call errquit('sym',0, UNKNOWN_ERR)
      if (iro(ce).ne.iro(cb)) call errquit('sym',0, UNKNOWN_ERR)
c
      sbk = mult(iro(bb),iro(kb))
      se = mult(iro(ab),sbk)
      sf = mult(iro(cb),sbk)
      eb = s_occ_virt(se,spine,occe)
      ee = s_occ_virt(se+1,spine,occe)-1
      fb = s_occ_virt(sf,spinf,occf)
      fe = s_occ_virt(sf+1,spinf,occf)-1
      la = ae-ab+1
      lc = ce-cb+1
      lea = (ee-eb+1)*la
      lcf = (fe-fb+1)*lc
      lbk = lenxy(sbk,spinb,spink,occb+occk)
      llbk = lbk
      len = lbk*(la+lc)*(ee-eb+1)*(fe-fb+1)
      nbatch = (size-1)/len + 1
      len = (lbk-1)/nbatch + 1
c
c...  get global integral array's
c
      if (.not.ga_create(MT_DBL,len,lea,'matA(a,e)',0,0,g_MATA))
     1   call errquit('matA',0, GA_ERR)
      if (.not.symmetric) then
       if (.not.ga_create(MT_DBL,len,lcf,'matA(c,f)',0,0,g_MATB))
     1     call errquit('matB',0, GA_ERR)
      end if
c
c
      bke = 0
      do batch=1,nbatch
       bkb = bke + 1
       bke = min(bke+len,llbk)
       lbk = (bke-bkb+1)
       call get_cct(g_MATA,typeA,bkb,bke,spinb,spink,
     1              eb,ee,spine,ab,ae,spina)
       if (.not.symmetric)
     1 call get_cct(g_MATB,typeB,bkb,bke,spinb,spink,
     2              fb,fe,spinf,ab,ae,spina)
c
c...  weigth mata (with sqrt(fac) if symmetric)
c
******       call cct_weigth(g_MATA,bkb,bke,sbk,spinb,spink,symmetric)
c   
       call ga_dgemm('T','N',lea,lcf,lbk,1.0d0,g_MATA,g_MATB,
     1               1.0d0,g_RESULT)
c
      end do
c
      if (symmetric) then
       if (.not.ga_destroy(g_MATB)) call errquit('destroy MATB',0,
     &       GA_ERR)
      end if
      if (.not.ga_destroy(g_MATA)) call errquit('destroy MATA',0,
     &       GA_ERR)
c
      return
      end
