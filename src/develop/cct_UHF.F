      subroutine CCT_UHF 
c 
c....  development version UHF CCSD(T)
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
      implicit none
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
c
      double precision r_a2a,r_a2b,r_b2b,r_a3a,r_a3b,r_b3b
      double precision r_a4a,r_a4b,r_b4b
      double precision cct_uhf_D
c
      integer ncase_a2a,ncase_a2b,ncase_b2b
      integer ncase_a3a,ncase_a3b,ncase_b3b
      integer ncase_a4a,ncase_a4b,ncase_b4b
      parameter (ncase_a2a=5,ncase_a2b=6,ncase_b2b=5)
      parameter (ncase_a3a=5,ncase_a3b=6,ncase_b3b=5)
      parameter (ncase_a4a=4,ncase_a4b=6,ncase_b4b=4)
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
c...  and (to be save) the inices to be weighted in the intermediate)
c...  type contains types for intgrals/t's in first and second mult
c                     1 2 3 4  5 6  7 8
c                     a e c f  b k  i j, fac 
      data cases_a3a/ 1,1,1,1, 0,0, 1,1,  1,
     2                2,2,2,2, 0,0, 2,2,  1,
     3                1,2,1,2, 2,1, 1,2,  2,
     4                2,1,2,1, 1,2, 1,2,  2,
     5                1,1,2,2, 0,0, 1,2,  4/
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
      r_a2a = cct_uhf_D(a2a,occ_a2a,cases_a2a,ncase_a2a,type_a2a)
      r_a2b = cct_uhf_D(a2b,occ_a2b,cases_a2b,ncase_a2b,type_a2b)
      r_b2b = cct_uhf_D(b2b,occ_b2b,cases_b2b,ncase_b2b,type_b2b)
c
      r_a3a = cct_uhf_D(a3a,occ_a3a,cases_a3a,ncase_a3a,type_a3a)
      r_a3b = cct_uhf_D(a3b,occ_a3b,cases_a3b,ncase_a3b,type_a3b)
      r_b3b = cct_uhf_D(b3b,occ_b3b,cases_b3b,ncase_b3b,type_b3b)
c
      r_a4a = cct_uhf_D(a4a,occ_a4a,cases_a4a,ncase_a4a,type_a4a)
      r_a4b = cct_uhf_D(a4b,occ_a4b,cases_a4b,ncase_a4b,type_a4b)
      r_b4b = cct_uhf_D(b4b,occ_b4b,cases_b4b,ncase_b4b,type_b4b)
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
        ispin(i) = spin
       end do
       do i=1,virt(spin)
        k = k + 1
        iocc(k) = 2
        ispin(i) = spin
       end do
      end do
c
      do i=1,norb
       iro(i) = iro_in(i)
       iro(i+norb) = iro_in(i+norb)
       do spin=1,2
        k = i+(spin-1)*norb
        if (i.lt.occ(spin)) then
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
      return
      end
      
      subroutine cct_gr_ia(rint,nia,nbc,cb,ce,db,de,type,
     1                     spc,spd,spi,spa)
c
c...  get ia integrals (type determines which ones
c
      double precision rint(nia,*)
      integer  nbc,cb,ce,db,de,spc,spd,spi,spa
      character*(*) type
c
#include "cct_table_uhf.fh"
c
      integer c,d
c
c...  checks (e.g. within one symmetry)
c
      if (iro(cb).ne.iro(ce).or.iro(db).ne.iro(de)) 
     1   call errquit('symmetry in cct_gr_ia',0)
      scd = mult(iro(cb),iro(db))
c
      if (type.eq.'(ia|cd)') then
       ncd = 0
       do c=cb,ce
        do d=db,min0(de,c)
         ncd = ncd + 1
         call errquit('gather ia integrals in rint(*,ncd)',-999)
        end do
       end do
      else if (type.eq.'<id||ac>') then
       nbc = 0
       do b=bb,be
        do c=cb,ce
         nbc = nbc + 1
         call errquit('gather ia integrals in rint(*,ncd)',-999)
        end do
       end do
      else
       call errquit('unknown type in cct_gr_ia',0)
      end if
c
      return
      end
      subroutine cct_gt_ia(t2,nia,jb,je,bb,be,spj,spb,spi)
c
c...  get ia t2's 
c
      double precision t2(nia,*)
      integer  jb,je,bb,be
c
#include "cct_table_uhf.fh"
c
      integer j,b,spj,spb,spi
c
c...  checks (e.g. within one symmetry)
c...  (spa=deducable from spj,spb,spi)
c
      if (iro(jb).ne.iro(je).or.iro(bb).ne.iro(be)) 
     1   call errquit('symmetry in cct_gt_ia',0)
      scd = mult(iro(cb),iro(db))
c
       njb = 0
       do j=jb,je
        do b=bb,be
         njb = njb + 1
         call errquit('gather ia t2s',-999)
        end do
       end do
c
      return
      end 
      double precision function cct_uhf_D(diagram,occ,cases,ncase,type)
c
c...  CCT diagram evaluator; modelled for A3A
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
c...           occ,cases,ncase,type : diagram definition (see below)
c
c     integer ncase,case
c     parameter (ncase=5)
c     integer cases(9,ncase),occ(9)
c     integer c,e,a,f,i,j,b,k
c...    (occupations of these)
c     character*20 type(2,2)
c...  cases contains the spin cases and the factor to multiply by
c...  occ contains the occupations (1=occ,2=virt) and the divide factor
c...  type contains types for intgrals/t's in first and second mult
c                 1 2 3 4  5 6  7 8
c                 a e c f  b k  i j  fac 
c     data cases/ 1,1,1,1, 0,0, 1,1,  1,
c    1            2,2,2,2, 0,0, 2,2,  1,
c    2            1,2,1,2, 2,1, 1,2,  2,
c    3            2,1,2,1, 1,2, 1,2,  2,
c    4            2,1,1,2, 0,0, 1,2,  4/
c     data occ  / 2,2,2,2, 1,2, 1,1,  2/
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
      integer occ(9),ncase,cases(9,ncase)
      character*(*) diagram,type(2,2)
c
#include "mafdecls.fh"
c
      logical ga_create,ga_destroy
      integer spinc,spine,spina,spinf,spini,spinj,spinb,spink
      double precision factor
c
c
c...   symbolic matrix element
c
      integer case
      integer c,e,a,f,i,j,b,k
c...    (occupations of these)
c
c...   get available global memory
c
       avail_ga = ga_memory_avail()/8
       call ga_igop(13,avail_ga,1,'min')
c...    remember that there is overhead (in space) as well
       ga_memory = ga_nnodes()*avail_ga*0.9
c     
c                 1 2 3 4  5 6  7 8
c                 a e c f  b k  i j  fac 
      a = occ(1)
      e = occ(2)
      c = occ(3)
      f = occ(4)
      b = occ(5)
      k = occ(6)
      i = occ(7)
      j = occ(8)
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
       factor = cases(9,case)/(occ(9)*1.0D0)

c
       bbe = s_occ_virt(1,spine,e)
       bbf = s_occ_virt(1,spinf,f)
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
c...   prepare batches over c and a
c
          symmetry = (se.eq.sf).and.spin.eq.spin 
          size = ne*nf*n_occ_virt(sa,a)*n_occ_virt(sc,c) 
c         if (.not.symmetry) size = size * 2 ????????
          nbatch = ga_memory/size 
c...   'squares' for now 
          na = dsqrt((nbatch-1)*1.0d0) + 1.00001d0 
          nc = (nbatch-1)/nbatcha + 1 
          na = (n_occ_virt(sa,a)-1)/na+1 
          nc = (n_occ_virt(sc,c)-1)/nc+1 
c
          do ab=s_occ_virt(sa,nspin,a),s_occ_virt(sa+1,nspin,a)-1,na
           ae = min0(ab+na,s_occ_virt(sa+1),a)-1
           do cb=s_occ_virt(sc,nspin,c),s_occ_virt(sc+1,nspin,c)-1,nc
            ce = min0(cb+nc,s_occ_virt(sc+1),c)-1
c
            lae = (ae-ab+1)*ne
            lcf = (ce-cb+1)*nf
            if (.not.ga_create(MT_DBL,lae,lcf,'Result',0,0,g_RESULT))
     1          call errquit('RESULT',0)
            call ga_zero(g_RESULT)
c
c...    loop over bk spin if sbk = 0 
c
            spin = spinb
            if (spin.eq.0.and.spink.ne.0) 
     1          call errquit('whats this spin',0)
            spinb = max(spinb,1)
            spink = max(spink,1)
5           bb = s_occ_virt(sbk,spinb,b)
            kb = s_occ_virt(sbk,spink,k)
c...   * I(aa+bb)(aaaa)
            call mult_cct(bb,kb,bbe,bbf,
     1                    type(1,1),ab,ae,
     1                    type(2,1),cb,ce,
     1                    g_RESULT)
            if (spin.eq.0) then
             spin = 1
             spinb = 2
             spink = 2
             go to 5
            end if
c      * T(aa)(aaaa)
            if (.not.ga_create(MT_DBL,laf,lce,'Tesult',0,0,g_TESULT))
     1          call errquit('TESULt',0)
            call ga_zero(g_TESULT)
c
            spin = spini
            if (spin.eq.0.and.spinj.ne.0) 
     1          call errqui('whats this spini',0)
            spini = max(spini,1)
            spinj = max(spinj,1)
6           ib = s_occ_virt(sbk,spini,i)
            jb = s_occ_virt(sbk,spinj,j)
            call mult_cct(ib,jb,bbe,bbf,
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
            sum = sum + cct_com(g_RESULT,g_TESULT,ab,ae,cb,ce,
     1                          bbe,bbf,'N')
c
            if (cb.eq.ab) then
             sum = sum + cct_com(g_RESULT,g_TESULT,ab,ae,cb,ce,
     1                           bbe,bbf,'T')
            else
             call ga_zero(g_TESULT)
             spin = spini
             spini = max(spini,1)
             spinj = max(spinj,1)
7            ib = s_occ_virt(sbk,spini,i)
             jb = s_occ_virt(sbk,spinj,j)
             call mult_cct(ib,jb,bbe,bbf,
     1                     type(1,2),cb,ce,
     1                     type(2,2),ab,ae,
     1                      g_TESULT)
             if (spin.eq.0) then
              spin = 1
              spini = 2
              spinj = 2
              go to 7 
             end if
c
             sum = sum + cct_com(g_RESULT,g_TESULT,ab,ae,cb,ce,
     1                           bbe,bbf,'N')
            end if
c
            if (.not.ga_destroy(g_TESULT)) 
     1          call errquit('destroy TESULT',0)
            if (.not.ga_destroy(g_RESULT)) 
     1          call errquit('destroy RESULT',0)
c
           end do
          end do
10       end do
        end do
       end do
      end do
c
      return
      end
      subroutine mult_cct(bb,kb,eb,fb,
     1                    typeA,ab,ae, typeB,cb,ce,
     1                    g_RESULT)
c
      implicit none
c
c...   general (:)) routine to do t*t or int*int (asymmetric)
c...   notation taken from A3A Integral*Integral
c...                       (bk)<ba||ek><bc||fk>
c...                        integrals as (bk,ea) 
c
c...   the first bb and bk are used to denote the full range
c...   of that spin/space/occupation "symmetry" => bbb,bbk
c...   the first be and bf are used to denote full range of
c...   that  spin/occupation "symmetry"  => bbbe,bbbf
c
c...   IS SYMMETRY FULLY EXPLOITED ????????????
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
      integer ga_memory_avail,ga_nnodes,avail_ga,avail_ga
      integer g_MATA,g_MATB,size
      integer se,sf,sbk,eb,ee,fb,fe,la,lc,lea,lcf,lbk,llbk,len
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
      if (iro(ae).ne.iro(ab)) call errquit('sym',0)
      if (iro(ce).ne.iro(cb)) call errquit('sym',0)
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
     1   call errquit('matA',0)
      if (.not.symmetric) then
       if (.not.ga_create(MT_DBL,len,lcf,'matA(c,f)',0,0,g_MATB))
     1     call errquit('matB',0)
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
       call cct_weigth(g_MATA,bkb,bke,sbk,spinb,spink,symmetric)
c   
       call ga_dgemm('T','N',lea,lcf,lbk,1.0d0,g_MATA,g_MATB,
     1               1.0d0,g_RESULT)
c
      end do
c
      if (symmetric) then
       if (.not.ga_destroy(g_MATB)) call errquit('destroy MATB',0)
      end if
      if (.not.ga_destroy(g_MATA)) call errquit('destroy MATA',0)
c
      return
      end
      double precision function cct_uhf_a4a(rint,t2,res1,res2,niabc)
c
c     niabc has to be calculated so that the corresponding
c     t2 and rint fit also
c
c     terms A4A 
c       t(ij,ce) <ab||ek> t(kj,af) <cb||fi>            raa_aa*raa_aa
c       t(ij,ce) (ae|bk)-ab t(jk,af)-ab (ci|bf)-ab     raa_ab*rab_ab
c       t(ij,ce) (ak|be)-ba t(jk,fa)-ab <cb||fi>       rab_aa*rab_ba
c      -t(ij,ec)-ab <ab||ek> t(kj,af)-ab (cf|bi)-ba    rab_ba*rab_aa
c       t(ij,ce)-ab (ak|be)-ab t(kj,af) (ci|bf)-ab     rab_aa*raa_ab
c      -t(ij,ce)-ab (ae|bk)-ba t(kj,fa)-ba <cb||fi>    rab_aa*rba_ba
c      -t(ji,ec)-ab <ab||ek> t(kj,af) (ci|bf)-ba       rab_ba*raa_aa
c       t(ji,ce)-ab (ak|be)-ab t(kj,af) (cf|bi)-ab     rab_ab*raa_ab
c
c       and all b(eta),a(lpha) interchanged
c
c       reorder and relabel (running index always ia)
c       (order 1-4,2-3 ; in 1-4 ic=>ia in 2-3 ka=>ia)
c        e=b,b=c,f=d summation is always (ia)
c        (note left ia is not *really* right ia)
c
c   t(ij,ab)    <ac||di>   t(ij,ad)    <ac||bi>      raa_aa*raa_aa
c   t(ij,ab)    (ai|cd)-ab t(ji,ad)-ab (ab|ci)-ab    raa_ab*rab_ab
c   t(ij,ab)    <ac||di>   t(ji,da)-ab (ai|cb)-ba    rab_aa*rab_ba
c  -t(ij,ba)-ab (ad|ci)-ba t(ij,ad)-ab <ac||bi>      rab_ba*rab_aa
c   t(ij,ab)-ab (ai|cd)-ab t(ij,ad)    (ai|cb)-ab    rab_aa*raa_ab
c  -t(ij,ab)-ab <ac||di>   t(ij,da)-ba (ab|ci)-ba    rab_aa*rba_ba
c  -t(ji,ba)-ab (ai|cd)-ba t(ij,ad)    <ac||bi>      rab_ba*raa_aa
c   t(ji,ab)-ab (ad|ci)-ab t(ij,ad)    (ai|cb)-ab    rab_ab*raa_ab
c+beta
c   t(ij,ab)    <ac||di>   t(ij,ad)    <ac||bi>      raa_aa*raa_aa
c   t(ij,ab)    (ai|cd)-ba t(ji,ad)-ba (ab|ci)-ba    raa_ab*rab_ab
c   t(ij,ab)    <ac||di>   t(ji,da)-ba (ai|cb)-ab    rab_aa*rab_ba
c  -t(ij,ba)-ba (ad|ci)-ab t(ij,ad)-ba <ac||bi>      rab_ba*rab_aa
c   t(ij,ab)-ba (ai|cd)-ba t(ij,ad)    (ai|cb)-ba    rab_aa*raa_ab
c  -t(ij,ab)-ba <ac||di>   t(ij,da)-ab (ab|ci)-ab    rab_aa*rba_ba
c  -t(ji,ba)-ba (ai|cd)-ab t(ij,ad)    <ac||bi>      rab_ba*raa_aa
c   t(ji,ab)-ba (ad|ci)-ba t(ij,ad)    (ai|cb)-ba    rab_ab*raa_ab
c ====   = bb
c   t(ij,ab)    <ac||di>   t(ij,ad)    <ac||bi>      raa_aa*raa_aa
c   t(ij,ab)    (ai|cd)-ba t(ij,da)-ab (ab|ci)-ba    raa_ab*rab_ab
c   t(ij,ab)    <ac||di>   t(ij,ad)-ab (ai|cb)-ab    rab_aa*rab_ba
c  -t(ji,ab)-ab (ad|ci)-ab t(ji,da)-ab <ac||bi>      rab_ba*rab_aa
c   t(ji,ba)-ab (ai|cd)-ba t(ij,ad)    (ai|cb)-ba    rab_aa*raa_ab
c  -t(ji,ba)-ab <ac||di>   t(ij,da)-ab (ab|ci)-ab    rab_aa*rba_ba
c  -t(ij,ab)-ab (ai|cd)-ab t(ij,ad)    <ac||bi>      rab_ba*raa_aa
c   t(ij,ba)-ab (ad|ci)-ba t(ij,ad)    (ai|cb)-ba    rab_ab*raa_ab
c...    group for ia (internal indices)
c
c      pure spin case : coul-exch 
c       t(ij,ce) (ab||ek) t(kj,af) (cb||fi)
c     = t(ij,ce) (cb||fi) t(kj,af) (ab||ek)
c        (ic)=>(ia)  / (ka)=>(ia)
c     = (ia) t(ij,ae) (ab||fi) * (ia) t(ij,af) (ab||ei)
c        (e=b,b=c,f=d)
c     = (ia) t(ij,ab) (ac||di) * (ia) t(ij,ad) (ac||bi)
c     sym : ia=jb=cd and b=d
      implicit none
c
      integer niabc
      double precision rint(*),t2(*),res1(niabc),res2(*)
c
#include "cct_table_uhf.fh"
      double precision sum,jbcd_jdcb
      integer spin
      integer c,sc,sbd,sb,sd,scd,sj,jb,je,bb,db
c
c...  loop over symmetries
c
c...  t2 * integral (ic) t(ij,ce) (cf,bi) => (j,e,b,f)
c...  t2 * integral (ia) t(ij,ab) (ia|cd) => (j,b,c,d) (coul-exch)
c     (ia) t(ij,ab) (ac||di) * (ia) t(ij,ad) (ac||bi)
c         => (jbcd) * (jdcb) so keep all b,d
c
c...   could batch over c
c
       spin = 1
       do c=occ(spin)+1,virt(spin)
        sc = iro(c)
        do sbd =1,nirr
         do sb=1,nirr
          bb = svirt(sb,spin)
          sd = mult(sbd,sb)
          db = svirt(sd,spin)
          scd = mult(sc,sd)
          if (sd.lt.sb) go to 100
          sj = mult(scd,sb)
          jb = socc(sj,spin)
          je = jb-1 + nocc(sj,spin)
c
c     type = '<id||ca>'
c
c     = (ia) t(ij,ae) (ab||fi) * (ia) t(ij,af) (ab||ei)
c     sym : ia=jb=cd
c         => (jebf) * (jfbe) so keep all e,f
c            (jbcd) * (jdcb)
c...  ia aa
          call cct_mult_ia(t2,rint,jb,je,bb,c,c,db,
     1                     res1, '<id||ca>','ia',
     2                     1,1,1,1,1,1)
c...  ia bb
          call cct_mult_ia(t2,rint,jb,je,bb,c,c,db,
     1                     res1, '<id||ca>','ia',
     2                     1,1,1,1,2,2)
          if (sb.ne.sd) then
c...  ia aa
           call cct_mult_ia(t2,rint,jb,je,db,c,c,bb,
     1                      res2, '<id||ca>','ia',
     2                      1,1,1,1,1,1)
c...  ia bb
           call cct_mult_ia(t2,rint,jb,je,db,c,c,bb,
     1                      res2, '<id||ca>','ia',
     2                      1,1,1,1,2,2)
          end if
          if (sb.eq.sd) then
            sum = sum + jbcd_jdcb(res1,res1,
     1               svirt(sb,1),nvirt(sb,1),jb,je-jb+1,nvirt(sd,1),1)
          else
            sum = sum + 2.0d0*jbcd_jdcb(res1,res2,
     1               svirt(sb,1),nvirt(sb,1),jb,je-jb+1,nvirt(sd,1),1)
          end if
c
        end do
100    end do
      end do
c
      cct_uhf_a4a = sum
c
      return
      end 
      double precision function jbcd_jdcb(res1,res2,bb,nb,nd,jb,nc,nj)
c
c...  result is nb(rows) * nd(columns) ordered as result(b,d,c,j)
c...  do part of jbcd + jbdc summation
c...  this is a place-holder (but n**4)
c
      double precision res1(nb,nd,nc,nj),res2(nd,nb,nc,nj)
      integer bj,bb,nj,nb,nc,nd
#include "cct_table_uhf.fh"
c
      integer j,b,c,d
c
      jbcd_jdcb = 0.0d0
      do j=1,nj
       do c=1,nc
        do d=1,nd
         do b=1,nb
          jbcd_jdcb = jbcd_jdcb + res1(b,j,d,c)*res2(d,j,b,c)*
     1                            fac(j+bj-1)*fac(b+bb-1)
         end do
        end do
       end do
      end do
c
      return
      end
c
