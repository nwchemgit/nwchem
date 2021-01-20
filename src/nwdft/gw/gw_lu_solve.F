      subroutine gw_lu_solve(matrix,vector,n)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"

      integer,intent(in) :: n
      double precision,intent(inout) :: matrix(n,n),vector(n)

      logical ok
      integer adrA, adrB, adrI
      integer dimA1, dimA2, dimB1, dimB2
      integer mpA, nqA, mpB, nqB, info, alen, blen
      integer ha, hb, hi, na
      integer me, i
      integer lda, ldb, elemA, elemB
      integer lipiv, block_dims_A(2), block_dims_B(2), blocks(2)
      integer gridA(2), gridB(2)
      integer descA(10),descB(10)

      integer,parameter :: nb = 64
      integer,external :: numroc

      me = ga_nodeid()

      call SLinit()

      mpA = numroc(n, nb, myrow, 0, nprow)
      nqA = numroc(n, nb, mycol, 0, npcol)
      mpB = numroc(n, nb, myrow, 0, nprow)
      nqB = numroc(1, nb, mycol, 0, npcol)
      lda = max(1,mpA)
      ldb = max(1,mpB)

      elemA = mpA*nqA
      if (elemA.ne.0) ok = ma_push_get(mt_dbl,elemA,'a',ha,adra)
      
      elemB = mpB*nqB
      if (elemB.ne.0) ok = ma_push_get(mt_dbl,elemB,'b',hb,adrb)

      lipiv = mpa+nb
      ok = ma_push_get(mt_int,lipiv,'ipiv',hi,adri)

      call gw_to_SL(matrix,n,n,nb,nb,dbl_mb(adrA),lda,mpA,nqA)
      call gw_to_SL(vector,n,1,nb,nb,dbl_mb(adrB),ldb,mpB,nqB)

      call descinit(descA, n, n, nb, nb, 0, 0, iSLctxt, lda, info)
      if (info.ne.0) then
        call ga_error('gw_lusolve: descinit A failed',-info)
      endif

      call descinit(descB, n, 1, nb, nb, 0, 0, iSLctxt, ldb, info)
      if (info.ne.0) then
        call ga_error('gw_lusolve: descinit B failed',-info)
      endif

      call pdgetrf(n, n, dbl_mb(adrA), 1, 1, descA, int_mb(adrI), info)
      if (info.ne.0) then
        call ga_error('gw_lusolve: pdgetrf failed',-info)
      endif

      call pdgetrs('n',n,1,dbl_mb(adrA),1,1,descA,int_mb(adrI),
     $              dbl_mb(adrB),1,1,descB,info)
      if (info.ne.0) then
        call ga_error('gw_lusolve: pdgetrs failed',-info)
      endif

      vector = 0d0
      call gw_from_SL(vector,n,1,nb,nb,dbl_mb(adrB),ldb,mpB,nqB)
      call ga_dgop(dbl_mb,vector,n,'+')

      ok = ma_pop_stack(hi)
      if (elemB.ne.0) ok = ma_pop_stack(hb)
      if (elemA.ne.0) ok = ma_pop_stack(ha)

      call slexit()

      return
      end

      subroutine gw_to_SL(array, dim1, dim2, nbr, nbc, s_a, lda, np, nq)
      implicit none
#include "scalapack.fh"
      integer,intent(in) :: dim1, dim2, nbr, nbc, lda, np, nq
      double precision, intent(in) :: array(dim1,dim2)
      double precision, intent(out) ::  s_a(lda,*)
      integer row, col, tcol, trow, rbase, cbase, rbasel, cbasel
      integer pcol, prow
      integer rowl, coll
      integer marg1, marg2

      rbase = 1
      cbase = 1
      trow = nbr*nprow
      tcol = nbc*npcol
      do col=1,dim2,nbc
        pcol = mod(col,tcol)/nbc
        if(mycol.ne.pcol) cycle
        do row=1,dim1,nbr
          prow = mod(row,trow)/nbr
          if (myrow.ne.prow) cycle
          marg1 = row + nbr - 1
          rowl = min(marg1,dim1)
          marg1 = col + nbc - 1
          coll = min(marg1,dim2)
          rbasel = rowl - row + 1
          cbasel = coll - col + 1
          s_a(rbase:rbasel,cbase:cbasel) = array(row:rowl,col:coll)
          rbase = rbase + nbr
        enddo
        rbase = 1
        cbase = cbase + nbc
      enddo
      return
      end


      subroutine gw_from_SL(array,dim1,dim2,nbr,nbc,s_a,lda,np,nq)
      implicit none
#include "scalapack.fh"
      integer,intent(in) :: dim1, dim2, nbr, nbc, lda, np, nq
      double precision, intent(out) :: array(dim1,dim2)
      double precision, intent(in) ::  s_a(lda,*)
      integer row, col, tcol, trow, rbase, cbase, rbasel, cbasel
      integer pcol, prow
      integer rowl, coll
      integer marg1, marg2

      rbase = 1
      cbase = 1
      trow = nbr*nprow
      tcol = nbc*npcol
      do col=1,dim2,nbc
        pcol = mod(col,tcol)/nbc
        if(mycol.ne.pcol) cycle
        do row=1,dim1,nbr
          prow = mod(row,trow)/nbr
          if (myrow.ne.prow) cycle
          marg1 = row + nbr - 1
          rowl = min(marg1,dim1)
          marg1 = col + nbc - 1
          coll = min(marg1,dim2)
          rbasel = rowl - row + 1
          cbasel = coll - col + 1
          array(row:rowl,col:coll) = s_a(rbase:rbasel,cbase:cbasel)
          rbase = rbase + nbr
        enddo
        rbase = 1
        cbase = cbase + nbc
      enddo
      return
      end
