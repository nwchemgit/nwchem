      subroutine gw_lu_solve(matrix,vector,n)
      implicit none
#include "mafdecls.fh"
#include "global.fh"

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

      integer maxproc,nprow,npcol,myrow,mycol,islctxt,ga_comm,nnodes
      logical oactive

      integer,parameter :: nb = 32
      integer,external :: numroc,slgetmxproc

      me = ga_nodeid()
      nnodes = ga_nnodes()
c
c     initialize scalapack context
c
      maxproc = slgetmxproc(n,nnodes)
      call FindGrid(maxproc, nprow, npcol)
      call ga_mpi_comm_pgroup_default(ga_comm)
      islctxt = ga_comm
      call blacs_gridinit(islctxt,'R',nprow,npcol)
      if (me.lt.maxproc) then
        call blacs_gridinfo(islctxt,nprow,npcol,myrow,mycol)
      else
        nprow = 0
        npcol = 0
        myrow = 0
        mycol = 0
      endif

      oactive = me.lt.maxproc
      if (.not.oactive) goto 100

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

      call gw_to_SL(matrix,n,n,nb,nb,dbl_mb(adrA),lda,mpA,nqA,nprow,
     $              npcol,myrow,mycol)
      call gw_to_SL(vector,n,1,nb,nb,dbl_mb(adrB),ldb,mpB,nqB,nprow,
     $              npcol,myrow,mycol)

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
      call gw_from_SL(vector,n,1,nb,nb,dbl_mb(adrB),ldb,mpB,nqB,nprow,
     $                npcol,myrow,mycol)

      ok = ma_pop_stack(hi)
      if (elemB.ne.0) ok = ma_pop_stack(hb)
      if (elemA.ne.0) ok = ma_pop_stack(ha)

      call blacs_gridexit(islctxt)

  100 continue
      if (.not.oactive) vector = 0d0
      call ga_sync()
      call ga_dgop(dbl_mb,vector,n,'+')

      return
      end

      subroutine gw_to_SL(array, dim1, dim2, nbr, nbc, s_a, lda, np, nq,
     $                    nprow, npcol, myrow, mycol)
      implicit none
      integer,intent(in) :: dim1, dim2, nbr, nbc, lda, np, nq
      integer,intent(in) :: nprow,npcol,myrow,mycol
      double precision, intent(in) :: array(dim1,dim2)
      double precision, intent(out) ::  s_a(lda,*)
      integer row, col, tcol, trow, rbase, cbase, rbasel, cbasel
      integer pcol, prow
      integer rowl, coll
      integer marg1, marg2
      integer r0i,r1i,r0im1,r1im1
      logical putpending


      rbase = 1
      cbase = 1
      trow = nbr*nprow
      tcol = nbc*npcol
      do col=1,dim2,nbc
        pcol = mod(col,tcol)/nbc
        if(mycol.ne.pcol) cycle
        r0im1=-9999
        r1im1=0
        putpending=.false.
        marg1 = col + nbc - 1
        coll = min(marg1,dim2)
        cbasel = cbase + coll - col 

        do row=1,dim1,nbr
          prow = mod(row,trow)/nbr
          if (myrow.ne.prow) cycle
          if(.not.putpending) then
            r0im1=row
            marg1=row+nbr-1
            marg2=dim1
            r1im1=min(marg1,marg2)
          endif
          r0i = row
          marg1 = row + nbr - 1
          r1i = min(marg1,dim1)
          if(r0i.eq.(r1im1+1).and.r1i.ne.dim1) then
            r1im1=r1i
            putpending=.true.
          else
            rowl = r1i
            rbasel = rbase + rowl - r0im1
            s_a(rbase:rbasel,cbase:cbasel) = array(r0im1:rowl,col:coll)
            rbase = rbase + rowl - r0im1 + 1
            r0im1 = -1
            r1im1 = -1
            putpending = .false.
          endif
        enddo
        rbase = 1
        cbase = cbase + nbc
      enddo
      return
      end


      subroutine gw_from_SL(array,dim1,dim2,nbr,nbc,s_a,lda,np,nq,
     $                      nprow,npcol,myrow,mycol)
      implicit none
      integer,intent(in) :: dim1, dim2, nbr, nbc, lda, np, nq
      integer,intent(in) :: nprow,npcol,myrow,mycol
      double precision, intent(out) :: array(dim1,dim2)
      double precision, intent(in) ::  s_a(lda,*)
      integer row, col, tcol, trow, rbase, cbase, rbasel, cbasel
      integer pcol, prow
      integer rowl, coll
      integer marg1, marg2
      integer r0i,r1i,r0im1,r1im1
      logical putpending


      rbase = 1
      cbase = 1
      trow = nbr*nprow
      tcol = nbc*npcol
      do col=1,dim2,nbc
        pcol = mod(col,tcol)/nbc
        if(mycol.ne.pcol) cycle
        r0im1=-9999
        r1im1=0
        putpending=.false.
        marg1 = col + nbc - 1
        coll = min(marg1,dim2)
        cbasel = cbase + coll - col
        do row=1,dim1,nbr
          prow = mod(row,trow)/nbr
          if (myrow.ne.prow) cycle
                  if(.not.putpending) then
                     r0im1=row
                     marg1=row+nbr-1
                     marg2=dim1
                     r1im1=min(marg1,marg2)
                  endif
                  r0i=row
          marg1 = row + nbr - 1
          r1i = min(marg1,dim1)
          if(r0i.eq.(r1im1+1).and.r1i.ne.dim1) then
            r1im1=r1i
            putpending=.true.
          else
            rowl = r1i
            rbasel = rbase + rowl - r0im1
            array(r0im1:rowl,col:coll) = s_a(rbase:rbasel,cbase:cbasel)
            rbase = rbasel + 1
            r0im1 = -1
            r1im1 = -1
            putpending = .false.
          endif
        enddo
        rbase = 1
        cbase = cbase + nbc
      enddo
      return
      end
