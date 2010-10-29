      subroutine jacod(f,v,nb,nb1,nb2,nmin,nmax,big,jbig,maxao)
C$Id$
      implicit none
c.......................................................................
c      f is the matrix to be diagonalized. f is stored triangular
c      v is the array of eigenvectors. quadratic array, dimension nb*nb
c      big and jbig are temporary scratch areas of dimension nb
c      the rotations among the first nmin basis functions are not
c      accounted for.
c      the rotations among the last nb-nmax basis functions are not
c      accounted for.
c.......................................................................
      integer nb, nb1, nb2, nmin, nmax, maxao
      integer jbig(nb)
      real *8 big(nb),f(nb1),v(nb2)
c
      integer i, k, ir
      integer ii, j, jj, i1, ib, ia, iaa, ibb, jaa, jbb
      integer iear, iebr, ieaa, ieab, it, kq, ir1
      real *8 sd, t, dif, t2x2, t2x25, cx, sx, tt
c
      real *8 root2, c2, c3, c4, c5, c6
      data root2 /0.707106781186548d0/
      data c2,c3,c4,c5,c6/1.0d-12,4.0d-16,2.0d-16,
     +                        8.0d-9,3.0d-6/
c
      ia = 0
      ib = 0 
      iaa = 0
      ibb = 0
      ieaa = 1
      ieab = 1
c
      if (nb.eq.1) then
         v(1) = 1.0d0
         return
      end if
      ii = 0
c.......................................................................
c      loop over rows (i) of triangular matrix
c.......................................................................
      do i = 1 , nb
         big(i) = 0.00d0
         jbig(i) = 0
         if (i.ge.nmin .and. i.ne.1) then
            j = min0(i-1,nmax)
c.......................................................................
c      loop over columns (k) of triangular matrix to determine
c      largest off-diagonal elements in row(i).
c.......................................................................
            do k = 1 , j
               if (abs(big(i)).lt.abs(f(ii+k))) then
                  big(i) = f(ii+k)
                  jbig(i) = k
               end if
            enddo
         end if
         ii = ii + i
      enddo
 40   sd = 1.050d0
c.......................................................................
c      find smallest diagonal element and corresponding largest
c      off-diagonal element.
c.......................................................................
      jj = 0
      do j = 1 , nb
         jj = jj + j
         sd = min(sd,abs(f(jj)))
      enddo
      sd = max(sd,c6)*c2
      i1 = max0(2,nmin)
      t = abs(big(i1))
      ib = i1
      do i = i1 , nb
         if (t.lt.abs(big(i))) then
            t = abs(big(i))
            ib = i
         end if
      enddo
c.......................................................................
c      test for convergence, then determine rotation.
c.......................................................................
      ia = jbig(ib)
      if (t.lt.sd) then
         return
      else
         iaa = ia*(ia-1)/2
         ibb = ib*(ib-1)/2
         jaa = (ia-1)*maxao
         jbb = (ib-1)*maxao
         dif = f(iaa+ia) - f(ibb+ib)
         if (abs(dif).gt.c3*t) then
            t2x2 = big(ib)/dif
            t2x25 = t2x2*t2x2
            if (t2x25.le.c4) then
               cx = 1.0d0
               sx = t2x2
            else if (t2x25.le.c5) then
               sx = t2x2*(1.0d0-1.5d0*t2x25)
               cx = 1.0d0 - 0.5d0*t2x25
            else if (t2x25.gt.c6) then
               t = 0.25d0/sqrt(0.25d0+t2x25)
               cx = sqrt(0.5d0+t)
               sx = dsign(sqrt(0.5d0-t),t2x2)
            else
               cx = 1.0d0 + t2x25*(t2x25*1.375d0-0.5d0)
               sx = t2x2*(1.00d0+t2x25*(t2x25*3.875d0-1.5d0))
            end if
         else
            sx = root2
            cx = root2
         end if
         iear = iaa + 1
         iebr = ibb + 1
         do 90 ir = 1 , nb
            t = f(iear)*sx
            f(iear) = f(iear)*cx + f(iebr)*sx
            f(iebr) = t - f(iebr)*cx
            if (ir.lt.ia) then
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else if (ir.eq.ia) then
               tt = f(iebr)
               ieaa = iear
               ieab = iebr
               f(iebr) = big(ib)
               iear = iear + ir - 1
               if (jbig(ir).ne.0) go to 70
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else
               t = f(iear)
               it = ia
               iear = iear + ir - 1
               if (ir.lt.ib) then
               else if (ir.eq.ib) then
                  f(ieaa) = f(ieaa)*cx + f(ieab)*sx
                  f(ieab) = tt*cx + f(iebr)*sx
                  f(iebr) = tt*sx - f(iebr)*cx
                  iebr = iebr + ir - 1
                  go to 70
               else
                  if (abs(t).lt.abs(f(iebr))) then
                     if (ib.le.nmax) then
                        t = f(iebr)
                        it = ib
                     end if
                  end if
                  iebr = iebr + ir - 1
               end if
            end if
            if (abs(t).ge.abs(big(ir))) then
               big(ir) = t
               jbig(ir) = it
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else if (ia.ne.jbig(ir) .and. ib.ne.jbig(ir)) then
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            end if
 70         kq = iear - ir - ia + 1
            big(ir) = 0.00d0
            ir1 = min0(ir-1,nmax)
            do i = 1 , ir1
               k = kq + i
               if (abs(big(ir)).lt.abs(f(k))) then
                  big(ir) = f(k)
                  jbig(ir) = i
               end if
            enddo
            iear = iear + 1
            iebr = iebr + 1
 90      continue
         do i = 1 , maxao
            t = v(jbb+i)*sx
            v(jbb+i) = v(jaa+i)*sx - v(jbb+i)*cx
            v(jaa+i) = v(jaa+i)*cx + t
         enddo
         go to 40
      end if
      end
