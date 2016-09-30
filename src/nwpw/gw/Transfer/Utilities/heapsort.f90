      subroutine hpsort(n,nm,ra)
! Numerical Recipes in Fortran 77, 2nd ed.  Press et al.  pg 329
      integer n,nm
      double precision ra(nm)
      integer i,ir,j,l
      double precision rra
      if (n.lt.2) return
      l=n/2+1
      ir=n
 100  continue
        if (l.gt.1) then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if (ir.eq.1) then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
 200    if (j.le.ir) then
          if (j.lt.ir) then
            if (ra(j).lt.ra(j+1)) j=j+1
          endif
          if (rra.lt.ra(j)) then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 200
        endif
        ra(i)=rra
      goto 100
      end

      subroutine indxhpsort(n,nm,ra,indx)
      integer n,nm
      double precision ra(nm)
      integer i,ir,j,l,indx(nm),indxt
      double precision rra
      do j=1,n
        indx(j)=j
      enddo
      if (n.lt.2) return
      l=n/2+1
      ir=n
 100  continue
        if (l.gt.1) then
          l=l-1
          indxt=indx(l)
          rra=ra(indxt)
        else
          indxt=indx(ir)
          indx(ir)=indx(1)
          rra=ra(indxt)
          ir=ir-1
          if (ir.eq.1) then
            indx(1)=indxt
            rra=ra(indxt)
            return
          endif
        endif
        i=l
        j=l+l
 200    if (j.le.ir) then
          if (j.lt.ir) then
            if (ra(indx(j)).lt.ra(indx(j+1))) j=j+1
          endif
          if (rra.lt.ra(indx(j))) then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 200
        endif
        indx(i)=indxt
      goto 100
      end
