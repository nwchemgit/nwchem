      subroutine txs_uniq_pair(nq, isbl_q, ijpres2, nij_uniq, ij_uniq,
     $     itmp)
      implicit none
c
      integer nq                ! [input] Number of quartets
      integer isbl_q(*)         ! [input] Mapping to request quartets
      integer ijpres2(*)        ! [input] List of request ij's
      integer nij_uniq          ! [output] No. of unique ij's
      integer ij_uniq(*)        ! [output] List of unique ij's
      integer itmp(*)           ! [scratch] Work space dimension nq
c
c     Return the number and list of unique pairs in the current
c     request.  The number of such pairs is bounded by min(nq,npar(ijbl)).
c
c     Use an NlogN sort followed by linear pass.
c
c     Gather list of IJ's in the request and sort them
c
      integer iq, ijprev
c      
      do iq = 1, nq
         itmp(iq) = ijpres2(isbl_q(iq))
      enddo
      call txs_integer_sort(nq, itmp)
c
c     Pass thru list and determine the unique ones
c
      ijprev = -1
      nij_uniq = 0
      do iq = 1, nq
         if (itmp(iq) .ne. ijprev) then
            ijprev = itmp(iq)
            nij_uniq = nij_uniq + 1
            ij_uniq(nij_uniq) = ijprev
         endif
      enddo
c
c     Debug
c
      write(6,*) ' nij_uniq = ', nij_uniq
      write(6,*) (ij_uniq(iq),iq=1,nij_uniq)
c
      end
      subroutine txs_integer_sort(n,ra)
      implicit integer (a-z)
      dimension ra(n)
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end

         
