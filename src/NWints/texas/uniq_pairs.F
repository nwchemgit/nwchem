* $Id$
      subroutine uniq_pairs_1(from,ibl,kbl,bl,ijpres2,klpres2, 
     *                        isbl_size,isbl_q,isbl_point)
      implicit real*8 (a-h,o-z)
      character*8 from
c this is maximum pairs allowed ; must be: maxpar.le.limpair
      common /max_pairs/ maxpar   
      common /time5/ tcheck,tuniq2,tmap4u
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,    ! output
     *                  nkl_uniqe, kl_uniqe_p, map_kl     ! output
      dimension bl(*)
      dimension ijpres2(*),klpres2(*) ! dim=nqrtpnl
      dimension isbl_q(*) ! dim=nqrtpnl
c----------------------------------------------------------------------
c This routine is called ONLY :
c
c from spec_calcint.f (blockint)    from='blockint'
c               and 
c from c texas_hf.F (supblks_split) from='ks_split'
c
c for a given super-block ISBL,made out of Ibl & Kbl pair-blocks 
c find uniqe (present) pairs in pair-blocks
c----------------------------------------------------------------------
      call txs_second(timb)
c----------------------------------------------------------------------
c memory allocations :
c
      if(from.eq.'blockint') then
         call memo1_int( maxpar, ij_uniqe)  
         call memo1_int( maxpar, kl_uniqe) 
      endif
      if(from.eq.'ks_split') then
         call memo1_int(isbl_size, ij_uniqe)  
         call memo1_int(isbl_size, kl_uniqe) 
      endif
c----------------------------------------------------------------------
      ipoint=isbl_point + 1
c
            call find_uniq_2(bl,isbl_size,isbl_q(ipoint),ijpres2,
     *                       bl(ij_uniqe) ,    nij_uniqe )
            call find_uniq_2(bl,isbl_size,isbl_q(ipoint),klpres2,
     *                       bl(kl_uniqe) ,    nkl_uniqe )
c
      ij_uniqe_p=ij_uniqe
      kl_uniqe_p=kl_uniqe
c
c on return :
c number of pairs and list of pairs for Ibl block (nij_uniqe & ij_uniqe)
c number of pairs and list of pairs for Kbl block (nkl_uniqe & kl_uniqe)
c
c kept in 
c     common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
c    *                  nkl_uniqe, kl_uniqe_p, map_kl
c----------------------------------------------------------------------
c 2 calls to getmem
c--------------------------------------------------------
      call txs_second(time)
      tuniq2=tuniq2+time-timb
c--------------------------------------------------------
      end
c===============================================================
      subroutine find_uniq_2(bl,isbl_size,isbl_q,ijpres2,
     *                       ibl_uniqe ,nib_uniqe)
      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension ijpres2(*)    ! dim=nqrtpnl
      dimension isbl_q(*)  
      dimension ibl_uniqe(*)  ! input-output
c
cjun22call getmem(isbl_size,ij_temp)
      call memo1_int(isbl_size,ij_temp)
c
         call do_tmp1(isbl_size,isbl_q,ijpres2,bl(ij_temp),iqu )
         call do_uniq(iqu,bl(ij_temp),   ibl_uniqe ,ij_u)
c
      call retmem(1)  ! release ij_temp
c
      nib_uniqe=ij_u
c
      end
c===============================================================
      subroutine do_tmp1(isbl_size,isbl_q,ijpres2,ij_temp,iqu )
      dimension isbl_q(*),ijpres2(*), ij_temp(*)
c
         iqu=0
         do 100 iqp=1,isbl_size
         iq=isbl_q(iqp)
         if(iq.gt.0) then
            iqu=iqu+1
            ij_temp(iqu)=ijpres2(iq)
         endif
  100    continue
c
      end
c===============================================================
      subroutine do_uniq(nq, itmp,ij_uniq,nij_uniq)
      implicit none
c
      integer nq                ! [input] Number of quartets
      integer itmp(*)           ! [input] select unique from it
      integer nij_uniq          ! [output] No. of unique ij's
      integer ij_uniq(*)        ! [output] List of unique ij's
c
c     Return the number and list of unique pairs in the current
c     request.The number of such pairs is bounded by min(nq,maxpar).
c
c     Use an NlogN sort followed by linear pass.
c
c     Gather list of IJ's in the request and sort them
c
      integer iq, ijprev
c-----------------------------------------------------
c
cccc  write(6,*) nq,' original pairs:',(itmp(iq),iq=1,nq)
c
c      
c     do iq = 1, nq
c        itmp(iq) = ijpres2(isbl_q(iq))
c     enddo
c
      IF(nq.gt.1) THEN
         call txs_integer_sort(nq, itmp)
c
c        Pass thru list and determine the unique ones
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
      ELSE
         ij_uniq(1)=itmp(1)
         nij_uniq=1
      ENDIF
c
c     Debug
cccc  write(6,*) nij_uniq,' unique   pairs:',(ij_uniq(iq),iq=1,nij_uniq)
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
c===============================================================
