      subroutine uniq_pairs(bl,nbl2,npar,
     *                      ijblock,klblock,ijpres2,klpres2, 
     *                      isbl_s,isbl_q)
      implicit real*8 (a-h,o-z)
      common /index/ maxsh,ifp,inx(1)
      common /pnl007/ nij_uniqe, ij_uniqe ! uniqe pairs for pair-blocks
      dimension bl(*)
      dimension npar(nbl2) 
      dimension ijblock(nbl2),klblock(nbl2)
      dimension ijpres2(*),klpres2(*) ! dim=nqrtpnl
      dimension isbl_q(*)  ! dim=nqrtpnl
      dimension isbl_s(*) ! dimension nbl2*(nbl2+1)/2 num of super-blks
c----------------------------------------------------------------------
c for memory allocation : 
c
      ncs_pairs=0
      do 10 ibl=1,nbl2
      ijbl=ijblock(ibl)
      klbl=klblock(ibl)
      if(ijbl.gt.0 .or. klbl.gt.0) ncs_pairs=ncs_pairs+npar(ibl)
   10 continue
c----------------------------------------------------------------------
c memory allocations :
c
c for final list of unique pairs in pair-blocks :
c
      call getmem(ncs_pairs, ibl_uniqe_final)
c----------------------------------------------------------------------
      call getmem(nbl2, nij_uniqe)
      call getmem(nbl2, nkl_uniqe)
      call getmem(ncs_pairs, ij_uniqe) 
      call getmem(ncs_pairs, kl_uniqe) 
c----------------------------------------------------------------------
c zero out number of uniqe pairs in IBL, KBL:
c
      call zero_out_uniq(bl(nij_uniqe),bl(nkl_uniqe),nbl2)
c----------------------------------------------------------------------
c
      ipoint=1
      isbl=0
      ibl_uniqe=ij_uniqe     ! starting address for Ibl pairs
      do 100 ibl=1,nbl2
      ijbl=ijblock(ibl)
      if(ijbl.eq.0) go to 100
c     get number of unique pairs in IBL pair-block :
      call get_nu2(bl(nij_uniqe),ibl,nib_uniqe)
c
         kbl_uniqe=kl_uniqe     ! starting address for Kbl pairs
         do 150 kbl=1,ibl
         klbl=klblock(kbl)
         if(klbl.eq.0) go to 150
c        get number of unique pairs in KBL pair-block :
         call get_nu2(bl(nkl_uniqe),kbl,nkb_uniqe)
         isbl=ibl*(ibl-1)/2+kbl
c
         isbl_size=isbl_s(isbl)
         if(isbl_size.gt.0) then
            call find_uniq_2(bl,isbl_size,isbl_q(ipoint),ijpres2,
     *                       bl(ibl_uniqe) ,    nib_uniqe )
c           set number of unique pairs in IBL pair-block :
            call set_nu2(bl(nij_uniqe),ibl,nib_uniqe)
            call find_uniq_2(bl,isbl_size,isbl_q(ipoint),klpres2,
     *                       bl(kbl_uniqe) ,    nkb_uniqe )
c           set number of unique pairs in KBL pair-block :
            call set_nu2(bl(nkl_uniqe),kbl,nkb_uniqe)
         endif
c
c        on return : lists of pairs  in bl(ibl_uniqe) and bl(kbl_uniqe)
c
         ipoint=ipoint+isbl_size
         kbl_uniqe=kbl_uniqe+npar(kbl)
  150    continue
      ibl_uniqe=ibl_uniqe+npar(ibl)
  100 continue
c----------------------------------------------------------------------
c put together ibl_uniqe and kbl_unique :
c
      ij_total=0
      ibl_u=ij_uniqe     ! starting address for Ibl pairs
      kbl_u=kl_uniqe     ! starting address for Kbl pairs
      ibl_uniqe=ibl_uniqe_final ! starting address for Final IBL pairs
      do 200 ibl=1,nbl2
c
      call get_nu2(bl(nij_uniqe),ibl,n_ibl_u)
      call get_nu2(bl(nkl_uniqe),ibl,n_kbl_u)
c
      call ibl_plus_kbl(    n_ibl_u , bl(ibl_u),    n_kbl_u , bl(kbl_u),
     *                  n_ibl_uniqe, bl(ibl_uniqe),bl )
c
c put final number of unique pairs in IBL pair-block (n_ibl_uniqe) :
c
      call set_nu2(bl(nij_uniqe),ibl,n_ibl_uniqe)
c
      ijbl=ijblock(ibl)
      if(ijbl.gt.0) ibl_u=ibl_u+npar(ibl)
      klbl=klblock(ibl)
      if(klbl.gt.0) kbl_u=kbl_u+npar(ibl)
c
      ibl_uniqe=ibl_uniqe+n_ibl_uniqe
      ij_total=ij_total+n_ibl_uniqe
  200 continue
c
c--------------------------------------------------------
      call retmem(2) ! release ij_uniqe and kl_uniqe
      call retmem(1) ! release nkl_uniqe
c--------------------------------------------------------
c final allocation in BL() for list of unique pairs in pair-blocks :
c
      ij_uniqe=ibl_uniqe_final
c
c final allocation in BL()  for number of unique pairs is 
c
c     nij_uniqe 
c
c These two addresses are transfered outside thru the common 
c block : common /pnl007/ nij_uniqe, ij_uniqe ! unique pairs 
c 2 calls to getmem
c--------------------------------------------------------
cTEST
c
c     write(6,*)'------from uniq_pairs-------------beg     -------'
c
c     ij_total=0
c     ibl_uniqe=ij_uniqe
c     do 300 ibl=1,nbl2
c     call get_nu2(bl(nij_uniqe),ibl,n_ibl_u)
c     call prt_uniqe(ibl,bl(ibl_uniqe),n_ibl_u,npar(ibl))
c
c     ibl_uniqe=ibl_uniqe+n_ibl_u
c     ij_total=ij_total+n_ibl_u
c 300 continue
c     write(6,*) 'Total number of UNIQE pairs in request=',ij_total
c     write(6,*)'------from uniq_pairs-------------end     -------'
c--------------------------------------------------------
c ij_total should be the TOTAL number of uniqe pairs in whole request !
c
      end
c===============================================================
c TEST 
      subroutine prt_uniqe(ibl,ibl_uniqe,n_ibl_uniqe,ibl_par)
      dimension ibl_uniqe(*)
c----------------------------
      write(6,*)
     *' PAIR-BLOCK=',ibl,' UNIQE PAIRS =',n_ibl_uniqe,'NPAR=',ibl_par
c
      write(6,66)(ibl_uniqe(i),i=1,n_ibl_uniqe)
   66 format(20(i3,1x))
      end
c===============================================================
      subroutine find_uniq_2(bl,isbl_size,isbl_q,ijpres2,
     *                       ibl_uniqe ,nib_uniqe)
      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension ijpres2(*) ! dim=nqrtpnl
      dimension isbl_q(*)  ! from ipoint 
      dimension ibl_uniqe(*)
c
         call getmem(isbl_size,i_uniqe)
         call getmem(isbl_size,ij_temp)
c
         call do_tmp1(isbl_size,isbl_q,ijpres2,bl(ij_temp),iqu )
         call do_uniq(iqu,bl(ij_temp),bl(i_uniqe),ij_u)
c
         call retmem(1)  ! release ij_temp
         call getmem(nib_uniqe + ij_u, ij_temp)
c
         call do_tmp2(ij_u,bl(i_uniqe),nib_uniqe,ibl_uniqe,bl(ij_temp))
         call do_uniq(ij_u + nib_uniqe, bl(ij_temp),ibl_uniqe,ij_u)
c
         nib_uniqe=ij_u
c
         call retmem(2) ! release ij_temp and i_uniqe
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
      subroutine do_tmp2(ij_u,i_uniqe, nij_u,ibl_uniqe,ij_temp)
      dimension i_uniqe(*), ibl_uniqe(*), ij_temp(*)
c
      do 100 ii=1,ij_u
      ij_temp(ii)=i_uniqe(ii)
  100 continue
c
      do 200 ii=1,nij_u
      jj=ii+ij_u
      ij_temp(jj)=ibl_uniqe(ii)
  200 continue
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
c     request.The number of such pairs is bounded by min(nq,npar(ijbl)).
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
      subroutine ibl_plus_kbl(n_ibl_u,ibl_u, n_kbl_u,kbl_u,
     *                        n_ibl_uniqe,ibl_uniqe, bl)
      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension ibl_u(*),kbl_u(*), ibl_uniqe(*)
c-----------------------------------------------
c Input :
c     n_ibl_u= number of uniqe ij in IBL
c     n_kbl_u= number of uniqe kl in KBL
c     ibl_u()- list of unique pairs in IBL
c     kbl_u()- list of unique pairs in KBL
c OUTPUT:
c     n_ibl_uniqe = final number of unique pairs in IBL
c     ibl_uniqe   - final list of unique pairs in IBL
c-----------------------------------------------
c
c     write(6,*)'from ibl_plus_kbl:  ij_u=',n_ibl_u,'  kl_u=',n_kbl_u
c
      if(n_ibl_u.ne.0 .and. n_kbl_u.ne.0) then
         call getmem(n_ibl_u + n_kbl_u, ik_tmp)
         call do_tmp2(n_ibl_u,ibl_u, n_kbl_u,kbl_u,  bl(ik_tmp) )
         call do_uniq(n_ibl_u + n_kbl_u, bl(ik_tmp),kbl_u,  n_ibl_uniqe)
         call retmem(1) ! release ik_tmp
c
         call tfer_i(kbl_u,ibl_uniqe,n_ibl_uniqe)
      endif
      if(n_ibl_u.ne.0 .and. n_kbl_u.eq.0) then
         n_ibl_uniqe=n_ibl_u
         call tfer_i(ibl_u,ibl_uniqe,n_ibl_uniqe)
      endif
      if(n_ibl_u.eq.0 .and. n_kbl_u.ne.0) then
         n_ibl_uniqe=n_kbl_u
         call tfer_i(kbl_u,ibl_uniqe,n_ibl_uniqe)
      endif
      if(n_ibl_u.eq.0 .and. n_kbl_u.eq.0) then
         n_ibl_uniqe=0
      endif
c
      end
c===============================================================
      subroutine get_uniq(ibl,nij_u,ij_u_point,  bl)
      implicit real*8 (a-h,o-z)
      common /pnl007/ nij_uniqe, ij_uniqe ! uniqe pairs for pair-blocks
      dimension bl(*)
c---------------------------------------------------------
c OUTPUT :
c nij_u      - number of unique pairs in IBL
c ij_u_point - pointer to the list of unique pairs (in bl)
c---------------------------------------------------------
      call get_nu2(bl(nij_uniqe),ibl,nij_u)
c
      ipoint=0
      do 10 iblock=1,ibl-1
      call get_nu2(bl(nij_uniqe),iblock,n_iblock_u)
      ipoint=ipoint+    n_iblock_u
   10 continue
      ij_u_point=ipoint + ij_uniqe
c
      end
c===============================================================
      subroutine zero_out_uniq(nij_uniqe,nkl_uniqe,nbl2)
      dimension nij_uniqe(*),nkl_uniqe(*)
      do 20 ibl=1,nbl2
      nij_uniqe(ibl)=0
      nkl_uniqe(ibl)=0
   20 continue
      end
c===============================================================
      subroutine get_nu2(nij_uniqe,ibl,nij_u)
c get number of unique pairs in ibl pair-block :
      dimension nij_uniqe(*)
      nij_u= nij_uniqe(ibl)
      end
c===============================================================
      subroutine set_nu2(nij_uniqe,ibl,nij_u)
c set number of unique pairs in ibl pair-block :
      dimension nij_uniqe(*)
      nij_uniqe(ibl)=nij_u
      end
c===============================================================
