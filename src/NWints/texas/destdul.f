* $Id$
      subroutine destdul(ikbl,nbls,nblok1,ncs,inx,buf,
     *     buffer, icfx,jcfx,kcfx,lcfx, q4, use_q4,
     *     icfg,jcfg,kcfg,lcfg,ngcd,lnijkl,indxp,ipres,iqorder,
     *     map_txs_pnl)
c----------------------------------------------------------------
c     gradient derivatives 
c     
c     This is called for PNL-requested set of contracted shell quartets.
c     Only non-zero Integrals return WITH labels and they do not have 
c     to be in PNL-requested order.
c     
c     buf           - in-comming integrals
c     
c     buffer        - outgoing integrals
c     icfx()-lcfx() - corresponding labels (PNL)
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer map_txs_pnl(*)    ! txs to pnl basis map = ncfunct
      logical use_q4
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /neglect/ eps,eps1,epsr
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
      common /intgop/ ncache,maxprice,iprint,iblock
c----------------------------------------------------------------------
      double precision savezerotol
      common /csavezerotol/ savezerotol ! Used in detbul,set in texas_hf
c----------------------------------------------------------------------
c     
      dimension icfx(*),jcfx(*),kcfx(*),lcfx(*)
      dimension nblok1(2,*)
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension inx(12,*)
c     
cccc  dimension buffer(9,*)
      dimension buffer(12,*)
c     
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension ipres(*), iqorder(*)
      dimension indxp(*)
      dimension q4(*)
      dimension lder(12)        ! to re-order derivativs according to atoms
      dimension iix(4)
c
      double precision xtmp(12)
c     
      double precision threshold ! For screening output integrals
c--------------------------------
c     do not zero out integ_n0 here
c----------------------------
c     loop over quartets belonging to the block IKBL :
c     
c     
      do 10  ijklp=1,nbls
         ijkl=indxp(ijklp)
         if(ijkl.eq.0) go to 10
         iqreq=ipres(ijkl)
         if(iqreq.eq.0) go to 10
         iorder=iqorder(iqreq)
c     test
c     write(6,*)'destDul iorder=',iorder
c     test
         call reorder_der1(iorder,lder)
         if(use_q4) THEN
            symfact=q4(iqreq)
         else
            symfact = 1.0d0
         endif
c     
         threshold = savezerotol/symfact
c     
c---------------------------------------
c     write(6 ,1230)  ijkl,iqreq,iorder
c     1230 format('quart=',i5,' req-quart=,i5,'  iorder=',i4 )
c---------------------------------------
         ijcs=nblok1(1,ijkl)   
         klcs=nblok1(2,ijkl)  
         call get_ij_half(ijcs,ics,jcs)
         call get_ij_half(klcs,kcs,lcs)
         if(ngcd.eq.1) then
            ngcq=1
            icfg(1)=inx(11,ics)
            jcfg(1)=inx(11,jcs)
            kcfg(1)=inx(11,kcs)
            lcfg(1)=inx(11,lcs)
         else
            call indexg(inx,ics,jcs,kcs,lcs,ijcs,klcs,
     *           ilen,jlen,klen,llen, icfg,jcfg,kcfg,lcfg,ngcq)
         endif
c     
         do iqu=1,ngcq
            icff=icfg(iqu)
            jcff=jcfg(iqu)
            kcff=kcfg(iqu)
            lcff=lcfg(iqu)
            icff=map_txs_pnl(icff+1)-1 ! Relies on txs order = pnl order
            jcff=map_txs_pnl(jcff+1)-1
            kcff=map_txs_pnl(kcff+1)-1
            lcff=map_txs_pnl(lcff+1)-1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            integ=0
            do icf=icff+1,icff+ilen
               do jcf=jcff+1,jcff+jlen
                  do kcf=kcff+1,kcff+klen
                     do lcf=lcff+1,lcff+llen
                        integ=integ+1
c------>                   xint0=buf(integ)
                        xtmp( 1)=buf(1,ijklp,integ,iqu) ! xinta
                        xtmp( 2)=buf(4,ijklp,integ,iqu) ! yinta
                        xtmp( 3)=buf(7,ijklp,integ,iqu) ! zinta
                        xtmp( 4)=buf(2,ijklp,integ,iqu) ! xintb
                        xtmp( 5)=buf(5,ijklp,integ,iqu) ! yintb
                        xtmp( 6)=buf(8,ijklp,integ,iqu) ! zintb
                        xtmp( 7)=buf(3,ijklp,integ,iqu) ! xintc
                        xtmp( 8)=buf(6,ijklp,integ,iqu) ! yintc
                        xtmp( 9)=buf(9,ijklp,integ,iqu) ! zintc
                        xnorm = 0.0d0
                        do i = 1, 9
                           xnorm = xnorm + xtmp(i)*xtmp(i)
                        enddo
                        if (xnorm .gt. threshold*threshold) then
                           xtmp(10)=-(xtmp(1)+xtmp(4)+xtmp(7))
                           xtmp(11)=-(xtmp(2)+xtmp(5)+xtmp(8))
                           xtmp(12)=-(xtmp(3)+xtmp(6)+xtmp(9))
                           integ_n0=integ_n0+1
                           do i = 1, 12
                              buffer(lder(i),integ_n0) = xtmp(i)*symfact
                           enddo
                           call lab_req(iorder,icf,jcf,kcf,lcf,iix)
c     
c---------------------------> icfx(integ_n0)=icf
c     jcfx(integ_n0)=jcf
c     kcfx(integ_n0)=kcf
c---------------------------> lcfx(integ_n0)=lcf
                           icfx(integ_n0)=iix(1)
                           jcfx(integ_n0)=iix(2)
                           kcfx(integ_n0)=iix(3)
                           lcfx(integ_n0)=iix(4)
c     
                           if(iprint.ge.2) then
                              call print_der1(ics,jcs,kcs,lcs,inx,
     *                             buf(1,ijklp,integ,iqu),
     *                             icf,jcf,kcf,lcf)
                           endif
                        endif   !   threshold
                     enddo
                  enddo
               enddo
            enddo
         enddo
c     
 10   continue
c--------------------------------------------------------
      end

