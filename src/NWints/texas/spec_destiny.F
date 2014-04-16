* $Id$
c==============================================================
c K.Wolinski : May,1997 : the use of 
c       ics=iis(ijcs)    
c       jcs=jjs(ijcs)   
c       kcs=iis(klcs)  
c       lcs=jjs(klcs) 
c has been eliminated. 
c Now, these shell numbers are generated on the fly
c       call get_ij_half(ijcs,ics,jcs)
c       call get_ij_half(klcs,kcs,lcs)
c==============================================================
      subroutine desthuf(ikbl,nbls,nblok1,ncs,inx,buf,
     *                   buffer, itxspnl, q4,use_q4,
     *                   icfg,jcfg,kcfg,lcfg,ngcd,lnijkl,
     *                   indxp,ipres,iqorder)
c----------------------------------------------------------------
c ONLY hessian  derivatives 
c
c This is called for PNL-requested ONE contracted shell quartet.
c All Integrals (including zeros) return WITHOUT labels but they 
c have to be in PNL-requested order.
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      logical use_q4
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /neglect/ eps,eps1,epsr
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
      common /intgop/ ncache,maxprice,iprint,iblock
c
      dimension nblok1(2,*)
      dimension buf(45,nbls,lnijkl,ngcd)
      dimension inx(12,*)
c
      dimension buffer(78,*)
cccc  dimension buffer(lnijkl*ngcd,78)    ! output buffer of sec.der.
c can not be like this for spherical harmonics !
c
      dimension itxspnl(*)
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension indxp(*),ipres(*)
      dimension q4(*)
      dimension iqorder(*)
      dimension lder(78) ! derivative's order
      dimension worker(78) ! local working scratch
c--------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c--------------------------------------------------------
      IF(use_q4) THEN
c-----symmetry is used---------------------------
c
        integral=0
        do 10  ijklp=1,nbls
        ijkl=indxp(ijklp)
        if(ijkl.eq.0) go to 10
        iqreq=ipres(ijkl)
        if(iqreq.eq.0) go to 10
        iorder=iqorder(iqreq)
c--
        call reorder_der2(iorder,lder)
c--
        symfact=q4(iqreq)
c
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
     *                 ilen,jlen,klen,llen, icfg,jcfg,kcfg,lcfg,ngcq)
        endif
c
          do 15  iqu=1,ngcq
          icff=icfg(iqu)
          jcff=jcfg(iqu)
          kcff=kcfg(iqu)
          lcff=lcfg(iqu)
c
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 20  iii=1,ilen
             icf=icff+iii
             do 20  jjj=1,jlen
             jcf=jcff+jjj
             do 20  kkk=1,klen
             kcf=kcff+kkk
             do 20  lll=1,llen
             lcf=lcff+lll
c---
             integ=integ+1
             integral=integral+1
c---
             ipnl=itxspnl(integral)
c---
c------------------------------------------------------------
c construct all 10 blocks of sec.der. (output) from 6 blocks:
c
c          AA AB AC AD                AA AB AC 
c             BB BC BD      from         BB BC
c                CC CD                      CC
c                   DD  
c      1-6, 7-15,16-24,25-33         1-6, 7-15,16-24
c          34-39,40-48,49-57             25-30,31-39
c                58-63,64-72                   40-45
c                      73-78
c------------------------------------------------------------
c
             call make_78_from_45(buf(1,ijklp,integ,iqu),worker)
c
c------------------------------------------------------------
             do m=1,78
ccccc          buffer(m,ipnl)=buf(m,ijklp,integ,iqu)*symfact
               buffer(lder(m),ipnl)=worker(m)*symfact
             enddo
c------------------------------------------------------------
               if(iprint.ge.2) then
                  call print_der2(ics,jcs,kcs,lcs,inx,
     *                            buf(1,ijklp,integ,iqu),
     *                            icf,jcf,kcf,lcf)
               endif
c---
  20         continue
  15      continue
  10    continue
c--------------------------------------------------------
      ELSE
c-----symmetry is not used---------------------------
c
        integral=0
        do 100 ijklp=1,nbls
        ijkl=indxp(ijklp)
        if(ijkl.eq.0) go to 100
        iqreq=ipres(ijkl)
        if(iqreq.eq.0) go to 100
        iorder=iqorder(iqreq)
c--
        call reorder_der2(iorder,lder)
c--
c
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
     *                 ilen,jlen,klen,llen, icfg,jcfg,kcfg,lcfg,ngcq)
        endif
c
          do 150 iqu=1,ngcq
          icff=icfg(iqu)
          jcff=jcfg(iqu)
          kcff=kcfg(iqu)
          lcff=lcfg(iqu)
c
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 200 iii=1,ilen
             icf=icff+iii
             do 200 jjj=1,jlen
             jcf=jcff+jjj
             do 200 kkk=1,klen
             kcf=kcff+kkk
             do 200 lll=1,llen
             lcf=lcff+lll
c---
             integ=integ+1
             integral=integral+1
c---
             ipnl=itxspnl(integral)
             call make_78_from_45(buf(1,ijklp,integ,iqu),worker)
c---
             do m=1,78
cno no no      buffer(ipnl,lder(m))=worker(m)
               buffer(lder(m),ipnl)=worker(m)
             enddo
c---
             if(iprint.ge.2) then
                call print_der2(ics,jcs,kcs,lcs,inx,
     *                          buf(1,ijklp,integ,iqu),
     *                          icf,jcf,kcf,lcf)
                call print_pnl2(ics,jcs,kcs,lcs,inx,
     *                          worker,
     *                          icf,jcf,kcf,lcf)
             endif
c---
  200        continue
  150     continue
  100   continue
c
c--------------------------------------------------------
      ENDIF
c--------------------------------------------------------
      integ_n0=integ_n0+integral
c--------------------------------------------------------
      end
c==============================================================
c   #if defined(IBM)
c   *IBM COMPILER OPTIONS JUST FOR DESTBUL
c   @PROCESS OPT(2)
c   #endif
      subroutine desthul(ikbl,nbls,nblok1,ncs,inx,buf,
     *     buffer, icfx,jcfx,kcfx,lcfx, q4, use_q4,
     *     icfg,jcfg,kcfg,lcfg,ngcd,lnijkl,indxp,ipres,iqorder,
     *     map_txs_pnl)
c----------------------------------------------------------------
c ONLY hessian  derivatives 
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
      integer map_txs_pnl(*)        ! txs to pnl basis map = ncfunct
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
      dimension buf(45,nbls,lnijkl,ngcd)
      dimension inx(12,*)
c     
ccc   dimension buffer(45,*)
      dimension buffer(78,*)
c     
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension ipres(*), iqorder(*)
      dimension indxp(*)
      dimension q4(*)
      dimension lder(78) ! derivative's order
      dimension iix(4)
      dimension worker(78)
c
      double precision threshold ! For screening output integrals
c----------------------------
c     do not zero out integ_n0 here
c----------------------------
c     loop over quartets belonging to the block IKBL :
c     
c     
      do 10  ijklp=1,nbls
         ijkl=indxp(ijklp)
c     
         if(ijkl.eq.0) go to 10
         iqreq=ipres(ijkl)
c     
         if(iqreq.eq.0) go to 10
         iorder=iqorder(iqreq)
         call reorder_der2(iorder,lder)
c
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
                              absmax=0.d0
                              do m=1,45
                                der2=buf(m,ijklp,integ,iqu)
                                absder2=abs(der2)
                                if(absder2.gt.absmax) absmax=absder2
                              enddo
                           if(absmax.gt. threshold ) then
                              integ_n0=integ_n0+1
c
             call make_78_from_45(buf(1,ijklp,integ,iqu),worker)
             do m=1,78
               buffer(lder(m),integ_n0)=worker(m)*symfact
             enddo
c
                           call lab_req(iorder,icf,jcf,kcf,lcf,iix)
c                             icfx(integ_n0)=icf
c                             jcfx(integ_n0)=jcf
c                             kcfx(integ_n0)=kcf
c                             lcfx(integ_n0)=lcf
                              icfx(integ_n0)=iix(1)
                              jcfx(integ_n0)=iix(2)
                              kcfx(integ_n0)=iix(3)
                              lcfx(integ_n0)=iix(4)
                              if(iprint.ge.2) then
                                 call print_der2(ics,jcs,kcs,lcs,inx,
     *                                           buf(1,ijklp,integ,iqu),
     *                                           icf,jcf,kcf,lcf)
                                 call print_pnl2(ics,jcs,kcs,lcs,inx,
     *                                           worker,
     *                                           icf,jcf,kcf,lcf)
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
c     
 10      continue
c--------------------------------------------------------
         end
c==============================================================
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
      subroutine print_der2(ics,jcs,kcs,lcs,inx,
     *                      der2, icf,jcf,kcf,lcf)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension iix(4)
      dimension der2(45)
c$$$      dimension atom_pairs(10,3 , 10,3)
c
         iat=inx(2,ics)
         jat=inx(2,jcs)
         kat=inx(2,kcs)
         lat=inx(2,lcs)
c
         if(iat.gt.10 .or. jat.gt.10) then
            write(6,*) 'too many atoms; derivatives can not be printed'
            return
         endif
         if(kat.gt.10 .or. lat.gt.10) then
            write(6,*) 'too many atoms; derivatives can not be printed'
            return
         endif
c
c-----test print only : begining -----------------------------
c
c block aa:
       axax=der2(1)
       axay=der2(2)
       axaz=der2(3)
       ayay=der2(4)
       ayaz=der2(5)
       azaz=der2(6)
c
c block ab:
       axbx=der2(7)
       axby=der2(8)
       axbz=der2(9)
       aybx=der2(10)
       ayby=der2(11)
       aybz=der2(12)
       azbx=der2(13)
       azby=der2(14)
       azbz=der2(15)
c block ac:
       axcx=der2(16)
       axcy=der2(17)
       axcz=der2(18)
       aycx=der2(19)
       aycy=der2(20)
       aycz=der2(21)
       azcx=der2(22)
       azcy=der2(23)
       azcz=der2(24)
c block bb:
       bxbx=der2(25)
       bxby=der2(26)
       bxbz=der2(27)
       byby=der2(28)
       bybz=der2(29)
       bzbz=der2(30)
c block bc:
       bxcx=der2(31)
       bxcy=der2(32)
       bxcz=der2(33)
       bycx=der2(34)
       bycy=der2(35)
       bycz=der2(36)
       bzcx=der2(37)
       bzcy=der2(38)
       bzcz=der2(39)
c block cc:
       cxcx=der2(40)
       cxcy=der2(41)
       cxcz=der2(42)
       cycy=der2(43)
       cycz=der2(44)
       czcz=der2(45)
c block ad: from transl. inv.
       axdx=-(axax+axbx+axcx)
       axdy=-(axay+axby+axcy)
       axdz=-(axaz+axbz+axcz)
       aydx=-(aXaY+aybx+aycx)
       aydy=-(ayay+ayby+aycy)
       aydz=-(ayaz+aybz+aycz)
       azdx=-(aXaZ+azbx+azcx)
       azdy=-(aYaZ+azby+azcy)
       azdz=-(azaz+azbz+azcz)
c block bd: from transl. inv.
       bxdx=-(AxBx+bxbx+bxcx)
       bxdy=-(AYBX+bxby+bxcy)
       bxdz=-(AZBX+bxbz+bxcz)
       bydx=-(AXBY+bXbY+bycx)
       bydy=-(AYBY+byby+bycy)
       bydz=-(AZBY+bybz+bycz)
       bzdx=-(AXBZ+bXbZ+bzcx)
       bzdy=-(AyBz+bybz+bzcy)
       bzdz=-(AZBZ+bzbz+bzcz)
c block cd: from transl. inv.
       cxdx=-(AXCX+BXCX+cxcx)
       cxdy=-(AYCX+BYCX+cxcy)
       cxdz=-(AZCX+BZCX+cxcz)
       cydx=-(AXCY+BXCY+cXcY)
       cydy=-(AYCY+BYCY+cycy)
       cydz=-(AZCY+BZCY+cycz)
       czdx=-(AXCZ+BxCZ+cXcZ)
       czdy=-(AYCZ+BYCZ+cYcZ)
       czdz=-(AZCZ+BZCZ+czcz)
c block dd:
       dxdx=-(AXDX+BXDX+CXDX)
cccc   dxdy=-(dxay+dxby+dxcy)
       dxdy=-(AYDX+BYDX+CYDX)
cccc   dxdz=-(dxaz+dxbz+dxcz)
       dxdz=-(AZDX+BZDX+CZDX)
cccc   dydy=-(dyay+dyby+dycy)
       dydy=-(AYDY+BYDY+CYDY)
ccccc  dydz=-(dyaz+dybz+dycz)
       dydz=-(AZDY+BZDY+CZDY)
cccc   dzdz=-(dzaz+dzbz+dzcz)
       dzdz=-(AZDZ+BZDZ+CZDZ)
c--------------------------------------------------------------------
canonical order :
c
                   call descend(icf,jcf,kcf,lcf,iix)
c--------------------------------------------------------------------
c block AA:
                    a_min=min( axax,axay,axaz, ayay,ayaz,azaz )
                    a_max=max( axax,axay,axaz, ayay,ayaz,azaz )
                    absmax=max( abs(a_min),abs(a_max) )
                    if(absmax.gt.1.d-9) then
c     write(6,60) ics,jcs,kcs,lcs, icf,jcf,kcf,lcf,iat,jat,kat,lat
c 60  format('shells=',4i2,' functions=',4i2,' centers=',4i2)
cccc                   write(6,61) iix(1),iix(2),iix(3),iix(4),
                       write(6,61) icf,jcf,kcf,lcf, axax,axay,axaz
                       write(6,61) icf,jcf,kcf,lcf, ayay,ayaz,azaz
                    endif
  61  format('d2/dAidAj: ijkl=',4i2,1x,3(f12.4,2x))
c
c block BB:
                    b_min=min( bxbx,bxby,bxbz, byby,bybz,bzbz )
                    b_max=max( bxbx,bxby,bxbz, byby,bybz,bzbz )
                    absmax=max( abs(b_min),abs(b_max) )
                    if(absmax.gt.1.d-9) then
ccc   write(6,60) ics,jcs,kcs,lcs, icf,jcf,kcf,lcf,iat,jat,kat,lat
cccc                   write(6,62) iix(1),iix(2),iix(3),iix(4),
                       write(6,62) icf,jcf,kcf,lcf, bxbx,bxby,bxbz
                       write(6,62) icf,jcf,kcf,lcf, byby,bybz,bzbz
                    endif
  62  format('d2/dBidBj: ijkl=',4i2,1x,3(f12.4,2x))
c
c block CC:
                    c_min=min( cxcx,cxcy,cxcz, cycy,cycz,czcz )
                    c_max=max( cxcx,cxcy,cxcz, cycy,cycz,czcz )
                    absmax=max( abs(c_min),abs(c_max) )
                    if(absmax.gt.1.d-9) then
ccc   write(6,60) ics,jcs,kcs,lcs, icf,jcf,kcf,lcf,iat,jat,kat,lat
cccc                   write(6,63) iix(1),iix(2),iix(3),iix(4),
                       write(6,63) icf,jcf,kcf,lcf, cxcx,cxcy,cxcz
                       write(6,63) icf,jcf,kcf,lcf, cycy,cycz,czcz
                    endif
  63  format('d2/dCidCj: ijkl=',4i2,1x,3(f12.4,2x))
c block AB:
        ab_min=min(axbx,axby,axbz,aybx,ayby,aybz,azbx,azby,azbz)
        ab_max=max(axbx,axby,axbz,aybx,ayby,aybz,azbx,azby,azbz)
           absmax=max( abs(ab_min),abs(ab_max) )
           if(absmax.gt.1.d-9) then
ccccc         write(6,64) iix(1),iix(2),iix(3),iix(4),
              write(6,64) icf,jcf,kcf,lcf, axbx,axby,axbz
              write(6,64) icf,jcf,kcf,lcf, aybx,ayby,aybz
              write(6,64) icf,jcf,kcf,lcf, azbx,azby,azbz
  64  format('d2/dAidBj: ijkl=',4i2,1x,3(f12.4,2x))
           endif
c block AC:
        ac_min=min(axcx,axcy,axcz,aycx,aycy,aycz,azcx,azcy,azcz)
        ac_max=max(axcx,axcy,axcz,aycx,aycy,aycz,azcx,azcy,azcz)
           absmax=max( abs(ac_min),abs(ac_max) )
           if(absmax.gt.1.d-9) then
              write(6,65) icf,jcf,kcf,lcf, axcx,axcy,axcz
              write(6,65) icf,jcf,kcf,lcf, aycx,aycy,aycz
              write(6,65) icf,jcf,kcf,lcf, azcx,azcy,azcz
  65  format('d2/dAidCj: ijkl=',4i2,1x,3(f12.4,2x))
           endif
c block AD:
        ad_min=min(axdx,axdy,axdz,aydx,aydy,aydz,azdx,azdy,azdz)
        ad_max=max(axdx,axdy,axdz,aydx,aydy,aydz,azdx,azdy,azdz)
           absmax=max( abs(ad_min),abs(ad_max) )
           if(absmax.gt.1.d-9) then
ccc           write(6,66) iix(1),iix(2),iix(3),iix(4),
              write(6,66) icf,jcf,kcf,lcf, axdx,axdy,axdz
              write(6,66) icf,jcf,kcf,lcf, aydx,aydy,aydz
              write(6,66) icf,jcf,kcf,lcf, azdx,azdy,azdz
  66  format('d2/dAidDj: ijkl=',4i2,1x,3(f12.4,2x))
           endif
c block BC:
        bc_min=min(bxcx,bxcy,bxcz,bycx,bycy,bycz,bzcx,bzcy,bzcz)
        bc_max=max(bxcx,bxcy,bxcz,bycx,bycy,bycz,bzcx,bzcy,bzcz)
           absmax=max( abs(bc_min),abs(bc_max) )
           if(absmax.gt.1.d-9) then
ccc           write(6,67) iix(1),iix(2),iix(3),iix(4),
              write(6,67) icf,jcf,kcf,lcf, bxcx,bxcy,bxcz
              write(6,67) icf,jcf,kcf,lcf, bycx,bycy,bycz
              write(6,67) icf,jcf,kcf,lcf, bzcx,bzcy,bzcz
  67  format('d2/dBidCj: ijkl=',4i2,1x,3(f12.4,2x))
           endif
c block BD:
        bd_min=min(bxdx,bxdy,bxdz,bydx,bydy,bydz,bzdx,bzdy,bzdz)
        bd_max=max(bxdx,bxdy,bxdz,bydx,bydy,bydz,bzdx,bzdy,bzdz)
           absmax=max( abs(bd_min),abs(bd_max) )
           if(absmax.gt.1.d-9) then
cccc          write(6,68) iix(1),iix(2),iix(3),iix(4),
              write(6,68) icf,jcf,kcf,lcf, bxdx,bxdy,bxdz
              write(6,68) icf,jcf,kcf,lcf, bydx,bydy,bydz
              write(6,68) icf,jcf,kcf,lcf, bzdx,bzdy,bzdz
  68  format('d2/dBidDj: ijkl=',4i2,1x,3(f12.4,2x))
           endif
c block CD:
        cd_min=min(cxdx,cxdy,cxdz,cydx,cydy,cydz,czdx,czdy,czdz)
        cd_max=max(cxdx,cxdy,cxdz,cydx,cydy,cydz,czdx,czdy,czdz)
           absmax=max( abs(cd_min),abs(cd_max) )
           if(absmax.gt.1.d-9) then
ccc           write(6,69) iix(1),iix(2),iix(3),iix(4),
              write(6,69) icf,jcf,kcf,lcf, cxdx,cxdy,cxdz
              write(6,69) icf,jcf,kcf,lcf, cydx,cydy,cydz
              write(6,69) icf,jcf,kcf,lcf, czdx,czdy,czdz
  69  format('d2/dCidDj: ijkl=',4i2,1x,3(f12.4,2x))
           endif
c blodk DD:
                    d_min=min( dxdx,dxdy,dxdz, dydy,dydz,dzdz )
                    d_max=max( dxdx,dxdy,dxdz, dydy,dydz,dzdz )
                    absmax=max( abs(d_min),abs(d_max) )
                    if(absmax.gt.1.d-9) then
ccc   write(6,60) ics,jcs,kcs,lcs, icf,jcf,kcf,lcf,iat,jat,kat,lat
cccc                   write(6,70) iix(1),iix(2),iix(3),iix(4),
                       write(6,70) icf,jcf,kcf,lcf, dxdx,dxdy,dxdz
                       write(6,70) icf,jcf,kcf,lcf, dydy,dydz,dzdz
                    endif
  70  format('d2/dDidDj: ijkl=',4i2,1x,3(f12.4,2x))
c--------------------------------------------------------------------
          RETURN
c--------------------------------------------------------------------
c
c$$$       do 1234 iii=1,10
c$$$       do 1234 i=1,3
c$$$       do 1234 jjj=1,10
c$$$       do 1234 j=1,3
c$$$       atom_pairs(iii,i ,jjj,j)=0.d0
c$$$ 1234  continue
c$$$c
c$$$c atoms iat,iat
c$$$       atom_pairs(iat,1,iat,1)=atom_pairs(iat,1,iat,1)+axax
c$$$       atom_pairs(iat,1,iat,2)=atom_pairs(iat,1,iat,2)+axay
c$$$       atom_pairs(iat,1,iat,3)=atom_pairs(iat,1,iat,3)+axaz
c$$$c
c$$$       atom_pairs(iat,2,iat,1)=atom_pairs(iat,2,iat,1)+axay
c$$$       atom_pairs(iat,2,iat,2)=atom_pairs(iat,2,iat,2)+ayay
c$$$       atom_pairs(iat,2,iat,3)=atom_pairs(iat,2,iat,3)+ayaz
c$$$c
c$$$       atom_pairs(iat,3,iat,1)=atom_pairs(iat,3,iat,1)+axaz
c$$$       atom_pairs(iat,3,iat,2)=atom_pairs(iat,3,iat,2)+ayaz
c$$$       atom_pairs(iat,3,iat,3)=atom_pairs(iat,3,iat,3)+azaz
c$$$c atoms jat,jat
c$$$       atom_pairs(jat,1,jat,1)=atom_pairs(jat,1,jat,1)+bxbx
c$$$       atom_pairs(jat,1,jat,2)=atom_pairs(jat,1,jat,2)+bxby
c$$$       atom_pairs(jat,1,jat,3)=atom_pairs(jat,1,jat,3)+bxbz
c$$$c
c$$$       atom_pairs(jat,2,jat,1)=atom_pairs(jat,2,jat,1)+bxby
c$$$       atom_pairs(jat,2,jat,2)=atom_pairs(jat,2,jat,2)+byby
c$$$       atom_pairs(jat,2,jat,3)=atom_pairs(jat,2,jat,3)+bybz
c$$$c
c$$$       atom_pairs(jat,3,jat,1)=atom_pairs(jat,3,jat,1)+bxbz
c$$$       atom_pairs(jat,3,jat,2)=atom_pairs(jat,3,jat,2)+bybz
c$$$       atom_pairs(jat,3,jat,3)=atom_pairs(jat,3,jat,3)+bzbz
c$$$c atoms kat,kat
c$$$       atom_pairs(kat,1,kat,1)=atom_pairs(kat,1,kat,1)+cxcx
c$$$       atom_pairs(kat,1,kat,2)=atom_pairs(kat,1,kat,2)+cxcy
c$$$       atom_pairs(kat,1,kat,3)=atom_pairs(kat,1,kat,3)+cxcz
c$$$c
c$$$       atom_pairs(kat,2,kat,1)=atom_pairs(kat,2,kat,1)+cxcy
c$$$       atom_pairs(kat,2,kat,2)=atom_pairs(kat,2,kat,2)+cycy
c$$$       atom_pairs(kat,2,kat,3)=atom_pairs(kat,2,kat,3)+cycz
c$$$c
c$$$       atom_pairs(kat,3,kat,1)=atom_pairs(kat,3,kat,1)+cxcz
c$$$       atom_pairs(kat,3,kat,2)=atom_pairs(kat,3,kat,2)+cycz
c$$$       atom_pairs(kat,3,kat,3)=atom_pairs(kat,3,kat,3)+czcz
c$$$c atoms lat,lat
c$$$       atom_pairs(lat,1,lat,1)=atom_pairs(lat,1,lat,1)+dxdx
c$$$       atom_pairs(lat,1,lat,2)=atom_pairs(lat,1,lat,2)+dxdy
c$$$       atom_pairs(lat,1,lat,3)=atom_pairs(lat,1,lat,3)+dxdz
c$$$c
c$$$       atom_pairs(lat,2,lat,1)=atom_pairs(lat,2,lat,1)+dxdy
c$$$       atom_pairs(lat,2,lat,2)=atom_pairs(lat,2,lat,2)+dydy
c$$$       atom_pairs(lat,2,lat,3)=atom_pairs(lat,2,lat,3)+dydz
c$$$c
c$$$       atom_pairs(lat,3,lat,1)=atom_pairs(lat,3,lat,1)+dxdz
c$$$       atom_pairs(lat,3,lat,2)=atom_pairs(lat,3,lat,2)+dydz
c$$$       atom_pairs(lat,3,lat,3)=atom_pairs(lat,3,lat,3)+dzdz
c$$$c
c$$$c atoms iat,jat:
c$$$       atom_pairs(iat,1,jat,1)=atom_pairs(iat,1,jat,1)+axbx
c$$$       atom_pairs(iat,1,jat,2)=atom_pairs(iat,1,jat,2)+axby
c$$$       atom_pairs(iat,1,jat,3)=atom_pairs(iat,1,jat,3)+axbz
c$$$c
c$$$       atom_pairs(iat,2,jat,1)=atom_pairs(iat,2,jat,1)+aybx
c$$$       atom_pairs(iat,2,jat,2)=atom_pairs(iat,2,jat,2)+ayby
c$$$       atom_pairs(iat,2,jat,3)=atom_pairs(iat,2,jat,3)+aybz
c$$$c
c$$$       atom_pairs(iat,3,jat,1)=atom_pairs(iat,3,jat,1)+azbx
c$$$       atom_pairs(iat,3,jat,2)=atom_pairs(iat,3,jat,2)+azby
c$$$       atom_pairs(iat,3,jat,3)=atom_pairs(iat,3,jat,3)+azbz
c$$$                   iatom=iat
c$$$cccc               write(6,*)' atoms: iat,iat=',iat,iat
c$$$                   do icart=1,3
c$$$                   do jcart=1,3
c$$$                      deriv=atom_pairs(iatom,icart,iatom,jcart)
c$$$                      if(abs(deriv).gt.1d-9) then
c$$$                              write(6,111) 
c$$$     $                             iatom,icart, iatom,jcart,
c$$$     $                             iix(1),iix(2),iix(3),iix(4),
c$$$     $                             deriv
c$$$                      endif
c$$$                   enddo
c$$$                   enddo
c$$$c------------------
c$$$c               if(jat.ne.iat) then
c$$$c                  iatom=jat
c$$$c                  write(6,*)' atoms: jat,jat=',jat,jat
c$$$c                  do icart=1,3
c$$$c                  do jcart=1,3
c$$$c                     deriv=atom_pairs(iatom,icart,iatom,jcart)
c$$$c                     if(abs(deriv).gt.1d-9) then
c$$$c                             write(6,111) 
c$$$c    $                             iatom,icart, iatom,jcart,
c$$$c    $                             iix(1),iix(2),iix(3),iix(4),
c$$$c    $                             deriv
c$$$c                     endif
c$$$c                  enddo
c$$$c                  enddo
c$$$c               endif
c$$$c------------------
c$$$c               if(kat.ne.iat .and. kat.ne.jat) then
c$$$c                  iatom=kat
c$$$c                  write(6,*)' atoms: kat,kat=',kat,kat
c$$$c                  do icart=1,3
c$$$c                  do jcart=1,3
c$$$c                     deriv=atom_pairs(iatom,icart,iatom,jcart)
c$$$c                     if(abs(deriv).gt.1d-9) then
c$$$c                             write(6,111) 
c$$$c    $                             iatom,icart, iatom,jcart,
c$$$c    $                             iix(1),iix(2),iix(3),iix(4),
c$$$c    $                             deriv
c$$$c                     endif
c$$$c                  enddo
c$$$c                  enddo
c$$$c               endif
c$$$c------------------
c$$$c               if(lat.ne.iat .and. lat.ne.jat .and. lat.ne.iat) then
c$$$c                  iatom=lat
c$$$c                  write(6,*)' atoms: lat,lat=',lat,lat
c$$$c                  do icart=1,3
c$$$c                  do jcart=1,3
c$$$c                     deriv=atom_pairs(iatom,icart,iatom,jcart)
c$$$c                     if(abs(deriv).gt.1d-9) then
c$$$c                             write(6,111) 
c$$$c    $                             iatom,icart, iatom,jcart,
c$$$c    $                             iix(1),iix(2),iix(3),iix(4),
c$$$c    $                             deriv
c$$$c                     endif
c$$$c                  enddo
c$$$c                  enddo
c$$$c               endif
c$$$c------------------
c$$$ 111  format('atom=',i2,', xyz=',i2,' atom=',i2,' xyz=',i2,
c$$$     $        ' ijkl=',4i3,', d2(ij|kl)/dpdq=',f12.7   )
c
c-----test print only : end      -----------------------------
c
      end
c==============================================================
      subroutine destduf(ikbl,nbls,nblok1,ncs,inx,buf,
     *                   buffer, itxspnl, q4,use_q4,
     *                   icfg,jcfg,kcfg,lcfg,ngcd,lnijkl,
     *                   indxp,ipres,iqorder)
c----------------------------------------------------------------
c gradient derivatives 
c
c This is called for PNL-requested ONE contracted shell quartet.
c All Integrals (including zeros) return WITHOUT labels but they 
c have to be in PNL-requested order.
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      logical use_q4
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /neglect/ eps,eps1,epsr
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
      common /intgop/ ncache,maxprice,iprint,iblock
c
      dimension nblok1(2,*)
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension inx(12,*)
c
      dimension buffer(12,*)
cnono dimension buffer(lnijkl*ngcd,12) !for one quartet at the time ONLY
c can not be like this for spherical harmonics !
c
      dimension itxspnl(*)
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension indxp(*),ipres(*)
      dimension q4(*)
      dimension iqorder(*)
      dimension lder(12) ! derivative's order
c--------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c--------------------------------------------------------
      IF(use_q4) THEN
c-----symmetry is used---------------------------
        integral=0
        do 10  ijklp=1,nbls
        ijkl=indxp(ijklp)
        if(ijkl.eq.0) go to 10
        iqreq=ipres(ijkl)
        if(iqreq.eq.0) go to 10
        iorder=iqorder(iqreq)
ctest
c        write(6,*)'destDuf iorder=',iorder
ctest
c--
        call reorder_der1(iorder,lder)
c--
        symfact=q4(iqreq)
c
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
     *                 ilen,jlen,klen,llen, icfg,jcfg,kcfg,lcfg,ngcq)
        endif
c
          do 15  iqu=1,ngcq
          icff=icfg(iqu)
          jcff=jcfg(iqu)
          kcff=kcfg(iqu)
          lcff=lcfg(iqu)
c
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 20  iii=1,ilen
             icf=icff+iii
             do 20  jjj=1,jlen
             jcf=jcff+jjj
             do 20  kkk=1,klen
             kcf=kcff+kkk
             do 20  lll=1,llen
             lcf=lcff+lll
c---
             integ=integ+1
             integral=integral+1
ccccc        xint0=buf(ijklp,integ,iqu)
c--
             xinta=buf(1,ijklp,integ,iqu)
             xintb=buf(2,ijklp,integ,iqu)
             xintc=buf(3,ijklp,integ,iqu)
             xintd=-(xinta+xintb+xintc) ! trans. inv.
c
             yinta=buf(4,ijklp,integ,iqu)
             yintb=buf(5,ijklp,integ,iqu)
             yintc=buf(6,ijklp,integ,iqu)
             yintd=-(yinta+yintb+yintc) ! trans. inv.
c
             zinta=buf(7,ijklp,integ,iqu)
             zintb=buf(8,ijklp,integ,iqu)
             zintc=buf(9,ijklp,integ,iqu)
             zintd=-(zinta+zintb+zintc) ! trans. inv.
c---
             ipnl=itxspnl(integral)
cccccc       buffer(ipnl    )=xint0*symfact
c--
             buffer(lder(1),ipnl)=xinta*symfact
             buffer(lder(2),ipnl)=yinta*symfact
             buffer(lder(3),ipnl)=zinta*symfact
c--
             buffer(lder(4),ipnl)=xintb*symfact
             buffer(lder(5),ipnl)=yintb*symfact
             buffer(lder(6),ipnl)=zintb*symfact
c--
             buffer(lder(7),ipnl)=xintc*symfact
             buffer(lder(8),ipnl)=yintc*symfact
             buffer(lder(9),ipnl)=zintc*symfact
c--
             buffer(lder(10),ipnl)=xintd*symfact
             buffer(lder(11),ipnl)=yintd*symfact
             buffer(lder(12),ipnl)=zintd*symfact
c--
               if(iprint.ge.2) then
                  call print_der1(ics,jcs,kcs,lcs,inx,
     *                 buf(1,ijklp,integ,iqu),
     *                 icf,jcf,kcf,lcf)
               endif
c---
  20         continue
  15      continue
  10    continue
c--------------------------------------------------------
      ELSE
c-----symmetry is not used---------------------------
c
        integral=0
        do 100 ijklp=1,nbls
        ijkl=indxp(ijklp)
        if(ijkl.eq.0) go to 100
        iqreq=ipres(ijkl)
        if(iqreq.eq.0) go to 100
        iorder=iqorder(iqreq)
ctest
c        write(6,*)'destDuf iorder=',iorder
ctest
c--
        call reorder_der1(iorder,lder)
c--
c
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
     *                 ilen,jlen,klen,llen, icfg,jcfg,kcfg,lcfg,ngcq)
        endif
c
          do 150 iqu=1,ngcq
          icff=icfg(iqu)
          jcff=jcfg(iqu)
          kcff=kcfg(iqu)
          lcff=lcfg(iqu)
c
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 200 iii=1,ilen
             icf=icff+iii
             do 200 jjj=1,jlen
             jcf=jcff+jjj
             do 200 kkk=1,klen
             kcf=kcff+kkk
             do 200 lll=1,llen
             lcf=lcff+lll
c---
             integ=integ+1
             integral=integral+1
cccccccc     xint0=buf(ijklp,integ,iqu)
c---
             xinta=buf(1,ijklp,integ,iqu)
             xintb=buf(2,ijklp,integ,iqu)
             xintc=buf(3,ijklp,integ,iqu)
             xintd=-(xinta+xintb+xintc) ! trans. inv.
c
             yinta=buf(4,ijklp,integ,iqu)
             yintb=buf(5,ijklp,integ,iqu)
             yintc=buf(6,ijklp,integ,iqu)
             yintd=-(yinta+yintb+yintc) ! trans. inv.
c
             zinta=buf(7,ijklp,integ,iqu)
             zintb=buf(8,ijklp,integ,iqu)
             zintc=buf(9,ijklp,integ,iqu)
             zintd=-(zinta+zintb+zintc) ! trans. inv.
c---
             ipnl=itxspnl(integral)
ccccccc      buffer(ipnl    )=xint0
c--
             buffer(lder(1),ipnl)=xinta
             buffer(lder(2),ipnl)=yinta
             buffer(lder(3),ipnl)=zinta
c--
             buffer(lder(4),ipnl)=xintb
             buffer(lder(5),ipnl)=yintb
             buffer(lder(6),ipnl)=zintb
c--
             buffer(lder(7),ipnl)=xintc
             buffer(lder(8),ipnl)=yintc
             buffer(lder(9),ipnl)=zintc
c--
             buffer(lder(10),ipnl)=xintd
             buffer(lder(11),ipnl)=yintd
             buffer(lder(12),ipnl)=zintd
c--
               if(iprint.ge.2) then
                  call print_der1(ics,jcs,kcs,lcs,inx,
     *                 buf(1,ijklp,integ,iqu),
     *                 icf,jcf,kcf,lcf)
               endif
c---
  200        continue
  150     continue
  100   continue
c
c--------------------------------------------------------
      ENDIF
c--------------------------------------------------------
      integ_n0=integ_n0+integral
c--------------------------------------------------------
      end
c==============================================================
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
c----------------------------
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

c==============================================================
      subroutine print_der1(ics,jcs,kcs,lcs,inx,
     *                      der1,
     *                      icf,jcf,kcf,lcf)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension iix(4)
      dimension der1(9)
      dimension at_der(3,100)
c
         iat=inx(2,ics)
         jat=inx(2,jcs)
         kat=inx(2,kcs)
         lat=inx(2,lcs)
c
         if(iat.gt.100 .or. jat.gt.100) then
            write(6,*) 'too many atoms; derivatives can not be printed'
            return
         endif
         if(kat.gt.100 .or. lat.gt.100) then
            write(6,*) 'too many atoms; derivatives can not be printed'
            return
         endif
c-----test print only : begining -----------------------------
c
                   xinta=der1(1)
                   xintb=der1(2)
                   xintc=der1(3)
                   yinta=der1(4)
                   yintb=der1(5)
                   yintc=der1(6)
                   zinta=der1(7)
                   zintb=der1(8)
                   zintc=der1(9)
                   xintd=-(xinta+xintb+xintc)
                   yintd=-(yinta+yintb+yintc)
                   zintd=-(zinta+zintb+zintc)
c
                   at_der(1,iat)=0.d0
                   at_der(2,iat)=0.d0
                   at_der(3,iat)=0.d0
                   at_der(1,jat)=0.d0
                   at_der(2,jat)=0.d0
                   at_der(3,jat)=0.d0
                   at_der(1,kat)=0.d0
                   at_der(2,kat)=0.d0
                   at_der(3,kat)=0.d0
                   at_der(1,lat)=0.d0
                   at_der(2,lat)=0.d0
                   at_der(3,lat)=0.d0
c
                   at_der(1,iat)=at_der(1,iat)+xinta
                   at_der(2,iat)=at_der(2,iat)+yinta
                   at_der(3,iat)=at_der(3,iat)+zinta
c
                   at_der(1,jat)=at_der(1,jat)+xintb
                   at_der(2,jat)=at_der(2,jat)+yintb
                   at_der(3,jat)=at_der(3,jat)+zintb
c
                   at_der(1,kat)=at_der(1,kat)+xintc
                   at_der(2,kat)=at_der(2,kat)+yintc
                   at_der(3,kat)=at_der(3,kat)+zintc
c
                   at_der(1,lat)=at_der(1,lat)+xintd
                   at_der(2,lat)=at_der(2,lat)+yintd
                   at_der(3,lat)=at_der(3,lat)+zintd
c
canonical order :
c
                   call descend(icf,jcf,kcf,lcf,iix)
c--------------------------------------------------------------------
c Only for testing second deriv. by numerical differentiation:
c
                    a_min=min( xinta,yinta,zinta )
                    a_max=max( xinta,yinta,zinta )
                    absmax=max( abs(a_min),abs(a_max) )
                    if(absmax.gt.1.d-9) then
c     write(6,60) ics,jcs,kcs,lcs, icf,jcf,kcf,lcf,iat,jat,kat,lat
  60  format('shells=',4i2,' functions=',4i2,' centers=',4i2)
ccccc                  write(6,61) iix(1),iix(2),iix(3),iix(4),
                       write(6,61) icf,jcf,kcf,lcf,
     *                             xinta,yinta,zinta 
                    endif
  61  format('d /dAi   : ijkl=',4i2,1x, 3(f12.7,2x))
c
                    b_min=min( xintb,yintb,zintb )
                    b_max=max( xintb,yintb,zintb )
                    absmax=max( abs(b_min),abs(b_max) )
                    if(absmax.gt.1.d-9) then
c     write(6,60) ics,jcs,kcs,lcs, icf,jcf,kcf,lcf,iat,jat,kat,lat
ccccc                  write(6,62) iix(1),iix(2),iix(3),iix(4),
                       write(6,62) icf,jcf,kcf,lcf,
     *                             xintb,yintb,zintb 
                    endif
  62  format('d /dBi   : ijkl=',4i2,1x, 3(f12.7,2x))
c
                    c_min=min( xintc,yintc,zintc )
                    c_max=max( xintc,yintc,zintc )
                    absmax=max( abs(c_min),abs(c_max) )
                    if(absmax.gt.1.d-9) then
c     write(6,60) ics,jcs,kcs,lcs, icf,jcf,kcf,lcf,iat,jat,kat,lat
cccc                   write(6,63) iix(1),iix(2),iix(3),iix(4),
                       write(6,63) icf,jcf,kcf,lcf,
     *                             xintc,yintc,zintc 
                    endif
  63  format('d /dCi   : ijkl=',4i2,1x, 3(f12.7,2x))
c
                    d_min=min( xintd,yintd,zintd )
                    d_max=max( xintd,yintd,zintd )
                    absmax=max( abs(d_min),abs(d_max) )
                    if(absmax.gt.1.d-9) then
c     write(6,60) ics,jcs,kcs,lcs, icf,jcf,kcf,lcf,iat,jat,kat,lat
cccc                   write(6,64) iix(1),iix(2),iix(3),iix(4),
                       write(6,64) icf,jcf,kcf,lcf,
     *                             xintd,yintd,zintd 
                    endif
  64  format('d /dDi   : ijkl=',4i2,1x, 3(f12.7,2x))
c--------------------------------------------------------------------
         RETURN
c--------------------------------------------------------------------
c
c$$$                   do icart=1,3
c$$$                      if(abs( at_der(icart,iat)) .gt.1d-9) then
c$$$                              write(6,111) 
c$$$     $                             iat       , icart, 
c$$$     $                             iix(1),iix(2),iix(3),iix(4),
c$$$     $                             at_der(icart,iat)
c$$$                      endif
c$$$                   enddo
c$$$                if(jat.ne.iat) then
c$$$                   do icart=1,3
c$$$                      if(abs( at_der(icart,jat)) .gt.1d-9) then
c$$$                              write(6,111) 
c$$$     $                             jat       , icart, 
c$$$cccc $                             icf,jcf,kcf,lcf,
c$$$     $                             iix(1),iix(2),iix(3),iix(4),
c$$$     $                             at_der(icart,jat)
c$$$                      endif
c$$$                   enddo
c$$$                endif
c$$$                if(kat.ne.iat .and. kat.ne.jat) then
c$$$                   do icart=1,3
c$$$                      if(abs( at_der(icart,kat)) .gt.1d-9) then
c$$$                              write(6,111) 
c$$$     $                             kat       , icart, 
c$$$ccc  $                             icf,jcf,kcf,lcf,
c$$$     $                             iix(1),iix(2),iix(3),iix(4),
c$$$     $                             at_der(icart,kat)
c$$$                      endif
c$$$                   enddo
c$$$                endif
c$$$                if(lat.ne.iat .and. lat.ne.jat .and. lat.ne.kat) then
c$$$                   do icart=1,3
c$$$                      if(abs( at_der(icart,lat)) .gt.1d-9) then
c$$$                              write(6,111) 
c$$$     $                             lat       , icart, 
c$$$ccc  $                             icf,jcf,kcf,lcf,
c$$$     $                             iix(1),iix(2),iix(3),iix(4),
c$$$     $                             at_der(icart,lat)
c$$$                      endif
c$$$                   enddo
c$$$                endif
c$$$ 111                          format(' atom=',i2,', xyz=',i2,', ijkl=',
c$$$     $                             4i3,', d(ij|kl)/dq=',f12.7   )
c
c-----test print only : end      -----------------------------
c
      end
c==============================================================
      subroutine destbuf(ikbl,nbls,nblok1,ncs,inx,buf,
     *                   buffer, itxspnl, q4,use_q4,
     *                   icfg,jcfg,kcfg,lcfg,ngcd,lnijkl,indxp,ipres)
c----------------------------------------------------------------
c This is called for PNL-requested ONE contracted shell quartet.
c All Integrals (including zeros) return WITHOUT labels but they 
c have to be in PNL-requested order.
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      logical use_q4
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /neglect/ eps,eps1,epsr
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
c
      dimension nblok1(2,*)
      dimension buf(nbls,lnijkl,ngcd)
      dimension inx(12,*)
c
      dimension buffer(*)
      dimension itxspnl(*)
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension indxp(*),ipres(*)
      dimension q4(*)
c--------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c--------------------------------------------------------
c      ijkl_size=ilen*jlen*klen*llen
c      write(6,*)' IJKL_SIZE=',ijkl_size
      IF(use_q4) THEN
c-----symmetry is used---------------------------
        integral=0
        do 10  ijklp=1,nbls
        ijkl=indxp(ijklp)
        if(ijkl.eq.0) go to 10
        iqreq=ipres(ijkl)
        symfact=q4(iqreq)
c
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
     *                 ilen,jlen,klen,llen, icfg,jcfg,kcfg,lcfg,ngcq)
        endif
c
          do 15  iqu=1,ngcq
          icff=icfg(iqu)
          jcff=jcfg(iqu)
          kcff=kcfg(iqu)
          lcff=lcfg(iqu)
c
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 20  iii=1,ilen
             icf=icff+iii
             do 20  jjj=1,jlen
             jcf=jcff+jjj
             do 20  kkk=1,klen
             kcf=kcff+kkk
             do 20  lll=1,llen
             lcf=lcff+lll
c---
             integ=integ+1
             integral=integral+1
             xint0=buf(ijklp,integ,iqu)
c---
             ipnl=itxspnl(integral)
             buffer(ipnl    )=xint0*symfact
c---
  20         continue
  15      continue
  10    continue
c--------------------------------------------------------
      ELSE
c-----symmetry is not used---------------------------
c
        integral=0
        do 100 ijklp=1,nbls
        ijkl=indxp(ijklp)
        if(ijkl.eq.0) go to 100
c
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
     *                 ilen,jlen,klen,llen, icfg,jcfg,kcfg,lcfg,ngcq)
        endif
c
          do 150 iqu=1,ngcq
          icff=icfg(iqu)
          jcff=jcfg(iqu)
          kcff=kcfg(iqu)
          lcff=lcfg(iqu)
c
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 200 iii=1,ilen
             icf=icff+iii
             do 200 jjj=1,jlen
             jcf=jcff+jjj
             do 200 kkk=1,klen
             kcf=kcff+kkk
             do 200 lll=1,llen
             lcf=lcff+lll
c---
             integ=integ+1
             integral=integral+1
             xint0=buf(ijklp,integ,iqu)
c---
             ipnl=itxspnl(integral)
             buffer(ipnl    )=xint0
c                   write(6,66) icf,jcf,kcf,lcf,xint0
c  66               format(4(i3,1x),2x,'INT=',f12.8)
c---
  200        continue
  150     continue
  100   continue
c
c--------------------------------------------------------
      ENDIF
c--------------------------------------------------------
      integ_n0=integ_n0+integral
c--------------------------------------------------------
      end
c==============================================================
#if defined(IBM)
*IBM COMPILER OPTIONS JUST FOR DESTBUL
@PROCESS OPT(2)
#endif
      subroutine destbul(ikbl,nbls,nblok1,ncs,inx,buf,
     *     buffer, icfx,jcfx,kcfx,lcfx, q4, use_q4,
     *     icfg,jcfg,kcfg,lcfg,ngcd,lnijkl,indxp,ipres,iqorder,
     *     map_txs_pnl)
c----------------------------------------------------------------
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
      integer map_txs_pnl(*)        ! txs to pnl basis map = ncfunct
      logical use_q4
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /neglect/ eps,eps1,epsr
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
c----------------------------------------------------------------------
      double precision savezerotol
      common /csavezerotol/ savezerotol ! Used in detbul, set in texas_hf
c----------------------------------------------------------------------
c     
      dimension icfx(*),jcfx(*),kcfx(*),lcfx(*)
      dimension nblok1(2,*)
      dimension buf(*)          ! buf(nbls,lnijkl,ngcd)
      dimension inx(12,*)
c     
      dimension buffer(*)
c     
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension ipres(*), iqorder(*)
      dimension indxp(*)
      dimension q4(*)
c
      double precision threshold ! For screening output integrals
c
c TEST --------------------------
c
c      IF( max(ilen,jlen,klen,llen).ge.5) then
c      write(6,*) ' shell-size=',ilen,jlen,klen,llen
c      ENDIF
c
c TEST --------------------------
c     do not zero out integ_n0 here
c----------------------------
c     loop over quartets belonging to the block IKBL :
c     
c      ijkl_size=ilen*jlen*klen*llen
c     
      do 10  ijklp=1,nbls
         ijkl=indxp(ijklp)
c     
         if(ijkl.eq.0) go to 10
         iqreq=ipres(ijkl)
c?   
         if(iqreq.eq.0) go to 10
         iorder=iqorder(iqreq)
ctest
c        write(6,*)'destBul iorder=',iorder
ctest
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
ctest
c        if(use_q4) THEN
c           write(6,666) ijklp,ijkl,iqreq, ics,jcs,kcs,lcs,symfact
c666        format('ijklp,ijkl,iqreq,=',3i4,' shells=',4i3,2x,f12.7)
c        endif
ctest
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
            icff=map_txs_pnl(icff+1)-1 ! Relies on txs order == pnl order
            jcff=map_txs_pnl(jcff+1)-1
            kcff=map_txs_pnl(kcff+1)-1
            lcff=map_txs_pnl(lcff+1)-1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CTEST
c-----test print only : begining -----------------------------
c      IF( max(ilen,jlen,klen,llen).ge.5) then
c              integ=ijklp + (iqu-1)*nbls*lnijkl
c              do icf=icff+1,icff+ilen
c                 do jcf=jcff+1,jcff+jlen
c                    do kcf=kcff+1,kcff+klen
c                       do lcf=lcff+1,lcff+llen
c                          xint0=buf(integ)
cccccccccccc-----------    if(abs(xint0).gt. threshold ) then
cccccccccccc-----------       integ_n0=integ_n0+1
cccccccccccc-----------       buffer(integ_n0)=xint0*symfact
cccccccccccc-----------    endif
c                          write(6,66) icf,jcf,kcf,lcf,xint0
c  66 format(4(i3,1x),2(f12.7,1x))
c                          integ=integ+nbls
c                       enddo
c                    enddo
c                 enddo
c              enddo
c      ENDIF
c-----test print only : end      -----------------------------
CTEST
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Indices and integrals in the quartet ijkl :
c     
            integ=ijklp + (iqu-1)*nbls*lnijkl
            if(iorder.eq.1234) then
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. threshold ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=icf
                              jcfx(integ_n0)=jcf
                              kcfx(integ_n0)=kcf
                              lcfx(integ_n0)=lcf
c                   write(6,66) icf,jcf,kcf,lcf,xint0
c  66               format(4(i3,1x),2x,'INT=',f12.8)
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            else if(iorder.eq.1243) then
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. threshold ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=icf
                              jcfx(integ_n0)=jcf
                              kcfx(integ_n0)=lcf
                              lcfx(integ_n0)=kcf
c                   write(6,66) icf,jcf,lcf,kcf,xint0
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            else if(iorder.eq.2134) then 
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. threshold ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=jcf
                              jcfx(integ_n0)=icf
                              kcfx(integ_n0)=kcf
                              lcfx(integ_n0)=lcf
c                   write(6,66) jcf,icf,kcf,lcf,xint0
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            else if(iorder.eq.2143) then 
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. threshold ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=jcf
                              jcfx(integ_n0)=icf
                              kcfx(integ_n0)=lcf
                              lcfx(integ_n0)=kcf
c                   write(6,66) jcf,icf,lcf,kcf,xint0
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            else if(iorder.eq.3412) then
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. threshold ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=kcf
                              jcfx(integ_n0)=lcf
                              kcfx(integ_n0)=icf
                              lcfx(integ_n0)=jcf
c                   write(6,66) kcf,lcf,icf,jcf,xint0
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            else if(iorder.eq.4312) then 
c                   means:    ijkl 
c                   i.e. i->l(4), j->k(3), k->i(1), l->j(2)
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. threshold ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=kcf
                              jcfx(integ_n0)=lcf
                              kcfx(integ_n0)=jcf
                              lcfx(integ_n0)=icf
c                   write(6,66) kcf,lcf,jcf,icf,xint0
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            else if(iorder.eq.3421) then 
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. threshold ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=lcf
                              jcfx(integ_n0)=kcf
                              kcfx(integ_n0)=icf
                              lcfx(integ_n0)=jcf
c                   write(6,66) lcf,kcf,icf,jcf,xint0
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            else if(iorder.eq.4321) then
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. threshold ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=lcf
                              jcfx(integ_n0)=kcf
                              kcfx(integ_n0)=jcf
                              lcfx(integ_n0)=icf
c                   write(6,66) lcf,kcf,jcf,icf,xint0
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            endif
         enddo
 10      continue
c--------------------------------------------------------
         end
c==============================================================
      subroutine lab_req(iorder,icf,jcf,kcf,lcf,iix)
      dimension iix(4)
         if(iorder.eq.1234) then
            iix(1)=icf
            iix(2)=jcf
            iix(3)=kcf
            iix(4)=lcf
            return
         endif
         if(iorder.eq.1243) then
            iix(1)=icf
            iix(2)=jcf
            iix(3)=lcf
            iix(4)=kcf
            return
         endif
         if(iorder.eq.2134) then 
            iix(1)=jcf
            iix(2)=icf
            iix(3)=kcf
            iix(4)=lcf
            return
         endif
         if(iorder.eq.2143) then 
            iix(1)=jcf
            iix(2)=icf
            iix(3)=lcf
            iix(4)=kcf
            return
         endif
         if(iorder.eq.3412) then
            iix(1)=kcf
            iix(2)=lcf
            iix(3)=icf
            iix(4)=jcf
            return
         endif
         if(iorder.eq.4312) then 
cccc        i->4, j->3, k->1, l->2
            iix(1)=kcf
            iix(2)=lcf
            iix(3)=jcf
            iix(4)=icf
            return
         endif
         if(iorder.eq.3421) then 
cccc        i->3, j->4, k->2, l->1
            iix(1)=lcf
            iix(2)=kcf
            iix(3)=icf
            iix(4)=jcf
            return
         endif
         if(iorder.eq.4321) then
            iix(1)=lcf
            iix(2)=kcf
            iix(3)=jcf
            iix(4)=icf
            return
         endif
c
      end
c==============================================================
      subroutine indexg(inx,ics,jcs,kcs,lcs,ijcs,klcs,
     *               ilen,jlen,klen,llen,
     *               icfg,jcfg,kcfg,lcfg,ngcq)
#include "errquit.fh"
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      dimension inx(12,*)
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      integer mxdim 
      parameter (mxdim=200)
      dimension iix(mxdim),jjx(mxdim),kkx(mxdim),llx(mxdim)
c dim. 100 should be enough since the max. ge.con is 9, so 81
c is actually the max. for iix,jjx,kkx,llx
c
c HvD: The above comment is incorrect. Someone must have changed 
c the code without adjusting this parameter so we better check.
ccccccccccc
c
         icff=inx(11,ics)
         jcff=inx(11,jcs)
         kcff=inx(11,kcs)
         lcff=inx(11,lcs)
c
         if ((ngci+1)*(ngcj+1).gt.mxdim) 
     +     call errquit("indexg: ijpg will exceed mxdim",
     +                  (ngci+1)*(ngcj+1),UERR)
         if ((ngck+1)*(ngcl+1).gt.mxdim) 
     +     call errquit("indexg: klpg will exceed mxdim",
     +                  (ngck+1)*(ngcl+1),UERR)
             ijpg=0
             icf=icff
             do 2041 igc=0,ngci
c            ngcjx=ngcj
ctry         if(jcs.eq.ics) ngcjx=igc
               jcf=jcff
               do 2042 jgc=0,ngcj
               ijpg=ijpg+1
               iix(ijpg)=icf
               jjx(ijpg)=jcf
c
               jcf=jcf+jlen
 2042          continue
             icf=icf+ilen
 2041        continue
c
             klpg=0
             kcf=kcff
             do 2043 kgc=0,ngck
c            ngclx=ngcl
ctry         if(lcs.eq.kcs) ngclx=kgc
               lcf=lcff
               do 2044 lgc=0,ngcl
               klpg=klpg+1
               kkx(klpg)=kcf
               llx(klpg)=lcf
               lcf=lcf+llen
 2044          continue
             kcf=kcf+klen
 2043        continue
c
             ijklg=0
             do 2045 ijp1=1,ijpg
             klpx=klpg
ckw          if(klcs.eq.ijcs) klpx=ijp1
             do 2045 klp1=1,klpx
             ijklg=ijklg+1
c
             icfg(ijklg)=iix(ijp1)
             jcfg(ijklg)=jjx(ijp1)
             kcfg(ijklg)=kkx(klp1)
             lcfg(ijklg)=llx(klp1)
c
 2045        continue
c
      ngcq=ijklg
c
      return
      end
c==============================================================
      subroutine descend(i,j,k,l,iix)
      dimension iix(4) 
c     this routine orders the 4 indices in canonical order
c    i>=j; k>=l; ij>=kl
      ij1=max0(i,j)
      ij0=min0(i,j)
      kl1=max0(k,l)
      kl0=min0(k,l)
      if(ij1.gt.kl1.or.(ij1.eq.kl1.and.ij0.ge.kl0)) then
        iix(1)=ij1
        iix(2)=ij0
        iix(3)=kl1
        iix(4)=kl0
      else
        iix(1)=kl1
        iix(2)=kl0
        iix(3)=ij1
        iix(4)=ij0
      end if
      end        
c==============================================================
      subroutine reorder_der1(iorder,lder)
      dimension lder(12)
c
c reordering for centers is DIFFERENT than for labels
c
c pqrs(=iorder)
c ABCD   p->A(123), q->B(456), r->C(789), s->D(10 11 12)
c
c for instance 
c iorder 3412 means that
c 3( 7  8  9) goes to position 1 i.e. lder(1)=7, lder(2)=8, lder(3)=9
c 4(10 11 12) goes to position 2 i.e. lder(4)=10 lder(5)=11 lder(6)=12
c 1( 1  2  3) goes to position 3 i.e. lder(7)=1  lder(8)=2  lder(9)=3
c 2( 4  5  6) goes to position 4 i.e. lder(10=4  lder(11=5  lder(12=6
c---------|-------------------------------------------
c lder(i) |                   iorder
c---------|-------------------------------------------
c         | 1234  1243 2134 2143  3412 3421 4312 4321
c---------|-------------------------------------------
c lder( 1)|  1     1    4    4     7    7   10   10
c lder( 2)|  2     2    5    5     8    8   11   11
c lder( 3)|  3     3    6    6     9    9   12   12
c
c lder( 4)|  4     4    1    1    10   10    7    7
c lder( 5)|  5     5    2    2    11   11    8    8
c lder( 6)|  6     6    3    3    12   12    9    9
c
c lder( 7)|  7    10    7   10     1    4    1    4
c lder( 8)|  8    11    8   11     2    5    2    5
c lder( 9)|  9    12    9   12     3    6    3    6
c
c lder(10)| 10     7   10    7     4    1    4    1
c lder(11)| 11     8   11    8     5    2    5    2
c lder(12)| 12     9   12    9     6    3    6    3
c-----------------------------------------------------------------------
c
      do ll=1,12
      lder(ll)=ll
      enddo
c
         if(iorder.eq.1234) return
c
         if(iorder.eq.1243) then
            lder(7)=10
            lder(8)=11
            lder(9)=12
            lder(10)=7
            lder(11)=8
            lder(12)=9
            return
         endif
         if(iorder.eq.2134) then 
            lder(1)=4
            lder(2)=5
            lder(3)=6
            lder(4)=1
            lder(5)=2
            lder(6)=3
            return
         endif
         if(iorder.eq.2143) then 
            lder(1)=4
            lder(2)=5
            lder(3)=6
            lder(4)=1
            lder(5)=2
            lder(6)=3
            lder(7)=10
            lder(8)=11
            lder(9)=12
            lder(10)=7
            lder(11)=8
            lder(12)=9
            return
         endif
         if(iorder.eq.3412) then
            lder(1)=7
            lder(2)=8
            lder(3)=9
            lder(4)=10
            lder(5)=11
            lder(6)=12
            lder(7)=1
            lder(8)=2
            lder(9)=3
            lder(10)=4
            lder(11)=5
            lder(12)=6
            return
         endif
         if(iorder.eq.4312) then 
            lder(1)=10
            lder(2)=11
            lder(3)=12
            lder(4)= 7
            lder(5)= 8
            lder(6)= 9
            lder(7)= 1
            lder(8)= 2
            lder(9)= 3
            lder(10)=4
            lder(11)=5
            lder(12)=6
            return
         endif
         if(iorder.eq.3421) then 
            lder(1)= 7
            lder(2)= 8
            lder(3)= 9
            lder(4)=10
            lder(5)=11
            lder(6)=12
            lder(7)= 4
            lder(8)= 5
            lder(9)= 6
            lder(10)=1
            lder(11)=2
            lder(12)=3
            return
         endif
         if(iorder.eq.4321) then
            lder(1)=10
            lder(2)=11
            lder(3)=12
            lder(4)=7
            lder(5)=8
            lder(6)=9
            lder(7)=4
            lder(8)=5
            lder(9)=6
            lder(10)=1
            lder(11)=2
            lder(12)=3
            return
         endif
c
      end
c==============================================================
      subroutine make_78_from_45(der2,worker)
      implicit real*8 (a-h,o-z)
      dimension der2(45),worker(78)
c
c block aa:
       axax=der2(1)
       axay=der2(2)
       axaz=der2(3)
       ayay=der2(4)
       ayaz=der2(5)
       azaz=der2(6)
c
c block ab:
       axbx=der2(7)
       axby=der2(8)
       axbz=der2(9)
       aybx=der2(10)
       ayby=der2(11)
       aybz=der2(12)
       azbx=der2(13)
       azby=der2(14)
       azbz=der2(15)
c block ac:
       axcx=der2(16)
       axcy=der2(17)
       axcz=der2(18)
       aycx=der2(19)
       aycy=der2(20)
       aycz=der2(21)
       azcx=der2(22)
       azcy=der2(23)
       azcz=der2(24)
c block bb:
       bxbx=der2(25)
       bxby=der2(26)
       bxbz=der2(27)
       byby=der2(28)
       bybz=der2(29)
       bzbz=der2(30)
c block bc:
       bxcx=der2(31)
       bxcy=der2(32)
       bxcz=der2(33)
       bycx=der2(34)
       bycy=der2(35)
       bycz=der2(36)
       bzcx=der2(37)
       bzcy=der2(38)
       bzcz=der2(39)
c block cc:
       cxcx=der2(40)
       cxcy=der2(41)
       cxcz=der2(42)
       cycy=der2(43)
       cycz=der2(44)
       czcz=der2(45)
c block ad: from transl. inv.
       axdx=-(axax+axbx+axcx)
       axdy=-(axay+axby+axcy)
       axdz=-(axaz+axbz+axcz)
       aydx=-(aXaY+aybx+aycx)
       aydy=-(ayay+ayby+aycy)
       aydz=-(ayaz+aybz+aycz)
       azdx=-(aXaZ+azbx+azcx)
       azdy=-(aYaZ+azby+azcy)
       azdz=-(azaz+azbz+azcz)
c block bd: from transl. inv.
       bxdx=-(AxBx+bxbx+bxcx)
       bxdy=-(AYBX+bxby+bxcy)
       bxdz=-(AZBX+bxbz+bxcz)
       bydx=-(AXBY+bXbY+bycx)
       bydy=-(AYBY+byby+bycy)
       bydz=-(AZBY+bybz+bycz)
       bzdx=-(AXBZ+bXbZ+bzcx)
       bzdy=-(AyBz+bybz+bzcy)
       bzdz=-(AZBZ+bzbz+bzcz)
c block cd: from transl. inv.
       cxdx=-(AXCX+BXCX+cxcx)
       cxdy=-(AYCX+BYCX+cxcy)
       cxdz=-(AZCX+BZCX+cxcz)
       cydx=-(AXCY+BXCY+cXcY)
       cydy=-(AYCY+BYCY+cycy)
       cydz=-(AZCY+BZCY+cycz)
       czdx=-(AXCZ+BxCZ+cXcZ)
       czdy=-(AYCZ+BYCZ+cYcZ)
       czdz=-(AZCZ+BZCZ+czcz)
c block dd:
       dxdx=-(AXDX+BXDX+CXDX)
cccc   dxdy=-(dxay+dxby+dxcy)
       dxdy=-(AYDX+BYDX+CYDX)
cccc   dxdz=-(dxaz+dxbz+dxcz)
       dxdz=-(AZDX+BZDX+CZDX)
cccc   dydy=-(dyay+dyby+dycy)
       dydy=-(AYDY+BYDY+CYDY)
ccccc  dydz=-(dyaz+dybz+dycz)
       dydz=-(AZDY+BZDY+CZDY)
cccc   dzdz=-(dzaz+dzbz+dzcz)
       dzdz=-(AZDZ+BZDZ+CZDZ)
c------------------------------------------------------------
c construct all 10 blocks of sec.der. (output) from 6 blocks:
c
c          AA AB AC AD                AA AB AC 
c             BB BC BD      from         BB BC
c                CC CD                      CC
c                   DD  
c      1-6, 7-15,16-24,25-33         1-6, 7-15,16-24
c          34-39,40-48,49-57             25-30,31-39
c                58-63,64-72                   40-45
c                      73-78
c------------------------------------------------------------
c first 24 derivatives are in right order :
c
c blocks AA,AB,AC :
       do m=1,24
          worker(m)=der2(m)
       enddo
c
c block AD :
       worker(25)=axdx    !  =-(axax+axbx+axcx)
       worker(26)=axdy    !  =-(axay+axby+axcy)
       worker(27)=axdz    !  =-(axaz+axbz+axcz)
       worker(28)=aydx    !  =-(aXaY+aybx+aycx)
       worker(29)=aydy    !  =-(ayay+ayby+aycy)
       worker(30)=aydz    !  =-(ayaz+aybz+aycz)
       worker(31)=azdx    !  =-(aXaZ+azbx+azcx)
       worker(32)=azdy    !  =-(aYaZ+azby+azcy)
       worker(33)=azdz    !  =-(azaz+azbz+azcz)
c
c block BB :
       worker(34)=bxbx    !  =der2(25)
       worker(35)=bxby    !  =der2(26)
       worker(36)=bxbz    !  =der2(27)
       worker(37)=byby    !  =der2(28)
       worker(38)=bybz    !  =der2(29)
       worker(39)=bzbz    !  =der2(30)
c
c block BC :
       worker(40)=bxcx    !  =der2(31)
       worker(41)=bxcy    !  =der2(32)
       worker(42)=bxcz    !  =der2(33)
       worker(43)=bycx    !  =der2(34)
       worker(44)=bycy    !  =der2(35)
       worker(45)=bycz    !  =der2(36)
       worker(46)=bzcx    !  =der2(37)
       worker(47)=bzcy    !  =der2(38)
       worker(48)=bzcz    !  =der2(39)
c
c block BD :
       worker(49)=bxdx    !  =-(AxBx+bxbx+bxcx)
       worker(50)=bxdy    !  =-(AYBX+bxby+bxcy)
       worker(51)=bxdz    !  =-(AZBX+bxbz+bxcz)
       worker(52)=bydx    !  =-(AXBY+bXbY+bycx)
       worker(53)=bydy    !  =-(AYBY+byby+bycy)
       worker(54)=bydz    !  =-(AZBY+bybz+bycz)
       worker(55)=bzdx    !  =-(AXBZ+bXbZ+bzcx)
       worker(56)=bzdy    !  =-(AyBz+bybz+bzcy)
       worker(57)=bzdz    !  =-(AZBZ+bzbz+bzcz)
c
c block CC :
       worker(58)=cxcx    !  =der2(40)
       worker(59)=cxcy    !  =der2(41)
       worker(60)=cxcz    !  =der2(42)
       worker(61)=cycy    !  =der2(43)
       worker(62)=cycz    !  =der2(44)
       worker(63)=czcz    !  =der2(45)
c
c block CD :
       worker(64)=cxdx    !  =-(AXCX+BXCX+cxcx)
       worker(65)=cxdy    !  =-(AYCX+BYCX+cxcy)
       worker(66)=cxdz    !  =-(AZCX+BZCX+cxcz)
       worker(67)=cydx    !  =-(AXCY+BXCY+cXcY)
       worker(68)=cydy    !  =-(AYCY+BYCY+cycy)
       worker(69)=cydz    !  =-(AZCY+BZCY+cycz)
       worker(70)=czdx    !  =-(AXCZ+BxCZ+cXcZ)
       worker(71)=czdy    !  =-(AYCZ+BYCZ+cYcZ)
       worker(72)=czdz    !  =-(AZCZ+BZCZ+czcz)
c
c block DD :
       worker(73)=dxdx    !  =-(AXDX+BXDX+CXDX)
       worker(74)=dxdy    !  =-(AYDX+BYDX+CYDX)
       worker(75)=dxdz    !  =-(AZDX+BZDX+CZDX)
       worker(76)=dydy    !  =-(AYDY+BYDY+CYDY)
       worker(77)=dydz    !  =-(AZDY+BZDY+CZDY)
       worker(78)=dzdz    !  =-(AZDZ+BZDZ+CZDZ)
c--------------------------------------------------------------------
      end
c==============================================================
      subroutine reorder_der2(iorder,lder)
      dimension lder(78)
      dimension i_1234(78), i_1243(78), i_2134(78), i_2143(78)
      dimension i_3412(78), i_4312(78), i_3421(78), i_4321(78)
c-----------------------------------------------------
c                  ABCD
c        AA AB AC AD  
c           BB BC BD 
c              CC CD
c                 DD 
      data i_1234/ 1, 2, 3, 4, 5, 6,              ! AA
     *             7, 8, 9,10,11,12,13,14,15,     ! AB
     *            16,17,18,19,20,21,22,23,24,     ! AC
     *            25,26,27,28,29,30,31,32,33,     ! AD
     *            34,35,36,37,38,39,              ! BB
     *            40,41,42,43,44,45,46,47,48,     ! BC
     *            49,50,51,52,53,54,55,56,57,     ! BD
     *            58,59,60,61,62,63,              ! CC
     *            64,65,66,67,68,69,70,71,72,     ! CD
     *            73,74,75,76,77,78          /    ! DD
c-----------------------------------------------------
c                  ABDC
c        AA AB AC AD   <====   AA AB AD AC
c           BB BC BD              BB BD BC
c              CC CD                 DD DC
c                 DD                    CC
      data i_1243/ 1, 2, 3, 4, 5, 6,              ! AA
     *             7, 8, 9,10,11,12,13,14,15,     ! AB
     *            25,26,27,28,29,30,31,32,33,     ! AC<=AD
     *            16,17,18,19,20,21,22,23,24,     ! AD<=AC
     *            34,35,36,37,38,39,              ! BB
     *            49,50,51,52,53,54,55,56,57,     ! BC<=BD
     *            40,41,42,43,44,45,46,47,48,     ! BD<=BC
     *            73,74,75,76,77,78,              ! CC<=DD
     *            64,67,70,65,68,71,66,69,72,     ! CD<=DC
     *            58,59,60,61,62,63          /    ! DD<=CC
c-----------------------------------------------------
c                  BACD
c        AA AB AC AD   <====   BB BA BC BD
c           BB BC BD              AA AC AD
c              CC CD                 CC CD
c                 DD                    DD
      data i_2134/34,35,36,37,38,39,              ! AA<=BB
     *             7,10,13, 8,11,14, 9,12,15,     ! AB<=BA
     *            40,41,42,43,44,45,46,47,48,     ! AC<=BC
     *            49,50,51,52,53,54,55,56,57,     ! AD<=BD
     *             1, 2, 3, 4, 5, 6,              ! BB<=AA
     *            16,17,18,19,20,21,22,23,24,     ! BC<=AC
     *            25,26,27,28,29,30,31,32,33,     ! BD<=AD
     *            58,59,60,61,62,63,              ! CC
     *            64,65,66,67,68,69,70,71,72,     ! CD
     *            73,74,75,76,77,78          /    ! DD

c-----------------------------------------------------
c                  BADC
c        AA AB AC AD   <====   BB BA BD BC
c           BB BC BD              AA AD AC
c              CC CD                 DD DC
c                 DD                    CC
      data i_2143/34,35,36,37,38,39,              ! AA<=BB
     *             7,10,13, 8,11,14, 9,12,15,     ! AB<=BA
     *            49,50,51,52,53,54,55,56,57,     ! AC<=BD
     *            40,41,42,43,44,45,46,47,48,     ! AD<=BC
     *             1, 2, 3, 4, 5, 6,              ! BB<=AA
     *            25,26,27,28,29,30,31,32,33,     ! BC<=AD
     *            16,17,18,19,20,21,22,23,24,     ! BD<=AC
     *            73,74,75,76,77,78,              ! CC<=DD
     *            64,67,70,65,68,71,66,69,72,     ! CD<=DC
     *            58,59,60,61,62,63/              ! DD<=CC
c
c-----------------------------------------------------
c     if(iorder.eq.3412) then
c                  CDAB
c        AA AB AC AD   <====   CC CD CA CB
c           BB BC BD              DD BA DB
c              CC CD                 AA AB
c                 DD                    BB
      data i_3412/58,59,60,61,62,63,              ! AA<=CC
     *            64,65,66,67,68,69,70,71,72,     ! AB<=CD
     *            16,19,22,17,20,23,18,21,24,     ! AC<=CA
     *            40,43,46,41,44,47,42,45,48,     ! AD<=CB
     *            73,74,75,76,77,78,              ! BB<=DD
     *            25,28,31,26,29,32,27,30,33,     ! BC<=BA
     *            49,52,55,50,53,56,51,54,57,     ! BD<=DB
     *             1, 2, 3, 4, 5, 6,              ! CC<=AA
     *             7, 8, 9,10,11,12,13,14,15,     ! CD<=AB
     *            34,35,36,37,38,39          /    ! DD<=BB
c
c-----------------------------------------------------
c     if(iorder.eq.4312) then 
c                  DCAB
c        AA AB AC AD   <====   DD DC DA DB
c           BB BC BD              CC CA CB
c              CC CD                 AA AB
c                 DD                    BB
      data i_4312/73,74,75,76,77,78,              ! AA<=DD
     *            64,67,70,65,68,71,66,69,72,     ! AB<=DC
     *            25,28,31,26,29,32,27,30,33,     ! AC<=DA
     *            49,52,55,50,53,56,51,54,57,     ! AD<=DB
     *            58,59,60,61,62,63,              ! BB<=CC
     *            16,19,22,17,20,23,18,21,24,     ! BC<=CA
     *            40,43,46,41,44,47,42,45,48,     ! BD<=CB
     *             1, 2, 3, 4, 5, 6,              ! CC<=AA
     *             7, 8, 9,10,11,12,13,14,15,     ! CD<=AB
     *            34,35,36,37,38,39          /    ! DD<=BB
c-----------------------------------------------------
c     if(iorder.eq.3421) then 
c                  CDBA
c        AA AB AC AD   <====   CC CD CB CA
c           BB BC BD              DD DB DA
c              CC CD                 BB BA
c                 DD                    AA
      data i_3421/58,59,60,61,62,63,              ! AA<=CC
     *            64,65,66,67,68,69,70,71,72,     ! AB<=CD
     *            40,43,46,41,44,47,42,45,48,     ! AC<=CB
     *            16,19,22,17,20,23,18,21,24,     ! AD<=CA
     *            73,74,75,76,77,78,              ! BB<=DD
     *            49,52,55,50,53,56,51,54,57,     ! BC<=DB
     *            25,28,31,26,29,32,27,30,33,     ! BD<=DA
     *            34,35,36,37,38,39,              ! CC<=BB
     *             7,10,13, 8,11,14, 9,12,15,     ! CD<=BA
     *             1, 2, 3, 4, 5, 6          /    ! DD<=AA
c-----------------------------------------------------
c     if(iorder.eq.4321) then
c                  DCBA
c        AA AB AC AD   <====   DD DC DB DA
c           BB BC BD              CC CB CA
c              CC CD                 BB BA
c                 DD                    AA
      data i_4321/73,74,75,76,77,78,              ! AA<=DD
     *            64,67,70,65,68,71,66,69,72,     ! AB<=DC
     *            49,52,55,50,53,56,51,54,57,     ! AC<=DB
     *            25,28,31,26,29,32,27,30,33,     ! AD<=DA
     *            58,59,60,61,62,63,              ! BB<=CC
     *            40,43,46,41,44,47,42,45,48,     ! BC<=CB
     *            16,19,22,17,20,23,18,21,24,     ! BD<=CA
     *            34,35,36,37,38,39,              ! CC<=BB
     *             7,10,13, 8,11,14, 9,12,15,     ! CD<=BA
     *             1, 2, 3, 4, 5, 6          /    ! DD<=AA
c-----------------------------------------------------
*      write(6,*)'reorder_der2:iorder:',iorder
      if(iorder.eq.1234) call tfer_i(i_1234,lder,78)
      if(iorder.eq.1243) call tfer_i(i_1243,lder,78)
      if(iorder.eq.2134) call tfer_i(i_2134,lder,78)
      if(iorder.eq.2143) call tfer_i(i_2143,lder,78)
      if(iorder.eq.3412) call tfer_i(i_3412,lder,78)
      if(iorder.eq.4312) call tfer_i(i_4312,lder,78)
      if(iorder.eq.3421) call tfer_i(i_3421,lder,78)
      if(iorder.eq.4321) call tfer_i(i_4321,lder,78)
c
      end
c==============================================================
      subroutine print_pnl2(ics,jcs,kcs,lcs,inx,
     *                      der2, icf,jcf,kcf,lcf)
c PRINT ONLY 
      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension iix(4)
      dimension der2(78)
c
         iat=inx(2,ics)
         jat=inx(2,jcs)
         kat=inx(2,kcs)
         lat=inx(2,lcs)
c
c--------------------------------------------------------------------
canonical order :
c
                   call descend(icf,jcf,kcf,lcf,iix)
c--------------------------------------------------------------------
c block AA:
                   write(8,61)'d2/dAidAj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=1,6)
  61  format(a16,4i2,1x,3(f12.4,2x)/39x,2(f12.4,2x)/53x,1(f12.4,2x))
c
c block AB:
                   write(8,62)'d2/dAidBj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=7,15)
  62  format(a16,4i2,1x,3(f12.4,2x)/25x,3(f12.4,2x)/25x,3(f12.4,2x))
c
c block AC:
                   write(8,62)'d2/dAidCj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=16,24)
c
c block AD:
                   write(8,62)'d2/dAidDj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=25,33)
c
c block BB:
                   write(8,61)'d2/dBidBj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=34,39)
c
c block BC:
                   write(8,62)'d2/dBidCj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=40,48)
c
c block BD:
                   write(8,62)'d2/dBidDj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=49,57)
c
c block CC:
                   write(8,61)'d2/dCidCj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=58,63)
c
c block CD:
                   write(8,62)'d2/dCidDj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=64,72)
c
c block DD:
                   write(8,61)'d2/dDidDj: ijkl=', icf,jcf,kcf,lcf, 
     *                 (der2(ii),ii=73,78)
c--------------------------------------------------------------------
c
      end
c==============================================================
