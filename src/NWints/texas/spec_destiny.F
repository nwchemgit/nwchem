c==============================================================
      subroutine destbuf(ikbl,nbls,nblok1,iis,jjs,ncs,inx,buf,
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
      dimension iis(*),jjs(*)
      dimension buf(nbls,lnijkl,ngcd)
      dimension inx(12,*)
      dimension iix(4)
c
      dimension buffer(*)
      dimension itxspnl(*)
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension indxp(*),ipres(*)
      dimension q4(*)
c--------------------------------
      data zero,half /0.d0,0.5d0/
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
        symfact=q4(iqreq)
c
        ijcs=nblok1(1,ijkl)   
        klcs=nblok1(2,ijkl)  
        ics=iis(ijcs)    
        jcs=jjs(ijcs)   
        kcs=iis(klcs)  
        lcs=jjs(klcs) 
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
             ii=icf*(icf-1)/2
             do 20  jjj=1,jlen
             jcf=jcff+jjj
             jj=jcf*(jcf-1)/2
             do 20  kkk=1,klen
             kcf=kcff+kkk
             kk=kcf*(kcf-1)/2
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
        ics=iis(ijcs)    
        jcs=jjs(ijcs)   
        kcs=iis(klcs)  
        lcs=jjs(klcs) 
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
             ii=icf*(icf-1)/2
             do 200 jjj=1,jlen
             jcf=jcff+jjj
             jj=jcf*(jcf-1)/2
             do 200 kkk=1,klen
             kcf=kcff+kkk
             kk=kcf*(kcf-1)/2
             do 200 lll=1,llen
             lcf=lcff+lll
c---
             integ=integ+1
             integral=integral+1
             xint0=buf(ijklp,integ,iqu)
c---
             ipnl=itxspnl(integral)
             buffer(ipnl    )=xint0
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
      subroutine destbul(ikbl,nbls,nblok1,iis,jjs,ncs,inx,buf,
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
c     
      dimension icfx(*),jcfx(*),kcfx(*),lcfx(*)
      dimension nblok1(2,*)
      dimension iis(*),jjs(*)
      dimension buf(*)          ! buf(nbls,lnijkl,ngcd)
      dimension inx(12,*)
      dimension iix(4)
c     
      dimension buffer(*)
c     
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension ipres(*), iqorder(*)
      dimension indxp(*)
      dimension q4(*)
c TEST --------------------------
c
c      IF( max(ilen,jlen,klen,llen).ge.5) then
c      write(6,*) ' BLOCK=',ikbl,' shell-size=',ilen,jlen,klen,llen
c      ENDIF
c
c TEST --------------------------
c--------------------------------
      data zero,half /0.d0,0.5d0/
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
c??   
         if(iqreq.eq.0) go to 10
         iorder=iqorder(iqreq)
         if(use_q4) THEN
            symfact=q4(iqreq)
         else
            symfact = 1.0d0
         endif
c---------------------------------------
c     write(6 ,1230)  ijkl,iqreq,iorder
c     1230 format('quart=',i5,' req-quart=,i5,'  iorder=',i4 )
c---------------------------------------
         ijcs=nblok1(1,ijkl)   
         klcs=nblok1(2,ijkl)  
         ics=iis(ijcs)    
         jcs=jjs(ijcs)   
         kcs=iis(klcs)  
         lcs=jjs(klcs) 
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
cccccccccccc-----------    if(abs(xint0).gt. eps ) then
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
                           if(abs(xint0).gt. eps ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=icf
                              jcfx(integ_n0)=jcf
                              kcfx(integ_n0)=kcf
                              lcfx(integ_n0)=lcf
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
                           if(abs(xint0).gt. eps ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=icf
                              jcfx(integ_n0)=jcf
                              kcfx(integ_n0)=lcf
                              lcfx(integ_n0)=kcf
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
                           if(abs(xint0).gt. eps ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=jcf
                              jcfx(integ_n0)=icf
                              kcfx(integ_n0)=kcf
                              lcfx(integ_n0)=lcf
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
                           if(abs(xint0).gt. eps ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=jcf
                              jcfx(integ_n0)=icf
                              kcfx(integ_n0)=lcf
                              lcfx(integ_n0)=kcf
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
                           if(abs(xint0).gt. eps ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=kcf
                              jcfx(integ_n0)=lcf
                              kcfx(integ_n0)=icf
                              lcfx(integ_n0)=jcf
                           endif
                           integ=integ+nbls
                        enddo
                     enddo
                  enddo
               enddo
            else if(iorder.eq.4312) then 
               do icf=icff+1,icff+ilen
                  do jcf=jcff+1,jcff+jlen
                     do kcf=kcff+1,kcff+klen
                        do lcf=lcff+1,lcff+llen
                           xint0=buf(integ)
                           if(abs(xint0).gt. eps ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=kcf
                              jcfx(integ_n0)=lcf
                              kcfx(integ_n0)=jcf
                              lcfx(integ_n0)=icf
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
                           if(abs(xint0).gt. eps ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=lcf
                              jcfx(integ_n0)=kcf
                              kcfx(integ_n0)=icf
                              lcfx(integ_n0)=jcf
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
                           if(abs(xint0).gt. eps ) then
                              integ_n0=integ_n0+1
                              buffer(integ_n0)=xint0*symfact
                              icfx(integ_n0)=lcf
                              jcfx(integ_n0)=kcf
                              kcfx(integ_n0)=jcf
                              lcfx(integ_n0)=icf
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
            iix(1)=kcf
            iix(2)=lcf
            iix(3)=jcf
            iix(4)=icf
            return
         endif
         if(iorder.eq.3421) then 
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
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      dimension inx(12,*)
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension iix(100),jjx(100),kkx(100),llx(100)
c dim. 100 should be enough since the max. ge.con is 9, so 81
c is actually the max. for iix,jjx,kkx,llx
ccccccccccc
c
         icff=inx(11,ics)
         jcff=inx(11,jcs)
         kcff=inx(11,kcs)
         lcff=inx(11,lcs)
c
             ijpg=0
             icf=icff
             do 2041 igc=0,ngci
             ngcjx=ngcj
             if(jcs.eq.ics) ngcjx=igc
               jcf=jcff
               do 2042 jgc=0,ngcjx
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
             ngclx=ngcl
             if(lcs.eq.kcs) ngclx=kgc
               lcf=lcff
               do 2044 lgc=0,ngclx
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
             if(klcs.eq.ijcs) klpx=ijp1
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
