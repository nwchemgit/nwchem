c=================================================================
c Only for test calculations :
c a set of contracted shell quartets is requested at the time :
c
c
      subroutine test_calc2(nshells, l_blsize)
      implicit real*8 (a-h,o-z)
      logical more_int
c-----------------------------------
c     parameter (leri= 2 000 000) 
c     parameter (nquart=1000)
c     parameter(l_blscr=3 000 000)
      parameter (leri= 100 000) 
      parameter (nquart=1000)
      parameter(l_blscr=500 000)
c-----------------------------------
      common /check_int/ integ_check
c-----------------------------------
      dimension blscr(l_blscr)
c-----------------------------------
c returning integrals and indeces :
      dimension eri(leri) 
      dimension icf(leri),jcf(leri),kcf(leri),lcf(leri)
      dimension ra(3),rb(3),rc(3),rd(3)
c requested shells :
      dimension ics(nquart),jcs(nquart),kcs(nquart),lcs(nquart)
      dimension q4(nquart)
c-------------------------
      dimension stsum(5,3)
c--------------------------------------------------------------
c
c     write(6,*) '***   ENTERING TEST_CALC2   ***'
c     write(6,*) ' estimated BL() scratch =',l_blsize
c     write(6,*) ' used      BL() scratch =',l_blscr 
c     write(6,*) '-------------------------------'
c
c--------------------------------------------------------------
      call txs_second(tcal1)
c-------------------------
      ij_basis=566
      kl_basis=566
c-------------------------
      integ_check=0
c-------------------------
      do 10 ii=1,3 
      do 10 jj=1,5
      stsum(jj,ii)=0.d0
   10 continue
c-------------------------
      more_int=.false.
c-------------------------
c total number of cont.shell quartets :
c
      nq1=nshells*(nshells+1)/2
      nquartets=nq1*(nq1+1)/2
c
      integrals=0
      ncalls=0
      nqrt=0
      ijkl=0
      ijsh=0
      do 100 ish=1,nshells
        do 200 jsh=1,ish
          ijsh=ijsh+1
          klsh=0
          do 300 ksh=1,nshells
            do 400 lsh=1,ksh
            klsh=klsh+1
            if(ijsh.ge.klsh) then
               ijkl=ijkl+1
               nqrt=nqrt+1
  410          continue
               if(nqrt.lt.nquart) then
                  ics(nqrt)=ish
                  jcs(nqrt)=jsh
                  kcs(nqrt)=ksh
                  lcs(nqrt)=lsh
                  if(ijkl.eq.nquartets) then
                     ncalls=ncalls+1
                     write(6,*)'* texas_hf called ',ncalls,' time *'
  451                continue
                     call texas_hf2_m(ij_basis,ics,jcs,kl_basis,kcs,lcs,
CCC  *                     nqrt,q4,use_q4,
     *                     nqrt,q4,.false.,
     *                     ra,rb,rc,rd,.false., eri,leri,
     *                     icf,jcf,kcf,lcf,integ_n0,.true.,
     *                     more_int,blscr,l_blscr)
c
        call print_int1(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
                     call check_sums(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
                     integrals=integrals+integ_n0
                     if(more_int) go to 451
                   endif
                else if(nqrt.eq.nquart) then
                  ics(nqrt)=ish
                  jcs(nqrt)=jsh
                  kcs(nqrt)=ksh
                  lcs(nqrt)=lsh
                     ncalls=ncalls+1
                     write(6,*)'* texas_hf called ',ncalls,' time *'
  452                continue
                     call texas_hf2_m(ij_basis,ics,jcs,kl_basis,kcs,lcs,
CCC  *                     nqrt,q4,use_q4,
     *                     nqrt,q4,.false.,
     *                     ra,rb,rc,rd,.false., eri,leri,
     *                     icf,jcf,kcf,lcf,integ_n0,.true.,
     *                     more_int,blscr,l_blscr)
        call print_int1(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
                     call check_sums(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
                     integrals=integrals+integ_n0
                     if(more_int) go to 452
                  nqrt=0
                else
                  nqrt=1
                  go to 410
                endif
            endif
  400       continue
  300     continue
  200   continue
  100 continue
c
c-------------------------
      write(6,*)'-------------------------------------------------'
      write(6,*)'total number of integrals processed =',integrals
      write(6,*)'total number of integrals checked in=',integ_check
      write(6,*)'-------------------------------------------------'
      write(6,*)'-------------------------------------------------'
      write(6,*)'total number of integrals processed =',integrals
      write(6,*)'total number of integrals checked in=',integ_check
      write(6,*)'-------------------------------------------------'
      call print_check(stsum)
c-------------------------
      call txs_second(tcal2)
      write(6 ,66) tcal2-tcal1
      write(6,66) tcal2-tcal1
   66 format('Total CPU check-sum time =',f10.2)
      write(6,*)'-------------------------------------------------'
c--------------------------------------------------------------
c
      write(6,*) '***   LEAVING TEST_CALC2   ***'
c
c-------------------------
      call texas_terminate()
c-------------------------
c--------------------------------------------------------------
      stop 'stopped after test calculations by kwol in test_calc2'
c-------------------------
c
      end
c=================================================================
      subroutine check_sums(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
      implicit real*8 (a-h,o-z)
      common /check_int/ integ_check
      dimension stsum(5,3)
      dimension icf(*),jcf(*),kcf(*),lcf(*)
      dimension eri(*)
c
      eps=10.d0**(-7.5d0)
c
      int=0
      do 100 ii=1,integ_n0
      x=eri(ii)
      if(abs(x).gt.eps) then
        int=int+1
        x2=x*x
        i=icf(ii)
        j=jcf(ii)
        k=kcf(ii)
        l=lcf(ii)
        denom=1.d0/(1.d0 + x2)
c
        s1=dble(i+j+k+l)
        xij=dble(i-j)
        xkl=dble(k-l)
        xij=abs(xij)
        xkl=abs(xkl)
        d1=abs(xij-xkl)
c
        stsum(1,1)=stsum(1,1)+x
        stsum(2,1)=stsum(2,1)+abs(x)
        stsum(3,1)=stsum(3,1)+x2
        stsum(4,1)=stsum(4,1)+x*denom
        stsum(5,1)=stsum(5,1)+x2*denom
c
        stsum(1,2)=stsum(1,2)+x        *s1
        stsum(2,2)=stsum(2,2)+abs(x)   *s1
        stsum(3,2)=stsum(3,2)+x2       *s1
        stsum(4,2)=stsum(4,2)+x*denom  *s1
        stsum(5,2)=stsum(5,2)+x2*denom *s1
c
        stsum(1,3)=stsum(1,3)+x        *d1
        stsum(2,3)=stsum(2,3)+abs(x)   *d1
        stsum(3,3)=stsum(3,3)+x2       *d1
        stsum(4,3)=stsum(4,3)+x*denom  *d1
        stsum(5,3)=stsum(5,3)+x2*denom *d1
      endif
c
  100 continue
c
      integ_check=integ_check+int
c
      end
c=================================================================
      subroutine print_check(stsum)
      implicit real*8 (a-h,o-z)
      dimension stsum(5,3)
c
      write(6,*)'-------------------------------------------------'
      write(6,*)'check_sums '
      write(6,*)' '
      do ii=1,3
      write(6,66) (stsum(jj,ii),jj=1,5)
      enddo
      write(6,*)'-------------------------------------------------'
c
c
      write(6,*)'-------------------------------------------------'
      write(6,*)'check_sums '
      write(6,*)' '
      do ii=1,3
      write(6,66) (stsum(jj,ii),jj=1,5)
      enddo
  66  format(5(1pd15.8,3x))
      write(6,*)'-------------------------------------------------'
c
      end
c=================================================================
      subroutine blprint(number,ibas,inuc,na,ncs,nsh,ncf,inx,
     *                   datbas,datnuc)
      implicit real*8 (a-h,o-z)
      dimension datbas(13,nsh),datnuc(5,*) 
      dimension inx(12,*)
c
      write(8,*)' from blprint call number=',number
      write(8,11) ncs,nsh,ncf, ibas,inuc
   11 format('ncs=',i3,' nsh=',i3,' ncf=',i3,' ibas=',i5,' inuc=',i5)
c
      write(8,*)' datbas :'
      do 100 i=1,13
      write(8,88) (datbas(i,ii),ii=1,nsh)
   88 format(5(f12.7,2x))
  100 continue
c     
      write(8,*)' datnuc :'
      do 200 i=1,5
      write(8,88) (datnuc(i,ii),ii=1,na)
  200 continue
c     
      write(8,*)' inx    :'
      do 300 i=1,12
      write(8,77) (inx(i,ii),ii=1,ncs)
   77 format(5(i5,2x))
  300 continue
c
      end
c=====================================================================
      subroutine print_int1(ncall, eri,leri,icf,jcf,kcf,lcf,integ_n0 )
      implicit real*8 (a-h,o-z)
      dimension eri(leri) 
      dimension icf(leri),jcf(leri),kcf(leri),lcf(leri)
c
        write(6,*)'call no=',ncall,' ; ',integ_n0,' non-zero integ.'
c
        do ii=1,integ_n0
        i1=icf(ii)
        j1=jcf(ii)
        k1=kcf(ii)
        l1=lcf(ii)
          write(6,44) i1,j1,k1,l1,eri(ii)
        enddo
c
cc 77 format('pnl_lab: ',4(i4,1x),1x,f20.10 )
   33 format(3(i4,1x),1x,f30.20,6x,a5 )
   44 format(4(i4,1x),1x,f30.20 )
      end
c=====================================================================
      subroutine print_int2(ncall,eri,leri,ics,jcs,kcs,lcs,integ_n0,
     *                      nquart)
c no labels 
      implicit real*8 (a-h,o-z)
      dimension eri(leri) 
      dimension ics(nquart),jcs(nquart),kcs(nquart),lcs(nquart)
c
        write(6,*)'call no=',ncall,' ; ',integ_n0,' integrals (all)'
        write(6,*)' no of quartets=',nquart
c
        do iis=1,nquart
          is=ics(iis)
          js=jcs(iis)
          ks=kcs(iis)
          ls=lcs(iis)
          write(6,22) is,js,ks,ls
   22     format('shells =',4(i4,1x))
        enddo
c
        do ii=1,integ_n0
          write(6,44) eri(ii)
        enddo
c
   33 format(3(i4,1x),1x,f30.20,6x,a5 )
   44 format(20x     ,1x,f30.20 )
      end
c=====================================================================
