* $Id$
c=================================================================
c Only for test calculations :
c
cccc  subroutine test_cent4(num_bas_ij,num_bas_kl,l_blsize)
      subroutine test_cent4(num_bas_ij,num_bas_kl,l_blscr )
      implicit real*8 (a-h,o-z)
      logical more_int
c-----------------------------------
c ordinary four center two-electron integrals
c-----------------------------------
      parameter (leri= 129 600) 
      parameter (nquart=10 000)
      parameter(l_blsize=3 000 000)
c-----------------------------------
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
c-----------------------------------
      common /check_int/ integ_check
c-----------------------------------
      dimension blscr(l_blsize)
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
      write(6,*) '-------------------------------'
      write(6,*) '***   ENTERING TEST_CENT4   ***'
      write(6,*) 'fixed   l_blsize=',l_blsize,' avail.l_blscr=',l_blscr
      if(l_blscr.gt.l_blsize) stop 'fixed size .lt. needed'    
      write(6,*) '-------------------------------'
c--------------------------------------------------------------
      call txs_second(tcal1)
c-------------------------
      ij_basis=num_bas_ij
      kl_basis=num_bas_kl
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
c ij 
c     if(num_bas_ij.eq.num_bas_1) then  ! take care of compiler warnings
         ish_b=1
         ish_e=ncs_bas_1
         jsh_b=1
         jsh_e=ncs_bas_1
         ij_shell=ncs_bas_1*(ncs_bas_1 + 1)/2
c     endif
      if(num_bas_ij.eq.num_bas_2) then
         ish_b=ncs_bas_1+1
         ish_e=ncs_bas_2
         jsh_b=ncs_bas_1+1
         jsh_e=ncs_bas_2
         n_shell=ncs_bas_2 - ncs_bas_1
         ij_shell=n_shell*(n_shell+1)/2
      endif
c kl 
c     if(num_bas_kl.eq.num_bas_1) then  ! take care of compiler warnings
         ksh_b=1
         ksh_e=ncs_bas_1
         lsh_b=1
         lsh_e=ncs_bas_1
         kl_shell=ncs_bas_1*(ncs_bas_1 + 1)/2
c     endif
      if(num_bas_kl.eq.num_bas_2) then
         ksh_b=ncs_bas_1+1
         ksh_e=ncs_bas_2
         lsh_b=ncs_bas_1+1
         lsh_e=ncs_bas_2
         n_shell=ncs_bas_2 - ncs_bas_1
         kl_shell=n_shell*(n_shell+1)/2
      endif
c
      if(ij_basis.eq.kl_basis) then
         nquartets=ij_shell*(ij_shell+1)/2
      else
         nquartets=ij_shell*kl_shell
      endif
c
c-------------------------
c
      integrals=0
      ncalls=0
      nqrt=0
      ijkl=0
      ijsh=0
      ish=0
c-------------------------
      nqrt=0
      do 100 ish1=ish_b,ish_e
      ish=ish+1
        jsh=0
        do 200 jsh1=jsh_b,ish1
        jsh=jsh+1
          ijsh=ijsh+1
          klsh=0
          ksh=0
          do 300 ksh1=ksh_b,ksh_e
          ksh=ksh+1
            lsh=0
            do 400 lsh1=lsh_b,ksh1
            lsh=lsh+1
            klsh=klsh+1
c>>>>       if(ijsh.ge.klsh) then
            if(ijsh.ge.klsh .OR. ij_basis.ne.kl_basis) then
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
ccccccccccccccc      write(6,*)'* texas_hf called ',ncalls,' time *'
  451                continue
                     call texas_hf2_m(ij_basis,ics,jcs,kl_basis,kcs,lcs,
CCC  *                     nqrt,q4,use_q4,
     *                     nqrt,q4,.false.,
     *                     ra,rb,rc,rd,.false., eri,leri,
     *                     icf,jcf,kcf,lcf,integ_n0,.true.,
     *                     more_int,blscr,l_blscr,0.0d0,'scfd_int')
ccc  *                     more_int,blscr,l_blscr,0.0d0,'der1_int')
ccc  *                     more_int,blscr,l_blscr,0.0d0,'der2_int')
c
ccccc   call print_int1(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
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
ccccccccccccccc      write(6,*)'* texas_hf called ',ncalls,' time *'
  452                continue
                     call texas_hf2_m(ij_basis,ics,jcs,kl_basis,kcs,lcs,
CCC  *                     nqrt,q4,use_q4,
     *                     nqrt,q4,.false.,
     *                     ra,rb,rc,rd,.false., eri,leri,
     *                     icf,jcf,kcf,lcf,integ_n0,.true.,
     *                     more_int,blscr,l_blscr,0.0d0,'scfd_int')
ccc  *                     more_int,blscr,l_blscr,0.0d0,'der1_int')
ccc  *                     more_int,blscr,l_blscr,0.0d0,'der2_int')
ccccc   call print_int1(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
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
      write(6,*)'       Four-center two-electron integrals    '
      write(6,*)'total number of integrals processed =',integrals
      write(6,*)'total number of integrals checked in=',integ_check
      write(6,*)'-------------------------------------------------'
      call print_check(stsum)
c-------------------------
      call txs_second(tcal2)
      write(6,66) tcal2-tcal1
   66 format('Total CPU check-sum time =',f10.2)
      write(6,*)'-------------------------------------------------'
c--------------------------------------------------------------
c
      write(6,*) '***   LEAVING TEST_CENT4   ***'
c
c--------------------------------------------------------------
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
ccc     write(6,*)'call no=',ncall,' ; ',integ_n0,' non-zero integ.'
c
        do ii=1,integ_n0
        i1=icf(ii)
        j1=jcf(ii)
        k1=kcf(ii)
        l1=lcf(ii)
             write(6,44) i1,j1,k1,l1,eri(ii)
        enddo
c
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
      subroutine test_cent2(num_bas_ij,num_bas_kl,l_blsize)
c-----------------------------------
c ONLY 2-center integrals (I0 | K0)
c nshells includes an additional S0-shell (the last one)
c-----------------------------------
      implicit real*8 (a-h,o-z)
      logical more_int
c-----------------------------------
c     parameter (leri= 2 000 000) 
c     parameter (nquart=1000)
c     parameter(l_blscr=3 000 000)
      parameter (leri= 100 000) 
      parameter (nquart=1000)
      parameter(l_blscr=100 000)
c-----------------------------------
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
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
      write(6,*) '-------------------------------'
      write(6,*) '***   ENTERING TEST_CENT2   ***'
      write(6,*) '-------------------------------'
c
c--------------------------------------------------------------
      call txs_second(tcal1)
c-------------------------
      ij_basis=num_bas_ij
      kl_basis=num_bas_kl
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
c (I 0 | K 0)
c
c ij 
c     if(num_bas_ij.eq.num_bas_1) then   ! take care of compiler warnings
         ish_b=1
         ish_e=ncs_bas_1
         jsh_b=1
         jsh_e=ncs_bas_1
         ij_shell=ncs_bas_1
c     endif
      if(num_bas_ij.eq.num_bas_2) then
         ish_b=ncs_bas_1+1
         ish_e=ncs_bas_2
         jsh_b=ncs_bas_1+1
         jsh_e=ncs_bas_2
         n_shell=ncs_bas_2 - ncs_bas_1
         ij_shell=n_shell
      endif
      if(num_bas_ij.eq.num_bas_3) then
         ish_b=ncs_bas_2+1
         ish_e=ncs_bas_3
         jsh_b=ncs_bas_2+1
         jsh_e=ncs_bas_3
         n_shell=ncs_bas_3 - ncs_bas_2
         ij_shell=n_shell
      endif
c kl 
c     if(num_bas_kl.eq.num_bas_1) then   ! take care of compiler warnings
         ksh_b=1
         ksh_e=ncs_bas_1
         lsh_b=1
         lsh_e=ncs_bas_1
         kl_shell=ncs_bas_1
c     endif
      if(num_bas_kl.eq.num_bas_2) then
         ksh_b=ncs_bas_1+1
         ksh_e=ncs_bas_2
         lsh_b=ncs_bas_1+1
         lsh_e=ncs_bas_2
         n_shell=ncs_bas_2 - ncs_bas_1
         kl_shell=n_shell
      endif
      if(num_bas_kl.eq.num_bas_3) then
         ksh_b=ncs_bas_2+1
         ksh_e=ncs_bas_3
         lsh_b=ncs_bas_2+1
         lsh_e=ncs_bas_3
         n_shell=ncs_bas_3 - ncs_bas_2
         kl_shell=n_shell
      endif
c
      if(ij_basis.eq.kl_basis) then
         nquartets=ij_shell*(ij_shell+1)/2
      else
         nquartets=ij_shell*kl_shell
      endif
c-------------------------
c
      integrals=0
      ncalls=0
      nqrt=0
      ijkl=0
      ijsh=0
c
      ish=0
      do 100 ish1=ish_b,ish_e
      ish=ish+1
        do 200 jsh=0,0
          ijsh=ijsh+1
          klsh=0
          ksh=0
          do 300 ksh1=ksh_b,ksh_e
          ksh=ksh+1
            do 400 lsh=0,0
            klsh=klsh+1
            if(ijsh.ge.klsh .OR. ij_basis.ne.kl_basis) then
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
     *                     more_int,blscr,l_blscr,0.0d0,'der2_int')
cccc *                     more_int,blscr,l_blscr,0.0d0,'der1_int')
c
c       call print_in2c(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
                     call check_sum2(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
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
     *                     more_int,blscr,l_blscr,0.0d0,'der2_int')
cccc *                     more_int,blscr,l_blscr,0.0d0,'der1_int')
c       call print_in2c(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
                     call check_sum2(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
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
      write(6,*)'       Two-center two-electron integrals    '
      write(6,*)'total number of integrals processed =',integrals
      write(6,*)'total number of integrals checked in=',integ_check
      write(6,*)'-------------------------------------------------'
      call print_check(stsum)
c-------------------------
      call txs_second(tcal2)
      write(6,66) tcal2-tcal1
   66 format('Total CPU check-sum time =',f10.2)
      write(6,*)'-------------------------------------------------'
c--------------------------------------------------------------
c
      write(6,*) '***   LEAVING TEST_CENT2   ***'
c
c--------------------------------------------------------------
c
      end
c=================================================================
      subroutine check_sum2(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
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
c>>>>   j=jcf(ii)
        j=0       
        k=kcf(ii)
c>>>>   l=lcf(ii)
        l=0
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
      subroutine print_in2c(ncall, eri,leri,icf,jcf,kcf,lcf,integ_n0 )
      implicit real*8 (a-h,o-z)
      dimension eri(leri) 
      dimension icf(leri),jcf(leri),kcf(leri),lcf(leri)
c
        write(6,*)'call no=',ncall,' ; ',integ_n0,' non-zero integ.'
c
        do ii=1,integ_n0
        i1=icf(ii)
c>>>>   j1=jcf(ii)
        k1=kcf(ii)
c>>>>   l1=lcf(ii)
c 2-center
             write(6,44) i1,k1,eri(ii)
        enddo
c
   44 format(2(i4,1x),1x,f30.20 )
      end
c=====================================================================
      subroutine test_cent3i_kl(num_bas_ij,num_bas_kl,l_blsize)
c-----------------------------------
c 3-center integrals (I 0 | K L)
c-----------------------------------
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
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
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
      write(6,*) '-------------------------------'
      write(6,*) '*** ENTERING TEST_CENT3i_kl ***'
c     write(6,*) '-------------------------------'
c
c--------------------------------------------------------------
      call txs_second(tcal1)
c-------------------------
      ij_basis=num_bas_ij
      kl_basis=num_bas_kl
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
c (I 0 | K L)
c
c ij 
c     if(num_bas_ij.eq.num_bas_1) then   ! take care of compiler warnings
         ish_b=1
         ish_e=ncs_bas_1
         jsh_b=1
         jsh_e=ncs_bas_1
         ij_shell=ncs_bas_1
c     endif
      if(num_bas_ij.eq.num_bas_2) then
         ish_b=ncs_bas_1+1
         ish_e=ncs_bas_2
         jsh_b=ncs_bas_1+1
         jsh_e=ncs_bas_2
         n_shell=ncs_bas_2 - ncs_bas_1
         ij_shell=n_shell
      endif
c kl 
c     if(num_bas_kl.eq.num_bas_1) then   ! take care of compiler warnings
         ksh_b=1
         ksh_e=ncs_bas_1
         lsh_b=1
         lsh_e=ncs_bas_1
         kl_shell=ncs_bas_1*(ncs_bas_1 + 1)/2
c     endif
      if(num_bas_kl.eq.num_bas_2) then
         ksh_b=ncs_bas_1+1
         ksh_e=ncs_bas_2
         lsh_b=ncs_bas_1+1
         lsh_e=ncs_bas_2
         n_shell=ncs_bas_2 - ncs_bas_1
         kl_shell=n_shell*(n_shell+1)/2
      endif
c
      nquartets=ij_shell*kl_shell
c-------------------------
c
      integrals=0
      ncalls=0
      nqrt=0
      ijkl=0
      ish=0
      do 100 ish1=ish_b,ish_e
      ish=ish+1
        do 200 jsh=0,0
          ksh=0
          do 300 ksh1=ksh_b,ksh_e
          ksh=ksh+1
            lsh=0
            do 400 lsh1=lsh_b,ksh1
            lsh=lsh+1
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
     *                     more_int,blscr,l_blscr,0.0d0,'der2_int')
ccc  *                     more_int,blscr,l_blscr,0.0d0,'der1_int')
c
c       call print_i_kl(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
                     call check_sikl(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
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
     *                     more_int,blscr,l_blscr,0.0d0,'der2_int')
cc   *                     more_int,blscr,l_blscr,0.0d0,'der1_int')
c       call print_i_kl(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
                     call check_sikl(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
                     integrals=integrals+integ_n0
                     if(more_int) go to 452
                  nqrt=0
               else
                  nqrt=1
                  go to 410
               endif
  400       continue
  300     continue
  200   continue
  100 continue
c
c-------------------------
      write(6,*)'-------------------------------------------------'
      write(6,*)'       Three-center two-electron integrals    '
      write(6,*)'total number of integrals processed =',integrals
      write(6,*)'total number of integrals checked in=',integ_check
      write(6,*)'-------------------------------------------------'
      call print_check(stsum)
c-------------------------
      call txs_second(tcal2)
      write(6,66) tcal2-tcal1
   66 format('Total CPU check-sum time =',f10.2)
      write(6,*)'-------------------------------------------------'
c--------------------------------------------------------------
c
      write(6,*) '*** LEAVING TEST_CENT3 i_kl ***'
c
c--------------------------------------------------------------
c
      end
c=================================================================
      subroutine check_sikl(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
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
c>>>>   j=jcf(ii)
        j=0       
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
      subroutine print_i_kl(ncall, eri,leri,icf,jcf,kcf,lcf,integ_n0 )
      implicit real*8 (a-h,o-z)
      dimension eri(leri) 
      dimension icf(leri),jcf(leri),kcf(leri),lcf(leri)
c
        write(6,*)'call no=',ncall,' ; ',integ_n0,' non-zero integ.'
c
        do ii=1,integ_n0
        i1=icf(ii)
c>>>>   j1=jcf(ii)
        k1=kcf(ii)
        l1=lcf(ii)
             write(6,44) i1,k1,l1,eri(ii)
        enddo
c
   44 format('< ',i3,' |',i3,1x,i3,' >',3x,f30.20 )
      end
c=====================================================================
      subroutine test_cent3ij_k(num_bas_ij,num_bas_kl,l_blsize)
c-----------------------------------
c ONLY 3-center integrals (I J | K 0)
c nshells includes an additional S0-shell (the last one)
c-----------------------------------
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
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
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
      write(6,*) '-------------------------------'
      write(6,*) '*** ENTERING TEST_CENT3ij_k ***'
      write(6,*) '-------------------------------'
c--------------------------------------------------------------
      call txs_second(tcal1)
c-------------------------
      ij_basis=num_bas_ij
      kl_basis=num_bas_kl
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
c (I J | K 0)
c
c ij 
c     if(num_bas_ij.eq.num_bas_1) then   ! take care of compiler warnings
         ish_b=1
         ish_e=ncs_bas_1
         jsh_b=1
         jsh_e=ncs_bas_1
         ij_shell=ncs_bas_1*(ncs_bas_1+1)/2
c     endif
      if(num_bas_ij.eq.num_bas_2) then
         ish_b=ncs_bas_1+1
         ish_e=ncs_bas_2
         jsh_b=ncs_bas_1+1
         jsh_e=ncs_bas_2
         n_shell=ncs_bas_2 - ncs_bas_1
         ij_shell=n_shell*(n_shell+1)/2
      endif
c kl 
c     if(num_bas_kl.eq.num_bas_1) then   ! take care of compiler warnings
         ksh_b=1
         ksh_e=ncs_bas_1
         lsh_b=1
         lsh_e=ncs_bas_1
         kl_shell=ncs_bas_1
c     endif
      if(num_bas_kl.eq.num_bas_2) then
         ksh_b=ncs_bas_1+1
         ksh_e=ncs_bas_2
         lsh_b=ncs_bas_1+1
         lsh_e=ncs_bas_2
         n_shell=ncs_bas_2 - ncs_bas_1
         kl_shell=n_shell
      endif
c
      nquartets=ij_shell*kl_shell
c-------------------------
c
      integrals=0
      ncalls=0
      nqrt=0
      ijkl=0
      ish=0
      do 100 ish1=ish_b,ish_e
      ish=ish+1
        jsh=0
        do 200 jsh1=jsh_b,ish1
        jsh=jsh+1
          ksh=0
          do 300 ksh1=ksh_b,ksh_e
          ksh=ksh+1
            do 400 lsh=0,0
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
     *                     more_int,blscr,l_blscr,0.0d0,'der2_int')
cc   *                     more_int,blscr,l_blscr,0.0d0,'der1_int')
c
c       call print_ij_k(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
                     call check_sijk(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
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
     *                     more_int,blscr,l_blscr,0.0d0,'der2_int')
ccc  *                     more_int,blscr,l_blscr,0.0d0,'der1_int')
c       call print_ij_k(ncalls,eri,leri,icf,jcf,kcf,lcf,integ_n0 )
                     call check_sijk(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
                     integrals=integrals+integ_n0
                     if(more_int) go to 452
                  nqrt=0
               else
                  nqrt=1
                  go to 410
               endif
  400       continue
  300     continue
  200   continue
  100 continue
c
c-------------------------
      write(6,*)'-------------------------------------------------'
      write(6,*)'       Three-center two-electron integrals    '
      write(6,*)'total number of integrals processed =',integrals
      write(6,*)'total number of integrals checked in=',integ_check
      write(6,*)'-------------------------------------------------'
      call print_check(stsum)
c-------------------------
      call txs_second(tcal2)
      write(6,66) tcal2-tcal1
   66 format('Total CPU check-sum time =',f10.2)
      write(6,*)'-------------------------------------------------'
c--------------------------------------------------------------
c
      write(6,*) '*** LEAVING TEST_CENT3 ij_k ***'
c
c--------------------------------------------------------------
c
      end
c=================================================================
      subroutine check_sijk(integ_n0,icf,jcf,kcf,lcf,eri,stsum)
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
c>>>>   l=lcf(ii)
        l=0
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
      subroutine print_ij_k(ncall, eri,leri,icf,jcf,kcf,lcf,integ_n0 )
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
c>>>>   l1=lcf(ii)
c 
             write(6,44) i1,j1,k1,eri(ii)
        enddo
   44 format('< ',i3,1x,i3,' |',i3,' >',3x,f30.20 )
c
      end
c=====================================================================
