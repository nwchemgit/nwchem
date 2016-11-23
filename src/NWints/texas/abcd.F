* $Id$                        
      subroutine abcd_1(nbls1,lci,lcj, lcij,lckl,
     *                  indxij,indxkl,index, nij_uniqe,nkl_uniqe,
     *                  abcd,apb,rcpd,cpd,rapb,
     *                  aaa,bbb,const,gcoef,ngcd)
c--------------------------------------------------------------------
c  OUTPUT
c  ABCD(ijkl)= apb*rcpd   OR  cpd*rapb abd logical stable(ijkl)
c--------------------------------------------------------------------
c because of possible numerical instablity in Tracy recursive 
c it is needed to check it out (see abcd_2) 
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
c numerical stability :
      logical stable_2
      logical stable_kind_1,stable_kind_2
      double precision avexp
      common /stability_2/ stable_2,stable_kind_1,stable_kind_2   ! use here and in num_s..
      logical stable 
      common /tracy_stability/ stable(5000) ! this is limit limblks
c
      common /allow_lost/ lost_allow
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      dimension indxij(*),indxkl(*),index(*)
c
      dimension apb(nij_uniqe,*),rapb(nij_uniqe,*)
      dimension cpd(nkl_uniqe,*),rcpd(nkl_uniqe,*)
      dimension aaa(nkl_uniqe,*),bbb(nij_uniqe,*)
      dimension const(nbls1),gcoef(ngcd,*)
      dimension abcd(nbls1)
c---------------------------------------------------------------
c For test only : 
c counter of all & unstable blocks of primitive quartets
c used here and in texas_face (init & terminate)
c
      common /count_primitive/ x_prim_all,
     *                         x_prim_wr1, x_prim_wr2, x_prim_wr3
c---------------------------------------------------------------
      x_prim_all=x_prim_all + nbls1
c---------------------------------------------------------------
ctest ONLY
c        stable_2=.false.
c        do i=1,nbls1
c           ijkl=index(i)
c           ijpar=indxij(ijkl)      
c           klpar=indxkl(ijkl)     
c           abcd(i)=rapb(ijpar,lcij)*cpd(klpar,lckl)
c           stable(i)=.false.
c        enddo
c
c      RETURN
c
ctest ONLY
c---------------------------------------------------------------
      if(stable_2) then
ctest
c       write(6,*)'abcd_1: stable_2=',stable_2,
c    *            ' kind_1,_2=',stable_kind_1,stable_kind_2,
c    *            ' a.m. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
ctest
         do i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)      
            klpar=indxkl(ijkl)     
            abcd(i)=apb(ijpar,lcij)*rcpd(klpar,lckl)
            stable(i)=.true.
ctest
c           avexp=aaa(ijpar,lci)
c           b_exp=bbb(ijpar,lcj)
c           write(6,*)' a,b=',0.5*avexp,0.5*b_exp
c           write(6,*)' a+b=',apb(ijpar,lcij),' c+d=',cpd(klpar,lckl)
c           write(6,*)' abcd=',abcd(i)
ctest
         enddo
      else
ctest
c       write(6,*)'abcd_1: stable_2=',stable_2,
c    *            ' kind_1,_2=',stable_kind_1,stable_kind_2,
c    *            ' a.m. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
ctest
         x_prim_wr1=x_prim_wr1 + nbls1              ! test 
c
         do i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)      
            klpar=indxkl(ijkl)     
c
c make it stable first :
c
            stable(i)=.true.
            abcd(i)=apb(ijpar,lcij)*rcpd(klpar,lckl)
c
            if(.not.stable_kind_1) then
               avexp=aaa(ijpar,lci)
cedo`               b_exp=bbb(ijpar,lcj)
cedo               boa=2.d0*b_exp/max(2.d0,a_exp) ! exp. are already *2 .
               if(avexp.gt.2.0d0) then
                 boa=2.0d0*bbb(ijpar,lcj)/avexp ! exp. are already *2 .
               else
                 boa=bbb(ijpar,lcj)
               endif
ctest
c              write(6,*)' a,b=',0.5*a_exp,0.5*b_exp
c              write(6,*)' a+b=',apb(ijpar,lcij),' c+d=',cpd(klpar,lckl)
ctest
               if(boa .ge. 10.d0) then
                  x_prim_wr2=x_prim_wr2 + 1.d0         !   test
                  xost_digit_1= log10(boa*rcpd(klpar,lckl))
                  xost_digit  = (nskl-1)*xost_digit_1 
                  lost_digit= nint( xost_digit )
                  if(ngcd.gt.1) then
                     call get_max_gcoef(gcoef(1,ijkl),ngcd,gc_max)
                     lost_const  =nint( log10( gc_max*abs(const(i)) ) )
                  else
                     lost_const  =nint( log10( abs( const(i) ) ) )
                  endif
                  lost_const = max( -14, lost_const )
                  lost_digit = lost_digit + lost_const
                  if(lost_digit.gt.lost_allow) then
                     x_prim_wr3=x_prim_wr3 + 1.d0                !   test
                     stable(i)=.false.
                     abcd(i)=rapb(ijpar,lcij)*cpd(klpar,lckl)
                  endif
ctest..................................................................
c                    write(6,*)'I st-kind stability: Quartet no=',i
c                    write(6,*)'const=',const(i) ,
c    *               ' a+b=',apb(ijpar,lcij),' c+d=',cpd(klpar,lckl)
c                    write(6,*)'Lost : in ',nskl-1,' steps =',
c    *                  int( xost_digit ),
c    *               '; due to const=',lost_const,
c    *               '; TOTAL=',lost_digit,'(allowed=',lost_allow,')'
c                    write(6,*)'==> stable=',stable(i),
c    *               ' a.m. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
ctest..................................................................
               endif       !    if(boa .ge. 10.d0) then
            endif          ! if(.not.stable_kind_1) then
c
            if(.not.stable_kind_2) then 
ctest
c              a_exp=aaa(ijpar,lci)
c              b_exp=bbb(ijpar,lcj)
c              write(6,*)' a,b=',0.5*a_exp,0.5*b_exp
c              write(6,*)'const=',const(i) ,
c    *         ' a+b=',apb(ijpar,lcij),' c+d=',cpd(klpar,lckl)
ctest
               if( stable(i)) then
                  abcd_max=abcd(i)
                  if(abcd_max.gt.  10.d0) then
                     x_prim_wr2=x_prim_wr2 + 1.d0         !   test
                     abcd_max=abcd(i)
                     xost_digit_1= log10(abcd_max) ! rjh was 0.5*
                     xost_digit  = abs(nskl-1)*xost_digit_1 
c                    xost_digit  = abs(nskl-3)*xost_digit_1 ! rjh - BAD
                     if(ngcd.gt.1) then
                        call get_max_gcoef(gcoef(1,ijkl),ngcd,gc_max)
                        lost_const  =nint(log10(gc_max*abs(const(i))))
                     else
                        lost_const  =nint(log10(abs(const(i))))
                     endif
                     lost_const = max( -14, lost_const )
                     lost_digit= nint( xost_digit ) + lost_const
                     if(lost_digit.gt.lost_allow) then
                        x_prim_wr3=x_prim_wr3 + 1.d0                !   test
                        stable(i)=.false.
                        abcd(i)=rapb(ijpar,lcij)*cpd(klpar,lckl)
                     endif
ctest..................................................................
c                       write(6,*)'IIed-kind stability: Quartet no=',i
c                       write(6,*)'const=',const(i) ,
c    *                  ' a+b=',apb(ijpar,lcij),' c+d=',cpd(klpar,lckl)
c                       write(6,*)'Lost : in ',nskl-1,' steps =',
c    *                     int( xost_digit ),
c    *                  '; due to const=',lost_const,
c    *                  '; TOTAL=',lost_digit,'(allowed=',lost_allow,')'
c                       write(6,*)'==> stable=',stable(i),
c    *                  ' a.m. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
ctest..................................................................
                  endif
               endif
            endif          ! if(.not.stable_kind_2) then
         enddo
      endif
c---------------------------------------------------------------
      end
c====================================================================
      subroutine abcd_2(nbls1, lcij,lckl,index,
     *                  abcd,apb,rcpd,cpd,rapb,
     *                  const,gcoef,ngcd)
c---------------------------------------------------------------
c const_max is maximum of abs values of const(nbls1) 
c OUTPUT : abcd & numerical stability info stable(1)
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension index(*)
c
      dimension apb(*),rapb(*)       ! apb(lc12),rapb(lc12)
      dimension cpd(*),rcpd(*)       ! cpd(lc34),rcpd(lc34)
      dimension const(nbls1), gcoef(ngcd,*)
c---------------------------------------------------------------
c accuracy :
      logical stable_2
      logical stable_kind_1,stable_kind_2
      common /stability_2/ stable_2,stable_kind_1,stable_kind_2 ! use here and in num_s..
      common /stability_3/ b_over_a(1000) 
c
c used only here and in trobsa.f :
c
      logical stable 
      common /tracy_stability/ stable(5000) ! this is limit limblks
      common /allow_lost/ lost_allow
c---------------------------------------------------------------
c For test only : counter of all & unstable blocks of primitive quartets
c used here and in texas_face (init & terminate)
c
      common /count_primitive/ x_prim_all,
     *                         x_prim_wr1, x_prim_wr2, x_prim_wr3
c---------------------------------------------------------------
c This routine tries to determine stability for Tracy's recursive
c for the WHOLE block , not for individual quartets (as it is done
c for iroute=1 in abcd_1 ).
c---------------------------------------------------------------
c for accuracy in Tracy's recursive : 
c final integrals are obtained by multiplying pruducts from Tracy's rec.
c by the CONST factor (involving contr.coeff.) which can be very small.
c In such cases the final integrals are also very small and we may
c allow ourselves to loose some digit. In the other words, small value 
c of const may let us "get back" some lost digit. However, we can not 
c get back more than 14 digit ( double precission limit ). Thus :
c---------------------------------------------------------------
      x_prim_all=x_prim_all + nbls1         ! test only
c---------------------------------------------------------------
ctest ONLY
c        stable_2=.false.
c        abcd=rapb(lcij)*cpd(lckl)
c        stable(1)=.false.
c
c      RETURN
c
ctest ONLY
c---------------------------------------------------------------
      if(stable_2) then
ctest
c        write(6,*)'abcd_2: stable_2=',stable_2,
c    *             ' kind_1,_2=',stable_kind_1,stable_kind_2,
c    *             ' a.m. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
c        write(6,*)' a+b=',apb(lcij),' c+d=',cpd(lckl)
ctest
         abcd=apb(lcij)*rcpd(lckl)
         stable(1)=.true.
      else
ctest
c        write(6,*)'abcd_2: stable_2=',stable_2,
c    *            ' kind_1,_2=',stable_kind_1,stable_kind_2,
c    *            ' a.m. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
c        write(6,*)' a+b=',apb(lcij),' c+d=',cpd(lckl)
ctest
         x_prim_wr1=x_prim_wr1 + nbls1              ! test 
c
c make it stable first :
c
         stable(1)=.true.
         abcd=apb(lcij)*rcpd(lckl)
c
         if(.not.stable_kind_1) then
            boa=b_over_a(lcij)
            if(boa .ge. 10.d0) then
               x_prim_wr2=x_prim_wr2 + nbls1        !   test
cedo               xost_digit_1= log10(boa*rcpd(lckl))
               xost_digit_1= log10(boa)+log10(rcpd(lckl))
               xost_digit  = (nskl-1)*xost_digit_1 
               lost_digit= nint( xost_digit )
          call get_max_const(nbls1,const,ngcd,index,gcoef,const_max)
               lost_const=nint( log10(const_max) )
               lost_const = max( -14, lost_const )
               lost_digit = lost_digit + lost_const
               if(lost_digit.gt.lost_allow) then
                  x_prim_wr3=x_prim_wr3 + nbls1               !   test
                  stable(1)=.false.
                  abcd=rapb(lcij)*cpd(lckl)
               endif
ctest..................................................................
c              write(6,*)'I st-kind stability: const_max=',const_max
c              write(6,*)'Lost : in ',nskl-1,' steps =',
c    *         int( xost_digit ),
c    *         '; due to const=',lost_const,
c    *         '; TOTAL=',lost_digit,'(allowed=',lost_allow,')'
c              write(6,*)'==> stable=',stable(1),
c    *         ' a.m. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
ctest..................................................................
            endif       !    if(boa .ge. 10.d0) then
         endif          ! if(.not.stable_kind_1) then
c
         if(.not.stable_kind_2) then 
            if( stable(1)) then
               abcd_max=abcd
               if(abcd_max.gt.  10.d0) then
                  x_prim_wr2=x_prim_wr2 + 1.d0         !   test
                  xost_digit_1= log10(abcd_max) ! rjh was 0.5*
                  xost_digit  = (nskl-1)*xost_digit_1 
c                  xost_digit  = (nskl-3)*xost_digit_1 ! RJH BAD
             call get_max_const(nbls1,const,ngcd,index,gcoef,const_max)
                  lost_const=nint( log10(const_max) )
                  lost_const = max( -14, lost_const )
                  lost_digit= nint( xost_digit ) + lost_const
                  if(lost_digit.gt.lost_allow) then
                     x_prim_wr3=x_prim_wr3 + nbls1               !   test
                     stable(1)=.false.
                     abcd=rapb(lcij)*cpd(lckl)
                  endif
ctest..................................................................
c                 write(6,*)'IIed-kind stability: const_max=',const_max 
c                 write(6,*)'Lost : in ',nskl-1,' steps =',
c    *            int( xost_digit ),
c    *            '; due to const=',lost_const,
c    *            '; TOTAL=',lost_digit,'(allowed=',lost_allow,')'
c                 write(6,*)'==> stable=',stable(1),
c    *            ' a.m. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
ctest..................................................................
               endif
            endif
         endif          ! if(.not.stable_kind_2) then
      endif
c---------------------------------------------------------------
c old fashion way :
c
c     if(nsij.ge.nskl) then
c        stable(1)=.true.
c        abcd=apb(lcij)*rcpd(lckl)
c     else
c        stable(1)=.false.
c        abcd=rapb(lcij)*cpd(lckl)
c     endif
c---------------------------------------------------------------
c
      end
c====================================================================
