      Subroutine xc_fitv( AO_bas_han, XC_bas_han,geom,
     &                   Tmat,nTrows,nTcols,
     &                   idirect,iga_vxc,iga_xcinv,IOLGC)

C$Id: xc_fitv.F,v 1.3 1995-04-27 03:16:56 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

c****
c**** nwchem handles
c****
      integer AO_bas_han, XC_bas_han,geom
c**** 

      Logical LResult,  IOLGC

      Parameter (itmp=98)

#include <bas.fh>

#include <mafdecls.h>
#include <global.fh>
#include <tcgmsg.fh>
#include <msgids.fh>

#if defined(DFT_TIMER)
#include <timer.fh>
#endif
c
c--> XC Energy & Potential
c
      dimension iga_vxc(4)
c
c--> T matrix
c
      Dimension Tmat(nTrows,nTcols)
c
c Compute the fitted XC energy and potential. The fitting coefficients are
c obtained by solving the normal equations for a standard linear least
c squares problem.
c
c                          _
c               ~   i     \          /      
c              (Vxc)   =  /_   C     | chi  chi  chi  dr
c                   pq       s  s,i  /    p    q    s
c                        _
c                       \     -1
c               C    =  /_   S    T
c                s,i       r  s,r  r,i
c
c                       /       i
c               T    =  | chi  V    dr
c                r,i    /    r  xc
c
c                                INDICES
c                               ---------
c                      pq --> product of AO basis functions (1,...,ntri_ao)
c                  r or s --> XC LSQ fitting functions (1,...,nbf_xc)
c                       i --> energy and potential index (1,...,nTcols)
c
c
c                      Energy and Potential Index Key
c
c              Value of     |     Definition of index "i"
c            ipol     nExc  |    1        2        3       4
c           --------------------------------------------------
c              1        1   |   Vxc      Exc   
c              2        1   |   Vxc^up   Vxc^dw   Exc
c              1        2   |   Vxc      Ex       Ec
c              2        2   |   Vxc^up   Vxc^dw   Ex       Ec    
c
c           nTcols = ipol + nExc
c
c******************************************************************************
#if defined(DFT_TIMER)
      call ga_sync
      call t_stat(10,-1)
#endif
c
      me=ga_nodeid()
      nproc = ga_nnodes()
c
c     Determine the characteristics of the AO and XC Gaussian basis sets.
c

      LResult = .true.
      LResult = LResult .and. bas_numcont( AO_bas_han, nshells_ao) 
      LResult = LResult .and. bas_numbf( AO_bas_han, nbf_ao)
      LResult = LResult .and. bas_nbf_cn_max( AO_bas_han, MXnshbf_ao)
      LResult = LResult .and. bas_numcont( XC_bas_han, nshells_XC) 
      LResult = LResult .and. bas_numbf( XC_bas_han, nbf_XC)
      LResult = LResult .and. bas_nbf_cn_max( XC_bas_han, MXnshbf_XC)
      if ( .not. LResult) then
        call errquit('Exit from fitvxc',1)
      endif 
c
c     Compute the XC fitting coefficients.
c
#if defined(DFT_TIMER)
      call dft_second(ts1)
#endif
c
      ndbl = nbf_xc*nTcols
      LResult = MA_Push_Get(MT_Dbl,ndbl,'xc_coef',lxc_coef,ixc_coef)
      LResult = MA_Push_Get(MT_Dbl,nbf_xc,'scratchinv',lscr,iscr)

      call dfill(ndbl,0.d0,DBL_MB(ixc_coef),1)

      if(IOLGC)             
     &     call gacrea(
     &     geom,XC_bas_han,nbf_xc,nbf_xc,'XCinv',iga_xcinv,'row')
      if(me.eq.0.and.IOLGC) then
        call dft_invio('XC',iga_xcinv,nbf_xc,'READ',lmiss)
      endif
      call ga_sync
c
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(1)='  XC READ '
      call tcllct(t_sub(1,1),ts2,ts1)
      call dft_second(ts1)
#endif
c
      next=nxtval(nproc)+1

      do k=1,nbf_xc
        if (k.eq.next) then
          next=nxtval(nproc)+1
c          call ga_get(iga_xcinv,
c     *         1,nbf_xc,k,k,Dbl_MB(iscr),1)
          call get_col(iga_xcinv,nbf_xc,k,Dbl_MB(iscr))
          istart=ixc_coef
          do j=1,nTcols
            call daxpy(nbf_xc,Tmat(k,j),DBL_MB(iscr),1,DBL_MB(istart),1)
            istart=istart+nbf_xc
          enddo
        endif
      enddo
      LResult = MA_pop_stack(lscr)
      next=nxtval(-nproc)
      call ga_dgop(Msg_XC_Coef,DBL_MB(ixc_coef),ndbl,'+')
      if(IOLGC)  call gadest(iga_xcinv)
c
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(2)='  xc_coef '
      call tcllct(t_sub(1,2),ts2,ts1)
      call dft_second(ts1)
#endif
c
c     Compute the fitted Vxc.
c
c
c     Allocate scratch.
c
      n3OI = (MXnshbf_ao**2)*MXnshbf_xc
      LResult = MA_Push_Get(MT_Dbl,n3OI,'3OI',l3OI,i3OI)
      LResult = MA_Push_Get(MT_Dbl,nTcols*n3OI,'vvxc',lvvxc,ivvxc)
      LResult = MA_Push_Get(MT_Dbl,n3OI,'vec2',lvec2,ivec2)
      call int_mem_3ov(mag,Max_scr)
      LResult = MA_Push_Get(MT_Dbl,Max_scr,'scr2',lscr2,iscr2)
c
c     Compute 3-ctr OIs between a given pair of AO shells and a set of XC shells.
c     
        call ga_sync
        n1 = nxtval(nproc)
        n2 = 0
        do 240 ishb = 1,nshells_ao
          if (.not. bas_cn2bfr( AO_bas_han,ishb,
     &         ifirstb,ilastb))
     &         call errquit('Exiting in fitvxc.',10)
          nshbfb=ilastb-ifirstb+1
          do 230 ishc = 1,ishb
            if( n1.eq.n2 )then
              if( abs(idirect).eq.1 )then
                if (.not. bas_cn2bfr( AO_bas_han,ishc,
     &               ifirstc,ilastc))
     &               call errquit('Exiting in fitvxc.',10)
                nshbfc=ilastc-ifirstc+1
                Nao2 = nshbfb*nshbfc
                iixc=0
                do ii=0,nTcols*Nao2-1
                  DBL_MB(ivvxc+ii)=0.d0
                enddo
                do 210 isha = 1,nshells_xc
                  if (.not. bas_cn2bfr( XC_bas_han,isha,
     &                 ifirsta,ilasta))
     &                 call errquit('Exiting in fitvxc.',12)
                  nshbfa=ilasta-ifirsta+1
                  
                  Nint = Nao2*nshbfa
                  nscr2=max_scr
                  call int_1e3ov( XC_bas_han, isha, AO_bas_han, ishb,
     &                 AO_bas_han, ishc,  
     &                 nscr2, Dbl_MB(iscr2),
     &                 Nint, Dbl_MB(i3OI) )
c      
c      Multiply the matrix of 3-ctr OIs by the XC expansion coefficients.
c      
                  do ii=0,nTcols-1
                    nxc=ixc_coef+ii*nbf_xc+iixc
                    npint=i3OI
                    istart=ivvxc+ii*Nao2
                    do ia=0,nshbfa-1
                      fact=DBL_MB(nxc+ia)
                      call daxpy(
     &                     Nao2,fact,DBL_MB(npint),1,DBL_MB(istart),1)
                      npint=npint+Nao2
                    enddo
                  enddo
                  iixc=iixc+nshbfa
 210            continue
c      
c      Scatter the computed values into Vxc.
c      
                do 220 k = 1,nTcols
                  call sca_lab_mat(AO_bas_han,
     &                 ishb,ishc,nshbfb,nshbfc,
     &                 iga_vxc(k),DBL_MB(ivvxc+(k-1)*Nao2),
     &                 dbl_mb(ivec2))
 220            continue
              endif
              
              n1 = nxtval(nproc)
              
            end if
            
            n2 = n2 + 1
            
 230      continue
 240    continue
        call ga_sync
c
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(3)='  Vxc     '
      call tcllct(t_sub(1,3),ts2,ts1)
      call dft_second(ts1)
#endif
c
      LResult = MA_Pop_Stack(lscr2)
      LResult = MA_Pop_Stack(lvec2)
      LResult = MA_Pop_Stack(lvvxc)
      LResult = MA_Pop_Stack(l3OI)
      LResult = MA_Pop_Stack(lxc_coef)

      n1 = nxtval(-nproc)
c
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(4)='  DGOP    '
      call tcllct(t_sub(1,4),ts2,ts1)
#endif
c
#if defined(DFT_TIMER)
      if( me.eq.0 )then
        ltime = 7
        write(ltime,*)
        write(ltime,9996) '      Analysis of FITVXC'
        write(ltime,*)
      end if
      subname=' Fitvxc   '
      call t_stat(4,1)
 9996 format(10x,a)
#endif
c
      return
      end
