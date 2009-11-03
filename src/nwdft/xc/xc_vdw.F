c     $Id: xc_vdw.F,v 1.6 2008-10-02 17:27:35 edo Exp $
      subroutine xc_vdw(rtdb,geom,exc,force,what)
      implicit none
      character *(*) what
      integer geom,rtdb
      double precision exc,force(*),s6
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "xc_vdw.fh"
c
      integer n
      integer itags,ltags,i_xyz,l_xyz,icharge,lcharge,
     I     l_fvdw,i_fvdw
      external xc_vdw_e
      double precision xc_vdw_e,evdw,scalea
      integer iz,lz,i
      logical xc_vdw_init
      external xc_vdw_init
      logical oprint,oprinth
c
      double precision get_scale6
      external get_scale6
c
      double precision delta,delta_default
c
      oprint = util_print('vdw', print_medium)
      oprinth = util_print('vdw high', print_high)
c
      if (.not. geom_ncent(geom, n)) 
     &   call errquit('xcvdw: geom_ncent failed',geom, GEOM_ERR)
      if (.not.MA_push_get(MT_Dbl,n*3,'xyz',l_xyz,i_xyz))
     &   call errquit('xcvdw: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_int,n,'atns',lz,iz))
     &   call errquit('xcvdw: cannot allocate atns',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,n,'charge',lcharge,icharge))
     &   call errquit('xcvdw: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,n*16,'tags',ltags,itags))
     &   call errquit('xcvdw: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, n, Byte_MB(itags),
     &                        Dbl_MB(i_xyz), Dbl_MB(icharge)))
     &   call errquit('xcvdw: geom_cart_get failed',74, GEOM_ERR)
      if (.not.ma_pop_stack(ltags))
     &   call errquit('xcvdw: cannot pop stack',2, MA_ERR)
c
c
      do i=1,n
          int_mb(iz+i-1)=int(dbl_mb(icharge+i-1))
      enddo
c
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('xcvdw: cannot pop stack',3, MA_ERR)

      if (.not.rtdb_get(rtdb,'dft:ivdw',MT_INT,1,ivdw))
     &      ivdw = 2
c
c     conversion factor angs 2 au
c
       if(.not.geom_get_ang2au(geom, scalea)) call
     S     errquit('xcvdw: gang2au failed',0,0)
c     write(6,*)'angs to au1', scalea
      if(.not.xc_vdw_init(scalea))
     &   call errquit('xcvdw: vwdinit failed',0, 0)
c
      if (.not.rtdb_get(rtdb,'dft:vdwalpha',MT_DBL,1,alpha))
     &      alpha = 20.0d0
c
      if (.not.rtdb_get(rtdb, 'dft:vdw', mt_dbl, 1, s6))
     &   s6=0.0d0
c
      if(s6.eq.0.0d0) s6 = get_scale6()
c
      if(what.eq.'energy') then
         if(oprinth.and.ga_nodeid().eq.0) then
           write(luout,*) ' s6 =',s6
           write(luout,*) ' alpha  =',alpha
           write(luout,*) ' ivdw  =', ivdw
         endif
         evdw=xc_vdw_e(s6,n,dbl_mb(i_xyz),int_mb(iz))
         if(oprint.and.ga_nodeid().eq.0) then
              write(luout,222) s6, evdw
 222  format(//
     &     '          s6 scale factor =', f22.12/
     &     '              vdW contrib =', f22.12//)
         endif
         Exc=Exc+evdw
      elseif(what.eq.'forces') then
      if (.not.MA_push_get(MT_Dbl,n*3,'xyz',l_fvdw,i_fvdw))
     &   call errquit('xcvdw: cannot allocate forcev',0, MA_ERR)
         call xc_vdw_der(s6,n,dbl_mb(i_xyz),int_mb(iz),
     D     dbl_mb(i_fvdw))
         if(oprinth.and.ga_nodeid().eq.0) then
             write(luout,*) ' vdW contrib for S6=',s6
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               dbl_mb(i_fvdw+3*(i-1)),
     Y               dbl_mb(i_fvdw+3*(i-1)+1),
     Z        dbl_mb(i_fvdw+3*(i-1)+2)
             enddo
             write(luout,*) ' before vdw contr @@@@@'
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               force(1+3*(i-1)),
     Y               force(1+3*(i-1)+1),
     Z               force(1+3*(i-1)+2)
             enddo

          endif
          call daxpy(3*n,1d0,dbl_mb(i_fvdw),1,force,1)
          if(oprinth.and.ga_nodeid().eq.0) then
             write(luout,*) ' after vdw contr @@@@@'
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               force(1+3*(i-1)),
     Y               force(1+3*(i-1)+1),
     Z               force(1+3*(i-1)+2)
             enddo
          endif
      elseif(what.eq.'hessian') then
c     get delta as used in a numerical hessian DFT calculation
      delta_default  =  0.01d0
      if (.not.rtdb_get(rtdb,'stpr_gen:delta',MT_DBL,1,delta))
     &      delta = delta_default 
         call xc_vdw_hess(delta,s6,n,dbl_mb(i_xyz),int_mb(iz))
         if(oprint.and.ga_nodeid().eq.0)  then
              write(luout,*) ' s6 = ',s6
              write(luout,*) ' vdw to hessian contribution is done'
         endif
      endif
      if (.not.ma_chop_stack(l_xyz)) 
     C   call errquit('xcvdw: cannot pop stack',4, MA_ERR)
      return
      end
      subroutine xc_vdw_der(s6,n,x,z,force)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2693 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
#include "xc_vdw.fh"
      double precision s6
      integer n
      double precision x(3,*),force(3,*)
      integer z(*)
c
      integer i,j,A
      double precision c6ij_sk,rAj
      external c6ij_sk
      double precision drajdxa,r0aj
      double precision ff1,rr,ff
      double precision fdmp,f1dmp
c
c     deriv
c
      do A=1,n
         force(1,A)=0d0
         force(2,A)=0d0
         force(3,A)=0d0
         do j=1,n
            if(A.ne.j) then
               rAj=sqrt(
     +              (x(1,A)-x(1,j))**2 +
     +              (x(2,A)-x(2,j))**2 +
     +              (x(3,A)-x(3,j))**2)
               r0aj=r0(z(A))+r0(z(j))
               ff= fdmp(rAj,r0aj)
               ff1= f1dmp(rAj,r0aj,ff)
               rr=c6ij_sk(A,j,z)/(rAj**6)*
c    *             (-6d0*ff/rAj+ff1)
     *             ((-6d0*ff/rAj)+ff1)
               do i=1,3
c                 drAjdxa=(1.0d0/rAj)*(x(i,A)-x(i,j))
                  drAjdxa=(x(i,A)-x(i,j))/rAj
                  force(i,A)=force(i,A)-rr*drAjdxa
               enddo
            endif
         enddo
      enddo
c
c     fudge factor if needed
c      
      if(abs(s6-1d0).gt.1d-9) 
     F     call dscal(3*n,s6,force,1)
#ifdef DEBUG
      write(6,*) ' gradient vdw called'
#endif
      return
      end
      double precision function c6ij_sk(i,j,z)
      implicit none
#include "xc_vdw.fh"
      integer i,j
      integer z(*)
c     
      integer zi,zj
      double precision ni,nj
      external xc_neff
      double precision xc_neff,ci,cj
c
      zi=z(i)
      zj=z(j)
      if(c6(zi).eq.-1d0.or.c6(zj).eq.-1d0) then
         write(6,*) ' i j',i,j,' c6 ',c6(zi),c6(zj)
         stop 'cra'
      endif
      ci=c6(zi)
      cj=c6(zj)
c
      c6ij_sk=0.0d0
      if(ivdw.eq.1)then
c version 1 
c     S. Grimme J Comp Chem 25, 1463 (2004)
        c6ij_sk =2.d0*ci*cj/(ci+cj)
        else if(ivdw.eq.2) then
c version 2
c     S. Grimme J Comp Chem 27, 1787 (2006)
          c6ij_sk =sqrt(ci*cj)
        else if(ivdw.eq.0) then
          ni=xc_neff(zi)
          nj=xc_neff(zj)
          c6ij_sk=2d0*(ci*ci*cj*cj*ni*nj)**(1d0/3d0)/
     /       ((nj*nj*ci)**(1d0/3d0)+(ni*ni*cj)**(1d0/3d0))
      endif
      return
      end
      double precision function xc_neff(z)
      implicit none
      integer z
c
      integer nv
c
      if(z.eq.1) then
         xc_neff=0.8d0
         return
      endif
      if(z.le.2) nv=z
      if(z.le.10) nv=z-2
      if(z.gt.10) stop 'neff stop '
      xc_neff=1.17d0+0.33d0*nv
      return
      end
      double precision function xc_vdw_e(s6,n,x,z)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2693 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
#include "xc_vdw.fh"
      double precision s6
      integer n
      double precision x(3,*)
      integer z(*)
c
      integer i,j
      double precision fdmp
      double precision c6ij_sk,rij
      external c6ij_sk
c
      xc_vdw_e=0d0
      do i=1,n-1
         do j=i+1,n
            rij=sqrt(
     +           (x(1,i)-x(1,j))**2 +
     +           (x(2,i)-x(2,j))**2 +
     +           (x(3,i)-x(3,j))**2)
            xc_vdw_e=xc_vdw_e-c6ij_sk(i,j,z)*
     *           fdmp(rij,r0(z(i))+r0(z(j)))*
     *           (rij)**(-6)
         enddo
      enddo
c
c     fudge factor if needed
c      
      xc_vdw_e=xc_vdw_e*s6
!      write(0,*) ' xc_vdw in kcal/mol ',xc_vdw_e
!      xc_vdw_e=xc_vdw_e/au2kcalmol
      return
      end
      logical function xc_vdw_init(angs2au)
      implicit none
#include "xc_vdw.fh"
      integer ii
      double precision angs2au
      double precision scalerad
      double precision r0_grimme1(54)
      double precision r0_grimme2(54)
c
c VDW radii are in ang. determined by
c atomic ROHF/TZV calculations and then taken as the radius 
c of the 0.01 density contour
c (scaling done below)
c data obatained using the 0.005 contour and scale=1.0 
c gave too large differences
c between the atoms i.e 1.1 1.5 1.42 1.36 1.3 for H,C,N,O,F
      data r0_grimme2
c H, He
     .     /0.91d0,0.92d0,
c Li-Ne
     .      0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0,
c Na-Ar
     .      1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0,
c K, Ca
     .      1.35d0,1.34d0,
c Sc-Zn
     .      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,
     .      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,
c Ga-Kr
     .      1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0,
c Rb, Sr
     .      1.48d0,1.46d0,
c Y-Cd
     .      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,
     .      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,
c In, Sn, Sb, Te, I, Xe
     .      1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0/
c    .      46*1.d+99/
c
      data r0_grimme1
c Original data radii from Grimme JCC(2004) now angstroms
c only for H,C,N,O,F,Ne
c H, He
     .     /1.11d0,0.0d0,
c Li-Ne
     .      0.00d0,0.00d0,0.00d0,1.61d0,1.55d0,1.49d0,1.43d0,1.38d0,
c Na-Ar
     .      0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
c K, Ca
     .      0.00d0,0.00d0,
c Sc-Zn
     .      0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
     .      0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
c Ga-Kr
     .      0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
c Rb, Sr
     .      0.00d0,0.00d0,
c Y-Cd
     .      0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
     .      0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
c In, Sn, Sb, Te, I, Xe
     .      0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0/
c    .      46*1.d+99/

c
      integer i
      double precision jnm6mol2au
c
      do i=1,128
         c6(i)=-1
         r0(i)=-1
      enddo
c
      if(ivdw.eq.1.or.ivdw.eq.2)then
c     if(ivdw.eq.2)then
        do i=1, 54
           r0(i)=r0_grimme2(i)
        enddo
c     else if(ivdw.eq.1.or.ivdw.eq.0)then
      else if(ivdw.eq.0)then
        do i=1, 54
          r0(i)=r0_grimme1(i)
        enddo
      endif
c
c to avoid divide by zero for dummy atoms (CP correction)
      call dfill(c6,1.0d-14,128,1)
c
      if(ivdw.eq.1)then
c scaled radii
c in J mol^-1 nm^6
        c6(1) =0.16d0
        c6(5) =1.65d0
        c6(6) =1.65d0
        c6(7) =1.11d0
        c6(8) =0.70d0
        c6(9) =0.57d0
        c6(10)=0.45d0
c not tested
        c6(14)=25.00d0
c quick and dirty fit for ph3 dimer (MP2/aug-cc-pVTZ)
        c6(15)=20.00d0
c quick and dirty fit for h2s dimer
        c6(16)=12.00d0
c quick and dirty fit for hcl dimer
        c6(17)=9.00d0
c quick and dirty fit for Ar2 (exp. value is 3.7)
        c6(18)=6.50d0
c not really tested
        c6(2) =0.15D0
        c6(24)=10.0D0
        c6(25)=10.0D0
        c6(26)=10.0D0
c estimate for bromine
        c6(35)=18.0D0
      else if(ivdw.eq.2)then
c scaled radii
c      c metals as average of preceeding rare gas and group 3 element
c      def2-QZVP PBE0 IP and alpha
        c6( 1)=    0.14
        c6( 2)=    0.08
        c6( 3)=1.61
        c6( 4)=1.61
        c6( 5)=    3.13
        c6( 6)=    1.75
        c6( 7)=    1.23
        c6( 8)=    0.70
        c6( 9)=    0.75
        c6(10)=    0.63
        c6(11)=5.71
        c6(12)=5.71
        c6(13)=   10.79
        c6(14)=    9.23
        c6(15)=    7.84
        c6(16)=    5.57
        c6(17)=    5.07
        c6(18)=    4.61
        do ii=19,30
           c6(ii)=10.8
        enddo
        c6(31)=   16.99
        c6(32)=   17.10
        c6(33)=   16.37
        c6(34)=   12.64
        c6(35)=   12.47
        c6(36)=   12.01
        do ii=37,48
           c6(ii)=24.67
        enddo
        c6(49)=   37.32
        c6(50)=   38.71
        c6(51)=   38.44
        c6(52)=   31.74
        c6(53)=   31.50
        c6(54)=   29.99
      else if(ivdw.eq.0)then
c scaled radii
        c6(1)=0.163d0
        c6(6)=1.881d0
        c6(7)=1.11d0
        c6(8)=0.684d0
        c6(9)=0.57d0
        c6(10)=0.45d0
      endif
c
c     if(ivdw.eq.1.or.ivdw.eq.0)then
      if (ivdw.eq.0) then
        scalerad=1.00d0
      else if (ivdw.eq.1) then
c true ivdw eq 1
        scalerad=1.220d0
      else if (ivdw.eq.2) then
        scalerad=1.10d0
      else
        scalerad=0.0d0
      end if
c     convert to au
c     do i=1,128
c     j nm6/mol to  kcal nm6/mol  to au nm6 to au
c        c6(i)=c6(i)*joule2kcal/au2kcalmol/((angs2nm/angs2au)**6)
c        c6(i)=c6(i)*(((10*angs2au)**6)/4186)/627.51
c        c6(i) = 0.001d0 * c6(i) / 2625.4999d0 * 0.052917726d0
c picom to angs to au
c        r0(i)=r0(i)*angs2au
c     enddo
c

c use scaled radii and units change
      call dscal(128,angs2au*scalerad,r0,1)
c use  units change j nm6/mol to au
      jnm6mol2au=(((10*angs2au)**6)/4186)/627.51
      call dscal(128,jnm6mol2au,c6,1)
      xc_vdw_init=.true.
      return
      end
      double precision function get_scale6()
      implicit none
#include "cdft.fh"
#include "xc_vdw.fh"
      double precision scalec6
c old DFT-D version
      if( ivdw.eq.1) then
c BLYP
c  Becke88
c  LYP
        if(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
          scalec6 = 1.40d0
c PBE
        else if(xfac(7).eq.1.0d0.and.cfac(12).eq.1.0d0) then
         scalec6 = 0.70d0
c BP86
         else if(xfac(3).eq.1.0d0.and.xfac(2).eq.1.0d0.and.
     &              cfac(1).eq.1.0d0.and.cfac(4).eq.1.d0) then
         scalec6 = 1.30d0
c default
        else 
         scalec6 = 1.00d0
         write(6,*)'WARNING:'
         write(6,*)'you are using HF or an untested DF together'
         write(6,*)'with the dispersion correction!'
         endif
       else if ( ivdw.eq.2) then
c new DFT-D version
c BLYP
         if(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
           scalec6 = 1.20d0
c PBE
        else if(xfac(7).eq.1.0d0.and.cfac(12).eq.1.0d0) then
           scalec6 = 0.75d0
c TPSS
         else if(xfac(21).eq.1.0d0.and.cfac(27).eq.1.0d0) then
           scalec6 = 1.00d0
c B3LYP
         else if (xccomb(2)) then
           scalec6 = 1.05d0
c B97-D
         else if (xccomb(37)) then
           scalec6 = 1.25d0
c SSB-D
         else if (xccomb(46)) then
           scalec6 = 0.847455d0
c BP86
         else if(xfac(3).eq.1.0d0.and.xfac(2).eq.1.0d0.and.
     &              cfac(1).eq.1.0d0.and.cfac(4).eq.1.d0) then
          scalec6 = 1.05d0
c B2-PLYP
c        scalec6 = 0.55d0
         else
c default
           scalec6 =1.20d0
           write(6,*)'WARNING:'
           write(6,*)'you are using HF or an untested DF together'
           write(6,*)'with the dispersion correction!'
        endif
      else
        scalec6 = 0.0
      endif
      get_scale6 = scalec6 
      return
      end
      subroutine xc_vdw_to_hessian(rtdb)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb
      integer  geom
      double precision dum
      integer dum1
      character*255 name
c
      if (.not. rtdb_cget(rtdb,'geometry', 1, name))
     $     name = 'geometry'
c
      if (.not. geom_create(geom, name))
     $ call errquit('xc_vdw_to_hessian: geom_create failed !', 
     $                                                0,GEOM_ERR)
c
      if (.not. geom_rtdb_load(rtdb, geom, name))
     $ call errquit('xc_vdw_to_hessian: no geometry load form rtdb', 0,
     $        GEOM_ERR)
c
      call xc_vdw(rtdb, geom, dum, dum1, 'hessian')
c
      if(.not. geom_destroy(geom))
     $ call errquit('xc_vdw_to_hessian: geom_create failed !', 
     $                                                0,GEOM_ERR)
c
      return
      end
      subroutine xc_vdw_hess(delta,s6,n,x,z)
      implicit none
c This function makes vdw empirical correction to the hessian
c must be called before thermochemical data and vibrational
c analysis are generated.
c#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      double precision s6     
      integer n                
      double precision x(3,n)
      double precision z(n)
      double precision l_force(3,n)
      double precision r_force(3,n)
      double precision hessvdw(3,n,3,n)
      double precision dval0
      double precision dval1
      double precision delta
      integer i, j, A, B
      integer n3xyz 
      integer nat2
      integer nhesst
      integer lenhess 
      integer l_exybs,k_exybs 
      integer l_exybt,k_exybt
c
      character*(nw_max_path_len) filehess
c
c  dispersion contribution to hessian 
c
c in principle this task is very fast,
c so only master node works, and read-write to disk.
c 
      if (ga_nodeid().eq.0) then
c
      call util_file_name('hess', .false., .false.,filehess)

      lenhess=inp_strlen(filehess)
      n3xyz=3*n
      nhesst=n3xyz*(n3xyz+1)/2
      nat2=n3xyz*n3xyz
c
      if (.not.ma_push_get(mt_dbl,nat2,'xcvdwhess exybs ',
     *  l_exybs,k_exybs))
     &   call errquit('xcvdwhess: cannot allocate exybs',0, MA_ERR)
      call dfill(nat2,0.0d00,dbl_mb(k_exybs),1)
c
      if (.not.ma_push_get(mt_dbl,nhesst,'xcvdwhess exybt ',
     *  l_exybt,k_exybt))
     &   call errquit('xcvdwhess: cannot allocate exybt',0, MA_ERR)
      call dfill(nhesst,0.0d00,dbl_mb(k_exybt),1)
c
      write(LuOut,* ) 'Read old hessian file : ', filehess 
c
      call dfill(dbl_mb(k_exybt),0.0d0,nhesst,1)
      !write(6,* ) 'leee '
      open(unit=69,file=filehess,form='formatted',status='old',
     &    err=99900,access='sequential')
      do j = 0,(nhesst-1)
        read(69,*,err=99901,end=99902) dval0
        dbl_mb(k_exybt+j) = dval0
      !write(6,* ) 'dval ', dval0
      enddo
      close(unit=69,status='keep')
c
c     complete the square matrix from triangle values
c
      call trin2squa(dbl_mb(k_exybs),dbl_mb(k_exybt),n3xyz)
c
      write(LuOut,* ) 'vdW contribution to hessian '
c
c     call output(dbl_mb(k_exybs), 1, n3xyz, 1, n3xyz, n3xyz, n3xyz, 1)
c
       do A = 1, n
         do i = 1, 3
           do B=1, n
            do j=1, 3
               r_force(j,B)=0.0d0
               l_force(j,B)=0.0d0
            enddo
           enddo
c right displacement
            x(i,A) = x(i,A) + delta
            call xc_vdw_der(s6,n,x,z,r_force)
c left displacement
            x(i,A) = x(i,A) - 2.0d0*delta
            call xc_vdw_der(s6,n,x,z,l_force)
c back to original position
            x(i,A) = x(i,A) + delta
c
            do B=1,n
              do  j=1, 3
                dval1 = (r_force(j,B)-l_force(j,B)) / (2.0d00*delta)
                hessvdw(i,A,j,B)=dval1
              enddo
            enddo
         enddo
      enddo
c
c         call output(hessvdw,1,n3xyz,1,n3xyz,
c    &        n3xyz,n3xyz,1)
c
      call daxpy(nat2,1d0,hessvdw,1,dbl_mb(k_exybs),1)
c
c:write the final hessian
c
          call output(dbl_mb(k_exybs),1,n3xyz,1,n3xyz,
     &        n3xyz,n3xyz,1)
      call  stpr_wrt_fd_from_sq(dbl_mb(k_exybs),n3xyz,filehess)
      write(LuOut,* ) 'New hessian file vdw corrected has been 
     .                      written:', filehess 
      if (.not.ma_chop_stack(l_exybs)) 
     C   call errquit('xcvdwhess: cannot pop stack exybs',4, MA_ERR)

      endif
c
      call ga_sync
c
      return
99900 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('hess_read: error opening file: "filehess"',911)
99901 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('hess_read: error reading file: "filehess"',911)
99902 continue
      write(luout,*)'hess_file => ',filehess
      call errquit
     & ('hess_read: unexpected EOF when reading file: "filehess"',911)
      end
c
c      fdmp(x1,x2)=1.0d0/(1.0d0+exp(-alpha*((x1/x2)-1.0d0)))
c
      double precision function fdmp(x1,x2)
      implicit none
#include "xc_vdw.fh"
      double precision x1,x2
      fdmp=1.0d0/(1.0d0+exp(-alpha*((x1/x2)-1.0d0)))
      return
      end
c
c      f1dmp(x1,x2,ff)=(alpha/x2)*ff*ff*exp(-alpha*((x1/x2)-1.0d0))
c
      double precision function f1dmp(x1,x2,ff)
      implicit none
#include "xc_vdw.fh"
      double precision x1,x2,ff
      f1dmp=(alpha/x2)*ff*ff*exp(-alpha*((x1/x2)-1.0d0))
      return
      end
