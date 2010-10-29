      logical function jvltest(rtdb)
*     
*     $Id$
*      calling routine for CCSD(T) trials
*
      implicit none
      integer rtdb
c...    rtdb is handle to database
c
      call triple_so(rtdb)
c
      jvltest = .true.
c
      end

      subroutine triple_so(rtdb)
c
c...  first stab at a triple code for CCSD(T)
c...  using spin-orbital integrals sitting in core 
c...  for now a 4-dimensional array of <ij||ab> 
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer nn_almlof,nn_laguer
      parameter (nn_almlof=8, nn_laguer=200)
      integer rtdb,n_almlof,np_almlof(nn_almlof),
     1             n_laguer,np_laguer(nn_laguer)
c
      integer occ,virt,nbasis
      logical oexist
      integer l_rint,k_rint,l_t1,k_t1,l_t2,k_t2,l_orb,k_orb
      integer l_orbn,k_orbn
      double precision delta,energy
      integer i
      character*26 date
c
c...  check existence of required files
c
      inquire(file='T2',exist=oexist)
      if (.not. oexist)
     $     call errquit('jvltest: no T2 ', 0, UNKNOWN_ERR)
      inquire(file='ASOINTS',exist=oexist)
      if (.not. oexist)
     $     call errquit('jvltest: no ASOINTS ', 0, UNKNOWN_ERR)
      inquire(file='EVALS',exist=oexist)
      if (.not. oexist)
     $     call errquit('jvltest: no EVALS ', 0, UNKNOWN_ERR)
c
c... get basis size from T2
c
      open(9,file='T2',form='formatted')
      read(9,*) occ,virt
      nbasis = occ + virt
      close(9)
c
c..     integrals
c     double precision rint(nbasis,nbasis,nbasis,nbasis)
      if (.not. ma_push_get(mt_dbl, nbasis**4,'rint',l_rint, k_rint))
     $    call errquit('rint',0, MA_ERR)
c..    amplitudes 
c     double precision t1(occ,virt), t2(occ,occ,virt,virt)
      if (.not. ma_push_get(mt_dbl, occ**2*virt**2,'t2',l_t2, k_t2))
     $    call errquit('t2',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, occ*virt,'t1',l_t1, k_t1))
     $    call errquit('t1',0, MA_ERR)
c..    the orb  (orbital energies)
c     double precision orb(nbasis)
      if (.not. ma_push_get(mt_dbl, nbasis,'orb',l_orb, k_orb))
     $    call errquit('orb',0, MA_ERR)
c
      call get_INTS(dbl_mb(k_rint),nbasis)
      call get_T1(dbl_mb(k_t1),occ,virt)
      call get_T2(dbl_mb(k_t2),occ,virt)
      call get_EVALS(dbl_mb(k_orb),nbasis,occ,virt)
c
*     call calc_robert(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
*    1                dbl_mb(k_orb),nbasis,occ,virt)
c
*     call calc_pople(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
*    1                dbl_mb(k_orb),nbasis,occ,virt)
c
c...  now make modified (homo-lumo) orbital energies and 
c...  call all + cullen
c
*     if (.not. ma_push_get(mt_dbl, nbasis,'orbn',l_orbn, k_orbn))
*    $    call errquit('orbn',0)
*     call mod_EVALS(dbl_mb(k_orb),dbl_mb(k_orbn),nbasis,occ,virt)
c
*     print *,' ORBITAL ENERGIES TO HOMO-LUMO '
*     delta = dbl_mb(k_orb+occ-1) - dbl_mb(k_orb+occ)
*     print *,' delta e used ',delta
*     delta = delta * 3.0d0
c
*     call calc_robert(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
*    1                dbl_mb(k_orbn),nbasis,occ,virt)
c
*     call calc_pople(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
*    1                dbl_mb(k_orbn),nbasis,occ,virt)
c
c...  set orbital energies to 1 to get almlof to do cullen
c...  delta contains de delta e (*3)
c
      call dfill(nbasis,1.0d0,dbl_mb(k_orbn),1)
      delta = 1.0d0
      call calc_almlof(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
     1                 delta,dbl_mb(k_orbn),
     2                 nbasis,occ,virt,energy,.true.)
c
c     try integrating
c
c     call mod_EVALS(dbl_mb(k_orb),dbl_mb(k_orbn),nbasis,occ,virt)
c     call dcopy(nbasis,dbl_mb(k_orbn),1,dbl_mb(k_orb),1)
c     print *,' cw keep ORBITAL ENERGIES as HOMO-LUMO '
c
*     if (.not.rtdb_get_info(rtdb,'cct_almlof',i,n_almlof,date))
*    1 then
*       n_almlof = 1
*       np_almlof(1) = 5
*     else
*      if (.not.rtdb_get(rtdb,'cct_almlof',MT_INT,nn_almlof,np_almlof))
*    1 call errquit(' cct_almlof ',0)
*     end if
c
*     do i=1,n_almlof
*     call int_almlof(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
*    1                dbl_mb(k_orb),dbl_mb(k_orbn),   
*    2                nbasis,occ,virt,energy,np_almlof(i))
*     print *,np_almlof(i),' point almlof integrated energy ',energy
*     end do
c
*     if (.not.rtdb_get_info(rtdb,'cct_laguer',i,n_laguer,date))
*    1 then
*       n_laguer = 1
*       np_laguer(1) = 5
*     else
*      if (.not.rtdb_get(rtdb,'cct_laguer',MT_INT,nn_laguer,np_laguer))
*    1 call errquit(' cct_laguer ',0)
*     end if
c
*     do i=1,n_laguer
*     call int_laguer(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
*    1                dbl_mb(k_orb),dbl_mb(k_orbn),   
*    2                nbasis,occ,virt,energy,np_laguer(i),rtdb)
*     print *,np_laguer(i),' point laguer integrated energy ',energy
*     end do
c
c...  free all core beyond rint
c
      if (.not. ma_chop_stack(l_rint)) call errquit(' ma chop?', 0,
     &       MA_ERR)
c
c      could also use ma_pop_stack(l_occ)  which is like gmem_free
c      (but using handle instead of pointer)
c
      end
      subroutine int_almlof(rint,t1,t2,orb,fac,nbasis,occ,virt,
     1                      energy,np)
c
c...   control routine for Almloef integration
c...   energy returns the finally computed energy
c...   np is thew # points to be used
c
      implicit none
#include "errquit.fh"
c
      integer nbasis,occ,virt,np
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt),t1(occ,virt)
      double precision orb(nbasis),fac(nbasis),energy
c
      integer ip,itw
      double precision enerp
      common/flop/ term1,term2,term3,term4
      double precision term1,term2,term3,term4,tt1,tt2,tt3,tt4
c
      integer npoint,nd
      parameter (npoint=8, nd=npoint*(npoint+1)/2)
      double precision t_alpha(nd), weight(nd)
c
      data t_alpha/0.2543066d0,
     2             0.1057374d0, 0.8263985d0,
     3             0.0664618d0, 0.4380785d0, 1.6101925d0,
     4             0.0135967d0, 0.1646374d0, 0.6569702d0, 2.0357815d0,
     5             0.0064785d0, 0.0669154d0, 0.3089156d0, 0.9478808d0,
     5             2.5786393d0,
     6             0.0060050d0, 0.0553016d0, 0.2346500d0, 0.6587781d0,
     6             1.5760606d0, 3.6502575d0,
     7             0.0035243d0, 0.0251939d0, 0.1193583d0, 0.3653936d0,
     7             0.8850526d0, 1.9454617d0, 4.2401468d0,
     8             0.0034250d0, 0.0238627d0, 0.1023430d0, 0.2947359d0,
     8             0.6736046d0, 1.3821837d0, 2.7133497d0, 5.4051090d0/
      data weight /0.7443778d0,
     2             0.2998529d0, 1.3750404d0,
     3             0.1793919d0, 0.6329602d0, 1.9996848d0,
     4             0.0527171d0, 0.2757647d0, 0.7858024d0, 2.2786679d0,
     5             0.0177806d0, 0.1278500d0, 0.3878375d0, 0.9741098d0,
     5             2.6125732d0,
     6             0.0164421d0, 0.0995694d0, 0.2759458d0, 0.6106806d0,
     6             1.3178664d0, 3.1857008d0,
     7             0.0095904d0, 0.0437062d0, 0.1566528d0, 0.3542937d0,
     7             0.7264528d0, 1.4928836d0, 3.4634348d0,
     8             0.0093035d0, 0.0390022d0, 0.1267465d0, 0.2692945d0,
     8             0.5112993d0, 0.9512673d0, 1.8163835d0, 3.9522146d0/
c
      tt1 = 0.0d0
      tt2 = 0.0d0
      tt3 = 0.0d0
      tt4 = 0.0d0
c
      if (np.gt.npoint) call errquit('almlof overflow',0, UNKNOWN_ERR)
      energy=0.0d0
      do ip=1,np
       energy=energy+dexp(-5.0d0*t_alpha(np*(np-1)/2 + ip))
     1                          *weight(np*(np-1)/2 + ip)
      end do
      print *,np,' points 0.2 is ', energy
c
      energy = 0.0d0
      do ip=1,np
       itw = np*(np-1)/2 + ip
       call cct_fac(orb,fac,nbasis,occ,virt,t_alpha(itw))
c
       call calc_almlof(rint,t1,t2,-1.0d0,fac,
     1                  nbasis,occ,virt,enerp,.false.)
c
       energy = energy + enerp * weight(itw)
c      print *,' point ',ip,' of ',np,' contr ',enerp,' energy ',energy
       tt1 = tt1 + term1*weight(itw)
       tt2 = tt2 + term2*weight(itw)
       tt3 = tt3 + term3*weight(itw)
       tt4 = tt4 + term4*weight(itw)
      end do
c      print *,' t ',tt1,tt2,tt3,tt4
c
      return
      end
      subroutine int_laguer(rint,t1,t2,orb,fac,nbasis,occ,virt,
     1                      energy,np,rtdb)
c
c...   control routine for Almloef integration
c...   using straight gauss-laguerre
c...   energy returns the finally computed energy
c...   np is thew # points to be used
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer nbasis,occ,virt,np,rtdb
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt),t1(occ,virt)
      double precision orb(nbasis),fac(nbasis),energy
c
      integer ip,i,n
      double precision enerp,sum,avocc,avirt,delta
c
      integer npoint
      parameter (npoint= 200)
      double precision t_alpha(npoint), weight(npoint) 
      double precision scra(npoint), dum(2),elimit(2)
      character*26 date
c
      if (np.gt.npoint) call errquit(' enlarge int_laguer ',0,
     &       UNKNOWN_ERR)
      call gaussq(6,np,0.0d0,0.0d0,0,dum,scra,t_alpha,weight)
c
c...  dermine integration scale required
c     default : average delta e (*3)
c               might exclude extremes or all
c
      if (.not.rtdb_get_info(rtdb,'cct_laguer_scale',i,n,date))
     1 then
c..    default (-5 to +5 for now)
        elimit(1) = -5.0d0
        elimit(2) = +5.0d0
      else
       if (n.eq.0) then
c..    next default (all)
        elimit(1) = -9.0d99
        elimit(2) = +9.0d99
       else 
        if (.not.rtdb_get(rtdb,'cct_laguer_scale',MT_DBL,2,elimit))
     1     call errquit(' cct_laguer_scale ',0, RTDB_ERR)
        if (n.eq.1) elimit(2) = elimit(1)
        elimit(1) = -elimit(1)
       end if
      end if
c
c...    compute scale
c
      sum = 0.0d0
      n = 0
      do i=1,occ
         if (orb(i).gt.elimit(1)) then
            n = n + 1
            sum = sum + orb(i)
         end if
      end do
      avocc = sum/max0(n,1)
      ip = n
      sum = 0.0d0
      n = 0
      do i=occ+1,nbasis
         if (orb(i).lt.elimit(2)) then
            n = n + 1
            sum = sum + orb(i)
         end if
      end do
      avirt = sum/max0(n,1)
c...    scale
      delta = 3.0d0*(avirt-avocc)
      if (elimit(1).eq.0.0d0.and.elimit(2).ne.0.0d0) delta = elimit(2)
      if (delta .lt.1.0d-2) delta = 1.0d0
      print *,' laguerre scaled by ',delta,' kept ',ip,'/',n
c
c...  add the ex(xi) to weights
c
      energy=0.0d0
      do ip=1,np
       energy=energy+dexp(-5.0d0*t_alpha(ip)/delta)*weight(ip)
     1                           *dexp(t_alpha(ip))/delta
      end do
      print *,np,' points 0.2 is ', energy
c
      energy = 0.0d0
      do ip=1,np
       call cct_fac(orb,fac,nbasis,occ,virt,t_alpha(ip)/delta)
c
       call calc_almlof(rint,t1,t2,-1.0d0,fac,
     1                  nbasis,occ,virt,enerp,.false.)
c
       energy = energy + enerp * weight(ip)*dexp(t_alpha(ip))/delta
      end do
c
      return
      end

      subroutine get_EVALS(orb,nbasis,occ,virt)
c
c...  read orbital energies from EVALS
c
      implicit none
#include "errquit.fh"
c
      integer nbasis,occ,virt
      double precision orb(nbasis)
c
      integer nbin,i
      logical oexist
c
      inquire(file='EVALS',exist=oexist)
      if (.not.oexist) call errquit(' no EVAL ',0, UNKNOWN_ERR)
      open(9,file='EVALS',form='formatted')
      read(9,*) nbin
      if (nbin.ne.nbasis) call errquit('nbasis wrong',1, BASIS_ERR)
      do i=1,nbasis
         read(9,'(5x,f20.10)') orb(i)
      end do
      close(9)
c
      print *,' orbital energies read '
c     print *, orb
      return
      end
      subroutine mod_EVALS(orb,orbn,nbasis,occ,virt)
c
c...  make modified orbital energies  (homo/lumo)
c
      implicit none
c
      integer nbasis,occ,virt
      double precision orb(nbasis),orbn(nbasis)
      integer i
c
       do i=1,occ
         orbn(i) = orb(occ)
       end do
       do i=occ+1,nbasis
         orbn(i) = orb(occ+1)
       end do
c
      return
      end
      subroutine get_T2(t2,occ,virt)
c
c...  real t2 amplitudes from T2
c
      implicit none
#include "errquit.fh"
c
      integer occ,virt
      double precision t2(occ,occ,virt,virt)
c
      integer ocin,vin,i,j,a,b
      logical oexist
c
      inquire(file='T2',exist=oexist)
      if (.not.oexist) call errquit(' no T2 ',0, UNKNOWN_ERR)
      open(9,file='T2',form='formatted')
      read(9,*) ocin,vin
      if (ocin.ne.occ.or.vin.ne.virt) call errquit(' occ,virt wrong',0,
     &       UNKNOWN_ERR)
c
      call dfill(occ*occ*virt*virt,0.0d0,t2,1)
c
10    read(9,'(4i5,f20.10)',end=20) i,j,a,b,t2(i,j,a,b)
      go to 10
c
20    close(9)
      print *,' T2 read '
c     print *,((((t2(i,j,a,b),i=1,3),j=1,3),a=1,3),b=1,3)
c
      return
      end
      subroutine get_T1(t1,occ,virt)
c
c...  real t1 amplitudes from T1
c
      implicit none
#include "errquit.fh"
c
      integer occ,virt
      double precision t1(occ,virt)
c
      integer ocin,vin,i,a
      logical oexist
c
      inquire(file='T1',exist=oexist)
      if (.not.oexist) then
         print *,' no T1 '
         call dfill(occ*virt,0.0d0,t1,1)
         return
      end if
      open(9,file='T1',form='formatted')
      read(9,*) ocin,vin
      if (ocin.ne.occ.or.vin.ne.virt) call errquit(' occ,virt wrong',1,
     &       UNKNOWN_ERR)
c
      call dfill(occ*virt,0.0d0,t1,1)
c
10    read(9,'(2i5,f20.10)',end=20) i,a,t1(i,a)
      go to 10
c
20    close(9)
      print *,' T1 read '
c
      return
      end
      subroutine get_INTS(rint,nbasis)
c
c...  real integrals from ASOINTS
c
      implicit none
#include "errquit.fh"
c
      integer nbasis
      double precision rint(nbasis,nbasis,nbasis,nbasis)
c
      integer in,i,j,k,l
      logical oexist
c
      inquire(file='ASOINTS',exist=oexist)
      if (.not.oexist) call errquit(' no INTS ',0, DISK_ERR)
      open(9,file='ASOINTS',form='formatted')
      read(9,*) in
      if (in.ne.nbasis) call errquit(' nbasis on INTS wrong',0,
     &       BASIS_ERR)
c
      call dfill(nbasis**4,0.0d0,rint,1)
c
10    read(9,'(4i5,f20.10)',end=20) i,j,k,l,rint(i,j,k,l)
      go to 10
c
20    close(9)
      print *,' integrals read '
c     print *,((((rint(i,j,k,l),i=1,3),j=1,3),k=1,3),l=1,3)
c
      return
      end
      subroutine cct_fac(orb,fac,nbasis,occ,virt,factor)
c
c...  generatewd factores for laplace integration
c...  calculated as exp(-orb)*factor for occ
c...                exp(orb)*factor for virt
c
      implicit none
#include "errquit.fh"
      integer nbasis,occ,virt
      double precision orb(nbasis),fac(nbasis),factor
c
      integer i
c
      do i=1,occ
       fac(i) = dexp(orb(i)*factor)
      end do
      do i=occ+1,nbasis
       fac(i) = dexp(-orb(i)*factor)
      end do
c
      if (orb(occ).gt.0.0d0.or.orb(occ+1).lt.0.0d0)
     1 call errquit(' OEPS orbital energies ',0, UNKNOWN_ERR)
c     print *,' factor ... ',orb(occ),orb(occ+1)
c     print *,fac(1),fac(nbasis)
c     print *,(fac(1)*fac(nbasis))**3
c
      return
      end 
      subroutine calc_robert(rint,t1,t2,orb,nbasis,occ,virt)
c
c...   use Robert's formulae
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t1(occ,virt), t2(occ,occ,virt,virt)
      double precision orb(nbasis)
c..     intermediates
      double precision w1,w2,w3,w4,v,w1w1,w1w2,w1w3,w1w4
      double precision vw1,vw2,vw3,vw4
c
      double precision e1,e2,energy,delta
      integer i,j,k,a,b,c,e,m
c
      w1w1 = 0.0d0
      w1w2 = 0.0d0
      w1w3 = 0.0d0
      w1w4 = 0.0d0
      vw1 = 0.0d0
      vw2 = 0.0d0
      vw3 = 0.0d0
      vw4 = 0.0d0
c
      do i=1,occ
       do j=1,occ
        do k=1,occ
         do a=1,virt
          do b=1,virt
           do c=1,virt
            w1 = 0.0d0
            w2 = 0.0d0
            w3 = 0.0d0
            w4 = 0.0d0
            do e=1,virt
             w1 = w1 + t2(i,j,c,e)*rint(a+occ,b+occ,e+occ,k)
             w2 = w2 + t2(k,j,c,e)*rint(a+occ,b+occ,e+occ,i)
             w3 = w3 + t2(i,j,a,e)*rint(c+occ,b+occ,e+occ,k)
             w4 = w4 + t2(k,j,a,e)*rint(c+occ,b+occ,e+occ,i)
            end do
            do m=1,occ
             w1 = w1 + t2(k,m,a,b)*rint(c+occ,m,i,j)
             w2 = w2 + t2(i,m,a,b)*rint(c+occ,m,k,j)
             w3 = w3 + t2(k,m,c,b)*rint(a+occ,m,i,j)
             w4 = w4 + t2(i,m,c,b)*rint(a+occ,m,k,j)
            end do
            v = t1(k,c)*rint(i,j,a,b)
            delta = (orb(i)+orb(j)+orb(k)
     1            - orb(a+occ)-orb(b+occ)-orb(c+occ))
            w1w1 = w1w1 + w1*w1 / delta
            w1w2 = w1w2 + w1*w2 / delta
            w1w3 = w1w3 + w1*w3 / delta
            w1w4 = w1w4 + w1*w4 / delta
            vw1 = vw1 + v*w1 / delta
            vw2 = vw2 + v*w2 / delta
            vw3 = vw3 + v*w3 / delta
            vw4 = vw4 + v*w4 / delta
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  calculate energies
c
      e1 = (vw1 - 2.0d0*vw2 -2.0d0*vw3 + 4.0d0*vw4) /4.0d0
      e2 = (w1w1 - 2.0d0*w1w2 -2.0d0*w1w3 + 4.0d0*w1w4) /4.0d0
      print *,' w1w1 ',w1w1
      print *,' w1w2 ',w1w2
      print *,' w1w3 ',w1w3
      print *,' w1w4 ',w1w4
      energy = e1 + e2
c
      print *,' Robert : e1= ',e1,' e2= ',e2,' e= ',energy
c
      end
      subroutine calc_pople(rint,t1,t2,orb,nbasis,occ,virt)
c
c...   use Pople's formulae (who can argue with a Nobel-prize)
c...   BUT ... typed in from Laplace triple document (jvl)
c...   ( + formulae(QCI ...))
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t1(occ,virt), t2(occ,occ,virt,virt)
      double precision orb(nbasis)
c..     intermediates
      double precision u1,u2,u3,u4,u5,u6,u7,u8,u9,v1
c
      double precision e1,e2,energy,delta
      integer i,j,k,a,b,c,e,m
c
      e1 = 0.0d0
      e2 = 0.0d0
c
      do i=1,occ
       do j=1,occ
        do k=1,occ
         do a=1,virt
          do b=1,virt
           do c=1,virt
            u1 = 0.0d0
            u2 = 0.0d0
            u3 = 0.0d0
            u4 = 0.0d0
            u5 = 0.0d0
            u6 = 0.0d0
            u7 = 0.0d0
            u8 = 0.0d0
            u9 = 0.0d0
            do e=1,virt
             u1 = u1 + t2(i,j,a,e)*rint(b+occ,c+occ,e+occ,k)
             u2 = u2 + t2(i,j,b,e)*rint(c+occ,a+occ,e+occ,k)
             u3 = u3 + t2(i,j,c,e)*rint(a+occ,b+occ,e+occ,k)
             u4 = u4 + t2(k,i,a,e)*rint(b+occ,c+occ,e+occ,j)
             u5 = u5 + t2(k,i,b,e)*rint(c+occ,a+occ,e+occ,j)
             u6 = u6 + t2(k,i,c,e)*rint(a+occ,b+occ,e+occ,j)
             u7 = u7 + t2(j,k,a,e)*rint(b+occ,c+occ,e+occ,i)
             u8 = u8 + t2(j,k,b,e)*rint(c+occ,a+occ,e+occ,i)
             u9 = u9 + t2(j,k,c,e)*rint(a+occ,b+occ,e+occ,i)
            end do
            do m=1,occ
             u1 = u1 + t2(i,m,a,b)*rint(c+occ,m,j,k)
             u2 = u2 + t2(i,m,b,c)*rint(a+occ,m,j,k)
             u3 = u3 + t2(i,m,c,a)*rint(b+occ,m,j,k)
             u4 = u4 + t2(j,m,a,b)*rint(c+occ,m,k,i)
             u5 = u5 + t2(j,m,b,c)*rint(a+occ,m,k,i)
             u6 = u6 + t2(j,m,c,a)*rint(b+occ,m,k,i)
             u7 = u7 + t2(k,m,a,b)*rint(c+occ,m,i,j)
             u8 = u8 + t2(k,m,b,c)*rint(a+occ,m,i,j)
             u9 = u9 + t2(k,m,c,a)*rint(b+occ,m,i,j)
            end do
            v1 = t1(i,a)*rint(j,k,b+occ,c+occ) 
     1         + t1(i,b)*rint(j,k,c+occ,a+occ) 
     2         + t1(i,c)*rint(j,k,a+occ,b+occ) 
     3         + t1(j,a)*rint(k,i,b+occ,c+occ) 
     4         + t1(j,b)*rint(k,i,c+occ,a+occ)
     5         + t1(j,c)*rint(k,i,a+occ,b+occ) 
     6         + t1(k,a)*rint(i,j,b+occ,c+occ) 
     7         + t1(k,b)*rint(i,j,c+occ,a+occ)
     8         + t1(k,c)*rint(i,j,a+occ,b+occ) 
            delta = (orb(i)+orb(j)+orb(k)
     1            - orb(a+occ)-orb(b+occ)-orb(c+occ))
            e1 = e1+ v1*(u1+u2+u3+u4+u5+u6+u7+u8+u9) / delta
            e2 = e2+ (u1+u2+u3+u4+u5+u6+u7+u8+u9)**2 / delta
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  calculate energies
c
      e1 = e1 / 36.0d0
      e2 = e2 / 36.0d0
      energy = e1 + e2
c
      print *,' Pople  : e1= ',e1,' e2= ',e2,' e= ',energy
c
      end
      subroutine calc_almlof(rint,t1,t2,delta,fac,nbasis,occ,virt,
     1                       energy,opr)
c
c...   this routine can do Cullen+Zerner's approach by
c...   setting all fac'a ,to 1 and delta to delta e
c...   for almlof delta is set to 1
c
c...   starting from "Robert's" formula
c     notation v(vv) : running index (virtual), o(oo)= running occupied)
c
c
      implicit none
#include "errquit.fh"
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t1(occ,virt), t2(occ,occ,virt,virt)
      double precision delta,fac(nbasis),energy
      logical opr
#include "global.fh"
#include "mafdecls.fh"
c
c...  functions
c
      double precision cct_so_a4a,cct_so_a4b,cct_so_b4b
      double precision cct_so_a3a,cct_so_a3b,cct_so_b3b
      double precision cct_so_a2a,cct_so_a2b,cct_so_b2b
      double precision cct_so_a1a,cct_so_a1b,cct_so_b1b
c
      double precision term4,term_a4a,term_a4b,term_b4b
      double precision term3,term_a3a,term_a3b,term_b3b
      double precision term2,term_a2a,term_a2b,term_b2b
      double precision term1,term_a1a,term_a1b,term_b1b
      double precision e1,e2
      common/flop/ term1,term2,term3,term4
c
      integer lenabcd,leniabc,lenijab,lenijka,lenijkl,lenab,lenia,lenij
      integer l_w,k_w
c
      lenabcd = virt*virt*virt*virt
      leniabc = occ*virt*virt*virt
      lenijab = occ*occ*virt*virt
      lenijka = occ*occ*occ*virt
      lenijkl = occ*occ*occ*occ
      lenab = virt*virt
      lenia = virt*occ
      lenij = occ*occ
      if (opr) print *,' triple delta e used ',delta
c
c term 1
c
      if (.not. ma_push_get(mt_dbl,lenab*2,'a2a', l_w, k_w))
     1    call errquit('push a1a',0, MA_ERR)
      term_a1a = cct_so_a1a(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenab))
      if (.not.ma_pop_stack(l_w)) call errquit('pop a2a',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lenia*2,'a1b', l_w, k_w))
     1    call errquit('push a1b',0, MA_ERR)
      term_a1b = cct_so_a1b(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenia))
      if (.not.ma_pop_stack(l_w)) call errquit('pop a1b',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lenij*2,'b1b', l_w, k_w))
     1    call errquit('push b1b',0, MA_ERR)
      term_b1b = cct_so_b1b(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenij))
      if (.not.ma_pop_stack(l_w)) call errquit('pop b1b',0, MA_ERR)
c
      term1 = term_a1a + 2.0d0*term_a1b + term_b1b
c
      if (opr) then
c        print *,' term a1a ',term_a1a/ delta
c        print *,' term a1b=b1a ',term_a1b/ delta
c        print *,' term b1b ',term_b1b/ delta
         print *,' term1 ',term1/ delta
      end if

c
c term 2
c
      if (.not. ma_push_get(mt_dbl,lenijab+lenijab,'a2a', l_w, k_w))
     1    call errquit('push a2a',0, MA_ERR)
      term_a2a = cct_so_a2a(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenijab))
      if (.not.ma_pop_stack(l_w)) call errquit('pop a2a',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lenijka*2,'a2a', l_w, k_w))
     1    call errquit('push a2b',0, MA_ERR)
      term_a2b = cct_so_a2b(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenijka))
      if (.not.ma_pop_stack(l_w)) call errquit('pop a2b',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lenijka+lenijkl,'b2b', l_w, k_w))
     1    call errquit('push b2b',0, MA_ERR)
      term_b2b = cct_so_b2b(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenijka))
      if (.not.ma_pop_stack(l_w)) call errquit('pop b2b',0, MA_ERR)
c
      term2 = term_a2a + 2.0d0*term_a2b + term_b2b
c
      if (opr) then
c        print *,' term a2a ',term_a2a/ delta
c        print *,' term a2b=b2a ',term_a2b/ delta
c        print *,' term b2b ',term_b2b/ delta
         print *,' term2 ',term2/ delta
      end if
c
c term 3
c
      if (.not. ma_push_get(mt_dbl,lenabcd+leniabc,'a3a', l_w, k_w))
     1    call errquit('push a3a',0, MA_ERR)
      term_a3a = cct_so_a3a(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenabcd))
      if (.not.ma_pop_stack(l_w)) call errquit('pop a3a',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl,leniabc+lenijka,'a3b', l_w, k_w))
     1    call errquit('push a3b',0, MA_ERR)
      term_a3b = cct_so_a3b(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+leniabc))
      if (.not.ma_pop_stack(l_w)) call errquit('pop a3b',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lenijab+lenijka,'b3b', l_w, k_w))
     1    call errquit('push b3b',0, MA_ERR)
      term_b3b = cct_so_b3b(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenijab))
      if (.not.ma_pop_stack(l_w)) call errquit('pop b3b',0, MA_ERR)
c
      print *,' term_a3a *.5 ',term_a3a*0.5d0
      print *,' term_a3b *1, ',term_a3b 
      print *,' term_b3b *.5 ',term_b3b*0.5d0 
      term3 = term_a3a + 2.0d0*term_a3b + term_b3b
c
      if (opr) then
c        print *,' term a3a ',term_a3a/ delta
c        print *,' term a3b=b3a ',term_a3b/ delta
c        print *,' term b3b ',term_b3b/ delta
         print *,' term3 ',term3/ delta
      end if

c
c term 4
c
      if (.not. ma_push_get(mt_dbl,leniabc+lenijab,'a4a', l_w, k_w))
     1    call errquit('push a4a',0, MA_ERR)
      term_a4a = cct_so_a4a(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+leniabc))
      if (.not.ma_pop_stack(l_w)) call errquit('pop a4a',0, MA_ERR)
c
      if (.not. ma_push_get(mt_dbl,lenijab+lenijka,'a4b', l_w, k_w))
     1    call errquit('push a4b',0, MA_ERR)
      term_a4b = cct_so_a4b(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenijab))
      if (.not.ma_pop_stack(l_w)) call errquit('pop a4b',0, MA_ERR)
c
      if (.not. ma_push_get(mt_dbl,lenijka*2,'b4b', l_w, k_w))
     1    call errquit('push b4b',0, MA_ERR)
      term_b4b = cct_so_b4b(rint,nbasis,t2,occ,virt,fac,
     1                      dbl_mb(k_w),dbl_mb(k_w+lenijka))
      if (.not.ma_pop_stack(l_w)) call errquit('pop b4b',0, MA_ERR)
c
      term4 = term_a4a + 2.0d0*term_a4b + term_b4b
c
      if (opr) then
c        print *,' term a4a ',term_a4a/ delta
c        print *,' term a4b=b4a ',term_a4b/ delta
c        print *,' term b4b ',term_b4b/ delta
         print *,' term4 ',term4/ delta
      end if
c
       e1 = 0.0d0
       e2 = (0.25d0*term1 - 0.5d0*term2 - 0.5d0*term3 + term4)/ delta
       energy = e1 + e2
      if (opr) print *,' Cullen  : e1= ',e1,' e2= ',e2,' e= ',energy
c
      end
      double precision function cct_so_a4a(rint,nbasis,t2,occ,virt,fac,
     1                                     wiabc,wijab)
c
c     term A4A 
c       t(ij,ce) (ab||ek) t(kj,af) (cb||fi)
c     = t(ij,ce) (cb||fi) (ab||ek) t(kj,af) 
c     weighting 1.
c
c...  fac has to be added to i,j,k,a,b,c (there is choice)
c...  early is done / late might be better ?
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt),fac(nbasis)
c
      integer j,e,b,f,i,c,k,a
      double precision wiabc(virt,virt,virt,occ),
     1                 wijab(virt,virt,occ,occ),sum
c
c...  t2 * integral (ic) t(ij,ce) (cb||fi) => (j,e,b,f)
c
      do j=1,occ
       do e=1,virt
        do b=1,virt
         do f=1,virt
          wiabc(f,b,e,j) = 0.0d0
          do i=1,occ
           do c=1,virt
            wiabc(f,b,e,j) = wiabc(f,b,e,j) + t2(i,j,c,e)
     1                                      * rint(c+occ,b+occ,f+occ,i)
     2                     * fac(j)*fac(b+occ)*fac(i)*fac(c+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (eb) (j,e,b,f) (ab||ek) => (j,k,f,a)  (choice)
c
      do j=1,occ
       do k=1,occ
        do f=1,virt
         do a=1,virt
          wijab(a,f,k,j) = 0.0d0
          do e=1,virt
           do b = 1,virt
            wijab(a,f,k,j) = wijab(a,f,k,j) + wiabc(f,b,e,j)
     1                                      * rint(a+occ,b+occ,e+occ,k)
     2                     * fac(k)*fac(a+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (jkfa)  (j,k,f,a) t(kj,af)
c
      sum = 0.0d0
      do j=1,occ
       do k=1,occ
        do f=1,virt
         do a=1,virt
          sum= sum + wijab(a,f,k,j) * t2(k,j,a,f)
         end do
        end do
       end do
      end do
c
      cct_so_a4a = sum
c
      return
      end 
      double precision function cct_so_a4b(rint,nbasis,t2,occ,virt,fac,
     1                                     wijab,wijka)
c
c     term A4B = B4A
c A4B = t(ij,ce) (ab||ek) t(in,cb) (an||kj)
c     = t(ij,ce) t(in,cb) (ab||ek) (an||kj)
c
c B4A = t(km,ab) (cm||ij) t(kj,af) (cb||fi)
c relabel in cb   an  kj    ij ce   ab  ek
c     weighting 1.
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt),fac(nbasis)
c
      integer j,n,e,b,i,c,k,a
      double precision wijab(virt,virt,occ,occ),
     1                 wijka(virt,occ,occ,occ),sum
c
c..   contract t's (ic) t(ij,ce) t(in,cb) => (j,n,e,b) choice
c
      do j=1,occ
       do n=1,occ
        do e=1,virt
         do b=1,virt
          wijab(b,e,n,j) = 0.0d0
          do i=1,occ
           do c=1,virt
            wijab(b,e,n,j) = wijab(b,e,n,j) + t2(i,j,c,e)
     1                                      * t2(i,n,c,b)
     2                     * fac(j)*fac(b+occ)*fac(i)*fac(c+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c..   (be)  (j,n,e,b) (ab||ek)  => (j,n,k,a) (choice)
c
      do j=1,occ
       do n=1,occ
        do k=1,occ
         do a=1,virt
          wijka(a,k,n,j) = 0.0d0
          do b=1,virt
           do e=1,virt
            wijka(a,k,n,j) = wijka(a,k,n,j) + wijab(b,e,n,j)
     1                                      * rint(a+occ,b+occ,e+occ,k)
     2                     * fac(k)*fac(a+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c     (jnka)  (j,n,k,a) ((an||kj) => sum
c
      sum = 0.0d0
      do j=1,occ
       do n=1,occ
        do k=1,occ
         do a=1,virt
          sum = sum + wijka(a,k,n,j) * rint(a+occ,n,k,j)
         end do
        end do
       end do
      end do
c
      cct_so_a4b = sum
c
      return
      end 
      double precision function cct_so_b4b(rint,nbasis,t2,occ,virt,fac,
     1                                     wijka,wijka2)
c
c     term B4B 
c       t(km,ab)(cm||ij) t(in,cb) (an||kj)
c     = t(km,ab)(an||kj) (cm||ij) t(in,cb)
c
c     weighting 1.
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer m,n,j,b,k,a,i,c
      double precision wijka(virt,occ,occ,occ),
     1                 wijka2(virt,occ,occ,occ),sum
c
c...  (ka) ..  t*rint (left) t(km,ab)(an||kj => (m,n,j,b)
c
      do m=1,occ
       do n=1,occ
        do j=1,occ
         do b=1,virt
          wijka(b,j,n,m) = 0.0d0
          do k=1,occ
           do a=1,virt
            wijka(b,j,n,m) = wijka(b,j,n,m) + t2(k,m,a,b)
     1                                      * rint(a+occ,n,k,j)
     2                     * fac(j)*fac(b+occ)*fac(k)*fac(a+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (nb) ..  (m,n,j,b) t(in,cb) => (m,j,i,c)  
c
      do m=1,occ
       do j=1,occ
        do i=1,occ
         do c=1,virt
          wijka2(c,i,j,m) = 0.0d0
          do n=1,occ
           do b=1,virt
            wijka2(c,i,j,m) = wijka2(c,i,j,m) + wijka(b,j,n,m)
     1                                        * t2(i,n,c,b)
     2                      * fac(i)*fac(c+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...    (mjic) ..  (m,j,i,c) ((cm||ij) => sum
c
      sum = 0.0d0
      do m=1,occ
       do j=1,occ
        do i=1,occ
         do c=1,virt      
          sum = sum + wijka2(c,i,j,m) * rint(c+occ,m,i,j)
         end do
        end do
       end do
      end do
c
      cct_so_b4b = sum
c
      return
      end 
      double precision function cct_so_a3a(rint,nbasis,t2,occ,virt,fac,
     1                                     wabcd,wiabc)
c
c     term A3A 
c       
c       t(ij,ce) (ab||ek) t(ij,af) (cb||fk)
c     = t(ij,ce) t(ij,af) (ab||ek) (cb||fk)
c
c     weighting 0.5
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer c,e,a,f,i,j,k,b
      double precision wabcd(virt,virt,virt,virt),
     1                 wiabc(virt,virt,virt,occ),sum
c
c...  (ij) t(ij,ce) t(ij,af) => (c,e,a,f)
c
      do c=1,virt
       do e=1,virt
        do a=1,virt
         do f=1,virt
          wabcd(f,a,e,c) = 0.0d0
          do i=1,occ
           do j=1,occ
            wabcd(f,a,e,c) = wabcd(f,a,e,c) + t2(i,j,c,e)*t2(i,j,a,f)
     2                     * fac(i)*fac(j)*fac(a+occ)*fac(c+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (ae) (c,e,a,f) (ab||ek) => (k,c,f,b)
c
      do k=1,occ
       do c=1,virt
        do f=1,virt
         do b=1,virt
          wiabc(b,f,c,k) = 0.0d0
          do a=1,virt
           do e=1,virt
            wiabc(b,f,c,k) = wiabc(b,f,c,k) + wabcd(f,a,e,c)
     1                                      * rint(a+occ,b+occ,e+occ,k)
     2                     * fac(k)*fac(b+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c..   (kcfb) (k,c,f,b) (cb||fk) => sum
c
      sum = 0.0d0
      do k=1,occ
       do c=1,virt
        do f=1,virt
         do b=1,virt
          sum = sum + wiabc(b,f,c,k)*rint(c+occ,b+occ,f+occ,k)
         end do
        end do
       end do
      end do
c
      cct_so_a3a = sum
c
      return 
      end
      double precision function cct_so_a3b(rint,nbasis,t2,occ,virt,fac,
     1                                     wiabc,wijka)
c
c     term A3B = B3A
c       
c       t(ij,ce) (ab||ek) t(kn,cb) (an||ij)
c     = (ab||ek) t(kn,cb) t(ij,ce) (an||ij)
c
c B3A = t(km,ab) (cm||ij) t(ij,af) (cb||fk)
c relabel kn cb   an  ij    ij ce   ab  ek
c     weighting 0.5
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer n,c,a,e,k,b,i,j
      double precision wiabc(virt,virt,virt,occ),
     1                 wijka(virt,occ,occ,occ),sum
c
c...  (kb) t(kn,cb) (ab||ek) => (n,c,a,e)
c
      do n=1,occ
       do c=1,virt
        do a=1,virt
         do e=1,virt
          wiabc(e,a,c,n) = 0.0d0
          do k=1,occ
           do b=1,virt
            wiabc(e,a,c,n) = wiabc(e,a,c,n) + t2(k,n,c,b)
     1                                      * rint(a+occ,b+occ,e+occ,k)
     2                     * fac(c+occ)*fac(a+occ)*fac(k)*fac(b+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (ce) (n,c,a,e) t(i,j,c,e) => (n,i,j,a)
c
      do n=1,occ
       do i=1,occ
        do j=1,occ
         do a=1,virt
          wijka(a,j,i,n) = 0.0d0
          do c=1,virt
           do e=1,virt
            wijka(a,j,i,n) = wijka(a,j,i,n) + wiabc(e,a,c,n)*t2(i,j,c,e)
     2                     * fac(i)*fac(j)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (nija) (n,i,j,a) (an||ij) => sum
c
      sum = 0.0d0
      do n=1,occ
       do i=1,occ
        do j=1,occ
         do a=1,virt
          sum = sum + wijka(a,j,i,n) * rint(a+occ,n,i,j)
         end do
        end do
       end do
      end do
c
      cct_so_a3b = sum
c
      return 
      end
      double precision function cct_so_b3b(rint,nbasis,t2,occ,virt,fac,
     1                                     wijab,wijka)
c
c     term B3B 
c       
c       t(km,ab) (cm||ij) t(kn,cb) (an||ij)
c     = t(km,ab) t(kn,cb) (cm||ij) (an||ij)
c
c     weighting 0.5
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer m,n,j,b,k,a,i,c
      double precision wijab(virt,virt,occ,occ),
     1                 wijka(virt,occ,occ,occ),sum
c
c...  (kb) t(km,ab) t(kn,cb) => (m,n,a,c)
c
      do m=1,occ
       do n=1,occ
        do a=1,virt
         do c=1,virt
          wijab(c,a,n,m) = 0.0d0
          do k=1,occ
           do b=1,virt
            wijab(c,a,n,m) = wijab(c,a,n,m) + t2(k,m,a,b) * t2(k,n,c,b)
     2                     * fac(a+occ)*fac(c+occ)*fac(k)*fac(b+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (mc) (m,n,a,c) (cm||ij) => (n,i,j,a)
c
      do n=1,occ
       do i=1,occ
        do j=1,occ
         do a=1,virt
          wijka(a,j,i,n) = 0.0d0
          do m=1,occ
           do c=1,virt
            wijka(a,j,i,n) = wijka(a,j,i,n) + wijab(c,a,n,m) 
     1                                      * rint(c+occ,m,i,j)
     2                     * fac(i)*fac(j)
           end do
          end do
         end do
        end do
       end do
      end do
c
c..   (nija) (n,i,j,a) (an||ij) => sum
c
      sum = 0.0d0
      do n=1,occ
       do i=1,occ
        do j=1,occ
         do a=1,virt
          sum = sum + wijka(a,j,i,n)*rint(a+occ,n,i,j)
         end do
        end do
       end do
      end do
c
      cct_so_b3b = sum
c
      return 
      end
      double precision function cct_so_a2a(rint,nbasis,t2,occ,virt,fac,
     1                                     wijab,wijab2)
c
c     term A2A 
c       
c       t(ij,ce) (ab||ek) t(kj,cf) (ab||fi)
c     = (ab||ek) (ab||fi) t(ij,ce) t(kj,cf)
c
c     weighting 0.5
c
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer k,i,e,f,a,b,j,c
      double precision wijab(virt,virt,occ,occ),
     1                 wijab2(virt,virt,occ,occ),sum
c
c...  (ab) (ab||ek) (ab||fi) => (k,i,e,f)
c
      do k=1,occ
       do i=1,occ
        do e=1,virt
         do f=1,virt
          wijab(f,e,i,k) = 0.0d0
          do a=1,virt
           do b=1,virt
             wijab(f,e,i,k) =  wijab(f,e,i,k)+rint(a+occ,b+occ,e+occ,k)
     1                                       *rint(a+occ,b+occ,f+occ,i)
     2                      * fac(k)*fac(i)*fac(a+occ)*fac(b+occ)
             if (.not. ma_verify_allocator_stuff()) call errquit('a',0,
     &       MA_ERR)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (ie) (k,i,e,f) t(ij,ce) => (k,j,c,f)
c
      do k=1,occ
       do j=1,occ
        do c=1,virt
         do f=1,virt
          wijab2(f,c,j,k) = 0.0d0
          do i=1,occ
           do e=1,virt
            wijab2(f,c,j,k) = wijab2(f,c,j,k) + wijab(f,e,i,k)
     1                                        * t2(i,j,c,e)
     2                      * fac(j)*fac(c+occ)
             if (.not. ma_verify_allocator_stuff()) call errquit('b',0,
     &       MA_ERR)
           end do
          end do
         end do
        end do
       end do
      end do
c
c..   (kjcf) (k,j,c,f) t(kj,cf) => sum
c
      sum = 0.0d0
      do k=1,occ
       do j=1,occ
        do c=1,virt
         do f=1,virt
          sum = sum + wijab2(f,c,j,k) * t2(k,j,c,f)
         end do
        end do
       end do
      end do
c
      cct_so_a2a = sum
c
      return 
      end
      double precision function cct_so_a2b(rint,nbasis,t2,occ,virt,fac,
     1                                     wijka,wijka2)
c
c     term A2B=B2A
c       
c       t(ij,ce) (ab||ek) t(in,ab) (cn||kj)
c     = (ab||ek) t(in,ab) t(ij,ce) (cn||kj)
c B2A = t(km,ab) (cm||ij) t(kj,cf) (ab||fi)
c rename  in ab   cn  kj    ij ce   ab  ek
c
c     weighting 0.5
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer i,n,k,e,a,b,j,c
      double precision wijka(virt,occ,occ,occ),
     1                 wijka2(virt,occ,occ,occ),sum
c
c...  (ab) t(in,ab) (ab||ek) => (i,n,k,e)
c
      do i=1,occ
       do n=1,occ
        do k=1,occ
         do e=1,virt
          wijka(e,k,n,i) = 0.0d0
          do a=1,virt
           do b=1,virt
            wijka(e,k,n,i) = wijka(e,k,n,i) + t2(i,n,a,b)
     1                                      * rint(a+occ,b+occ,e+occ,k)
     2                     * fac(i)*fac(k)*fac(a+occ)*fac(b+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (ie) (i,n,k,e) t(ij,ce) => (n,k,j,c)
c
      do n=1,occ
       do k=1,occ
        do j=1,occ
         do c=1,virt
          wijka2(c,j,k,n) = 0.0d0
          do i=1,occ
           do e=1,virt
            wijka2(c,j,k,n) = wijka2(c,j,k,n) + wijka(e,k,n,i) 
     1                                        * t2(i,j,c,e)
     2                      * fac(j)*fac(c+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c..   (nkjc) (n,k,j,c) (cn||kj) => sum
c
      sum = 0.0d0
      do n=1,occ
       do k=1,occ
        do j=1,occ
         do c=1,virt
          sum = sum + wijka2(c,j,k,n)*rint(c+occ,n,k,j)
         end do
        end do
       end do
      end do
c
      cct_so_a2b = sum
c
      return 
      end
      double precision function cct_so_b2b(rint,nbasis,t2,occ,virt,fac,
     1                                     wijka,wijkl)
c
c     term B2B 
c       
c       t(km,ab) (cm||ij) t(in,ab) (cn||kj)
c     = t(km,ab) t(in,ab) (cm||ij) (cn||kj)
c
c     weighting 0.5
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer k,m,i,n,a,b,j,c
      double precision wijka(virt,occ,occ,occ),
     1                 wijkl(occ,occ,occ,occ),sum
c
c...  (ab) t(km,ab) t(in,ab) => (k,m,i,n)
c
      do k=1,occ
       do m=1,occ
        do i=1,occ
         do n=1,occ
          wijkl(n,i,m,k) = 0.0d0
          do a=1,virt
           do b=1,virt
             wijkl(n,i,m,k) =  wijkl(n,i,m,k) + t2(k,m,a,b)*t2(i,n,a,b)
     2                      * fac(k)*fac(i)*fac(a+occ)*fac(b+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  (mi) (k,m,i,n) (cm||ij) => (k,n,j,c)
c
      do k=1,occ
       do n=1,occ
        do j=1,occ
         do c=1,virt
          wijka(c,j,n,k) = 0.0d0
          do m=1,occ
           do i=1,occ
            wijka(c,j,n,k) = wijka(c,j,n,k) + wijkl(n,i,m,k) 
     1                                      * rint(c+occ,m,i,j)
     2                     * fac(j)*fac(c+occ)
           end do
          end do
         end do
        end do
       end do
      end do
c
c..   (knjc) (k,n,j,c) (cn||kj) => sum
c
      sum = 0.0d0
      do k=1,occ
       do n=1,occ
        do j=1,occ
         do c=1,virt
          sum = sum + wijka(c,j,n,k)*rint(c+occ,n,k,j)
         end do
        end do
       end do
      end do
c
      cct_so_b2b = sum
c
      return 
      end
      double precision function cct_so_a1a(rint,nbasis,t2,occ,virt,fac,
     1                                     wab,wab2)
c
c     term A1A 
c     *seems might be combined with A2A*
c       
c       t(ij,ce) (ab||ek) t(ij,cf) (ab||fk)
c     = (ab||ek) (ab||fk) t(ij,ce) t(ij,cf)
c
c     weighting 0.25
c
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer e,f,k,a,b,i,j,c
      double precision wab(virt,virt),
     1                 wab2(virt,virt),sum,ddot
c
c...  (abk) (ab||ek) (ab||fk) => (e,f)
c
      do e=1,virt
       do f=1,virt
        wab(f,e) = 0.0d0
        do k=1,occ
         do a=1,virt
          do b=1,virt
           wab(f,e) =  wab(f,e) + rint(a+occ,b+occ,e+occ,k)
     1                          * rint(a+occ,b+occ,f+occ,k)
     2              * fac(k)*fac(a+occ)*fac(b+occ)
          end do
         end do
        end do
       end do
      end do
c
c...  (ijc) t(ij,ce) t(ijcf)  => (e,f)
c
      do e=1,virt
       do f=1,virt
        wab2(f,e) = 0.0d0
        do i=1,occ
         do j=1,occ
          do c=1,virt
            wab2(f,e) = wab2(f,e) + t2(i,j,c,e) * t2(i,j,c,f)
     2                * fac(i)*fac(j)*fac(c+occ)
          end do
         end do
        end do
       end do
      end do
c
      sum = ddot(virt*virt,wab,1,wab2,1)
c
      cct_so_a1a = sum
c
      return 
      end
      double precision function cct_so_a1b(rint,nbasis,t2,occ,virt,fac,
     1                                     wia,wia2)
c
c     term A1B=B1A
c       
c       t(ij,ce) (ab||ek) t(kn,ab) (cn||ij)
c     = (ab||ek) t(kn,ab) t(ij,ce) (cn||ij)
c
c     weighting 0.25
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer n,e,k,a,b,i,j,c
      double precision wia(virt,occ),
     1                 wia2(virt,occ),sum,ddot
c
c...  (abk) t(kn,ab) (ab||ek) => (n,e)
c
      do n=1,occ
       do e=1,virt
        wia(e,n) = 0.0d0
        do k=1,occ
         do a=1,virt
          do b=1,virt
           wia(e,n) = wia(e,n) + t2(k,n,a,b) * rint(a+occ,b+occ,e+occ,k)
     2              * fac(k)*fac(a+occ)*fac(b+occ)
          end do
         end do
        end do
       end do
      end do
c
c...  (ijc) (cn||ij) t(ij,ce) => (n,e)
c
      do n=1,occ
       do e=1,virt
        wia2(e,n) = 0.0d0
        do i=1,occ
         do j=1,occ
          do c=1,virt
           wia2(e,n) = wia2(e,n) + t2(i,j,c,e) * rint(c+occ,n,i,j)
     2               * fac(i)*fac(j)*fac(c+occ)
          end do
         end do
        end do
       end do
      end do
c
      sum = ddot(virt*occ,wia,1,wia2,1)
c
      cct_so_a1b = sum
c
      return 
      end
      double precision function cct_so_b1b(rint,nbasis,t2,occ,virt,fac,
     1                                     wij,wij2)
c
c     term B1B 
c       
c       t(km,ab) (cm||ij) t(kn,ab) (cn||ij)
c     = t(km,ab) t(kn,ab) (cm||ij) (cn||ij)
c
c     weighting 0.25
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision fac(nbasis)
c
      integer  m,n,k,a,b,j,i,c
      double precision wij(occ,occ),
     1                 wij2(occ,occ),sum,ddot
c
c...  (abk) t(km,ab) t(kn,ab) => (m,n)
c
      do m=1,occ
       do n=1,occ
        wij(n,m) = 0.0d0
        do k=1,occ
         do a=1,virt
          do b=1,virt
            wij(n,m) = wij(n,m) + t2(k,m,a,b) * t2(k,n,a,b)
     2               * fac(k)*fac(a+occ)*fac(b+occ)
          end do
         end do
        end do
       end do
      end do
c
c...  (cij) (cn||ij)  (cm||ij) => (m,n)
c
      do m=1,occ
       do n=1,occ
        wij2(n,m) = 0.0d0
        do j=1,occ
         do i=1,occ
          do c=1,virt
           wij2(n,m) = wij2(n,m) + rint(c+occ,n,i,j) * rint(c+occ,m,i,j)
     2               * fac(i)*fac(j)*fac(c+occ)
          end do
         end do
        end do
       end do
      end do
c
      sum = ddot(occ*occ,wij,1,wij2,1)
c
      cct_so_b1b = sum
c
      return 
      end
#if(0)
      subroutine w_ii(rint,nbasis,occ,virt,
     1                  w,con)
c
c...  calculated the integral only contracted intermediate
c...  (ab||ei)*(ab||fj) => (ij,ef)  ('ab')
c
      implicit none
#include "errquit.fh"
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      character*(*) con
c
      double precision w(occ,occ,virt,virt)
c
      if (con.eq.'ab') then
c
       do e=1,virt
        do f=1,virt
         do i=1,occ
          do j=1,occ
c
c..      wab(occ,occ,virt,virt)=sum(ab)(ab||ei)*(ab||fj)
c
           w(i,j,e,f) = 0.0d0
           do a=1,virt
            do b=1,virt
             w(i,j,e,f) = rint(a+occ,b+occ,e+occ,i)
     1                    * rint(a+occ,b+occ,f+occ,j)
            end do
           end do
c
          end do
         end do
        end do
       end do
c
      else if (con.eq.'ij') then
c
       do e=1,virt
        do f=1,virt
         do i=1,occ
          do j=1,occ
c
c..      w(occ,occ,virt,virt)=sum(ab)(ab||ei)*(ab||fj)
c
           w(i,j,e,f) = 0.0d0
           do a=1,occ
            do b=1,occ
             w(i,j,e,f) = rint(a,b,e+occ,i)
     1                    * rint(a,b,f+occ,j)
            end do
           end do
c
          end do
         end do
        end do
       end do
c
      else 
         call errquit(' what ?',0, UNKNOWN_ERR)
      end if
c
      end
      subroutine cct_tr_ij_vv_kc(rint,nbasis,t2,occ,virt,w)
c
c...  calculated t2(ij,vv) (vv//kc) contracted intermediate
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision w(i,j,k,c)
c
c
       do i=1,occ
        do j=1,occ
         do k=1,occ
          do c=1,virt
c
           w(i,j,k,c) = 0.0d0
           do a=1,virt
            do b=1,virt
             w(i,j,k,c) = rint(a+occ,b+occ,k,c+occ)
     1                  * t2(i,j,a,b)
            end do
           end do
c
          end do
         end do
        end do
       end do
c
      end
      subroutine cct_rt_avbo_oicv(rint,nbasis,t2,occ,virt,wiabc)
c
c...  calculated  intermediate
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t2(occ,occ,virt,virt)
      double precision w(i,j,k,c)
c
      integer i,a,b,c,v,o
c
      do i=1,occ
       do a=1,virt
        do b=1,virt
         do c=1,virt
c
           w(i,a,b,c) = 0.0d0
           do o=1,occ
            do v=1,virt
             w(i,a,b,c) = rint(a+occ,v+occ,b+occ,o)
     1                  * t2(o,i,c,v)
            end do
           end do
c
          end do
         end do
        end do
       end do
c
      end
#endif
