      subroutine grid_rtrunc(rtdb, zprim, 
     .     llind,exex,shiftm,ll,
     &     tol_reset, range_max)
c
C$Id:
c
      implicit none
c      
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "bas.fh"
#include "cdft.fh"
c
      integer rtdb, me,nproc
      logical tol_reset
      logical llind ! [in]
      integer ll(*) ! [out] for lindh grid
      double precision exex(*) ! [out] for lindh grid
      double precision shiftm(*) ! [out] for lindh grid
      double precision job_grid_acc
      double precision range_max
      double precision zprim(nbf_ao_mxprim)
c
      integer ll_at
      double precision exex_at,shiftm_2g
      integer g_s
      integer k_minexp,l_minexp,k_ldiff,l_ldiff,iat,
     .     lpmat,ipmat,l_ics,k_ics,iat_diff,k_coord,l_coord
c
c     determine node id
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     preliminary ops
c
      call grid_rtinit(rtdb,l_minexp, k_minexp,k_coord,l_coord,
     .     l_ldiff, k_ldiff,ipmat, lpmat, l_ics,k_ics,g_s)
      call ifill(2*ntypes, 0, rad_cutoff, 1)
      call ifill(ntypes, 0, ll, 1)
      call dfill(ntypes, 0d0, exex, 1)
c
c     get desired grid accuracy
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     &   call errquit('rad_trunc: rtdb_get failed', 911)
      call grid_rtrunc0(zprim,job_grid_acc,
     ,     int_mb(k_ics),int_mb(k_ldiff),dbl_mb(k_minexp))
      
c
c     loop over atoms again
c
      do iat=1+me,ncenters,nproc
         call grid_rtrunc1(iat,dbl_mb(k_minexp),int_mb(k_ldiff),
     .        g_s, dbl_mb(ipmat),job_grid_acc,iat_diff,range_max,
     .        dbl_mb(k_coord),
     .        ll_at,exex_at,shiftm_2g)
         if (tol_reset)then
            if((nint(range_max)).gt.rad_cutoff(1,iatype(iat))) then
               if(llind) then
                  ll(iatype(iat)) = ll_at
                  exex(iatype(iat)) = exex_at
                  shiftm(iatype(iat)) = shiftm_2g
               endif
               rad_cutoff(1,iatype(iat)) = nint(range_max) 
               rad_cutoff(2,iatype(iat)) = int_mb(k_ics+iat-1)
            endif
         endif
      enddo
c
c     check input flag to see if tolerance modifications desired:
c     accAOfunc, accqrad, radius
c
      call  grid_rtclea(l_minexp, l_coord,
     .     l_ldiff, lpmat,l_ics,g_s)

      call ga_igop(msg_igop_tru,rad_cutoff,2*ntypes,'absmax') 
      call ga_dgop(msg_igop_tru-1,exex,ntypes,'absmax') 
      call ga_igop(msg_igop_tru-2,ll,ntypes,'absmax') 
      if(me.eq.0) 
     .     call i_output(rad_cutoff, 1, 2, 1, ntypes, 2, ntypes, 1)
      if (tol_reset)then
         radius = range_max + 3.0d0        
         if (.not. rtdb_put(rtdb, 'dft:radius',
     &      mt_dbl, 1,radius))
     &      call errquit('rad_trunc: rtdb_put failed', 213)
      endif
c
      return
      end
      double precision function grid_rnegl(k, alpha, eps)
      implicit none
c
      integer k
      double precision alpha, eps
c
c     For a function f(r) = r^k*exp(-alpha*r^2) determine
c     the radial distance r such that the fraction of the 
c     function norm that is neglected if the 3D volume 
c     integration is terminated at a distance r is less
c     than or equal to eps.
c
      double precision r, test, grid_rnegl0, step
      external grid_rnegl0
c
      step = 0.5d0
      r = 1.0d0
 10   test = grid_rnegl0(k,alpha,r)
      if (test .gt. eps) then
         r = r + step
      else 
         r = r - step
         if (r .lt. 0.0d0) r = 0.0d0
         step = step*0.5d0
         r = r + step
      endif
      if (step .gt. 0.01d0) goto 10
c
      grid_rnegl = r
c
      end
      double precision function grid_rnegl0(k,alpha,r)
      implicit none
c
      integer k
      double precision alpha, r
c
c     For a function f(r) = r^k*exp(-alpha*r^2) determine
c     the fraction of the function norm that is neglected
c     if the 3D volume integration is terminated at a 
c     distance r.
c
c     neglected = int(t^2*f(t),t=r..infinity)/int(t^2*f(t),t=0..infinity)
c
      double precision grid_ik
      external grid_ik
c
      grid_rnegl0 = grid_ik(k+2,alpha,r)/grid_ik(k+2,alpha,0.0d0)
c
      end
      double precision function grid_ik(k,alpha,r)
      implicit none
c
      integer k
      double precision alpha, r
c
c     I(k) = int(t^k exp(-alpha*t^2), t=0..infinity)
c
c     I(k) = [(k-1)*I(k-2) + r^(k-1)*exp(-alpha*r^2)]/(2*alpha)
c
      integer i, ilo
      double precision value
#if defined(SGI) || defined(DECOSF) || defined(WIN32)
      double precision derfc
#else
      double precision erfc
#endif
c
      ilo = mod(k,2)
c
      if (ilo .eq. 0) then
#if defined(SGI) || defined(DECOSF) || defined(WIN32)
         value = 0.5d0*sqrt(4.0d0*atan(1.0d0)/alpha)*
     $        derfc(sqrt(alpha)*r)
#else
         value = 0.5d0*sqrt(4.0d0*atan(1.0d0)/alpha)*
     $        erfc(sqrt(alpha)*r)
      else
#endif
         value = exp(-alpha*r*r)/(2.0d0*alpha)
      endif
c
      do i = ilo+2,k,2
         value = ((i-1)*value + r**(i-1)*exp(-alpha*r*r))/(2.0d0*alpha)
      enddo
c
      grid_ik = value
c
      end
      subroutine grid_rtrunc1(iat,expo,l,g_s,pmat,my_acc,
     .     iat_diff,range_max,
     .     coord,
     .     ll_at,exex_at,shiftm_2g)
      implicit none
#include "cdft.fh"
      integer iat ! [in]
      integer iat_diff ! [out]
      integer g_s ! [in]
      double precision expo(*) ! [in]
      double precision coord(3,*) ! [in]
      integer l(*) ! [in]
      double precision pmat(*) ! [local]
      double precision my_acc ! [in]
      integer ll_at ! [out]
      double precision exex_at ! [out]
      double precision shiftm_2g ! [out]
      
c
      integer jat,dimi,dimj,l_iat
      double precision s_max,dabsmax,ov_trunc,exp_iat
      double precision range,range_max
      external dabsmax
      parameter(ov_trunc=1d-5)
      external grid_rnegl
      double precision grid_rnegl,shift_2g
c
      range_max=0d0
      l_iat=l(iat)
      exp_iat=expo(iat)
         do jat=1,ncenters
            s_max=10d0
            if(iat.ne.jat) then
               call get_atom_block(g_s,ao_bas_han,
     $              iat, jat, pmat, dimi, dimj)
               s_max=dabsmax(dimi*dimj,pmat)
            endif
            if(abs(s_max).gt.ov_trunc) then
               write(0,*) ' atom ',iat,' - neighbor  & S ',jat,s_max
               range=grid_rnegl(l_iat+l(jat)+2,exp_iat+expo(jat),
     .              my_acc) 
c
c     bit to add from gaussian product formula
c               
               shift_2g=sqrt(
     .              (coord(1,iat)-coord(1,jat))**2+
     .              (coord(2,iat)-coord(2,jat))**2+
     .              (coord(3,iat)-coord(3,jat))**2)*exp_iat/
     /              (exp_iat+expo(jat))
               range=range+shift_2g
cedo               write(0,*) ' expA expB shift_2g ',exp_iat,
cedo     .              expo(jat),shift_2g
               if(range.gt.range_max) then
                  shiftm_2g=shift_2g
                  range_max =  range
                  iat_diff=jat
                  ll_at=l_iat+l(jat)
                  exex_at=exp_iat+expo(jat)
               endif
            endif
         enddo
         write(0,1234)' zmin = ',exp_iat*expo(iat_diff),
     ,        ' lprod = ',l_iat*l(iat_diff),
     ,        ' range = ',range_max, ' acc ',my_acc
         return
 1234          format(A,f15.9,A,I3,A,f12.2,A,d10.2)
         end
      subroutine grid_rtrunc0(zprim,job_grid_acc,cs,ldiff,mexp)
      implicit none
#include "bas.fh"
#include "cdft.fh"
c
      double precision zprim(nbf_ao_mxprim)! [local]
      integer cs(ncenters)! [out]
      integer ldiff(ncenters)! [out]
      double precision mexp(ncenters)! [out]
      double precision job_grid_acc ! [in]
c
      integer iat,ic1,ic2,l,nprimo,ncontr,isphere,iprimo,
     .     icset,ncontrset
      double precision min_exp,range_max,range
      double precision grid_rnegl
      external grid_rnegl
c
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in rad_trunc',1)
c
c     loop over atoms to get min exp for each atom
c
      do iat=1,ncenters
c
c     get minimum exponent of all basis functions
c
         if (.not.bas_ce2cnr(ao_bas_han,iat,ic1,ic2))
     &        call errquit('Exiting in radtur.',11)
         min_exp=1d20
         range_max=0d0
         do icset = ic1,ic2
            if (.not.bas_continfo(AO_bas_han, icset,
     &           l, nprimo, ncontr, isphere))
     &           call errquit('Exiting in rad_trunc.',3)
            if(nprimo.gt.nbf_ao_mxprim) call errquit(
     .           'rtrunc: incorrect dimens for zprim',nprimo)
            if (.not.bas_get_exponent(AO_bas_han, icset, zprim))
     &           call errquit('Exiting in rad_trunc.',7)
            do  iprimo = 1,nprimo
               min_exp = min(zprim(iprimo),min_exp)
            enddo
            range = grid_rnegl(l+l+2, min_exp+min_exp, job_grid_acc)
            if(range.gt.range_max) then
               range_max =  range
               cs(iat)=icset
               mexp(iat)=min_exp
               ldiff(iat)=l
            endif
         enddo
      enddo
      return
      end
      subroutine grid_rtinit(rtdb,l_minexp, k_minexp,k_coord,l_coord,
     .     l_ldiff, k_ldiff,ipmat, lpmat,l_ics,k_ics,g_s)
      implicit none
#include "mafdecls.fh"
#include "geom.fh"
#include "cdft.fh"
c
      integer rtdb ! [in]
      integer k_coord,l_coord
      integer l_minexp, k_minexp, l_ics,k_ics,
     .     l_ldiff, k_ldiff,ipmat, lpmat,g_s ! [out]
c
      integer itags,ltags,k_charge,l_charge
      external ga_create_atom_blocked
      integer ga_create_atom_blocked
c
      if (.not.MA_Push_Get(MT_dbl, ncenters, 'minexp',
     &     l_minexp, k_minexp))
     &         call errquit('rdatru: cannot allocate degen minexp',0)
      if (.not.MA_Push_Get(MT_int, ncenters, 'diffuse l',
     &     l_ldiff, k_ldiff))
     &         call errquit('rdatru: cannot allocate degen minexp',0)
      if (.not.MA_Push_Get(MT_int, ncenters, 'diffuse ic',
     &     l_ics, k_ics))
     &         call errquit('rdatru: cannot allocate degen minexp',0)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce,
     ,     'matrix patch',lPmat,iPmat))
     &   call errquit('rdatru: cannot allocate DENS',0)
c
c     get overlap matrix
c
      call int_init(rtdb, 1, AO_bas_han) 
      g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      call int_1e_ga(ao_bas_han,ao_bas_han,g_s,'overlap',.false.)
      if (oskel)call sym_symmetrize(geom, AO_bas_han, .false., g_s)
      call int_terminate
c
c     get coord
c
      if (.not.MA_Push_Get(MT_dbl,3*ncenters,
     ,     'coord       ',l_coord,k_coord))
     &   call errquit('rdatru: cannot allocate DENS',1)
      if (.not.MA_Push_Get(MT_dbl,ncenters,
     ,     'coord       ',l_charge,k_charge))
     &   call errquit('rdatru: cannot allocate DENS',2)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit('gridtru: failed to 
     &   alloc center tags',0)

      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(k_coord), Dbl_MB(k_charge)))
     &     call errquit('grid_atom_type_info: geom_cart_get failed',74)
      if (.not.ma_pop_stack(ltags))
     &   call errquit('radtru: cannot pop stack',0)
      if (.not.ma_pop_stack(l_charge))
     &   call errquit('radtru: cannot pop stack',0)
      
      return
      end
      subroutine grid_rtclea(l_minexp, l_coord,
     .     l_ldiff, lpmat,l_ics,g_s)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer l_minexp,  l_ics, l_coord,
     .      l_ldiff, lpmat,g_s ! [in]
c
      if (.not. ga_destroy(g_s)) call errquit
     &     ('radtrun: could not destroy g_s', 0)
      if (.not.ma_pop_stack(l_coord))
     &   call errquit('radtru: cannot pop stack',0)
      if (.not.ma_pop_stack(lpmat))
     &   call errquit('radtru: cannot pop stack',1)
      if (.not.ma_pop_stack(l_ics))
     &   call errquit('radtru: cannot pop stack',2)
      if (.not.ma_pop_stack(l_ldiff))
     &   call errquit('radtru: cannot pop stack',3)
      if (.not.ma_pop_stack(l_minexp))
     &   call errquit('radtru: cannot pop stack',4)
      return
      end
