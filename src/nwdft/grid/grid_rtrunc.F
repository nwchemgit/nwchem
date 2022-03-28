      subroutine grid_rtrunc(rtdb, zprim, 
     .     llind,exex,shiftm,ll,
     &     tol_reset, range_max,check_s)
c 
c$Id$
c
      implicit none
#include "errquit.fh"
c      
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "bas.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "util.fh"
c
      integer rtdb, me,nproc
      logical tol_reset
      logical llind,check_s ! [in]
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
      integer iat,
     .     lpmat,ipmat,l_ics,k_ics,iat_diff,k_coord,l_coord
c
      double precision alphaa,alphab,xa(3),xb(3)
      integer lla,llb,itype,maxangg
      integer xc_leb_nrd2acc
      external xc_leb_nrd2acc
      logical oprint
      integer i
c
c     Set print options.
c
      oprint = util_print('griddebug', print_high)
      
c
c     determine node id
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     preliminary ops
c
      call grid_rtinit(rtdb,k_coord,l_coord,
     .     ipmat, lpmat, l_ics,k_ics,g_s)
      call ifill(2*ntypes, 0, rad_cutoff, 1)
      if(llind) then
         call ifill(ntypes, 0, ll, 1)
         call dfill(ntypes, 0d0, exex, 1)
         call dfill(ntypes, 0d0, shiftm, 1)
      endif
c
c     get desired grid accuracy
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     &   call errquit('rad_trunc: rtdb_get failed', 911, RTDB_ERR)
c
c   accuracy must be at least 10-7
c
c      job_grid_acc=min(job_grid_acc,1d-7)
      
      call grid_rtrunc0(zprim,job_grid_acc,int_mb(k_ics))
cc AJL/Begin/FDE
      if (frozemb_fde)
     &   call grid_rtrunc0_fde(zprim, job_grid_acc, int_mb(k_ics))
cc AJL/End
      range_max=0d0
      
c
c     loop over atoms again
c
      do iat=1+me,ncenters,nproc
         if(.not.iatype_pt_chg(iat))then

         call grid_rtrunc1(iat,
     .        g_s, dbl_mb(ipmat),job_grid_acc,iat_diff,range_max,
     .        dbl_mb(k_coord),
     .        ll_at,exex_at,shiftm_2g,oprint,LuOut)
         if (tol_reset)then
            if((nint(range_max)).gt.rad_cutoff(1,iatype(iat))) then
               if(llind) then
                  ll(iatype(iat)) = ll_at
                  exex(iatype(iat)) = exex_at
                  shiftm(iatype(iat)) = shiftm_2g
               endif
               rad_cutoff(1,iatype(iat)) = nint(range_max) 
               rad_cutoff(2,iatype(iat)) = iat_diff
c               rad_cutoff(2,iatype(iat)) = int_mb(k_ics+iat-1)
            endif
         endif
         else
               rad_cutoff(1,iatype(iat)) = 0
               rad_cutoff(2,iatype(iat)) = 0
         endif
      enddo
      call ga_igop(msg_igop_tru,rad_cutoff,2*ntypes,'max') 
      call ga_dgop(msg_igop_tru3,range_max,1,'max') 
      if(llind) then
         call ga_dgop(msg_igop_tru1,exex,ntypes,'max') 
         call ga_igop(msg_igop_tru2,ll,ntypes,'max') 
         call ga_dgop(msg_igop_tru4,shiftm,ntypes,'max') 
      endif
      if(check_s) then 
      call ga_igop(msg_igop_tru5,ldiff,ncenters,'max')
      call ga_dgop(msg_igop_tru6,minexp,ncenters,'max')

         do itype=1,ntypes
            iat=ityp2ctr(itype)
            iat_diff=rad_cutoff(2,itype)
            if(me.eq.0.and.oprint)
     .           write(LuOut,*) ' A and B are ',iat,iat_diff
            alphab=minexp(iat_diff)
            alphaa=exex(itype)-alphab
            llb=ldiff(iat_diff)
            lla=ll(iatype(iat))-llb
            xa(1)=dbl_mb(k_coord+3*(iat-1))
            xa(2)=dbl_mb(k_coord+3*(iat-1)+1)
            xa(3)=dbl_mb(k_coord+3*(iat-1)+2)
            xb(1)=dbl_mb(k_coord+3*(iat_diff-1))
            xb(2)=dbl_mb(k_coord+3*(iat_diff-1)+1)
            xb(3)=dbl_mb(k_coord+3*(iat_diff-1)+2)
            
            call grid_checks(job_grid_acc,iat,alphaa,lla,xa,
     .           iat_diff,alphab,llb,xb)
         enddo
c
c     recheck range_max
c
         do itype=1,ntypes
            range_max=max(range_max,dble(rad_cutoff(itype,1)))
         enddo
c
c     recheck iangquad
c
         if(leb) then
            maxangg=0
            do itype=1,ntypes
               maxangg=max(maxangg,nang(itype))
            enddo
            iangquad=xc_leb_nrd2acc(maxangg)
         endif
      endif
      
c
c     check input flag to see if tolerance modifications desired:
c     accAOfunc, accqrad, radius
c
      call  grid_rtclea(l_coord,
     .     lpmat,l_ics,g_s)

c      if(me.eq.0) 
c     .     call i_output(rad_cutoff, 1, 2, 1, ntypes, 2, ntypes, 1)
      if (tol_reset)then
         radius = range_max + 3.0d0        
         if (.not. rtdb_put(rtdb, 'dft:radius',
     &      mt_dbl, 1,radius))
     &      call errquit('rad_trunc: rtdb_put failed', 213, RTDB_ERR)
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
      double precision value,valexp,dexpo
#if defined(WIN32) ||defined(LINUX)
! DERFC is a GNU extension
      double precision derfc
#else
! ERFC is Fortran 2008
      double precision erfc
#endif
c
      ilo = mod(k,2)
c
         
      if (ilo .eq. 0) then
         dexpo=sqrt(alpha)*r
         if(dexpo.gt.40d0) then
            value=0d0
         else
#if defined(WIN32) ||defined(LINUX)
            value = 0.5d0*sqrt(4.0d0*atan(1.0d0)/alpha)*
     $           derfc(dexpo)
#else 
            value = 0.5d0*sqrt(4.0d0*atan(1.0d0)/alpha)*
     $           erfc(dexpo)
#endif
         endif
      endif
      if(alpha*r*r.gt.500d0) then
         valexp=0d0
         value =0d0
      else
         valexp=exp(-alpha*r*r)
         value =valexp/(2.0d0*alpha)
      endif
c
      do i = ilo+2,k,2
            value = ((i-1)*value + r**(i-1)*valexp)/
     /           (2.0d0*alpha)
      enddo
c
      grid_ik = value
c
      end
      subroutine grid_rtrunc1(iat,g_s,pmat,my_acc,
     .     iat_diff,range_max,
     .     coord,
     .     ll_at,exex_at,shiftm_2g,oprint,LuOut)
      implicit none
#include "cdft.fh"
      integer iat ! [in]
      integer iat_diff ! [out]
      integer g_s ! [in]
      double precision coord(3,*) ! [in]
      double precision pmat(*) ! [local]
      double precision my_acc ! [in]
      integer ll_at ! [out]
      double precision exex_at ! [out]
      double precision shiftm_2g ! [out]
      logical oprint ! [in]
      integer LuOut ! [in]
      
c
      integer jat,dimi,dimj,l_iat
      double precision s_max,dabsmax,ov_trunc,exp_iat
      double precision range,range_max
      external dabsmax
      parameter(ov_trunc=1d-5)
      external grid_rnegl
      double precision grid_rnegl,shift_2g,rrr
c
      range_max=0d0
      shiftm_2g=0d0
      l_iat=ldiff(iat)
      exp_iat=minexp(iat)
         do jat=1,ncenters
            s_max=10d0
            if(iat.ne.jat) then
               call get_atom_block(g_s,ao_bas_han,
     $              iat, jat, pmat, dimi, dimj)
               s_max=dabsmax(dimi*dimj,pmat)
            endif
            if(abs(s_max).gt.ov_trunc) then
              if(oprint) 
     .              write(LuOut,*) 
     .              ' atom ',iat,' - neighbor  & S ',jat,
     .              s_max,' exp AB ',exp_iat,minexp(jat)
              range=grid_rnegl(l_iat+ldiff(jat)+2,
     .             exp_iat+minexp(jat),    my_acc) 
c
c     bit to add from gaussian product formula
c               
              rrr=   sqrt(
     .             (coord(1,iat)-coord(1,jat))*
     .             (coord(1,iat)-coord(1,jat))+
     .              (coord(2,iat)-coord(2,jat))*
     .              (coord(2,iat)-coord(2,jat))+
     .              (coord(3,iat)-coord(3,jat))*
     .              (coord(3,iat)-coord(3,jat)))
               shift_2g=rrr*exp_iat/
     /              (exp_iat+minexp(jat))
c               if(shift_2g.gt.1e-12) shift_2g=shift_2g*2d0
               range=range+shift_2g
c
c     this hack gets <a|a> for nearest neighbors
c
c               if(iat.ne.jat) range=range+6d0/rrr
               if(iat.ne.jat) range=range+rrr
               if(range.gt.range_max) then
                  shiftm_2g=shift_2g
                  range_max =  range
                  iat_diff=jat
                  ll_at=l_iat+ldiff(jat)
                  exex_at=exp_iat+minexp(jat)
               endif
            endif
         enddo
         if(oprint)
     .    write(LuOut,1234)' zmins = ',exp_iat,minexp(iat_diff),
     ,        ' lprod = ',l_iat+ldiff(iat_diff),
     ,        ' range = ',range_max, ' acc ',my_acc,
     ,        ' sh ',shiftm_2g
         return
 1234          format(A,2f15.9,A,I3,A,f12.2,A,d10.2)
         end
      subroutine grid_rtrunc0(zprim,job_grid_acc,cs)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "cdft.fh"
c
      double precision zprim(nbf_ao_mxprim)! [local]
      integer cs(ncenters)! [out]
      double precision job_grid_acc ! [in]
c
      integer iat,ic1,ic2,l,nprimo,ncontr,isphere,iprimo,
     .     icset,ncontrset
      double precision min_exp,range_max,range
      double precision grid_rnegl
      character*16 tag
      double precision ictr_coord(3), ictr_chg
      external grid_rnegl
c
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in rad_trunc',1, BASIS_ERR)
c
c     loop over atoms to get min exp for each atom
c
      do iat=1,ncenters
         if(.not.iatype_pt_chg(iat))then

c
c     get minimum exponent of all basis functions
c
         if (.not.bas_ce2cnr(ao_bas_han,iat,ic1,ic2))
     &        call errquit('Exiting in radtur.',11, BASIS_ERR)
         if(geom_isbq(geom,iat).and.bqdontcare) then
            range_max=10d0
            minexp(iat)=.075d0
            ldiff(iat)=2
         else
            if(ic2.lt.0) then
               if (.not. geom_cent_get(geom, iat, tag, 
     &              ictr_coord, ictr_chg))call errquit
     &              ('grid_atom_type_info: geom_cent_get failed', 0,
     &              GEOM_ERR)
               write(0,*) ' atom   ',iat
               write(0,*) ' tag    ',tag
               write(0,'(A,3F12.6)') ' coord  ',ictr_coord
               write(0,*) ' charge ',ictr_chg
               write(0,*) ' ic1 ic2 ',ic1,ic2
               call errquit(
     .         'gridrtunc: no basis on ghost atom?',iat, BASIS_ERR)
            endif
         endif ! geom_isbq
         range_max=0d0
         do icset = ic1,ic2
            min_exp=1d20
            if (.not.bas_continfo(AO_bas_han, icset,
     &           l, nprimo, ncontr, isphere))
     &           call errquit('Exiting in rad_trunc.',3, BASIS_ERR)
            if(nprimo.gt.nbf_ao_mxprim) call errquit(
     .           'rtrunc: incorrect dimens for zprim',nprimo,
     &       UNKNOWN_ERR)
            if (.not.bas_get_exponent(AO_bas_han, icset, zprim))
     &           call errquit('Exiting in rad_trunc.',7, BASIS_ERR)
            do  iprimo = 1,nprimo
               min_exp = min(zprim(iprimo),min_exp)
            enddo
            range = grid_rnegl(l+l+2, min_exp+min_exp, job_grid_acc)
            if(range.gt.range_max) then
               range_max =  range
               cs(iat)=icset
               minexp(iat)=min_exp
               ldiff(iat)=l
            endif
         enddo
      endif
      enddo
      return
      end
cc AJL/Begin/FDE
      subroutine grid_rtrunc0_fde(zprim,job_grid_acc,cs)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "cdft.fh"
c
      double precision zprim(nbf_ao_mxprim)! [local]
      integer cs(ncenters)! [out]
      double precision job_grid_acc ! [in]
c
      integer iat,ic1,ic2,l,nprimo,ncontr,isphere,iprimo,
     .     icset,ncontrset
      double precision min_exp,range_max,range
      double precision grid_rnegl
      character*16 tag
      double precision ictr_coord(3), ictr_chg
      external grid_rnegl
c
      integer ncenters_fde
c
      if (.not.bas_numcont(FD_bas_han, ncontrset))
     &  call errquit('Exiting in grid_rtrunc0_fde',1, BASIS_ERR)
c
c Get the number of FDE centers
c
      if (.not. geom_ncent(geom_fde, ncenters_fde))
     &        call errquit('grid_rtrunc0_fde: 
     &                      geom_ncent failed',73, GEOM_ERR)
c
c     loop over atoms to get min exp for each atom
c
      do iat=1,ncenters_fde
         if(.not.iatype_pt_chg_fde(iat))then

c
c     get minimum exponent of all basis functions
c
         if (.not.bas_ce2cnr(FD_bas_han,iat,ic1,ic2))
     &        call errquit('Exiting in radtur.',11, BASIS_ERR)
         if(geom_isbq(geom_fde,iat).and.bqdontcare) then
            range_max=10d0
            minexp_fde(iat)=.075d0
            ldiff_fde(iat)=2
         else
           if(ic2.lt.0) then
              if (.not. geom_cent_get(geom_fde, iat, tag,
     &             ictr_coord, ictr_chg))call errquit
     &             ('grid_rtrunc0_fde: geom_cent_get failed', 0,
     &         GEOM_ERR)
              write(0,*) ' atom   ',iat
              write(0,*) ' tag    ',tag
              write(0,'(A,3F12.6)') ' coord  ',ictr_coord
              write(0,*) ' charge ',ictr_chg
              write(0,*) ' ic1 ic2 ',ic1,ic2
              call errquit(
     .        'grid_rtunc_fde: no basis on ghost atom?',iat, BASIS_ERR)
           endif
         endif ! geom_isbq
         range_max=0d0
         do icset = ic1,ic2
            min_exp=1d20
            if (.not.bas_continfo(FD_bas_han, icset,
     &           l, nprimo, ncontr, isphere))
     &           call errquit('Exiting in rad_trunc.',3, BASIS_ERR)
            if(nprimo.gt.nbf_fde_mxprim) call errquit(
     .           'grid_rtrunc0_fde: incorrect dimens for zprim',nprimo,
     &       UNKNOWN_ERR)
            if (.not.bas_get_exponent(FD_bas_han, icset, zprim))
     &           call errquit('Exiting in rad_trunc.',7, BASIS_ERR)
            do  iprimo = 1,nprimo
               min_exp = min(zprim(iprimo),min_exp)
            enddo
            range = grid_rnegl(l+l+2, min_exp+min_exp, job_grid_acc)
            if(range.gt.range_max) then
               range_max =  range
c               cs(iat)=icset
               minexp_fde(iat)=min_exp
               ldiff_fde(iat)=l
            endif
         enddo
      endif
      enddo
      return
      end
cc AJL/End
      subroutine grid_rtinit(rtdb,k_coord,l_coord,
     .     ipmat, lpmat,l_ics,k_ics,g_s)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "cdft.fh"
c
      integer rtdb ! [in]
      integer k_coord,l_coord
      integer l_ics,k_ics,
     .   ipmat, lpmat,g_s ! [out]
c
      integer itags,ltags,k_charge,l_charge
      external ga_create_atom_blocked
      integer ga_create_atom_blocked
c
      if (.not.MA_Push_Get(MT_int, ncenters, 'diffuse ic',
     &     l_ics, k_ics))
     &         call errquit('rdatru: cannot allocate degen minexp',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce,
     ,     'matrix patch',lPmat,iPmat))
     &   call errquit('rdatru: cannot allocate DENS',0, MA_ERR)
c
c     get overlap matrix
c
      g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      call ga_zero(g_s)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_s,'overlap',oskel)
      if (oskel)call sym_symmetrize(geom, AO_bas_han, .false., g_s)
c
c     get coord
c
      if (.not.MA_Push_Get(MT_dbl,3*ncenters,
     ,     'coord       ',l_coord,k_coord))
     &   call errquit('rdatru: cannot allocate DENS',1, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,ncenters,
     ,     'coord       ',l_charge,k_charge))
     &   call errquit('rdatru: cannot allocate DENS',2, MA_ERR)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit(
     .     'gridtru: failed to alloc center tags',0, MA_ERR)

      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(k_coord), Dbl_MB(k_charge)))
     &     call errquit('grid_atom_type_info: geom_cart_get failed',74,
     &       GEOM_ERR)
      if (.not.ma_pop_stack(ltags))
     &   call errquit('radtru: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_charge))
     &   call errquit('radtru: cannot pop stack',0, MA_ERR)
      
      return
      end
      subroutine grid_rtclea(l_coord,
     .     lpmat,l_ics,g_s)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer  l_ics, l_coord,
     .       lpmat,g_s ! [in]
c
      if (.not. ga_destroy(g_s)) call errquit
     &     ('radtrun: could not destroy g_s', 0, GA_ERR)
      if (.not.ma_pop_stack(l_coord))
     &   call errquit('radtru: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lpmat))
     &   call errquit('radtru: cannot pop stack',1, MA_ERR)
      if (.not.ma_pop_stack(l_ics))
     &   call errquit('radtru: cannot pop stack',2, MA_ERR)
      return
      end
