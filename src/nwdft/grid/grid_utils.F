      subroutine qshell_sort
c
c$Id: grid_utils.F,v 1.14 2000-03-20 21:11:50 edo Exp $
c
c
c     Sort quadrature shells, most work to least, based on radius 
c     (assuming largest radius has most angular pieces).
c
      implicit none
#include "cdft.fh"
      integer temp1, temp2, temp3, temp4
      integer i, hsize
c
c     build heap
c
      do i = nqshells/2, 1, -1
         call heapify(nqshells, i)
      enddo
c
c     main part of sort algorithm
c
      hsize = nqshells
      do i    = nqshells, 2, -1
c
c        swap element i and 1
c
         temp1 = iqshell(1,1)
         temp2 = iqshell(2,1)
         temp3 = iqshell(3,1)
         temp4 = iqshell(4,1)
c
         iqshell(1,1) = iqshell(1,i)
         iqshell(2,1) = iqshell(2,i)
         iqshell(3,1) = iqshell(3,i)
         iqshell(4,1) = iqshell(4,i)
c
         iqshell(1,i) = temp1
         iqshell(2,i) = temp2
         iqshell(3,i) = temp3
         iqshell(4,i) = temp4
c
c        maintain heap property from element 1 down
c
         hsize= hsize - 1
         call heapify(hsize, 1)
c
      enddo
      return
      end 
      subroutine heapify(n, elem)
c
c     establish heap property for a tree branch rooted at elem 
c
      implicit none
c
#include "cdft.fh"
c
      integer n, elem 
      integer left, right, smallest, i
      integer temp1, temp2, temp3, temp4
      integer ictr_left, irsh_left, ictr_i, irsh_i
      double precision rpts_left, rpts_i
      integer ictr_right, irsh_right, ictr_smallest, irsh_smallest
      double precision rpts_right, rpts_smallest
c
      i = elem
c
c     Main Loop
c
100   continue 
        left  = 2*i
        right = 2*i + 1
        if (left. gt. n .and. right .gt. n) return   !we traversed entire branch
c
c       check heap property among element i and its children
c
        ictr_left = iatype(iqshell(3,left))
        irsh_left = iqshell(1,left)
        rpts_left = rpts(irsh_left,ictr_left)
        ictr_i = iatype(iqshell(3,i))
        irsh_i = iqshell(1,i)
        rpts_i = rpts(irsh_i,ictr_i)
        if (left .le. n .and. rpts_left .lt. rpts_i) then
           smallest = left
        else
           smallest = i
        endif
        ictr_right = iatype(iqshell(3,right))
        irsh_right = iqshell(1,right)
        rpts_right = rpts(irsh_right,ictr_right)
        ictr_smallest = iatype(iqshell(3,smallest))
        irsh_smallest = iqshell(1,smallest)
        rpts_smallest = rpts(irsh_smallest,ictr_smallest)
        if (right .le. n .and. rpts_right .lt. rpts_smallest) 
     &     smallest = right
c
        if (smallest .ne. i) then
c
c          swap array elements if smallest is not i
c
           temp1 = iqshell(1,i)
           temp2 = iqshell(2,i)
           temp3 = iqshell(3,i)
           temp4 = iqshell(4,i)
c
           iqshell(1,i) = iqshell(1,smallest)
           iqshell(2,i) = iqshell(2,smallest)
           iqshell(3,i) = iqshell(3,smallest)
           iqshell(4,i) = iqshell(4,smallest)
c
           iqshell(1,smallest) = temp1
           iqshell(2,smallest) = temp2
           iqshell(3,smallest) = temp3
           iqshell(4,smallest) = temp4
c
c          traverse down the tree
c
           i = smallest
      goto 100
      endif
      return
      end
      Subroutine mbf_ao_max(rtdb, qxyz, rq0, zprim, coord)
c
C$Id: grid_utils.F,v 1.14 2000-03-20 21:11:50 edo Exp $
c
      implicit none
c      
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "bas.fh"
c
      integer rtdb
c
c     Cartesian Coordinates of Integration Center
c
      double precision coord(3,ncenters)
c
c     Cartesian Coordinates of Quadrature Points
c
      double precision qxyz(3,*)
c
c     Compute the quadrature points for a given 
c     set of radial shells.
c
      integer irsh, iang, iqsh,  l, ia_ictr
      double precision r
c
c     Distance Squared between Sampling Points and Centers
c
      double precision rq0(ncenters)
c
      double precision zprim(nbf_ao_mxprim)
      integer ncontrset, n1, icset, ictr, itype, nprimo, ncontr, 
     &        isphere, nbf_ang, iprimo
      double precision zmin
      double precision acc_AO_gauss
      integer mbf, nq, mbfnq
      integer m, me,nproc
      integer avail
      double precision x, y, z, r2
      integer nxtask,nt1,nt2
      external nxtask
c
c     determine node id
c
      me = ga_nodeid()
      nproc=ga_nnodes()
c
c
c     Loop over all the radial shells.
c
      acc_AO_gauss = dble(iAOacc)
      max_mbf = 0
      max_pr_mbf = 0
      max_pr_nq = 0
      max_pr_mbfnq = 0
      nt1 = 0
      nt2 = nxtask(nproc,1)
      do 70 iqsh = 1, nqshells
        if(nt1.eq.nt2) then
c
         irsh = iqshell(1,iqsh)
         ictr = iqshell(3,iqsh)
         iang = iqshell(4,iqsh)
c
         ia_ictr = iatype(ictr)
         r = rpts(irsh,ia_ictr)
c         
         nq = 0
         call dfill(ncenters,0.d0,rq0,1)
c
         if (leb.or.daul) then
c         
           nq=nq+ntheta(iang)
         else
           nq=nq+ntheta(iang)*nphi(iang)
c
         endif
         do 40 m = 1, ncenters
               x = coord(1,ictr) - coord(1,m)
               y = coord(2,ictr) - coord(2,m)
               z = coord(3,ictr) - coord(3,m)
               r2 = sqrt(x*x + y*y + z*z)
               rq0(m)=(r2-r)**2
   40    continue
c
         if (.not.bas_numcont(AO_bas_han, ncontrset))
     &      call errquit('Exiting in mbf_ao_max.',1)
c
         n1 = 0
c
         do 60 icset = 1,ncontrset
            if (.not.bas_cn2ce(AO_bas_han, icset, ictr))
     &         call errquit('Exiting in mbf_ao_max.',2)
c         
c           get info about current contraction set
c         
            if (.not.bas_continfo(AO_bas_han, icset,
     &         itype, nprimo, ncontr, isphere))
     &         call errquit('Exiting in mbf_ao_max.',3)
c         
c           angular momentum
c         
            l = 0
            if (itype .lt. 0)then
               call errquit('mbf_ao_max: sp-type orbital not coded', 5)
            else 
               l = itype
            endif 
c         
c           cartesian/spherical harmonic
c         
            nbf_ang = 0
            if (isphere .eq. 0)then !  cartesian set
               nbf_ang = (l+1)*(l+2)/2
            elseif (isphere .eq. 1)then !  spherical harmonic
               nbf_ang = 2*l+1
            else 
               call errquit('mbf_ao_max: illegal isphere value', 6)
            endif 
c         
c           get exponents and contraction coefficients for this contraction set
c         
            call dfill(nbf_ao_mxprim,0.d0,zprim,1)
            if (.not.bas_get_exponent(AO_bas_han, icset, zprim))
     &         call errquit('Exiting in mbf_ao_max.',7)
c   
c           Determine the minimum Gaussian exponent.
c   
            zmin = 1.D+06
            do 50 iprimo = 1,nprimo
               zmin = min(zprim(iprimo),zmin)
   50       continue
c
c           Only include those basis functions that are "non-zero" for at least one 
c           point in the sampling set.
c
            if (zmin*rq0(ictr).gt.acc_AO_gauss)goto 60
            if (l.eq.0)then
c
c              =============>  S Contractions  <============= 
c
               n1 = n1 + ncontr
            elseif (l.eq.1)then
c
c              =============>  P Contractions  <============= 
c
               n1 = n1+ncontr*3
            elseif (l.eq.2)then
c
c              =============>  D Contractions  <============= 
c
               n1 = n1 + ncontr*6
            else
c
c              =============>  General Case  <============= 
c
               n1 = n1 + ncontr*nbf_ang
            endif
c
  60     continue
c
         mbf = n1
c
c        need to determine max_mbf and max(mbf,nq) pair
c
         if (mbf.gt.max_mbf)then
            max_mbf = mbf
         endif
         mbfnq = mbf*nq
         if (mbfnq.gt.max_pr_mbfnq)then
            max_pr_mbfnq = mbfnq
            max_pr_mbf = mbf
            max_pr_nq = nq
         endif

c
            nt1 = nt1 + 1
            nt2 = nxtask(nproc,1)
         else
            nt1 = nt1 + 1
         endif
c
c
   70 continue
      nt2 = nxtask(-nproc,1)
      call ga_igop(dft_mxmbfnq, max_pr_mbfnq, 1, 'max')
      call ga_igop(dft_mxmbf, max_mbf, 1, 'max')
c
c     This can be further reduced by blocking the computed 
c     angular grid.
c
c     Assume we want to use no more than 1/3 of physical memory (stack + heap) for
c     the quadrature.  So, lets put chi(nq,mbf) and delchi(nq,3,mbf)
c     in roughly 1/3 - 8Mb by chunking up nq.
c
c     find - (minimum)amount local available memory on all nodes
c     
      call ga_sync
      avail = MA_inquire_avail(mt_dbl)
      call ga_igop(msg_min_stack_avail, avail, 1, 'min')

c      write(6,*)' avail = ',avail
      avail = avail/3
c      write(6,*)' avail/3 = ',avail
      avail = avail - 1024*1024
      if(avail.lt.0)
     &         call errquit('xc_setquad: out of memory',avail)
c      write(6,*)' amt to be used for xc = ',avail
c     
      if ( (4*max_pr_mbfnq) .gt. avail )then
         nq_chunk = avail/4/max_mbf
c
c        redefine max_pr_mbfnq
c
         max_pr_mbfnq = max_mbf*nq_chunk
c
c        reset store_wght to false (not working yet)
c
         store_wght=.false.
      else
c
c        everything fits no chunking necessary
c
         nq_chunk = 0
c
c        redefine nq_task
c
Cedo         if(nquad_task.eq.1) then
          if (.not. rtdb_get(rtdb, 'dft:nquad_task', mt_int, 1,
     &        nquad_task))then
          nquad_task = min(avail/2/(4*max_pr_mbfnq),6)
          if(nquad_task.lt.1) nquad_task=1
           nqmax=nqmax*nquad_task
           if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1, 
     &          nquad_task))
     &          call errquit('mbf_ao_max: rtdb_put failed', 911)
         endif
         
      endif
c      write(6,*)' nq_chunk = ',nq_chunk
      if (rtdb_get(rtdb, 'dft:nq_chunk', mt_int, 1, nq_chunk))then
         if (me.eq.0)write(LuOut,*)' nq_chunk input override= ',nq_chunk
c
c        redefine max_pr_mbfnq
c
         max_pr_mbfnq = max_mbf*nq_chunk
      endif
c
c     check nquad_task (if .ne. 1 makes no sense with chunking turned on ... reset)
c
      if (nq_chunk.ne.0)then
         if (nquad_task.gt.1)then
            nquad_task = 1
c
         endif
      endif
      if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1, 
     &     nquad_task))
     &     call errquit('mbf_ao_max: rtdb_put failed', 911)
c
c      write(6,*)' iAOacc, acc_AO_gauss: ',iAOacc, acc_AO_gauss
c      write(6,*)' max_mbf, max_pr_mbfnq, max_pr_mbf, max_pr_nq: ',
c     &            max_mbf, max_pr_mbfnq, max_pr_mbf, max_pr_nq
c
      return
      end
      subroutine rad_trunc(rtdb, zprim, shell_chg, 
     &     tol_reset, mul_reset,range_max)
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
      logical tol_reset, mul_reset
      double precision job_grid_acc, my_acc
      double precision range, range_max
      double precision r_neglected
      double precision zprim(nbf_ao_mxprim), shell_chg(nshells_ao)
      integer ncontrset, icset,  nprimo, ncontr, 
     &        isphere, iprimo, l,g_s
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      integer k_minexp,l_minexp,k_ldiff,l_ldiff,iat,ic1,ic2,jat,
     .     dimi,dimj,lpmat,ipmat,l_iat,l_jat,l_ics,k_ics
      double precision min_exp,dabsmax,ov_trunc,s_max,
     ,     exp_iat,exp_jat
      parameter(ov_trunc=1d-5)
      external dabsmax
c
c     determine node id
c
      me = ga_nodeid()
      nproc = ga_nnodes()

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
     &   call errquit('xc_quadv0: cannot allocate DENS',0)
c
c     get desired grid accuracy
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     &   call errquit('rad_trunc: rtdb_get failed', 911)
c
c     get overlap matrix
c
      call int_init(rtdb, 1, AO_bas_han) 
      g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      call int_1e_ga(ao_bas_han,ao_bas_han,g_s,'overlap',.false.)
      call int_terminate
      call ifill(2*ntypes, 0, rad_cutoff, 1)
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
            if (.not.bas_get_exponent(AO_bas_han, icset, zprim))
     &           call errquit('Exiting in rad_trunc.',7)
            do  iprimo = 1,nprimo
               min_exp = min(zprim(iprimo),min_exp)
            enddo
            range = r_neglected(2*l, 2d0*min_exp, job_grid_acc)
            if(range.gt.range_max) then
               range_max =  range
               int_mb(k_ics+iat-1)=icset
               dbl_mb(k_minexp+iat-1)=min_exp
               int_mb(k_ldiff+iat-1)=l
            endif
         enddo
      enddo
      
c
c     loop over atoms again
c
c      do iat=me+1,ncenters,nproc
      do iat=1,ncenters
         write(me,*) 
         write(me,*) ' me ',me,' atom ',iat
         exp_iat=dbl_mb(k_minexp+iat-1)
         l_iat=int_mb(k_ldiff+iat-1)
         range_max=0d0
         do jat=1,ncenters
            s_max=10d0
            if(iat.ne.jat) then
               call get_atom_block(g_s,ao_bas_han,
     $              iat, jat, dbl_mb(ipmat), dimi, dimj)
               s_max=dabsmax(dimi*dimj,dbl_mb(ipmat))
            endif
            if(abs(s_max).gt.ov_trunc) then
               write(me,*) '  - neighbor ',jat,s_max
               exp_jat=dbl_mb(k_minexp+jat-1)
               l_jat=int_mb(k_ldiff+jat-1)

c
               if(mul_reset)then
                  my_acc = job_grid_acc/
     /                 shell_chg(int_mb(k_ics+iat-1))
               else
                  my_acc = job_grid_acc
               endif
               range=r_neglected(l_iat+l_jat,exp_iat+exp_jat,my_acc)
               range_max = max(range_max, range)
               write(me,1234)' zmin = ',exp_iat*exp_jat,
     ,              ' lprod = ',l_iat*l_jat,
     ,              ' range = ',range, ' acc ',my_acc
 1234          format(A,f15.9,A,I3,A,f12.2,A,d10.2)
c     
            endif
         enddo
         if (tol_reset)then
            if((nint(range_max)).gt.rad_cutoff(1,iatype(iat))) then
               rad_cutoff(1,iatype(iat)) = nint(range_max) 
               rad_cutoff(2,iatype(iat)) = int_mb(k_ics+iat-1)
            endif
            write(me,*)' rad_cutoff(1-2,iatype(ictr)): ',
     .           iatype(iat),rad_cutoff(1,iatype(iat)),
     .           rad_cutoff(2,iatype(iat))
         endif
      enddo
cedo      call i_output(rad_cutoff, 1, 2, 1, ntypes, 2, ntypes, 1)
c
c     check input flag to see if tolerance modifications desired:
c     accAOfunc, accqrad, radius
c
      if (.not. ga_destroy(g_s)) call errquit
     &     ('radtrun: could not destroy g_s', 0)
      if (.not.ma_pop_stack(lpmat))
     &   call errquit('radtru: cannot pop stack',0)
      if (.not.ma_pop_stack(l_ics))
     &   call errquit('radtru: cannot pop stack',0)
      if (.not.ma_pop_stack(l_ldiff))
     &   call errquit('radtru: cannot pop stack',1)
      if (.not.ma_pop_stack(l_minexp))
     &   call errquit('radtru: cannot pop stack',2)
      write(me,*) ' me ',me, ' in dgop '
      call setdbg(1)
      call flush(me)
      write(me,*) ' me ',me,msg_igop_tru,2*ntypes
      call ga_igop(msg_igop_tru,2*ntypes,
     ,     rad_cutoff,'absmax') 
      call setdbg(0)
      write(me,*) ' me ',me, ' off dgop '
      if (tol_reset)then
         radius = range_max + 3.0d0        
         if (.not. rtdb_put(rtdb, 'dft:radius',
     &      mt_dbl, 1,radius))
     &      call errquit('rad_trunc: rtdb_put failed', 213)
      endif
c
      return
      end
      double precision function dft_gaussian_range(n, alpha, eps)
      implicit none
c
      integer n
      double precision alpha, eps
c
c     Return an approximation to the outer solution of 
c     .     r^n*exp(-ar^2) = eps
c     .     r = (n*ln(-ln(eps)) - n*ln(a) - 4*ln(eps)) /
c     .         4*sqrt(-alpha*ln(eps))
c
c     Accuracy improves with smaller eps.
c
      double precision logeps
c
      logeps = log(eps)
c      
      dft_gaussian_range = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      end

      double precision function r_neglected(k, alpha, eps)
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
      double precision r, test, neglected, step
c
      step = 0.5d0
      r = 1.0d0
 10   test = neglected(k,alpha,r)
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
      r_neglected = r
c
      end
      double precision function neglected(k,alpha,r)
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
      double precision ik
c
      neglected = ik(k+2,alpha,r)/ik(k+2,alpha,0.0d0)
c
      end
      double precision function ik(k,alpha,r)
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
      ik = value
c
      end

      Subroutine grid_repack(xyzw, qxyz, qwght, nq,rad,istep)
c
C$Id: grid_utils.F,v 1.14 2000-03-20 21:11:50 edo Exp $
c
      implicit none
c
      integer nq
      double precision xyzw(4,nq), qxyz(3,nq), qwght(nq),
     .     rad
c
      integer n,istep
c
      istep=istep+1
      nq=dble(xyzw(1,istep))
      rad=xyzw(2,istep)
      do 30 n = 1, nq
c
        qxyz(1,n) = xyzw(1,n+istep)
        qxyz(2,n) = xyzw(2,n+istep)
        qxyz(3,n) = xyzw(3,n+istep)
c
        qwght(n) = xyzw(4,n+istep)
c
!            write(0,'(A,i4,4F16.9)') 
!     .           ' RE ',n+istep,
!     ,       qxyz(1,n),qxyz(2,n),qxyz(3,n),qwght(n)
   30 continue
      istep=istep+nq
c      write(6,*)' repacked buffer '
      return
      end
