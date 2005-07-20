c
c $Id: qmmm_cons.F,v 1.5 2005-07-20 20:49:16 marat Exp $
c
c
      subroutine qmmm_cons_init(irtdb)
      implicit none
c
#include "mm_utils.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
      integer irtdb
c
      character*32 pname
      logical fixw
      logical fixs
      logical freeqmlink
      logical freeqm
      logical fixqm
      logical fixqmlink
c
      if (.not. rtdb_get(irtdb,'qmmm:fixsolvent',mt_log,1,fixw))
     $     fixw=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixsolute',mt_log,1,fixs))
     $     fixs=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqmlink',mt_log,1,freeqmlink))
     $     freeqmlink=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqm',mt_log,1,freeqm))
     $     freeqm=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixqm',mt_log,1,fixqm))
     $     fixqm=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixqmlink',mt_log,1,fixqmlink))
     $     fixqmlink=.false.

       
     
      if(fixw) call mm_fix_solvent()
      if(fixs) call mm_fix_solute(mm_classical)
      if(fixqm) then
         write(*,*) "fixing quantum atoms"
         call mm_fix_solute(mm_quant)
      end if
      if(fixqmlink) then
         write(*,*) "fixing quantum+link atoms"
         call mm_fix_solute(mm_quant_link)
      end if

      if(freeqmlink) call mm_free_solute(mm_quant_link)
      if(freeqm) call mm_free_solute(mm_quant)

      return
      end
c
      subroutine qmmm_cons_reaction(irtdb,do_grad,etot,g)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
      integer irtdb
      logical do_grad
      double precision g(3,*)
      double precision etot
c
      character*32 pname
      integer igeom 
      integer ncent
      integer i_c,h_c
      integer nfg
      integer i_g,i_ifg
      integer ns
      integer i_is,h_is
      integer i_as,h_as
      integer i_rs,h_rs
      double precision k0,r0
c
      if (.not. rtdb_get(irtdb,'qmmm:n_spring',mt_int,1,ns))
     $     ns = 0

      if(ns.eq.0) return

     
      pname = "qmmm_cons_reaction"

      if(qmmm_print_default()) 
     +  write(*,*) "in",pname

      nfg = qmmm_geom_get_nfg()


      if(.not.ma_alloc_get(mt_int,2*ns,'qmmm rcons ind',h_is,i_is))
     + call errquit(pname//'Failed to allocate memory for is',
     + 2*ns, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,ns,'qmmm rcons r',h_rs,i_rs))
     + call errquit(pname//'Failed to allocate memory for rs',
     + ns, MA_ERR)
        call dfill(ns,0.0d0,dbl_mb(i_rs),1)

      if(.not.ma_alloc_get(mt_dbl,ns,'qmmm rcons coeff',h_as,i_as))
     + call errquit(pname//'Failed to allocate memory for as',
     + ns, MA_ERR)


      if (.not. rtdb_get(irtdb,'qmmm:i_spring',mt_int,
     +                    2*ns,int_mb(i_is)))
     + call errquit(pname//'Failed to allocate memory for is',
     + 3*ns, MA_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:a_spring',mt_dbl,
     +                    ns,dbl_mb(i_as)))
     +  call dfill(ns,1.0d0,dbl_mb(i_as),1)

      if (.not. rtdb_get(irtdb,'qmmm:k_spring',mt_dbl,
     +                    1,k0))
     + call errquit(pname//'Failed to get k0',
     + 1, RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:r_spring',mt_dbl,
     +                    1,r0))
     + call errquit(pname//'Failed to get r0',
     + 1, RTDB_ERR)


c    --------------------------------------
c    get coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom, ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)

      if(nfg.ne.ncent) 
     + call errquit(pname//': geometry mismatch',0, GEOM_ERR)

      if(.not.ma_alloc_get(mt_dbl,3*ncent,'ctmp',h_c,i_c))
     + call errquit( pname//'Failed to allocate memory for c',
     + 3*ncent, MA_ERR)

      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_c)))
     + call errquit(pname//': Failed to get geometry',0, GEOM_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit(pname//'Failed to destroy geometry',0, GEOM_ERR)

      call qmmm_cons_spring_egrad(ns,nfg,
     >                       do_grad,
     >                       k0,
     >                       r0,
     >                       etot,
     >                       int_mb(i_is),
     >                       dbl_mb(i_as),
     >                       dbl_mb(i_rs),
     >                       dbl_mb(i_c),
     >                       g)
 
      if(.not.ma_free_heap(h_c))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_as))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_rs))
     + call errquit(pname//'Failed to deallocate memory for rs',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_is))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      return
      end
c
      subroutine qmmm_cons_spring_egrad(ns,
     >                       nat,
     >                       do_grad,
     >                       k0,
     >                       r0,
     >                       etot,
     >                       ind,
     >                       a,
     >                       r,
     >                       c,
     >                       g)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      integer ns,nat
      logical do_grad
      double precision k0,r0
      integer ind(2,ns)
      double precision etot
      double precision  a(ns)
      double precision  r(ns)
      double precision g(3,nat)
      double precision c(3,nat)

      integer is,k,i1,i2
      double precision r1(3),r2(3)
      double precision rr,e,g0,ehalf
      character*32 pname

      pname = "qmmm_cons_spring_egrad"
      do is=1,ns
        i1 = ind(1,is)
        i2 = ind(2,is)
        do k=1,3
          r1(k) = c(k,i1)
          r2(k) = c(k,i2)
        end do
        rr = 0.0d0
        do k=1,3
          rr = rr+(r1(k)-r2(k))**2
        end do
        r(is) = sqrt(rr)
      end do
      
      ehalf = -r0
      do is=1,ns
        ehalf = ehalf+a(is)*r(is)
      end do
      e = ehalf*ehalf*k0
      etot = etot + e

      if(.not.do_grad) return

      do is=1,ns
        i1 = ind(1,is)
        i2 = ind(2,is)
        do k=1,3
          r1(k) = c(k,i1)
          r2(k) = c(k,i2)
        end do 
        do k=1,3
          g0=2*k0*ehalf*a(is)*(r1(k)-r2(k))/r(is)
          g(k,i1) = g(k,i1)+g0
          g(k,i2) = g(k,i2)-g0
        end do 
      end do

      return 
      end



