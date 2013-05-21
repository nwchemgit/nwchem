      subroutine hnd_efgmap(rtdb,basis,geom)
c
c $Id$
c
c     This routine calculates the electric field gradient and 
c     the orientation of the EFG for a given density at the
c     atomic positions.
c
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cosmo.fh"
c
      integer rtdb      ! [Input] rtdb  
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
c
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, i, j, ij
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt, l_efgs, k_efgs
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      integer nefc, l_efcc, k_efcc, l_efcz, k_efcz
      character*3 scftyp
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision vec(3,3), eig(3), a(6)
      double precision pi, deg, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz
      double precision rr, rr5, eta
      logical status
c     bq variables (MV)
      logical dobq
      integer bq_ncent
      integer i_cbq
      integer i_qbq
      double precision elpotbq
c
      logical docosmo
      integer ncosbq
c
c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Get density matrix
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     ----- calculate electric field gradient -----
c
      if (ga_nodeid().eq.0) write(luout,9999)
      if (ga_nodeid().eq.0) write(luout,9994)
c
      pi  = acos(-1.0d0)
      deg = 180.0d0/pi
c
      call ecce_print_module_entry('efg')
c
c
c     ----- define points for calculation -----
c           1. nuclei
c
      status=geom_ncent(geom,nat)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_efgmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,6*nat,'efg pnt',l_efgs,k_efgs))
     &    call errquit('hnd_efgmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_efgmap: ma failed',911,MA_ERR)
c
      do 30 iat=1,nat
         status=geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &                        dbl_mb(k_zanpt+iat-1))
   30 continue
c
      call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(k_xyzpt),nat,
     &                dbl_mb(k_efgs),2)
c
c     get bq structures if any (MV)
c     -----------------------------
      dobq = .false.
      if(geom_extbq_on()) then
        dobq = .true.
        bq_ncent = geom_extbq_ncenter()
        i_cbq = geom_extbq_coord()
        i_qbq = geom_extbq_charge()
      end if
c
      docosmo = .false.
      ncosbq = 0
      if (rtdb_get(rtdb,'cosmo:nefc',mt_int,1,ncosbq).and.(ncosbq.gt.0))
     &   docosmo = .true.
c
c     ----- collect and output results of all points -----
c 
cc      if (ga_nodeid().gt.0) goto 300
c
      if (docosmo) then
         if (.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc))
     &         call errquit('hnd_efgmap: rtdb get failed for nefc ',911,
     &         RTDB_ERR)
         if (.not.ma_push_get(mt_dbl,nefc*3,'efcc',l_efcc,k_efcc))
     &         call errquit('hnd_efgmap: malloc k_efcc fail',911,ma_err)
         if (.not.ma_push_get(mt_dbl,nefc,'efcz',l_efcz,k_efcz))
     &         call errquit('hnd_efgmap: malloc k_efcz fail',911,ma_err)
         if (.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,
     &         dbl_mb(k_efcc))) call 
     &         errquit('hnd_efgmap: rtdb get failed efcc',912,rtdb_err)
         if (.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,nefc,
     &         dbl_mb(k_efcz))) call 
     &         errquit('hnd_efgmap: rtdb get failed efcz',913,rtdb_err)
      end if ! docosmo
c
      do 230  iat=1,nat
         xp = dbl_mb(k_xyzpt  +3*(iat-1))
         yp = dbl_mb(k_xyzpt+1+3*(iat-1))
         zp = dbl_mb(k_xyzpt+2+3*(iat-1))
c
c     ----- add nuclear contribution -----
c
         efgxx = dbl_mb(k_efgs  +6*(iat-1))/3.0d0
         efgyy = dbl_mb(k_efgs+1+6*(iat-1))/3.0d0
         efgzz = dbl_mb(k_efgs+2+6*(iat-1))/3.0d0
         efgxy = dbl_mb(k_efgs+3+6*(iat-1))/3.0d0
         efgxz = dbl_mb(k_efgs+4+6*(iat-1))/3.0d0
         efgyz = dbl_mb(k_efgs+5+6*(iat-1))/3.0d0
         do 210 i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr = sqrt(xn*xn + yn*yn + zn*zn)
            if (rr.lt.1.0d-3) go to 210
            rr5=rr*rr*rr*rr*rr
            efgxx = efgxx - zan*xn*xn/rr5
            efgyy = efgyy - zan*yn*yn/rr5
            efgzz = efgzz - zan*zn*zn/rr5
            efgxy = efgxy - zan*xn*yn/rr5
            efgxz = efgxz - zan*xn*zn/rr5
            efgyz = efgyz - zan*yn*zn/rr5
  210    continue
c
c     ----- form -efc- contribution -----
c           from cosmo point charges !!!!
c
         if (docosmo) then
            do i = 1,nefc
              xn = dbl_mb(k_efcc+3*(i-1)  ) - xp
              yn = dbl_mb(k_efcc+3*(i-1)+1) - yp
              zn = dbl_mb(k_efcc+3*(i-1)+2) - zp
              rr =  sqrt(xn*xn + yn*yn + zn*zn)
              if (rr.lt.1.0d-3) then
                if (ga_nodeid().eq.0) write(luout,9993) xp,yp,zp,i
              else
               rr5=rr*rr*rr*rr*rr
               efgxx = efgxx - dbl_mb(k_efcz+i-1)*xn*xn/rr5
               efgyy = efgyy - dbl_mb(k_efcz+i-1)*yn*yn/rr5
               efgzz = efgzz - dbl_mb(k_efcz+i-1)*zn*zn/rr5
               efgxy = efgxy - dbl_mb(k_efcz+i-1)*xn*yn/rr5
               efgxz = efgxz - dbl_mb(k_efcz+i-1)*xn*zn/rr5
               efgyz = efgyz - dbl_mb(k_efcz+i-1)*yn*zn/rr5
              endif
            enddo
         end if ! docosmo
c
c        adding external bq contributions(MV)
c        ----------------------------------
         if (dobq) then
            do i = 1,bq_ncent
               xn = dbl_mb(i_cbq+3*(i-1)  ) - xp
               yn = dbl_mb(i_cbq+3*(i-1)+1) - yp
               zn = dbl_mb(i_cbq+3*(i-1)+2) - zp
               rr =  sqrt(xn*xn + yn*yn + zn*zn)
               if (rr.lt.1.0d-3) then
                  write(luout,9993) xp,yp,zp,i
               else
               rr5=rr*rr*rr*rr*rr
               efgxx = efgxx - dbl_mb(i_qbq+i-1)*xn*xn/rr5
               efgyy = efgyy - dbl_mb(i_qbq+i-1)*yn*yn/rr5
               efgzz = efgzz - dbl_mb(i_qbq+i-1)*zn*zn/rr5
               efgxy = efgxy - dbl_mb(i_qbq+i-1)*xn*yn/rr5
               efgxz = efgxz - dbl_mb(i_qbq+i-1)*xn*zn/rr5
               efgyz = efgyz - dbl_mb(i_qbq+i-1)*yn*zn/rr5
               endif
            end do
         end if
c
         dbl_mb(k_efgs  +6*(iat-1)) = 2.0d0*efgxx - efgyy - efgzz
         dbl_mb(k_efgs+1+6*(iat-1)) = 2.0d0*efgyy - efgxx - efgzz
         dbl_mb(k_efgs+2+6*(iat-1)) = 2.0d0*efgzz - efgxx - efgyy
         dbl_mb(k_efgs+3+6*(iat-1)) = 3.0d0*efgxy
         dbl_mb(k_efgs+4+6*(iat-1)) = 3.0d0*efgxz
         dbl_mb(k_efgs+5+6*(iat-1)) = 3.0d0*efgyz
c
c        ----- reorder into a as xx xy yy xz yz zz to form matrix -----
c
         a(1) = dbl_mb(k_efgs  +6*(iat-1))        
         a(2) = dbl_mb(k_efgs+3+6*(iat-1))        
         a(3) = dbl_mb(k_efgs+1+6*(iat-1))        
         a(4) = dbl_mb(k_efgs+4+6*(iat-1))        
         a(5) = dbl_mb(k_efgs+5+6*(iat-1))        
         a(6) = dbl_mb(k_efgs+2+6*(iat-1))        
         ij=0
         do 241 i = 1, 3
         do 241 j = 1, i
            ij = ij + 1
            vec(i,j) = a(ij) 
            vec(j,i) = a(ij)
  241    continue
c
c        ----- store ecce data -----
c
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_efgmap: geom_cent_tag failed',0,GEOM_ERR)
c        geom_tag_to_element returns false for Bq elements (MV)
c        -----------------------------------------------------
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) then 
            if(symbol.ne."bq") call
     &      errquit('hnd_efgmap: geom_tag_to_element failed',0,GEOM_ERR)
         end if
c
c         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) call
c     &      errquit('hnd_efgmap: geom_tag_to_element failed',0,GEOM_ERR)
         call ecce_print1_char('atom name',symbol,1)
         call ecce_print2('EFG tensor',MT_DBL,vec,3,3,3)
c
c        ----- print tensor components -----
c
         if (ga_nodeid().eq.0) then
           write(luout,9998) iat,symbol,xp,yp,zp
           write(luout,9997) 
           write(luout,9995) (dbl_mb(k_efgs+6*(iat-1)+i),i=0,5)
         end if
c
c        ----- diagonalize to get principal components and vectors -----
c
         call hnd_diag(vec,eig,3,.true.,.false.)
         eta  = abs( (eig(3)-eig(2)) / eig(1) )
c
         call ecce_print1('EFG eigenvalues',MT_DBL,eig,3)
         call ecce_print2('EFG eigenvectors',MT_DBL,vec,3,3,3)
         call ecce_print1('EFG asymmetry',MT_DBL,eta,1)
c
         if (ga_nodeid().eq.0) then
           write(luout,9992)
           write(luout,9991) eig(1),eig(2),eig(3),eta
           write(luout,9988) ((vec(i,j),j=1,3),i=1,3)
           write(luout,*) ' '
         end if
c
  230 continue ! Assemblin and printing next atom
c
      if (docosmo) then
        if (.not.ma_chop_stack(l_efcc)) 
     &    call errquit('hnd_efgmap: chop stack l_efcc',913,ma_err)
      end if ! docosmo
c
      call ecce_print_module_exit('EFG','ok')
      call util_flush(luout)
c
c     ----- release memory block -----
c
  300 call ga_sync()
c
c     ------- Deallocate MA memory ------
c
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_efgmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_efgs)) call errquit
     &   ('hnd_efgmap, ma_pop_stack of l_efgs failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_efgmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('efgmap: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
 9999 format(/,10x,23(1h-),/,10x,'Electric field gradient',
     1       /,10x,23(1h-),/)
 9998 format(/,1x,60(1h-),/,3x,'Atom',6x,'X',9x,'Y',9x,'Z',/,1x,60(1h-),
     1       /,i5,1x,a2,3f10.5,/,1x,60(1h-),/)
 9997 format(1x,'Electric field gradient in molecular frame (a.u.)',/,
     2 9x,'XX',13x,'YY',13x,'ZZ',13x,'XY',13x,'XZ',13x,'YZ',/,
     3 1x,90(1h-))
 9996 format(' --- Warning - electric field gradient at ',
     1 3F10.5,' . contribution from nucleus ',i3,' ignored')
 9995 format(1x,6f15.6,/)
 9994 format(' 1 a.u. = 0.324123 10**(16) esu/cm**3 ',
     1       ' ( or statvolts/cm**2 )',' = 0.97174 10**(22) v/m**2 ',/)
 9993 format(' --- Warning - electric field gradient at ',
     1 3f10.5,' . contribution from  -efc-  ',i3,' ignored')
 9992 format(1x,'Principal components (a.u.) and orientation ',
     1       /,' of principal axis w.r.t. absolute frame',
     2       22x,'Asymmetry parameter eta',/,1x,86(1h-))
 9991 format(1x,3f15.6,16x,f15.6,/)
 9988 format(1X,3F15.6)
      end
