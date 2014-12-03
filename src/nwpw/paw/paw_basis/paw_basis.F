*
* $Id$
*

!**************************************************
!
!       Name: init_paw_basis
!
!       Purpose 
!
!       Created:        7/11/2002
!**************************************************
      subroutine init_paw_basis()
      implicit none
#include "errquit.fh"

      
#include "paw_basis.fh"
#include "bafdecls.fh"
  
*     !*** local variables ***
      integer it, ir, nt
      integer l_nb,l_ng
      integer i_nb,i_ng
      integer i_lo,i_hi,ii,mult_l
      integer  max_mult_l
      integer max_ic
      logical ok
      double precision log_amesh,amesh
      character*255 fname

c      !*** external functions ***
      integer ion_nkatm
      external ion_nkatm
      
*     *** check for existence of _basis files ***
      call paw_basis_exist_check()

      nt = ion_nkatm()

*     **** allocating memory from stack ****
      ok = BA_push_get(MT_INT,nt, "nb_tmp_array", l_nb,i_nb)
      ok = ok.and.
     >     BA_push_get(MT_INT,nt, "ng_tmp_array", l_ng,i_ng)
      if (.not.ok) 
     > call errquit('init_paw_basis: out of stack memory',0, MA_ERR)
  
      do it = 1, nt
       call ion_atom_plus_suffix(it,'_basis',fname)
       call util_file_name_resolve(fname, .false.)
       call read_paw_basis_dims(int_mb(i_nb-1 + it),
     >                          int_mb(i_ng-1 + it),
     >                          fname)
      end do

      call alloc_paw_basis_data(nt,int_mb(i_nb),int_mb(i_ng))


      do it = 1,nt 
         call ion_atom_plus_suffix(it,'_basis',fname)
         call util_file_name_resolve(fname, .false.)
         call read_paw_basis(  dbl_mb(paw_basis_i_r1(it)),
     >                         dbl_mb(paw_basis_i_rmax(it)),
     >                         int_mb(paw_basis_i_ngrid(it)),
     >                         int_mb(paw_basis_i_nbasis(it)),
     >                         int_mb(paw_basis_i_i_cut(it)),
     >                         int_mb(paw_basis_i_pr_n(it)),
     >                         int_mb(paw_basis_i_pr_n_ps(it)),
     >                         int_mb(paw_basis_i_orb_l(it)),
     >                         dbl_mb(paw_basis_i_sigma(it)),
     >                         dbl_mb(paw_basis_i_ion_charge(it)),
     >                         dbl_mb(paw_basis_i_core_kin_energy(it)),
     >                         dbl_mb(paw_basis_i_eig(it)),
     >                         dbl_mb(paw_basis_i_phi_ae(it)),
     >                         dbl_mb(paw_basis_i_phi_ae_prime(it)),
     >                         dbl_mb(paw_basis_i_phi_ps(it)),
     >                         dbl_mb(paw_basis_i_phi_ps_prime(it)),
     >                         dbl_mb(paw_basis_i_prj_ps(it)),
     >                         dbl_mb(paw_basis_i_prj_ps0(it)),
     >                         dbl_mb(paw_basis_i_core_density(it)),
     >                         dbl_mb(paw_basis_i_ps_core_density(it)),
     >                         dbl_mb(paw_basis_i_v_ps(it)),
     >                         fname)

c        **** compute radial derivative of core densities ****
  
         log_amesh = LOG(dbl_mb(paw_basis_i_rmax(it))
     >                  /dbl_mb(paw_basis_i_r1(it)))
     >                 /(int_mb(paw_basis_i_ngrid(it))-1)
         amesh     = EXP(log_amesh)
    
         dbl_mb(paw_basis_i_log_amesh(it)) = log_amesh
         dbl_mb(paw_basis_i_rgrid(it)) = dbl_mb(paw_basis_i_r1(it))
  
         do ir=1,int_mb(paw_basis_i_ngrid(it))-1
          dbl_mb(paw_basis_i_rgrid(it)+ir)  = 
     >           amesh*dbl_mb(paw_basis_i_rgrid(it)+ir-1)
         end do
         
c       !***** determine max orbital moment for this kind of atom ****
         i_lo =  paw_basis_i_orb_l(it)
         i_hi =  paw_basis_i_orb_l(it)  
     >        + int_mb(paw_basis_i_nbasis(it)) - 1
         mult_l = int_mb(i_lo)
         do ii=i_lo+1,i_hi
           if (int_mb(ii).gt.mult_l) mult_l = int_mb(ii)
         end do
         int_mb(paw_basis_i_mult_l(it)) = 2*mult_l
           
c        **** compute radial derivative of core densities - needed for gga's****
         call paw_basis_derivative_ngrid(
     >            int_mb(paw_basis_i_ngrid(it)),
     >            log_amesh,
     >            dbl_mb(paw_basis_i_rgrid(it)),  
     >            dbl_mb(paw_basis_i_core_density(it)),  
     >            dbl_mb(paw_basis_i_core_density_prime(it)))
         call paw_basis_derivative_ngrid(
     >            int_mb(paw_basis_i_ngrid(it)),
     >            log_amesh,
     >            dbl_mb(paw_basis_i_rgrid(it)),  
     >            dbl_mb(paw_basis_i_ps_core_density(it)),  
     >            dbl_mb(paw_basis_i_ps_core_density_prime(it)))
      end do


c     !***** determine largest max orbital moment ****
      max_mult_l = 0
      do it=1,nt
        if (int_mb(paw_basis_i_mult_l(it)).gt.max_mult_l) 
     >    max_mult_l = int_mb(paw_basis_i_mult_l(it))
      end do
      call set_max_mult_l(max_mult_l)

c     !***** determine largest sphere cutoff ****
      max_ic = 0
      do it=1,nt
        if (int_mb(paw_basis_i_i_cut(it)).gt.max_ic) 
     >    max_ic = int_mb(paw_basis_i_i_cut(it))
      end do
      call set_max_i_cut(max_ic)


c     !*** set core charges ***
      do it = 1,nt 
        call calc_paw_basis_core_charge(
     (        it,dbl_mb(paw_basis_i_core_charge(it)))
      end do


c     !*** set paw_radgrid_map ***
      do it = 1,nt 
        call paw_radgrid_map(int_mb(paw_basis_i_i_cut(it)),
     >                       int_mb(paw_basis_i_i_start(it)),
     >                       int_mb(paw_basis_i_i_end(it)))
        
        int_mb(paw_basis_i_npoints(it)) 
     >  = int_mb(paw_basis_i_i_end(it))
     >  - int_mb(paw_basis_i_i_start(it)) + 1

c        write(*,*) "nt:",it,int_mb(paw_basis_i_ngrid(it)),
c     >             int_mb(paw_basis_i_nbasis(it)),
c     >             int_mb(paw_basis_i_i_cut(it)),
c     >             paw_basis_i_start(it),
c     >             paw_basis_i_end(it),
c     >             paw_basis_npoints(it)
      end do


*     **** deallocating memory from stack ****
      ok =        BA_pop_stack(l_ng)
      ok = ok.and.BA_pop_stack(l_nb)
      if (.not.ok) 
     > call errquit('init_paw_basis:error popping stack',1, MA_ERR)

      return
      end 

!**************************************************
!
!       Name: read_paw_basis_dims
!
!       Purpose 
!
!       Created:        7/11/2002
!**************************************************
      subroutine read_paw_basis_dims(nb,ng,filename)
      implicit none
      integer          nb
      integer          ng
      character*(*)     filename

#include "errquit.fh"
#include "stdio.fh"
#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"

*     !*** local variables ***
      integer     fp
      double precision dummy_r1,zv
      double precision dummy_rmax
      integer MASTER,taskid
      parameter (MASTER=0)
      integer msglen,psp_type
      character*2 atom

*     !*** external functions ****
      logical  assert_file
      external assert_file

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then

      if( .not. assert_file(filename)) then
         write(luout,*)  ' cannot find ',filename
        call errquit("paswbasis: cannot find file",1, DISK_ERR)
      end if

      call open_file(filename,"r",fp)
        read(fp,*) psp_type
        read(fp,'(A)') atom
        read(fp,*) zv

        READ(fp,*) dummy_r1,dummy_rmax,ng,nb
      close(fp)
      endif
      msglen = 1
      call BRDCST(9+MSGINT,ng,             mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,nb,           mitob(msglen),MASTER)
      return
      end


!**************************************************
!
!       Name: read_paw_basis
!
!       Purpose 
!
!       Created:        7/11/2002
!************************************************** 
      subroutine read_paw_basis(r1,
     >                          rmax,
     >                          ngrid,
     >                          nbasis,
     >                          i_cut,
     >                          prin_n,
     >                          prin_n_ps,
     >                          orb_l,
     >                          sigma,
     >                          ion_charge,
     >                          core_kin_energy,
     >                          eig,
     >                          phi_ae,
     >                          phi_ae_prime,
     >                          phi_ps,
     >                          phi_ps_prime,
     >                          prj_ps,
     >                          prj_ps0,
     >                          core_density,
     >                          ps_core_density,
     >                          v_ps,
     >                          aname)
      implicit none
#include "errquit.fh"
      double precision r1
      double precision rmax
      integer          ngrid
      integer          nbasis
      integer          i_cut
      integer          prin_n(nbasis)
      integer          prin_n_ps(nbasis)
      integer          orb_l(nbasis)
      double precision sigma
      double precision ion_charge
      double precision core_kin_energy
      double precision eig(nbasis)
      double precision phi_ae(nbasis*ngrid)
      double precision phi_ae_prime(nbasis*ngrid)
      double precision phi_ps(nbasis*ngrid)
      double precision phi_ps_prime(nbasis*ngrid)
      double precision prj_ps(nbasis*ngrid)
      double precision prj_ps0(nbasis*ngrid)
      double precision core_density(ngrid)
      double precision ps_core_density(ngrid)
      double precision v_ps(ngrid)
      character*(*)    aname

#include "stdio.fh"
#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"


*     !*** local variables
      integer MASTER,taskid
      parameter (MASTER=0)
      integer i,j,msglen,idum
      integer fp
      !character*255 filename

      integer psp_type
      real*8  zv,rc(25)
      character*2  atom
      character*80 comment

*     !*** external functions ***
      logical  assert_file
      external assert_file

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then

      if( .not. assert_file(aname)) then
        write(luout,*)  ' cannot find ',aname
        call errquit('paswbasis: cannot find file',2, DISK_ERR)
      end if

      call open_file(aname,"r",fp)
        read(fp,*) psp_type
        read(fp,'(A)') atom
        read(fp,*) zv

        READ(fp,*) r1,rmax,ngrid,nbasis

        read(fp,*) (rc(i),i=1,nbasis)

        READ(fp,*) i_cut

        read(fp,'(A)') comment

        READ(fp,*) core_kin_energy
        READ(fp,*) (prin_n(i),eig(i),
     >              prin_n_ps(i),orb_l(i),i=1,nbasis)
        READ(fp,*) ((phi_ae((i-1)*ngrid+j), j=1,ngrid),
     >                                      i=1,nbasis)
        READ(fp,*) ((phi_ae_prime((i-1)*ngrid+j), j=1,ngrid),
     >                                            i=1,nbasis)
        READ(fp,*) ((phi_ps((i-1)*ngrid+j), j=1,ngrid),
     >                                      i=1,nbasis)
        READ(fp,*) ((phi_ps_prime((i-1)*ngrid+j), j=1,ngrid),
     >                                            i=1,nbasis)
        READ(fp,*) ((prj_ps((i-1)*ngrid+j), j=1,ngrid), 
     >                                      i=1,nbasis)
        READ(fp,*) (core_density(i),   i=1,ngrid)
        READ(fp,*) (ps_core_density(i),i=1,ngrid)
        READ(fp,*) (v_ps(i),           i=1,ngrid)
        READ(fp,*) sigma
        READ(fp,*) ion_charge
        READ(fp,*) ((prj_ps0((i-1)*ngrid+j), j=1,ngrid), 
     >                                       i=1,nbasis)
      close(fp)
      end if

*     **** send header data to all processors ****
      msglen = 1
      call BRDCST(9+MSGDBL,r1,             mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,rmax,           mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,core_kin_energy,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,sigma,          mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,ion_charge,     mdtob(msglen),MASTER)

      call BRDCST(9+MSGINT,nbasis,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,ngrid, mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,i_cut, mitob(msglen),MASTER)


      msglen = nbasis
      call BRDCST(9+MSGINT,prin_n,   mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,prin_n_ps,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,orb_l,    mitob(msglen),MASTER)

      call BRDCST(9+MSGDBL,eig,      mdtob(msglen),MASTER)


      msglen = ngrid
      call BRDCST(9+MSGDBL,core_density,   mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,ps_core_density,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,v_ps,           mdtob(msglen),MASTER)


      msglen = ngrid*nbasis
      call BRDCST(9+MSGDBL,phi_ae,      mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,phi_ps,      mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,phi_ae_prime,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,phi_ps_prime,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,prj_ps,      mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,prj_ps0,     mdtob(msglen),MASTER)


      return
      end 




!**************************************************
!
!       Name: paw_basis_exist_check
!
!       Purpose:
!
!       Created:        2/12/2002
!**************************************************
      subroutine paw_basis_exist_check()
      implicit none

#include "btdb.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical found
      integer ia,lmax,locp
      character*4 atom
      character*255 fname,libfile
      real*8        sum,rlocal

*     **** external functions ****
      logical     nwpw_filefind
      character*4 ion_atom
      integer     ion_nkatm
      external    nwpw_filefind
      external    ion_atom
      external    ion_nkatm

      do ia = 1, ion_nkatm()
        call ion_atom_plus_suffix(ia,'_basis',fname)

*       **** generate one-dimensional paw basis file ****
        if (.not.nwpw_filefind(fname)) then
          atom = ion_atom(ia)
          call psp_generator_auto(2,atom,lmax,locp,rlocal)
        end if

      end do

      call ga_sync()

      !*** sync - hack ***
      !sum = 1.0d0
      !call D3dB_SumAll(sum)

      return
      end


