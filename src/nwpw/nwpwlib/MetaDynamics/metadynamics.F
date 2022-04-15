c $Id$

c     **********************************************
c     *                                            *
c     *            meta_initialize                 *
c     *                                            *
c     **********************************************
      subroutine meta_initialize(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,cnmeta_read
      integer d,i,j,n1,n2,s1,s2,yshift,ntype
      real*8 x,y,z,r,rmax,dr,pi,da,maxcoord,a,b,n,m,sigma,w,r0,s0,v0,dv0
      real*8 n12
      real*8 boundary(2),temp,fac,sn2(12)
      character*80 rtdb_name,potential_filename
      character*60 gauss_filename
      character*4 celement1,celement2
      character*500 eqnstring

*     **** external functions ****
      logical     control_Nose
      real*8      lattice_unita,ion_Temperature,control_Nose_Tr
      character*7 c_index_name
      integer     control_it_out
      external    control_Nose
      external    lattice_unita,ion_Temperature,control_Nose_Tr
      external    c_index_name
      external    control_it_out

      call Parallel_taskid(taskid)
      metaprintcount = 0
      rtdb_name = 'metadynamics_print'
      if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,maxmetaprintcount)) 
     >   maxmetaprintcount=10

      rtdb_name = 'metadynamics_nmeta'
      if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nmeta)) nmeta = 0

      metafound = (nmeta.gt.0)

      if (metafound) then

         rtdb_name = 'metadynamics_update'
         if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,maxmetacount)) 
     >      maxmetacount = 1

         rtdb_name = 'metadynamics_metacount'
         if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,metacount)) 
     >      metacount = 0

         rtdb_name = 'metadynamics_print_shift'
         if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,metarayshift)) 
     >      metarayshift = maxmetaprintcount

         rtdb_name = 'metadynamics_metaraycount'
         if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,metaraycount)) 
     >      metaraycount = 0

         rtdb_name = 'metadynamics_tempered'
         if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,1,dTtempered)) 
     >      dTtempered = -1.0d0

         rtdb_name = 'metadynamics_boundary'
         if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,2,boundary)) then
            boundary(1) = 0.0d0
            boundary(2) = 0.0d0
         end if

         rtdb_name = 'metadynamics_sn2-surface'
         if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,12,sn2)) then
            call dcopy(12,0.0d0,0,sn2,1)
         end if

         value = .true.
         if (nmeta.gt.0) then
            rtdb_name = 'metadynamics_nindxmeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nindxmeta)) 
     >         nindxmeta = 0
            rtdb_name = 'metadynamics_nparammeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nparammeta)) 
     >         nparammeta = 0
            rtdb_name = 'metadynamics_nxmeta_all'
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nxmeta_all)) 
     >         nxmeta_all = 0

*           **** allocate and read nxmeta,sindxmeta,sparameta,and indxmeta ****
            value = value.and.
     >              BA_alloc_get(mt_int,nmeta,'nxmeta',
     >                           nxmeta(2),nxmeta(1))
            value = value.and.
     >              BA_alloc_get(mt_int,nmeta,'sindxmeta',
     >                           sindxmeta(2),sindxmeta(1))
            value = value.and.
     >              BA_alloc_get(mt_int,nmeta,'sparammeta',
     >                           sparammeta(2),sparammeta(1))
            value = value.and.
     >              BA_alloc_get(mt_int,nindxmeta,'indxmeta',
     >                           indxmeta(2),indxmeta(1))

            rtdb_name = 'metadynamics_nxmeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,nmeta,
     >                        int_mb(nxmeta(1)))) then
               do d=1,nmeta
                  int_mb(nxmeta(1)+d-1) = 0
               end do
            end if
            rtdb_name = 'metadynamics_sindxmeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,nmeta,
     >                        int_mb(sindxmeta(1)))) then
               do d=1,nmeta
                  int_mb(sindxmeta(1)+d-1) = 0
               end do
            end if
            rtdb_name = 'metadynamics_sparammeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,nmeta,
     >                        int_mb(sparammeta(1)))) then
               do d=1,nmeta
                  int_mb(sindxmeta(1)+d-1) = 0
               end do
            end if
            rtdb_name = 'metadynamics_indxmeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,nindxmeta,
     >                        int_mb(indxmeta(1)))) then
               do d=1,nindxmeta
                  int_mb(indxmeta(1)+d-1) = 0
               end do
            end if



*           ***** allocate and read pmeta,ameta,bmeta,parammeta, and ymeta ****
            value = value.and.
     >              BA_alloc_get(mt_int,nmeta,'pmeta',pmeta(2),pmeta(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nmeta,'ameta',ameta(2),ameta(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nmeta,'bmeta',bmeta(2),bmeta(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nparammeta,'parammeta',
     >                           parammeta(2),parammeta(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nxmeta_all,
     >                           'ymeta',ymeta(2),ymeta(1))

            rtdb_name = 'metadynamics_pmeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,nmeta,
     >                        int_mb(pmeta(1)))) then
               do d=1,nmeta
                  int_mb(pmeta(1)+d-1) = 0
               end do
            end if
            rtdb_name = 'metadynamics_ameta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,nmeta,
     >                        dbl_mb(ameta(1)))) then
               do d=1,nmeta
                  dbl_mb(ameta(1)+d-1) = 0.0d0
               end do
            end if
            rtdb_name = 'metadynamics_bmeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,nmeta,
     >                        dbl_mb(bmeta(1)))) then
               do d=1,nmeta
                  dbl_mb(bmeta(1)+d-1) = 0.0d0
               end do
            end if
            rtdb_name = 'metadynamics_parammeta'
            if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,nparammeta,
     >                        dbl_mb(parammeta(1)))) then
               do d=1,nparammeta
                  dbl_mb(parammeta(1)+d-1) = 0.0d0
               end do
            end if

*           **** set default limits ****
            x = lattice_unita(1,1)
            y = lattice_unita(2,1)
            z = lattice_unita(3,1)
            r = dsqrt(x*x+y*y+z*z)
            rmax = r
            x = lattice_unita(1,2)
            y = lattice_unita(2,2)
            z = lattice_unita(3,2)
            r = dsqrt(x*x+y*y+z*z)
            if (r.gt.rmax) rmax = r
            x = lattice_unita(1,3)
            y = lattice_unita(2,3)
            z = lattice_unita(3,3)
            r = dsqrt(x*x+y*y+z*z)
            if (r.gt.rmax) rmax = r
            do d=1,nmeta
               s1 = int_mb(sindxmeta(1) +d-1)
               if ((int_mb(indxmeta(1)+s1).eq.1).or.
     >             (int_mb(indxmeta(1)+s1).eq.5) .or.
     >             (int_mb(indxmeta(1)+s1).eq.13)) then
                  if (int_mb(pmeta(1)+d-1).lt.0)
     >               int_mb(pmeta(1)+d-1) = 0
                  if (dbl_mb(ameta(1)+d-1).lt.0.0d0) 
     >               dbl_mb(ameta(1)+d-1) = 0.0d0
                  if (dbl_mb(bmeta(1)+d-1).lt.0.0d0)
     >               dbl_mb(bmeta(1)+d-1) = rmax
               end if
               if (int_mb(indxmeta(1)+s1).eq.2) then
                  if (int_mb(pmeta(1)+d-1).lt.0)
     >               int_mb(pmeta(1)+d-1) = 0
                  if (dbl_mb(ameta(1)+d-1).lt.(-99.0d0))
     >               dbl_mb(ameta(1)+d-1) = 0.0d0
                  if (dbl_mb(bmeta(1)+d-1).lt.(-99.0d0))
     >               dbl_mb(bmeta(1)+d-1) = 4.0d0*datan(1.0d0)
               end if
               if (int_mb(indxmeta(1)+s1).eq.3) then
                  pi = 4.0d0*datan(1.0d0)
                  if (dbl_mb(ameta(1)+d-1).lt.(-pi))
     >               dbl_mb(ameta(1)+d-1) = -pi
                  if (dbl_mb(bmeta(1)+d-1).gt.(pi))
     >               dbl_mb(bmeta(1)+d-1) = pi
               end if
            end do
         end if

*        **** load old spline fits ****
         call dcopy(nxmeta_all,0.0d0,0,dbl_mb(ymeta(1)),1)
         rtdb_name = 'metadynamics_ymeta'
         if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,nxmeta_all,
     >                       dbl_mb(ymeta(1)))) then
            call dcopy(nxmeta_all,0.0d0,0,dbl_mb(ymeta(1)),1)
            if (boundary(1).gt.0.0d0) then
               if (taskid.eq.MASTER)
     >            write(luout,'(A,2F11.6)') 
     >            " Metadynamics: setting boundary, w,sigma=",
     >                 boundary(1),boundary(2)
               call meta_setboundary(boundary(1),boundary(2))
            end if
            if (dabs(sn2(1)).gt.0.0d0) then
               if (taskid.eq.MASTER)
     >            write(luout,'(A)')
     >            " Metadynamics: setting sn2 surface"
               call meta_setsn2surface(sn2)
            end if
         else
            if (taskid.eq.MASTER) 
     >        write(luout,*) 
     >        "Metadynamics: reading metadynamics data from rtdb"

            !*** re-scale tempered ***
            fac = 1.0d0
            if (dTtempered.gt.0.0d0) then
               if (taskid.eq.MASTER) 
     >            write(luout,*) 
     >            "Metadynamics: re-scaling to be consistent",
     >            " with tempering"
               if (control_Nose()) then
                  temp = control_Nose_Tr()
               else
                  temp = ion_Temperature()
               end if
               fac = dTtempered/(temp+dTtempered)
            end if
            do i=1,nxmeta_all
               dbl_mb(ymeta(1)+i-1) = fac*dbl_mb(ymeta(1)+i-1)
            end do
         end if


         call nwpw_interp_init(nmeta,8)

         if (.not.value)
     >    call errquit('cannot allocate heap memory for metadynamics',0,
     >       MA_ERR)


*        **** iniitialize meta_gauss_write ****
         rtdb_name = 'metadynamics_gaussian_filename'
         gauss_filename = 
     >   '                                                           '
         if(.not.btdb_cget(rtdb,rtdb_name,1,gauss_filename))
     >      call util_file_prefix('meta_gaussians',gauss_filename)
         call meta_gauss_write_init(gauss_filename)

*        **** start nwpw_expression ****
         call nwpw_expression_start(rtdb)


*     **** write out header info ****
      call Parallel_taskid(taskid)
      value = btdb_parallel(.false.)
      if (taskid.eq.MASTER) then
         write(luout,*) 
         write(luout,*)  "MetaDynamics parameters:"
         write(luout,'(A,5x,I6," * inner iterations")')
     >   "   - update      = ",maxmetacount
         write(luout,'(A,5x,I6," * inner iterations")')
     >   "   - print       = ",maxmetaprintcount
         write(luout,'(A,5x,I6)') "   - print_shift = ",metarayshift
         if (dTtempered.gt.0.0d0) then
            write(luout,'(A,5x,F11.1)') 
     >      "   - Tempered dT (K) = ",dTtempered
         end if
         do d=1,nmeta
            s1 = int_mb(sindxmeta(1) +d-1)
            s2 = int_mb(sparammeta(1)+d-1)
            ntype = int_mb(indxmeta(1)+s1)
            w     = dbl_mb(parammeta(1)+s2)
            sigma = dbl_mb(parammeta(1)+s2+1)

            if (ntype.eq.1) then
             write(luout,'(A,5x,2I4,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - Bond Parameters =  ",
     >           (int_mb(indxmeta(1)+s1+j),j=1,2),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic=',int_mb(pmeta(1)+d-1) 

            else if (ntype.eq.2) then
             write(luout,'(A,3I4,4x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - Angle Parameters =  ",
     >           (int_mb(indxmeta(1)+s1+j),j=1,3),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic=',int_mb(pmeta(1)+d-1) 

            else if (ntype.eq.3) then
             write(luout,'(A,4I4,4x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - Dihedral Parameters =  ",
     >           (int_mb(indxmeta(1)+s1+j),j=1,4),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic=',int_mb(pmeta(1)+d-1) 

            else if (ntype.eq.4) then
             n1 = int_mb(indxmeta(1)+s1+1)
             n2 = int_mb(indxmeta(1)+s1+2)
             n  =  dbl_mb(parammeta(1)+s2+2)
             m  =  dbl_mb(parammeta(1)+s2+3)
             r0 =  dbl_mb(parammeta(1)+s2+4)
             n12 = dbl_mb(parammeta(1)+s2+6)
             write(luout,1001) (int_mb(indxmeta(1)+s1+3+j-1),j=1,n1)
             write(luout,1002) (int_mb(indxmeta(1)+s1+3+n1+j-1),j=1,n2)

             write(luout,'(A,/6x,A,F11.6,4x,A,F11.6,4x,A,F11.6,
     >                   /6x,A,F11.6,4x,A,F11.6,
     >                   /6x,A,F11.6,
     >                   /6x,A,F11.6,F11.6,
     >                   /6x,A,I5,
     >                   /6x,A,I2)')
     >       "  - Coordination Number Parameters: ",
     >           'n= ',n,'m=       ',m,'n12= ',n12,
     >           'w= ',w,'sigma=   ',sigma,
     >           'r0=',r0,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'ngrid=  ',int_mb(nxmeta(1)+d-1),
     >           'periodic=  ',int_mb(pmeta(1)+d-1)
              if (dbl_mb(parammeta(1)+s2+5).lt.0) then
                 write(luout,'(6x,A)') '- LJ function form'
              else
                 write(luout,'(6x,A)') '- Sprik function form'
              end if

            else if (ntype.eq.5) then
             n1 = int_mb(indxmeta(1)+s1+2)
             write(luout,1001) int_mb(indxmeta(1)+s1+1)
             write(luout,1002) (int_mb(indxmeta(1)+s1+3+j-1),j=1,n1)
             write(luout,'(A,5x,
     >                /6x,A,F11.6,F11.6,F11.6,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "  - n-plane Parameters =  ",
     >           'normal=',
     >           dbl_mb(parammeta(1)+s2+2),
     >           dbl_mb(parammeta(1)+s2+3),
     >           dbl_mb(parammeta(1)+s2+4),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic= ',int_mb(pmeta(1)+d-1)
            else if (ntype.eq.6) then

             if (int_mb(indxmeta(1)+s1+1).eq.1) then
             write(luout,'(A,5x,I4,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - X Parameters =  ",
     >           int_mb(indxmeta(1)+s1+2),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic= ',int_mb(pmeta(1)+d-1)
             else if (int_mb(indxmeta(1)+s1+1).eq.2) then
             write(luout,'(A,5x,I4,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - Y Parameters =  ",
     >           int_mb(indxmeta(1)+s1+2),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic= ',int_mb(pmeta(1)+d-1)
             else 
             write(luout,'(A,5x,I4,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - Z Parameters =  ",
     >           int_mb(indxmeta(1)+s1+2),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic= ',int_mb(pmeta(1)+d-1)
             end if

            else if (ntype.eq.7) then
             write(luout,'(A,5x,A,I4,4x,A,I4,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - local_density_trace Parameters =  ",
     >           'atom index=',int_mb(indxmeta(1)+s1+1),
     >           'l=',int_mb(indxmeta(1)+s1+2),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic=',int_mb(pmeta(1)+d-1) 

            else if (ntype.eq.8) then
             write(luout,'(A,
     >                /6x,A,2I4,
     >                /6x,A,2I4,
     >                /6x,A,2F11.6,
     >                /6x,A,2F11.6,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - 2Bonds Parameters =  ",
     >           'bond 1 indexes=',(int_mb(indxmeta(1)+s1+j),j=1,2),
     >           'bond 2 indexes=',(int_mb(indxmeta(1)+s1+j),j=3,4),
     >           'bond center 1 =',(dbl_mb(parammeta(1)+s2+j),j=2,3),
     >           'bond center 2 =',(dbl_mb(parammeta(1)+s2+j),j=4,5),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic=',int_mb(pmeta(1)+d-1)

            else if (ntype.eq.9) then
             write(luout,'(A,5x,A,I4,4x,A,I4,4x,A,I4,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - local_density_trace2_diff Parameters =  ",
     >           'atom index1=',int_mb(indxmeta(1)+s1+1),
     >           'atom index2=',int_mb(indxmeta(1)+s1+2),
     >           'l=',int_mb(indxmeta(1)+s1+3),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic=',int_mb(pmeta(1)+d-1) 

            else if (ntype.eq.10) then
             write(luout,'(A,
     >                /6x,A,3I4,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >       "   - Bond_Difference2 Parameters =  ",
     >           'bond indexes=',(int_mb(indxmeta(1)+s1+j),j=1,3),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic=',int_mb(pmeta(1)+d-1)

            else if (ntype.eq.11) then
             n1 = int_mb(indxmeta(1)+s1+1)
             n2 = int_mb(indxmeta(1)+s1+2)
             n  =  dbl_mb(parammeta(1)+s2+2)
             m  =  dbl_mb(parammeta(1)+s2+3)
             r0 =  dbl_mb(parammeta(1)+s2+4)
             v0  =  dbl_mb(parammeta(1)+s2+6)
             dv0 =  dbl_mb(parammeta(1)+s2+7)
             write(luout,1003) (int_mb(indxmeta(1)+s1+3+j-1),j=1,n1)
             write(luout,1004) (int_mb(indxmeta(1)+s1+3+n1+j-1),j=1,n2)

             write(luout,'(A,/6x,A,F11.6,4x,A,F11.6,
     >                   /6x,A,F11.6,4x,A,F11.6,
     >                   /6x,A,F11.6,4x,A,F11.6,F11.6,
     >                   /6x,A,F15.6,F15.6,
     >                   /6x,A,I5,
     >                   /6x,A,I2)')
     >       "  - Bond_Difference_Shell Parameters: ",
     >           'n= ',n,'m=       ',m,
     >           'w= ',w,'sigma=   ',sigma,
     >           'r0=',r0,'V0,dV0= ',v0,dv0,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'ngrid=  ',int_mb(nxmeta(1)+d-1),
     >           'periodic=  ',int_mb(pmeta(1)+d-1)
              if (dbl_mb(parammeta(1)+s2+5).lt.0) then
                 write(luout,'(6x,A)') '- LJ function form'
              else
                 write(luout,'(6x,A)') '- Sprik function form'
              end if

            else if (ntype.eq.12) then
               call nwpw_expression_eqnstring(rtdb,
     >                                        int_mb(indxmeta(1)+s1+1),
     >                                        eqnstring)
               write(luout,'(A,/6x,A,A,
     >                /6x,A,F11.6,4x,A,F11.6,
     >                /6x,A,F11.6,F11.6,4x,A,I5,4x,A,I2)')
     >           '  - Equation Parameters: ',
     >           'equation= ',trim(eqnstring),
     >           'w=',w,'sigma=',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'nrange= ',int_mb(nxmeta(1)+d-1),
     >           'periodic=',int_mb(pmeta(1)+d-1) 

            else if (ntype.eq.13) then
             n1 = int_mb(indxmeta(1)+s1+1)
             n2 = int_mb(indxmeta(1)+s1+2)
             n  =  dbl_mb(parammeta(1)+s2+2)
             write(luout,1005) (int_mb(indxmeta(1)+s1+3+j-1),j=1,n1)
             write(luout,1006) (int_mb(indxmeta(1)+s1+3+n1+j-1),j=1,n2)

             write(luout,'(A,/6x,A,F11.6,
     >                   /6x,A,F11.6,4x,A,F11.6,
     >                   /6x,A,F11.6,F11.6,
     >                   /6x,A,I5,
     >                   /6x,A,I2)')
     >       "  - Group_Distance Parameters: ",
     >           'n= ',n,
     >           'w= ',w,'sigma=   ',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'ngrid=  ',int_mb(nxmeta(1)+d-1),
     >           'periodic=  ',int_mb(pmeta(1)+d-1)

            else if (ntype.eq.14) then
             n1 = int_mb(indxmeta(1)+s1+1)
             write(luout,'(A,
     >                   /6x,A,F11.6,4x,A,F11.6,
     >                   /6x,A,F11.6,F11.6,
     >                   /6x,A,I5,
     >                   /6x,A,I2)')
     >       "  - Bondings Parameters: ",
     >           'w= ',w,'sigma=   ',sigma,
     >           'range=  ',dbl_mb(ameta(1)+d-1),dbl_mb(bmeta(1)+d-1),
     >           'ngrid=  ',int_mb(nxmeta(1)+d-1),
     >           'periodic=  ',int_mb(pmeta(1)+d-1)

             write(luout,'(A)')" coefficient index1 index2    :"
             do j=1,n1
               write(luout,'(F12.6,2I7)') dbl_mb(parammeta(1)+s2+1+j),
     >                                 int_mb(indxmeta(1)+s1+2*j),
     >                                 int_mb(indxmeta(1)+s1+2*j+1)
             end do

            end if

         end do

 1001 FORMAT(2x,"- Coorination Number (Index1) :",10I5)
 1002 FORMAT(2x,"- Coorination Number (Index2) :",10I5)
 1003 FORMAT(2x,"- Bond Difference Indexes (Index1) :",10I5)
 1004 FORMAT(2x,"- Bath Indexes            (Index2) :",10I5)
 1005 FORMAT(2x,"- Group_Distance Indexes (Index1) :",10I5)
 1006 FORMAT(2x,"- Group_Distance Indexes (Index2) :",10I5)

      end if
      value = btdb_parallel(.true.)
      
      end if
      return
      end

c     **********************************************
c     *                                            *
c     *            meta_finalize                   *
c     *                                            *
c     **********************************************
      subroutine meta_finalize(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer i,j,yshift,k
      real*8 temp,fac
      character*80 filename,rtdb_name
      character*255 full_filename

*     **** external functions ****
      logical  control_Nose
      external control_Nose
      character*7 c_index_name
      external    c_index_name
      real*8   control_Nose_Tr,ion_Temperature
      external control_Nose_Tr,ion_Temperature

      call Parallel_taskid(taskid)

      if (metafound) then

*        **** print out potentials ******
         !metaraycount = metaraycount + metarayshift
         call meta_print_potentials(0)

         value=.true.
         rtdb_name='metadynamics_metacount'
         value= value.and.btdb_put(rtdb,rtdb_name,mt_int,1,metacount)
         rtdb_name='metadynamics_metaraycount'
         value= value.and.btdb_put(rtdb,rtdb_name,mt_int,1,metaraycount)

         if (nmeta.gt.0) then

            !*** un-scale tempered ***
            fac = 1.0d0
            if (dTtempered.gt.0.0d0) then
               if (control_Nose()) then
                  temp = control_Nose_Tr()
               else
                  temp = ion_Temperature()
               end if
               fac = (temp+dTtempered)/dTtempered
            end if
            do i=1,nxmeta_all
               dbl_mb(ymeta(1)+i-1) = fac*dbl_mb(ymeta(1)+i-1)
            end do 

*           **** write out current spline fits ****
            rtdb_name = 'metadynamics_ymeta'
            value = value.and.btdb_put(rtdb,rtdb_name,mt_dbl,
     >                        nxmeta_all,dbl_mb(ymeta(1)))

            value = value.and.BA_free_heap(pmeta(2))
            value = value.and.BA_free_heap(ameta(2))
            value = value.and.BA_free_heap(bmeta(2))
            value = value.and.BA_free_heap(nxmeta(2))
            value = value.and.BA_free_heap(indxmeta(2))
            value = value.and.BA_free_heap(parammeta(2))
            value = value.and.BA_free_heap(sindxmeta(2))
            value = value.and.BA_free_heap(sparammeta(2))
            value = value.and.BA_free_heap(ymeta(2))
            call nwpw_interp_end()
            call meta_gauss_write_end()

            call nwpw_expression_end(rtdb)

         end if

      end if

      return
      end


c     **********************************************
c     *                                            *
c     *            meta_collective                 *
c     *                                            *
c     **********************************************
      real*8 function meta_collective(indx,param,ispin,ne,psi1)
      implicit none
      integer indx(*)
      real*8  param(*)
      integer ispin,ne(2)
      real*8  psi1(*)

#include "bafdecls.fh"

*     **** local variables ****
      logical sprik
      integer i,j,ntype,n1,n2,isgn,nion
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz,r,r1,r2,r3
      real*8 dx1,dy1,dz1,dx3,dy3,dz3,theta,n,m,r0,f,la,lb
      real*8 nx,ny,nz,x4,y4,z4,dx2,dy2,dz2,d1a,d2a,d1b,d2b
      real*8 cf(100),v0,dv0,n12

*     **** external functions ****
      integer  ion_rion_ptr,ion_nion
      external ion_rion_ptr,ion_nion
      real*8   ion_rion,psp_ld_trace,metadynamics_coordspherediff
      external ion_rion,psp_ld_trace,metadynamics_coordspherediff
      real*8   nwpw_expression_f
      external nwpw_expression_f

      ntype = indx(1)

      f = 0.0d0
c     **** bond distance ***
      if (ntype.eq.1) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         dx = x1-x2
         dy = y1-y2
         dz = z1-z2
         call lattice_min_difference(dx,dy,dz)
         r  = dsqrt(dx**2 + dy**2 + dz**2)
         f  = r

c     **** bond angle ***
      else if (ntype.eq.2) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         x3 = ion_rion(1,indx(4))
         y3 = ion_rion(2,indx(4))
         z3 = ion_rion(3,indx(4))
         dx1 = x1-x2
         dy1 = y1-y2
         dz1 = z1-z2
         call lattice_min_difference(dx1,dy1,dz1)
         r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
         dx3 = x3-x2
         dy3 = y3-y2
         dz3 = z3-z2
         call lattice_min_difference(dx3,dy3,dz3)
         r3  = dsqrt(dx3**2 + dy3**2 + dz3**2)
         theta = (dx1*dx3 + dy1*dy3 + dz1*dz3)/(r1*r3)
         if (theta.gt.1.0d0)  theta = 1.0d0
         if (theta.lt.-1.0d0) theta = -1.0d0
         theta = dacos(theta)
         f = theta

c     **** bond dihedral ****
      else if (ntype.eq.3) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         x3 = ion_rion(1,indx(4))
         y3 = ion_rion(2,indx(4))
         z3 = ion_rion(3,indx(4))
         x4 = ion_rion(1,indx(5))
         y4 = ion_rion(2,indx(5))
         z4 = ion_rion(3,indx(5))
         dx1 = x2-x1
         dy1 = y2-y1
         dz1 = z2-z1
         call lattice_min_difference(dx1,dy1,dz1)
         dx2 = x3-x2
         dy2 = y3-y2
         dz2 = z3-z2
         call lattice_min_difference(dx2,dy2,dz2)
         dx3 = x4-x3
         dy3 = y4-y3
         dz3 = z4-z3
         call lattice_min_difference(dx3,dy3,dz3)

         dy = dx1*(dy2*dz3-dy3*dz2) 
     >      + dy1*(dz2*dx3-dz3*dx2) 
     >      + dz1*(dx2*dy3-dx3*dy2)
         dy = dy*dsqrt(dx2**2 + dy2**2 + dz2**2)

         dx = (dy2*dz3 - dy3*dz2)*(dy1*dz2 - dy2*dz1)
     >      + (dz2*dx3 - dz3*dx2)*(dz1*dx2 - dz2*dx1)
     >      + (dx2*dy3 - dx3*dy2)*(dx1*dy2 - dx2*dy1)
         f = datan2(dy,dx)

c     **** coordination number ****
      else if (ntype.eq.4) then
         n1=indx(2)
         n2=indx(3)
         if (param(6).lt.0.0d0) then
            f = 0.0d0
            n  = param(3)
            m  = param(4)
            r0 = param(5)
            n12 = param(7)
            do i=1,n1
               do j=1,n2
                  x1 = ion_rion(1,indx(3+i))
                  y1 = ion_rion(2,indx(3+i))
                  z1 = ion_rion(3,indx(3+i))
                  x2 = ion_rion(1,indx(3+n1+j))
                  y2 = ion_rion(2,indx(3+n1+j))
                  z2 = ion_rion(3,indx(3+n1+j))
                  dx1 = x1-x2
                  dy1 = y1-y2
                  dz1 = z1-z2
                  call lattice_min_difference(dx1,dy1,dz1)
                  r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
                  f = f + (1.0d0-(r1/r0)**n)/(1.0d0-(r1/r0)**m)
               end do
            end do
            f = f/n12
         else
            f = 0.0d0
            n  = param(3)
            m  = param(4)
            r0 = param(5)
            n12 = param(7)
            do i=1,n1
               do j=1,n2
                  x1 = ion_rion(1,indx(3+i))
                  y1 = ion_rion(2,indx(3+i))
                  z1 = ion_rion(3,indx(3+i))
                  x2 = ion_rion(1,indx(3+n1+j))
                  y2 = ion_rion(2,indx(3+n1+j))
                  z2 = ion_rion(3,indx(3+n1+j))
                  dx1 = x1-x2
                  dy1 = y1-y2
                  dz1 = z1-z2
                  call lattice_min_difference(dx1,dy1,dz1)
                  r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
c                 f = f + (1.0d0-(r1/r0)**n)/(1.0d0-(r1/r0)**m)
                  f = f + 1.d0/(1.d0+dexp(n*(r1-r0)))
               end do
            end do
            f = f/n12
         end if

c     **** n-plane ****
      else if (ntype.eq.5) then
         f = 0.0d0
         nx = param(3)
         ny = param(4)
         nz = param(5)
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         n1=indx(3)
         do i=1,n1
            x2 = ion_rion(1,indx(3+i))
            y2 = ion_rion(2,indx(3+i))
            z2 = ion_rion(3,indx(3+i))
            dx = x1-x2
            dy = y1-y2
            dz = z1-z2
            call lattice_min_difference(dx,dy,dz)
            r =  nx*dx + ny*dy + nz*dz
            f = f + dsqrt(r*r)
         end do
         f = f/dble(n1)

*     **** x,y,z ****
      else if (ntype.eq.6) then
         if (indx(2).eq.1)  then
            f = ion_rion(1,indx(3))
         else if (indx(2).eq.2)  then
            f = ion_rion(2,indx(3))
         else
            f = ion_rion(3,indx(3))
         end if

*     **** ld_trace ****
      else if (ntype.eq.7) then
         f =  psp_ld_trace(indx(2),indx(3),ispin,ne,psi1)

*     *** 2bonds ****
      else if (ntype.eq.8) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         x3 = ion_rion(1,indx(4))
         y3 = ion_rion(2,indx(4))
         z3 = ion_rion(3,indx(4))
         x4 = ion_rion(1,indx(5))
         y4 = ion_rion(2,indx(5))
         z4 = ion_rion(3,indx(5))
         dx1 = x2-x1
         dy1 = y2-y1
         dz1 = z2-z1
         call lattice_min_difference(dx1,dy1,dz1)
         dx2 = x4-x3
         dy2 = y4-y3
         dz2 = z4-z3
         call lattice_min_difference(dx2,dy2,dz2)
         d1a = param(3)
         d2a = param(4)
         d1b = param(5)
         d2b = param(6)
         r1 = dsqrt(dx1**2 + dy1**2 + dz1**2)
         r2 = dsqrt(dx2**2 + dy2**2 + dz2**2)
         la = dsqrt((r1-d1a)**2 + (r2-d2a)**2)
         lb = dsqrt((r1-d1b)**2 + (r2-d2b)**2)
         f = la/(la+lb)

*     **** ld_trace2diff ****
      else if (ntype.eq.9) then
         f =  psp_ld_trace(indx(2),indx(4),ispin,ne,psi1)
     >     -  psp_ld_trace(indx(3),indx(4),ispin,ne,psi1)

*     *** bond_difference2 ****
      else if (ntype.eq.10) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         x3 = ion_rion(1,indx(4))
         y3 = ion_rion(2,indx(4))
         z3 = ion_rion(3,indx(4))
         dx1 = x2-x1
         dy1 = y2-y1
         dz1 = z2-z1
         call lattice_min_difference(dx1,dy1,dz1)
         r1 = (dx1**2 + dy1**2 + dz1**2)
         dx2 = x3-x1
         dy2 = y3-y1
         dz2 = z3-z1
         call lattice_min_difference(dx2,dy2,dz2)
         r2 = (dx2**2 + dy2**2 + dz2**2)
         f = r1 - r2

c     **** bond_difference_shell ****
      else if (ntype.eq.11) then
         n1=indx(2)
         n2=indx(3)
         n   = param(3)
         m   = param(4)
         r0  = param(5)
         sprik = (param(6).ge.0.0d0)
         v0  = param(7)
         dv0 = param(8)
         f = metadynamics_coordspherediff(v0,dv0,sprik,n,m,r0,
     >                                    n1,indx(3   +1),
     >                                    n2,indx(3+n1+1),
     >                                    dbl_mb(ion_rion_ptr()))

      else if (ntype.eq.12) then
          nion = ion_nion()
          f = nwpw_expression_f(indx(2),nion,dbl_mb(ion_rion_ptr()))

c     **** group distance ****
      else if (ntype.eq.13) then
         n1=indx(2)
         n2=indx(3)
         f = 0.0d0
         n  = param(3)
         n12 = dble(n1*n2)
         do i=1,n1
            do j=1,n2
               x1 = ion_rion(1,indx(3+i))
               y1 = ion_rion(2,indx(3+i))
               z1 = ion_rion(3,indx(3+i))
               x2 = ion_rion(1,indx(3+n1+j))
               y2 = ion_rion(2,indx(3+n1+j))
               z2 = ion_rion(3,indx(3+n1+j))
               dx1 = x1-x2
               dy1 = y1-y2
               dz1 = z1-z2
               call lattice_min_difference(dx1,dy1,dz1)
               r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
               f = f + (1.0/r1)**n
            end do
         end do
         f = (f/n12)**(-1.0d0/n)
  
c     **** bondings ****
      else if (ntype.eq.14) then
         n1=indx(2)
         f = 0.0d0
         do i=1,n1
            x1 = ion_rion(1,indx(2*i+1))
            y1 = ion_rion(2,indx(2*i+1))
            z1 = ion_rion(3,indx(2*i+1))
            x2 = ion_rion(1,indx(2*i+2))
            y2 = ion_rion(2,indx(2*i+2))
            z2 = ion_rion(3,indx(2*i+2))
            dx1 = x1-x2
            dy1 = y1-y2
            dz1 = z1-z2
            call lattice_min_difference(dx1,dy1,dz1)
            r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
            f = f + param(i+2)*r1
         end do


      end if

      meta_collective = f
      return
      end

c     **********************************************
c     *                                            *
c     *           meta_collective_force            *
c     *                                            *
c     **********************************************
      subroutine meta_collective_force(indx,param,dv,
     >                                 ispin,ne,psi1,fmeta_psi1,
     >                                 move,fmeta)
      implicit none
      integer indx(*)
      real*8 param(*),dv
      integer ispin,ne(2)
      real*8 psi1(*),fmeta_psi1(*)
      logical move
      real*8 fmeta(3,*)

#include "bafdecls.fh"

*     **** local variables ****
      logical sprik
      integer i,j,ntype,n1,n2,nion
      real*8 x1,y1,z1,r1,vx1,vx2,vy1,vy2,vz1,vz2
      real*8 x2,y2,z2,aa,a11,a12,a22,r,df,f,n,m,r0,rn,rm
      real*8 x3,y3,z3,r3,ctheta,stheta,denom
      real*8 dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3,nx,ny,nz,la,lb
      real*8 x4,y4,z4,dx2,dy2,dz2,r2,df1,df2,d1a,d2a,d1b,d2b
      real*8 dfla,dflb,dlar1,dlar2,dlbr1,dlbr2
      real*8 v0,dv0,n12

*     **** external functions ****
      integer  ion_rion_ptr,ion_nion
      external ion_rion_ptr,ion_nion
      real*8   ion_rion
      external ion_rion

      ntype = indx(1)

      if (move) then

c     **** bond distance ***
      if (ntype.eq.1) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         dx = x1-x2
         dy = y1-y2
         dz = z1-z2
         call lattice_min_difference(dx,dy,dz)
         r  = dsqrt(dx**2 + dy**2 + dz**2)

         fmeta(1,indx(2)) = fmeta(1,indx(2)) - (dx/r)*dv
         fmeta(2,indx(2)) = fmeta(2,indx(2)) - (dy/r)*dv
         fmeta(3,indx(2)) = fmeta(3,indx(2)) - (dz/r)*dv
         fmeta(1,indx(3)) = fmeta(1,indx(3)) + (dx/r)*dv
         fmeta(2,indx(3)) = fmeta(2,indx(3)) + (dy/r)*dv
         fmeta(3,indx(3)) = fmeta(3,indx(3)) + (dz/r)*dv

c     **** bond angle ***
      else if (ntype.eq.2) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         x3 = ion_rion(1,indx(4))
         y3 = ion_rion(2,indx(4))
         z3 = ion_rion(3,indx(4))
         dx1 = x1-x2
         dy1 = y1-y2
         dz1 = z1-z2
         call lattice_min_difference(dx1,dy1,dz1)
         r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
         dx3 = x3-x2
         dy3 = y3-y2
         dz3 = z3-z2
         call lattice_min_difference(dx3,dy3,dz3)
         r3  = dsqrt(dx3**2 + dy3**2 + dz3**2)
         denom = r1*r3
         if (denom.gt.1.0d-11) then
            ctheta = (dx1*dx3 + dy1*dy3 + dz1*dz3)/(denom)
            if (ctheta.gt.1.0d0)  ctheta =  1.0d0
            if (ctheta.lt.-1.0d0) ctheta = -1.0d0
            stheta = dsqrt(1.0d0-ctheta*ctheta)
            if (stheta.lt.0.001d0) stheta = 0.001d0
            stheta = 1.0d0/stheta

            aa  =  dv*stheta
            a11 =  aa*ctheta/r1
            a12 = -aa/(denom)
            a22 =  aa*ctheta/r3

            vx1 = a11*dx1 + a12*dx3
            vx2 = a22*dx3 + a12*dx1

            vy1 = a11*dy1 + a12*dy3
            vy2 = a22*dy3 + a12*dy1

            vz1 = a11*dz1 + a12*dz3
            vz2 = a22*dz3 + a12*dz1

            fmeta(1,indx(2)) = fmeta(1,indx(2)) - vx1
            fmeta(2,indx(2)) = fmeta(2,indx(2)) - vy1
            fmeta(3,indx(2)) = fmeta(3,indx(2)) - vz1

            fmeta(1,indx(3)) = fmeta(1,indx(3)) + vx1 + vx2
            fmeta(2,indx(3)) = fmeta(2,indx(3)) + vy1 + vy2
            fmeta(3,indx(3)) = fmeta(3,indx(3)) + vz1 + vz2

            fmeta(1,indx(4)) = fmeta(1,indx(4)) - vx2
            fmeta(2,indx(4)) = fmeta(2,indx(4)) - vy2
            fmeta(3,indx(4)) = fmeta(3,indx(4)) - vz2
         end if

c     **** bond dihedral ****
      else if (ntype.eq.3) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         x3 = ion_rion(1,indx(4))
         y3 = ion_rion(2,indx(4))
         z3 = ion_rion(3,indx(4))
         x4 = ion_rion(1,indx(5))
         y4 = ion_rion(2,indx(5))
         z4 = ion_rion(3,indx(5))
         dx1 = x2-x1
         dy1 = y2-y1
         dz1 = z2-z1
         call lattice_min_difference(dx1,dy1,dz1)
         dx2 = x3-x2
         dy2 = y3-y2
         dz2 = z3-z2
         call lattice_min_difference(dx2,dy2,dz2)
         dx3 = x4-x3
         dy3 = y4-y3
         dz3 = z4-z3
         call lattice_min_difference(dx3,dy3,dz3)

         dy = dx1*(dy2*dz3-dy3*dz2)
     >      + dy1*(dz2*dx3-dz3*dx2)
     >      + dz1*(dx2*dy3-dx3*dy2)
         dy = dy*dsqrt(dx2**2 + dy2**2 + dz2**2)

         dx = (dy2*dz3 - dy3*dz2)*(dy1*dz2 - dy2*dz1)
     >      + (dz2*dx3 - dz3*dx2)*(dz1*dx2 - dz2*dx1)
     >      + (dx2*dy3 - dx3*dy2)*(dx1*dy2 - dx2*dy1)
         r2 = datan2(dy,dx)
         vx1 = -dy/(dx**2 + dy**2)
         vy1 =  dx/(dx**2 + dy**2)

c          dphi/dx  (dx/db1 db1/dr1 + dx/db2 db2/dr1 + dx/db3 db3/dr1)
c          dphi/dy  (dy/db1 db1/dr1 + dy/db2 db2/dr1 + dy/db3 db3/dr1)

            fmeta(1,indx(2)) = fmeta(1,indx(2)) 
            fmeta(2,indx(2)) = fmeta(2,indx(2)) 
            fmeta(3,indx(2)) = fmeta(3,indx(2)) 

            fmeta(1,indx(3)) = fmeta(1,indx(3)) 
            fmeta(2,indx(3)) = fmeta(2,indx(3)) 
            fmeta(3,indx(3)) = fmeta(3,indx(3)) 

            fmeta(1,indx(4)) = fmeta(1,indx(4)) 
            fmeta(2,indx(4)) = fmeta(2,indx(4)) 
            fmeta(3,indx(4)) = fmeta(3,indx(4)) 

            fmeta(1,indx(5)) = fmeta(1,indx(5)) 
            fmeta(2,indx(5)) = fmeta(2,indx(5)) 
            fmeta(3,indx(5)) = fmeta(3,indx(5)) 



c     **** coordination number ***
      else if (ntype.eq.4) then
         n1 = indx(2)
         n2 = indx(3)
         if (param(6).lt.0.0d0) then
            n  = param(3)
            m  = param(4)
            r0 = param(5)
            n12= param(7)
            do i=1,n1
               do j=1,n2
                  x1 = ion_rion(1,indx(3+i))
                  y1 = ion_rion(2,indx(3+i))
                  z1 = ion_rion(3,indx(3+i))
                  x2 = ion_rion(1,indx(3+n1+j))
                  y2 = ion_rion(2,indx(3+n1+j))
                  z2 = ion_rion(3,indx(3+n1+j))
                  dx = x1-x2
                  dy = y1-y2
                  dz = z1-z2
                  call lattice_min_difference(dx,dy,dz)
                  r  = dsqrt(dx**2 + dy**2 + dz**2)
                  rn  = (1.0d0-(r/r0)**n)
                  rm  = (1.0d0-(r/r0)**m)
                  df = (-n*rm/r0*(r/r0)**(n-1) + m*rn/r0*(r/r0)**(m-1))
     >                 / (rm**2)
                  df = df/n12
                  fmeta(1,indx(3+i)) = fmeta(1,indx(3+i)) - (dx/r)*df*dv
                  fmeta(2,indx(3+i)) = fmeta(2,indx(3+i)) - (dy/r)*df*dv
                  fmeta(3,indx(3+i)) = fmeta(3,indx(3+i)) - (dz/r)*df*dv
                  fmeta(1,indx(3+n1+j)) = fmeta(1,indx(3+n1+j)) 
     >                                  + (dx/r)*df*dv
                  fmeta(2,indx(3+n1+j)) = fmeta(2,indx(3+n1+j)) 
     >                                  + (dy/r)*df*dv
                  fmeta(3,indx(3+n1+j)) = fmeta(3,indx(3+n1+j)) 
     >                                  + (dz/r)*df*dv
               end do
            end do
         else
            n  = param(3)
            m  = param(4)
            r0 = param(5)
            n12= param(7)
            do i=1,n1
               do j=1,n2
                  x1 = ion_rion(1,indx(3+i))
                  y1 = ion_rion(2,indx(3+i))
                  z1 = ion_rion(3,indx(3+i))
                  x2 = ion_rion(1,indx(3+n1+j))
                  y2 = ion_rion(2,indx(3+n1+j))
                  z2 = ion_rion(3,indx(3+n1+j))
                  dx = x1-x2
                  dy = y1-y2
                  dz = z1-z2
                  call lattice_min_difference(dx,dy,dz)
                  r  = dsqrt(dx**2 + dy**2 + dz**2)
c                 rn  = (1.0d0-(r/r0)**n)
c                 rm  = (1.0d0-(r/r0)**m)
c                 df = (-n*rm/r0*(r/r0)**(n-1) + m*rn/r0*(r/r0)**(m-1))
c    >                 / (rm**2)
                  rn  = 1.d0+dexp(n*(r-r0))
                  df = -n*(rn-1.d0)/(rn*rn)
                  df = df/n12
                  fmeta(1,indx(3+i)) = fmeta(1,indx(3+i)) - (dx/r)*df*dv
                  fmeta(2,indx(3+i)) = fmeta(2,indx(3+i)) - (dy/r)*df*dv
                  fmeta(3,indx(3+i)) = fmeta(3,indx(3+i)) - (dz/r)*df*dv
                  fmeta(1,indx(3+n1+j)) = fmeta(1,indx(3+n1+j)) 
     >                                  + (dx/r)*df*dv
                  fmeta(2,indx(3+n1+j)) = fmeta(2,indx(3+n1+j)) 
     >                                  + (dy/r)*df*dv
                  fmeta(3,indx(3+n1+j)) = fmeta(3,indx(3+n1+j)) 
     >                                  + (dz/r)*df*dv
               end do
            end do
         endif

c     **** n-plane ****
      else if (ntype.eq.5) then
         nx = param(3)
         ny = param(4)
         nz = param(5)
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         n1=indx(3)
         do i=1,n1
            x2 = ion_rion(1,indx(3+i))
            y2 = ion_rion(2,indx(3+i))
            z2 = ion_rion(3,indx(3+i))
            dx = x1-x2
            dy = y1-y2
            dz = z1-z2
            call lattice_min_difference(dx,dy,dz)
            r =  nx*dx + ny*dy + nz*dz
            r0 =  dsqrt(r*r)
            df = (r/r0)*dv/dble(n1)
            fmeta(1,indx(2)) = fmeta(1,indx(2)) - nx*df
            fmeta(2,indx(2)) = fmeta(2,indx(2)) - ny*df
            fmeta(3,indx(2)) = fmeta(3,indx(2)) - nz*df

            fmeta(1,indx(3+i)) = fmeta(1,indx(3+i)) + nx*df
            fmeta(2,indx(3+i)) = fmeta(2,indx(3+i)) + ny*df
            fmeta(3,indx(3+i)) = fmeta(3,indx(3+i)) + nz*df
         end do

*     **** x,y,z ****
      else if (ntype.eq.6) then
         if (indx(2).eq.1)  then
            fmeta(1,indx(3)) = fmeta(1,indx(3)) - dv
         else if (indx(2).eq.2)  then
            fmeta(2,indx(3)) = fmeta(2,indx(3)) - dv
         else
            fmeta(3,indx(3)) = fmeta(3,indx(3)) - dv
         end if

c     **** ld_trace ***
      else if (ntype.eq.7) then
         call psp_ld_trace_gradient(indx(2),indx(3),ispin,ne,psi1,
     >                                 .false.,param,
     >                                 (-dv),fmeta_psi1,
     >                                 move,fmeta)

c     **** 2bonds ****
      else if (ntype.eq.8) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         x3 = ion_rion(1,indx(4))
         y3 = ion_rion(2,indx(4))
         z3 = ion_rion(3,indx(4))
         x4 = ion_rion(1,indx(5))
         y4 = ion_rion(2,indx(5))
         z4 = ion_rion(3,indx(5))
         dx1 = x2-x1
         dy1 = y2-y1
         dz1 = z2-z1
         call lattice_min_difference(dx1,dy1,dz1)
         r1 = dsqrt(dx1**2 + dy1**2 + dz1**2)
         dx2 = x4-x3
         dy2 = y4-y3
         dz2 = z4-z3
         call lattice_min_difference(dx2,dy2,dz2)
         r2 = dsqrt(dx2**2 + dy2**2 + dz2**2)
         d1a = param(3)
         d2a = param(4)
         d1b = param(5)
         d2b = param(6)
         la = dsqrt((r1-d1a)**2 + (r2-d2a)**2)
         lb = dsqrt((r1-d1b)**2 + (r2-d2b)**2)
         if (dabs(la).lt.1.0d-6) then
            dlar1 = dtanh((r1-d1a)/1.0d-7)
            dlar2 = dtanh((r2-d2a)/1.0d-7)
         else
            dlar1 = (r1-d1a)/la
            dlar2 = (r2-d2a)/la
         end if
         if (dabs(lb).lt.1.0d-6) then
            dlbr1 = dtanh((r1-d1b)/1.0d-7)
            dlbr2 = dtanh((r2-d2b)/1.0d-7)
         else
            dlbr1 = (r1-d1b)/lb
            dlbr2 = (r2-d2b)/lb
         end if
         !f = la/(la+lb)
         dfla = 1.0d0/(la+lb) - la/(la+lb)**2
         dflb = -la/(la+lb)**2
         df1 = dfla*dlar1 + dflb*dlbr1
         df2 = dfla*dlar2 + dflb*dlbr2

         fmeta(1,indx(2)) = fmeta(1,indx(2)) + (dx1/r1)*dv*df1
         fmeta(2,indx(2)) = fmeta(2,indx(2)) + (dy1/r1)*dv*df1
         fmeta(3,indx(2)) = fmeta(3,indx(2)) + (dz1/r1)*dv*df1

         fmeta(1,indx(3)) = fmeta(1,indx(3)) - (dx1/r1)*dv*df1
         fmeta(2,indx(3)) = fmeta(2,indx(3)) - (dy1/r1)*dv*df1
         fmeta(3,indx(3)) = fmeta(3,indx(3)) - (dz1/r1)*dv*df1

         fmeta(1,indx(4)) = fmeta(1,indx(4)) + (dx2/r2)*dv*df2
         fmeta(2,indx(4)) = fmeta(2,indx(4)) + (dy2/r2)*dv*df2
         fmeta(3,indx(4)) = fmeta(3,indx(4)) + (dz2/r2)*dv*df2

         fmeta(1,indx(5)) = fmeta(1,indx(5)) - (dx2/r2)*dv*df2
         fmeta(2,indx(5)) = fmeta(2,indx(5)) - (dy2/r2)*dv*df2
         fmeta(3,indx(5)) = fmeta(3,indx(5)) - (dz2/r2)*dv*df2

c     **** ld_trace2diff ***
      else if (ntype.eq.9) then
         call psp_ld_trace_gradient(indx(2),indx(4),ispin,ne,psi1,
     >                                 .false.,param,
     >                                 (-dv),fmeta_psi1,
     >                                 move,fmeta)
         call psp_ld_trace_gradient(indx(3),indx(4),ispin,ne,psi1,
     >                                 .false.,param,
     >                                 (dv),fmeta_psi1,
     >                                 move,fmeta)

c     **** bond_difference2 ****
      else if (ntype.eq.10) then
         x1 = ion_rion(1,indx(2))
         y1 = ion_rion(2,indx(2))
         z1 = ion_rion(3,indx(2))
         x2 = ion_rion(1,indx(3))
         y2 = ion_rion(2,indx(3))
         z2 = ion_rion(3,indx(3))
         x3 = ion_rion(1,indx(4))
         y3 = ion_rion(2,indx(4))
         z3 = ion_rion(3,indx(4))
         dx1 = x2-x1
         dy1 = y2-y1
         dz1 = z2-z1
         call lattice_min_difference(dx1,dy1,dz1)
         r1 = (dx1**2 + dy1**2 + dz1**2)
         dx2 = x3-x1
         dy2 = y3-y1
         dz2 = z3-z1
         call lattice_min_difference(dx2,dy2,dz2)
         r2 = (dx2**2 + dy2**2 + dz2**2)
         !f = r1 - r2
         df1 =  1
         df2 = -1
         fmeta(1,indx(2)) = fmeta(1,indx(2)) +2.0d0*(dx1*df1+dx2*df2)*dv
         fmeta(2,indx(2)) = fmeta(2,indx(2)) +2.0d0*(dy1*df1+dy2*df2)*dv
         fmeta(3,indx(2)) = fmeta(3,indx(2)) +2.0d0*(dz1*df1+dz2*df2)*dv

         fmeta(1,indx(3)) = fmeta(1,indx(3)) - (2.0d0*dx1*df1)*dv
         fmeta(2,indx(3)) = fmeta(2,indx(3)) - (2.0d0*dy1*df1)*dv
         fmeta(3,indx(3)) = fmeta(3,indx(3)) - (2.0d0*dz1*df1)*dv

         fmeta(1,indx(4)) = fmeta(1,indx(4)) - (2.0d0*dx2*df2)*dv
         fmeta(2,indx(4)) = fmeta(2,indx(4)) - (2.0d0*dy2*df2)*dv
         fmeta(3,indx(4)) = fmeta(3,indx(4)) - (2.0d0*dz2*df2)*dv

      else if (ntype.eq.11) then
        n1=indx(2)
        n2=indx(3)
        n   = param(3)
        m   = param(4)
        r0  = param(5)
        sprik = (param(6).ge.0.0d0)
        v0  = param(7)
        dv0 = param(8)
        call metadynamics_coordspherediff_force(dv,v0,dv0,
     >                                        sprik,n,m,r0,
     >                                        n1,indx(3   +1),
     >                                        n2,indx(3+n1+1),
     >                                        dbl_mb(ion_rion_ptr()),
     >                                        fmeta)
      else if (ntype.eq.12) then
         nion = ion_nion()
         call nwpw_expression_fion(dv,indx(2),
     >                        nion,dbl_mb(ion_rion_ptr()),
     >                        fmeta)

c     **** group distances ****
      else if (ntype.eq.13) then
         n1=indx(2)
         n2=indx(3)
         f = 0.0d0
         n  = param(3)
         n12 = dble(n1*n2)
         do i=1,n1
            do j=1,n2
               x1 = ion_rion(1,indx(3+i))
               y1 = ion_rion(2,indx(3+i))
               z1 = ion_rion(3,indx(3+i))
               x2 = ion_rion(1,indx(3+n1+j))
               y2 = ion_rion(2,indx(3+n1+j))
               z2 = ion_rion(3,indx(3+n1+j))
               dx1 = x1-x2
               dy1 = y1-y2
               dz1 = z1-z2
               call lattice_min_difference(dx1,dy1,dz1)
               r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
               f = f + (1.0/r1)**n
            end do
         end do
         f = (f/n12)**(-1.0d0/n)

         do i=1,n1
            do j=1,n2
               x1 = ion_rion(1,indx(3+i))
               y1 = ion_rion(2,indx(3+i))
               z1 = ion_rion(3,indx(3+i))
               x2 = ion_rion(1,indx(3+n1+j))
               y2 = ion_rion(2,indx(3+n1+j))
               z2 = ion_rion(3,indx(3+n1+j))
               dx = x1-x2
               dy = y1-y2
               dz = z1-z2
               call lattice_min_difference(dx,dy,dz)
               r1  = dsqrt(dx**2 + dy**2 + dz**2)
               df = ((f/r1)**(n+1))/n12

               fmeta(1,indx(3+i)) = fmeta(1,indx(3+i)) - (dx/r1)*df*dv
               fmeta(2,indx(3+i)) = fmeta(2,indx(3+i)) - (dy/r1)*df*dv
               fmeta(3,indx(3+i)) = fmeta(3,indx(3+i)) - (dz/r1)*df*dv
               fmeta(1,indx(3+n1+j)) = fmeta(1,indx(3+n1+j)) 
     >                               + (dx/r1)*df*dv
               fmeta(2,indx(3+n1+j)) = fmeta(2,indx(3+n1+j)) 
     >                               + (dy/r1)*df*dv
               fmeta(3,indx(3+n1+j)) = fmeta(3,indx(3+n1+j)) 
     >                               + (dz/r1)*df*dv


            end do
         end do

c     **** bondings ****
      else if (ntype.eq.14) then
         n1=indx(2)
         do i=1,n1
            x1 = ion_rion(1,indx(2*i+1))
            y1 = ion_rion(2,indx(2*i+1))
            z1 = ion_rion(3,indx(2*i+1))
            x2 = ion_rion(1,indx(2*i+2))
            y2 = ion_rion(2,indx(2*i+2))
            z2 = ion_rion(3,indx(2*i+2))
            dx1 = x1-x2
            dy1 = y1-y2
            dz1 = z1-z2
            call lattice_min_difference(dx1,dy1,dz1)
            r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
            fmeta(1,indx(2*i+1)) = fmeta(1,indx(2*i+1))
     >                           - param(i+2)*dx1/r1*dv
            fmeta(2,indx(2*i+1)) = fmeta(2,indx(2*i+1))
     >                           - param(i+2)*dy1/r1*dv
            fmeta(3,indx(2*i+1)) = fmeta(3,indx(2*i+1))
     >                           - param(i+2)*dz1/r1*dv
            fmeta(1,indx(2*i+2)) = fmeta(1,indx(2*i+2))
     >                           + param(i+2)*dx1/r1*dv
            fmeta(2,indx(2*i+2)) = fmeta(2,indx(2*i+2))
     >                           + param(i+2)*dy1/r1*dv
            fmeta(3,indx(2*i+2)) = fmeta(3,indx(2*i+2))
     >                           + param(i+2)*dz1/r1*dv
         end do

      end if

*     ***** move = .false. ******
      else
         if (ntype.eq.7) then
         call psp_ld_trace_gradient(indx(2),indx(3),ispin,ne,psi1,
     >                                 .false.,param,
     >                                 (-dv),fmeta_psi1,
     >                                 move,fmeta)
         else if (ntype.eq.9) then
         call psp_ld_trace_gradient(indx(2),indx(4),ispin,ne,psi1,
     >                                 .false.,param,
     >                                 (-dv),fmeta_psi1,
     >                                 move,fmeta)
         call psp_ld_trace_gradient(indx(3),indx(4),ispin,ne,psi1,
     >                                 .false.,param,
     >                                 (dv),fmeta_psi1,
     >                                 move,fmeta)
         end if
      end if

      return
      end



*     ***********************************************
*     *                                             *
*     *              meta_update                    *
*     *                                             *
*     ***********************************************
      subroutine meta_update(ispin,ne,psi1,E)
      implicit none
      integer ispin,ne(2)
      real*8 psi1(*)
      real*8 E(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      double precision kb
      parameter (kb=3.16679d-6)

      integer d,i,j(10),s1,s2
      real*8  r0(10),sigma0(10),r,x,w,ww,sigma,v,expv

*     **** external functions ****
      real*8   meta_collective,meta_energy
      external meta_collective,meta_energy

      if (metafound) then

         metacount = metacount + 1
         if (metacount.ge.maxmetacount) then
         metacount = 0

         if (dTtempered.gt.0.0d0) then
            v    = meta_energy(ispin,ne,psi1)
            expv = dexp(-v/(kb*dTtempered))
         else
            v    = 0.0d0
            expv = 1.0d0
         end if

         if (nmeta.gt.0) then
            w = 1.0d0
            do d=1,nmeta
               s1 = int_mb(sindxmeta(1) +d-1)
               s2 = int_mb(sparammeta(1)+d-1)
               r0(d) = meta_collective(int_mb(indxmeta(1)+s1),
     >                                 dbl_mb(parammeta(1)+s2),
     >                                 ispin,ne,psi1)
               sigma0(d) = dbl_mb(parammeta(1)+s2+1)
               w  = w*dbl_mb(parammeta(1)+s2)
               j(d) = 0
            end do
            w = w**(1.0d0/dble(nmeta))
           
            do i=1,nxmeta_all
               ww = 1.0d0
               do d=1,nmeta
                  s1 = int_mb(sindxmeta(1) +d-1)
                  s2 = int_mb(sparammeta(1)+d-1)
c                  r0    = meta_collective(int_mb(indxmeta(1)+s1),
c     >                                    dbl_mb(parammeta(1)+s2),
c     >                                    ispin,ne,psi1)
                  sigma = dbl_mb(parammeta(1)+s2+1)
                  if (int_mb(pmeta(1)+d-1).eq.0) then
                     x = dbl_mb(ameta(1)+d-1) 
     >                + j(d)*(dbl_mb(bmeta(1)+d-1)-dbl_mb(ameta(1)+d-1))
     >                    /dble(int_mb(nxmeta(1)+d-1)-1)
                  else
                     x = dbl_mb(ameta(1)+d-1) 
     >                + j(d)*(dbl_mb(bmeta(1)+d-1)-dbl_mb(ameta(1)+d-1))
     >                      /dble(int_mb(nxmeta(1)+d-1))
                  end if
                  r = x-r0(d)
                  ww = ww * dexp(-0.5d0*(r/sigma)**2)
               end do
            
               dbl_mb(ymeta(1)+i-1) = dbl_mb(ymeta(1)+i-1) + ww*w*expv

               j(1) = j(1) + 1
               do d=1,nmeta-1
                  if (j(d).ge.int_mb(nxmeta(1)+d-1)) then
                     j(d) = 0
                     j(d+1) = j(d+1)+1
                  end if
               end do
            end do
            call meta_gauss_write(nmeta,r0,sigma0,w*expv,E)
         end if

         end if

         metaprintcount = metaprintcount + 1
         if (metaprintcount.ge.maxmetaprintcount) then
            metaprintcount = 0
            metaraycount   = metaraycount + metarayshift
            call meta_print_potentials(metaraycount)
         end if

      end if
 
      return
      end

*     ***********************************************
*     *                                             *
*     *              meta_force                     *
*     *                                             *
*     ***********************************************
      subroutine meta_force(ispin,ne,psi1,fmeta_psi1,move,fmeta)
      implicit none
      integer ispin,ne(2)
      real*8 psi1(*),fmeta_psi1(*)
      logical move
      real*8 fmeta(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer d,s1,s2
      real*8 r(10),dv(10)

*     **** external functions ****
      integer  ion_nion
      external ion_nion
      real*8   meta_collective
      external meta_collective
      
      if (metafound) then

         if (nmeta.gt.0) then
            do d=1,nmeta
               s1 = int_mb(sindxmeta(1) +d-1)
               s2 = int_mb(sparammeta(1)+d-1)
               r(d) = meta_collective(int_mb(indxmeta(1)+s1),
     >                                dbl_mb(parammeta(1)+s2),
     >                                ispin,ne,psi1)
            end do

            call nwpw_dinterp(nmeta,int_mb(nxmeta(1)),
     >                        dbl_mb(ameta(1)),dbl_mb(bmeta(1)),
     >                        int_mb(pmeta(1)),
     >                        dbl_mb(ymeta(1)),r,dv)

            do d=1,nmeta
               s1 = int_mb(sindxmeta(1) +d-1)
               s2 = int_mb(sparammeta(1)+d-1)
               call meta_collective_force(int_mb(indxmeta(1) +s1),
     >                                    dbl_mb(parammeta(1)+s2),
     >                                    dv(d),
     .                                    ispin,ne,psi1,fmeta_psi1,
     >                                    move,fmeta)
            end do
      
         end if

      end if

      return
      end


*     ***********************************************
*     *                                             *
*     *              meta_energy                    *
*     *                                             *
*     ***********************************************
      real*8 function meta_energy(ispin,ne,psi1)
      implicit none
      integer ispin,ne(2)
      real*8 psi1(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer d,s1,s2
      real*8 r(10),v

*     **** external functions ****
      real*8   nwpw_interp,meta_collective
      external nwpw_interp,meta_collective

      v = 0.0d0
      if (metafound) then

         if (nmeta.gt.0) then
            do d=1,nmeta
               s1 = int_mb(sindxmeta(1) +d-1)
               s2 = int_mb(sparammeta(1)+d-1)
               r(d) = meta_collective(int_mb(indxmeta(1)+s1),
     >                                dbl_mb(parammeta(1)+s2),
     >                                ispin,ne,psi1)
            end do
            v = nwpw_interp(nmeta,int_mb(nxmeta(1)),
     >                      dbl_mb(ameta(1)),dbl_mb(bmeta(1)),
     >                      int_mb(pmeta(1)),
     >                      dbl_mb(ymeta(1)),r)

            call Parallel_Brdcst_value(0,v)  !*** syncing to master taskid ***
         end if

      end if

      meta_energy = v 
      return
      end





*     ***********************************************
*     *                                             *
*     *              meta_energypotential           *
*     *                                             *
*     ***********************************************
      subroutine meta_energypotential(ispin,ne,psi1,e,p)
      implicit none
      integer ispin,ne(2)
      real*8 psi1(*)
      real*8 e,p

#include "bafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer d,s1,s2,ntype
      real*8 r(10),dv(10),v

*     **** external functions ****
      real*8   nwpw_interp,meta_collective
      external nwpw_interp,meta_collective

      e = 0.0d0
      p = 0.0d0
      if (metafound) then

         if (nmeta.gt.0) then
            do d=1,nmeta
               s1 = int_mb(sindxmeta(1) +d-1)
               s2 = int_mb(sparammeta(1)+d-1)
               r(d) = meta_collective(int_mb(indxmeta(1)+s1),
     >                                dbl_mb(parammeta(1)+s2),
     >                                ispin,ne,psi1)
            end do
            e = nwpw_interp(nmeta,int_mb(nxmeta(1)),
     >                      dbl_mb(ameta(1)),dbl_mb(bmeta(1)),
     >                      int_mb(pmeta(1)),
     >                      dbl_mb(ymeta(1)),r)

            call nwpw_dinterp(nmeta,int_mb(nxmeta(1)),
     >                        dbl_mb(ameta(1)),dbl_mb(bmeta(1)),
     >                        int_mb(pmeta(1)),
     >                        dbl_mb(ymeta(1)),r,dv)

            do d=1,nmeta
               s1 = int_mb(sindxmeta(1) +d-1)
               ntype = int_mb(indxmeta(1)+s1)
               if ((ntype.eq.7).or.(ntype.eq.9)) then
                  p = p + r(d)*dv(d)
               end if
            end do

            call Parallel_Brdcst_value(0,e)  !*** syncing to master taskid ***
            call Parallel_Brdcst_value(0,p)
         end if

      end if

      return
      end


*     ***********************************************
*     *                                             *
*     *              meta_found                     *
*     *                                             *
*     ***********************************************
      logical function meta_found()
      implicit none

#include "metadynamics.fh"

      meta_found = metafound
      return
      end


c     **********************************************
c     *                                            *
c     *            meta_print_potentials           *
c     *                                            *
c     **********************************************
      subroutine meta_print_potentials(icount)
      implicit none
      integer icount

#include "bafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      integer d,i,j(10)
      real*8  x(10),r,v,dv(10),fac,temp
      character*80 filename
      character*255 full_filename

*     **** external functions ****
      logical  control_Nose
      external control_Nose
      character*7 c_index_name
      external    c_index_name
      real*8   nwpw_interp,control_Nose_Tr,ion_Temperature
      external nwpw_interp,control_Nose_Tr,ion_Temperature

      call Parallel_taskid(taskid)

      if (metafound) then

*        **** determine tempered metadynamics factor ****
         fac = 1.0d0
         if (dTtempered.gt.0.0d0) then
            if (control_Nose()) then
               temp = control_Nose_Tr()
            else
               temp = ion_Temperature()
            end if
            fac = (temp+dTtempered)/dTtempered
         end if


*        **** print out potentials ******
         if (nmeta.gt.0) then
            if (taskid.eq.MASTER) then
               if (icount.gt.0) then
                  filename = "meta"//c_index_name(icount)//".dat"
               else
                  filename = "meta"//".dat"
               end if

               call util_file_name_noprefix(filename,.false.,
     >                                      .false.,
     >                                      full_filename)
               open(unit=53,file=full_filename,form='formatted')

               write(53,'(A,A)') "#metadynamics filename=",filename
               write(53,'(A,I6)') "#nmeta ",nmeta
               write(53,'(A,10I6)') "#nxmeta ",
     >                              (int_mb(nxmeta(1)+i-1),i=1,nmeta) 
               write(53,'(A,10I6)') "#sindxmeta ",
     >                     (int_mb(sindxmeta(1)+i-1),i=1,nmeta) 
               write(53,'(A,10I6)') "#sparammeta ",
     >                     (int_mb(sparammeta(1)+i-1),i=1,nmeta) 
               write(53,'(A,10F18.12)') "#ameta ",
     >                                 (dbl_mb(ameta(1)+i-1),i=1,nmeta)
               write(53,'(A,10F18.12)') "#bmeta ",
     >                                (dbl_mb(bmeta(1)+i-1),i=1,nmeta)
               write(53,'(A,10I3)') "#pmeta ",
     >                              (int_mb(pmeta(1)+i-1),i=1,nmeta)
               write(53,'(A,I8)') "#nindxmeta ",nindxmeta
               write(53,'(A,500I6)') "#indxmeta ",
     >                     (int_mb(indxmeta(1)+i-1),i=1,nindxmeta) 
               write(53,'(A,I8)') "#nparammeta ",nparammeta
               write(53,'(A,500F18.12)') "#parammeta ",
     >                     (dbl_mb(parammeta(1)+i-1),i=1,nparammeta)
               write(53,'(A,I8)') "#maxmetacount ",maxmetacount
               write(53,'(A,I8)') "#metacount ",metacount
               write(53,'(A,I8)') "#metarayshift ",metarayshift
               write(53,'(A,I8)') "#metaraycount ",metaraycount
               write(53,'(A,I8)')"#maxmetaprintcount ",maxmetaprintcount
               write(53,'(A,F18.12)') "#dTempered ",dTtempered
               write(53,'(A,I8)') "#nxmeta_all ",nxmeta_all

               do d=1,nmeta
                  j(d) = 0
               end do
               do i=1,nxmeta_all
                  do d=1,nmeta
                     x(d) = dbl_mb(ameta(1)+d-1)
     >                + j(d)*(dbl_mb(bmeta(1)+d-1)-dbl_mb(ameta(1)+d-1))
     >                      /dble(int_mb(nxmeta(1)+d-1)-1)
                  end do
                  write(53,'(12F15.6)') 
     >                 (x(d),d=1,nmeta),dbl_mb(ymeta(1)+i-1)*fac
                  j(1) = j(1) + 1
                  if (j(1).ge.int_mb(nxmeta(1))) then
                     write(53,*)
                  end if
                  do d=1,nmeta-1
                     if (j(d).ge.int_mb(nxmeta(1)+d-1)) then
                        j(d) = 0
                        j(d+1) = j(d+1)+1
                     end if
                  end do
               end do

               close(53)
               

c               filename = "test.dat"
c               call util_file_name_noprefix(filename,.false.,
c     >                                      .false.,
c     >                                      full_filename)
c               open(unit=53,file=full_filename,form='formatted')
c               do i=1,501
c                  x(1) = 1.0d0 + (i-1)*8.0d0/501.0d0
c                
c                  v = nwpw_interp(nmeta,int_mb(nxmeta(1)),
c     >                      dbl_mb(ameta(1)),dbl_mb(bmeta(1)),
c     >                      dbl_mb(ymeta(1)),x)
c                  call nwpw_dinterp(nmeta,int_mb(nxmeta(1)),
c     >                        dbl_mb(ameta(1)),dbl_mb(bmeta(1)),
c     >                        dbl_mb(ymeta(1)),x,dv)
c
c                  write(53,'(12F15.6)') x(1),v,dv(1)
c
c               end do
c               close(53)

            end if
         end if

      end if

      return
      end



*     ***********************************************
*     *                                             *
*     *              meta_setboundary               *
*     *                                             *
*     ***********************************************
      subroutine meta_setboundary(bheight,bsigma)
      implicit none
      real*8 bheight,bsigma

#include "bafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer d,i,j(10)
      real*8  r0,r,x,w,ww,wa,wb,ra,rb

      if (metafound) then
         if (nmeta.gt.0) then
            do d=1,nmeta
               j(d) = 0
            end do
            do i=1,nxmeta_all
               ww = 0.0d0
               do d=1,nmeta
                  if (int_mb(pmeta(1)+d-1).eq.0) then
                     x = dbl_mb(ameta(1)+d-1) 
     >                + j(d)*(dbl_mb(bmeta(1)+d-1)-dbl_mb(ameta(1)+d-1))
     >                      /dble(int_mb(nxmeta(1)+d-1)-1)
                  else
                     x = dbl_mb(ameta(1)+d-1) 
     >                + j(d)*(dbl_mb(bmeta(1)+d-1)-dbl_mb(ameta(1)+d-1))
     >                      /dble(int_mb(nxmeta(1)+d-1))
                  end if
                  ra= dbl_mb(ameta(1)+d-1)
                  rb= dbl_mb(bmeta(1)+d-1)
                  wa = dexp(-0.5d0*((x-ra)/bsigma)**2)
                  wb = dexp(-0.5d0*((x-rb)/bsigma)**2)
                  ww = ww + wa + wb
               end do

               !*** make sure the maximum barrier height is bheight ***
               if (ww.gt.1.0d0) ww = 1.0d0 
            
               dbl_mb(ymeta(1)+i-1) = dbl_mb(ymeta(1)+i-1) + bheight*ww

               j(1) = j(1) + 1
               do d=1,nmeta-1
                  if (j(d).ge.int_mb(nxmeta(1)+d-1)) then
                     j(d) = 0
                     j(d+1) = j(d+1)+1
                  end if
               end do
            end do
         end if
      end if
      return
      end



*     ***********************************************
*     *                                             *
*     *          meta_gauss_write_init              *
*     *                                             *
*     ***********************************************
      subroutine meta_gauss_write_init(filename)
      implicit none
      character*(*) filename

#include "bafdecls.fh"
#include "stdio.fh"
#include "metadynamics.fh"

      integer   MASTER
      parameter (MASTER=0)

      logical found,found_bak
      integer taskid,l1,l2,i
      character*255 full_filename,full_bak

      call Parallel_taskid(taskid)


*     **** produce META-GAUSS FILE ****
      if (taskid.eq.MASTER) then

         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)

*        **** check for backup file ***
         call util_file_name_noprefix('META-GAUSS99-bak',.false.,
     >                                  .false.,
     >                                  full_bak)
         inquire(file=full_bak,exist=found_bak)
         if (found_bak) then
            write(luout,*)
            write(luout,*) "META-GAUSS99-bak exists:"
            l1=index(full_bak,' ')
            l2=index(full_filename,' ')
            write(luout,*) "   Copying ",full_bak(1:l2),
     >                 " to ",full_filename(1:l2)
            write(*,*)
            call util_file_copy(full_bak,full_filename)
         end if
         inquire(file=full_filename,exist=found)
         if (found) then

*           **** make a new backup file ***
            call util_file_copy(full_filename,full_bak)

           open(unit=54,file=full_filename,form='formatted',
     >          status='old')
           do while (found)
             read(54,*,end=100)
           end do
  100      continue
#if defined(FUJITSU) || defined(SOLARIS) || defined(__crayx1) || defined(GCC46)
           backspace 54
#endif
         else
           open(unit=54,file=full_filename,form='formatted',
     >           status='new')

           write(54,'(A,A40)') 
     >        "#metadynamics gaussian filename=",filename
           write(54,'(A,I6)') "#nmeta ",nmeta
           write(54,'(A,10I6)') "#nxmeta ",
     >                          (int_mb(nxmeta(1)+i-1),i=1,nmeta) 
           write(54,'(A,10I6)') "#sindxmeta ",
     >                 (int_mb(sindxmeta(1)+i-1),i=1,nmeta) 
           write(54,'(A,10I6)') "#sparammeta ",
     >                 (int_mb(sparammeta(1)+i-1),i=1,nmeta) 
           write(54,'(A,10F18.12)') "#ameta ",
     >                             (dbl_mb(ameta(1)+i-1),i=1,nmeta)
           write(54,'(A,10F18.12)') "#bmeta ",
     >                            (dbl_mb(bmeta(1)+i-1),i=1,nmeta)
           write(54,'(A,10I3)') "#pmeta ",
     >                          (int_mb(pmeta(1)+i-1),i=1,nmeta)
           write(54,'(A,I8)') "#nindxmeta ",nindxmeta
           write(54,'(A,500I6)') "#indxmeta ",
     >                 (int_mb(indxmeta(1)+i-1),i=1,nindxmeta) 
           write(54,'(A,I8)') "#nparammeta ",nparammeta
           write(54,'(A,500F18.12)') "#parammeta ",
     >                 (dbl_mb(parammeta(1)+i-1),i=1,nparammeta)
           write(54,'(A,I8)') "#maxmetacount ",maxmetacount
           write(54,'(A,I8)') "#metacount ",metacount
           write(54,'(A,I8)') "#metarayshift ",metarayshift
           write(54,'(A,I8)') "#metaraycount ",metaraycount
           write(54,'(A,I8)') "#maxmetaprintcount ",maxmetaprintcount
           write(54,'(2A)') "# (r(i),i=1,n),(sigma(i),i=1,n),w,",
     >                     "E(2)-E(31)+E(32),E(3),E(4),E(31),E(32)"
         end if
      end if

      return
      end




*     ***********************************************
*     *                                             *
*     *          meta_gauss_write_end               *
*     *                                             *
*     ***********************************************
      subroutine meta_gauss_write_end()
      implicit none

      integer   MASTER
      parameter (MASTER=0)

      integer taskid
      character*255 full_bak

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         close(unit=54)

*        **** remove backup file ***
         call util_file_name_noprefix('META-GAUSS99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         call util_file_unlink(full_bak)
      end if


      return
      end

*     ***********************************************
*     *                                             *
*     *          meta_gauss_write                   *
*     *                                             *
*     ***********************************************
      subroutine meta_gauss_write(n,r,sigma,w,E)
      implicit none
      integer n
      real*8 r(*),sigma(*),w,E(*)

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      integer i

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         write(54,111) (r(i),i=1,n),(sigma(i),i=1,n),w,
     >                 E(2)-E(31)+E(32),E(3),E(4),E(31),E(32)
         call util_flush(54)
      end if
  111 format(99e14.6)

      return
      end



*     ***********************************************
*     *                                             *
*     *              meta_setsn2surface             *
*     *                                             *
*     ***********************************************
      subroutine meta_setsn2surface(sn2)
      implicit none
      real*8 sn2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      logical sn2found,blaisbunker
      integer d1,d2,i,j(10),d
      real*8  xp,xr,ww
      real*8 er,ep,erp,betar,betap,betarp,xr0,xp0,xrp0
      real*8 wr,ar,br,wp,ap,bp,fr,fp,frp,gr,gp

      sn2found = dabs(sn2(1)).gt.1.0d-9
      if (sn2found) then

      blaisbunker = (sn2(1).lt.0.0d0)
    
      if (blaisbunker) then
      er  = sn2(2)
      ep  = sn2(3)
      erp = sn2(4)
      betar  = sn2(5)
      betap  = sn2(6)
      betarp = sn2(7)
      xr0  = sn2(8)
      xp0  = sn2(9)
      xrp0 = sn2(9)
      ap = sn2(10)
      bp = sn2(11)

      else
      xr0 = sn2(2)
      wr  = sn2(3)
      ar  = sn2(4)
      br  = sn2(5)

      xp0 = sn2(6)
      wp  = sn2(7)
      ap  = sn2(8)
      bp  = sn2(9)

      er  = sn2(10)
      ep  = sn2(11)
      erp = sn2(12)
      end if

      if (metafound) then
         if (nmeta.gt.0) then
            do d=1,nmeta
               j(d) = 0
            end do
            do i=1,nxmeta_all
               ww = 0.0d0
               d1 = 1
               d2 = 2
               xr = dbl_mb(ameta(1)+d1-1) 
     >           + j(d1)*(dbl_mb(bmeta(1)+d1-1)-dbl_mb(ameta(1)+d1-1))
     >                 /dble(int_mb(nxmeta(1)+d1-1)-1)
               xp = dbl_mb(ameta(1)+d2-1) 
     >           + j(d2)*(dbl_mb(bmeta(1)+d2-1)-dbl_mb(ameta(1)+d2-1))
     >                 /dble(int_mb(nxmeta(1)+d2-1)-1)

               if (blaisbunker) then
                  fr = er - er*(1.0d0-dexp(-betar*(xr-xr0)))**2
                  fp = ep - ep*(1.0d0-dexp(-betap*(xp-xp0)))**2
     >               + (1.0d0-dtanh(ap*xr+bp))
     >                *(ep - ep*dexp(-betap*(xp-xp0)))
                  frp = erp - erp*(1.0d0-dexp(-betarp*(xr+xp-xrp0)))**2
                  ww = fr+fp+frp
               else
                  fr = er + (erp-er)*(1.0d0+dtanh( (xp0-xp)/ar))
                  gr = dexp(-((xr-xr0)/wr)**2)
     >                *(0.5d0+0.5d0*dtanh((xp-xp0)/br))

                  fp = ep + (erp-ep)*(1.0d0+dtanh( (xr0-xr)/ap))
                  gp = dexp(-((xp-xp0)/wp)**2)
     >                 *(0.5d0+0.5d0*dtanh((xr-xr0)/bp))
                  ww = fr*gr + fp*gp
               end if

               if (ww.gt.0.0d0) ww = 0.0d0
               dbl_mb(ymeta(1)+i-1) = dbl_mb(ymeta(1)+i-1) - ww

               j(1) = j(1) + 1
               do d=1,nmeta-1
                  if (j(d).ge.int_mb(nxmeta(1)+d-1)) then
                     j(d) = 0
                     j(d+1) = j(d+1)+1
                  end if
               end do
            end do
         end if
      end if
      end if
      return
      end

