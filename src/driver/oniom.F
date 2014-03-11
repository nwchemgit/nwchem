      logical function oniom_energy(rtdb)
*
* $Id$
*
      implicit none
      integer rtdb
      logical oniom
      external oniom
      oniom_energy = oniom(rtdb,'energy')
      end
      logical function oniom_gradient(rtdb)
      implicit none
      integer rtdb
      logical oniom
      external oniom
      oniom_gradient = oniom(rtdb,'gradient')
      end
      logical function oniom(rtdb, task)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "geom.fh"
#include "inp.fh"
#include "nwc_const.fh"
      integer rtdb
      character*(*) task
c     
c     Implements a two- and three-layer ONIOM model
c     
c     E(High,Real) = E(High,Model) + E(Low,Real) - E(Low,Model)
c     
c     task can currently be energy or gradient
c     
      logical oprint, ignore, status, oprint_head
      character*256 title
      character*32 high_theory, low_theory, medium_theory, 
     $     gname, bname, ename
      character*32 high_basis, low_basis, medium_basis
      character*32 high_ecp, low_ecp, medium_ecp
c
      character*256 v_low_real, v_low_model, v_high_model,
     $     v_medium_model, v_medium_inter, v_low_inter
c
      integer real_nat, model_nat, nat, nattest, inter_nat
      integer geom
      integer i, j, k, bond, nlayer, dummy(2)
      integer maxbonds
      parameter (maxbonds = 50)
      integer model_bonds(2,maxbonds), model_nbonds
      double precision model_p(maxbonds), model_charge
      character*32 model_tags(maxbonds)
      integer inter_bonds(2,maxbonds), inter_nbonds
      double precision inter_p(maxbonds), inter_charge
      character*32 inter_tags(maxbonds)
      double precision p, real_charge
      character*1024 high_input, medium_input, low_input
c     
c     For manipulating the geometry
c     
      double precision coords(3,nw_max_atom), charge(nw_max_atom)
      character*16 tags(nw_max_atom)
c     
      double precision g(3,nw_max_atom), g_low_model(3,nw_max_atom),
     $     g_high_model(3,nw_max_atom), g_low_real(3,nw_max_atom),
     $     g_medium_model(3,nw_max_atom), g_medium_inter(3,nw_max_atom),
     $     g_low_inter(3,nw_max_atom)
c     
      double precision energy, e_low_real, e_low_model, e_high_model,
     $     e_medium_model, e_medium_inter, e_low_inter, e2, e3
      double precision gdummy
c
      character*2 symbol
      character*8 element
      integer atn
c
      integer k_actlist, l_actlist, ma_type, nactive
c     
      logical task_energy_doit, task_gradient_doit, geom_zmtmak
      logical usesymmol
      integer nops
      external task_energy_doit, task_gradient_doit, geom_zmtmak
c     
      status = rtdb_parallel(.true.)
         if(.not.rtdb_get(rtdb,'geom:symmol',mt_log,1,
     $        usesymmol)) usesymmol=.false.
      call util_print_push
      call util_print_rtdb_load(rtdb, 'oniom')
      call ecce_print_module_entry('oniom')
      oprint = util_print('information', print_low) .and. 
     $     ga_nodeid().eq.0
      oprint_head = util_print('headers', print_low) .and. 
     $     ga_nodeid().eq.0
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom,'geometry')) 
     $     call errquit('oniom: failed creating geometry?',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('oniom: no geometry ', 0, RTDB_ERR)
      if (.not. geom_rtdb_store(rtdb, geom, 'oniom_geometry'))
     $     call errquit('oniom: failed storing geometry',0, RTDB_ERR)
      if (.not. geom_ncent(geom,real_nat))
     $     call errquit('oniom: number of atoms from geometry?',0,
     &       GEOM_ERR)
      if (.not. geom_destroy(geom))
     $     call errquit('oniom: destroying geometry',0, GEOM_ERR)
c
c     Save the list active list of atoms for gradient calculations
c
      if (task .eq. 'gradient') then
         nactive = real_nat
         if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type, 
     $        nactive, l_actlist)) then
            if (.not. ma_get_index(l_actlist, k_actlist))
     $           call errquit('oniom: ma_get_index failed',l_actlist,
     &       MA_ERR)
            if (.not. rtdb_put(rtdb, 'oniom:actlist', mt_int, nactive, 
     $           int_mb(k_actlist))) call errquit
     $           ('oniom: saving actlist failed',0, RTDB_ERR)
            if (.not. ma_free_heap(l_actlist)) call errquit
     $           ('oniom: free of actlist?',0, MA_ERR)
         end if
      else
         ignore = rtdb_delete(rtdb,'oniom:actlist')
      end if
c     
c     Get ONIOM parameters from the database
c
      if (.not. rtdb_cget(rtdb, 'oniom:high:theory', 1, high_theory))
     $     call errquit('oniom: high_theory?', 0, RTDB_ERR)
      if (.not. rtdb_cget(rtdb, 'oniom:medium:theory', 1,medium_theory))
     $     medium_theory = ' '
      if (.not. rtdb_cget(rtdb, 'oniom:low:theory', 1, low_theory))
     $     call errquit('oniom: low_theory?', 0, RTDB_ERR)
c     
      if (.not. rtdb_cget(rtdb, 'oniom:high:basis', 1, high_basis))
     $     high_basis = 'ao basis'
      if (.not. rtdb_cget(rtdb, 'oniom:medium:basis', 1, medium_basis))
     $     medium_basis = high_basis
      if (.not. rtdb_cget(rtdb, 'oniom:low:basis', 1, low_basis))
     $     low_basis = medium_basis
c
      if (.not. rtdb_cget(rtdb, 'oniom:high:ecp', 1, high_ecp))
     $     high_ecp = ' '
      if (.not. rtdb_cget(rtdb, 'oniom:medium:ecp', 1, medium_ecp))
     $     medium_ecp = ' ' 
      if (.not. rtdb_cget(rtdb, 'oniom:low:ecp', 1, low_ecp))
     $     low_ecp = ' '
c
      if (.not. rtdb_cget(rtdb, 'oniom:high:input', 1, high_input))
     $     high_input = ' '
      if (.not. rtdb_cget(rtdb, 'oniom:medium:input', 1, medium_input))
     $     medium_input = ' ' 
      if (.not. rtdb_cget(rtdb, 'oniom:low:input', 1, low_input))
     $     low_input = ' '
c
      if (.not. rtdb_cget(rtdb, 'oniom:v_low_real', 1, 
     $     v_low_real)) 
     $     call util_file_name('lrmos', .false.,.false., v_low_real)
      if (.not. rtdb_cget(rtdb, 'oniom:v_low_model', 1, 
     $     v_low_model)) 
     $     call util_file_name('lmmos', .false.,.false., v_low_model)
      if (.not. rtdb_cget(rtdb, 'oniom:v_high_model', 1, 
     $     v_high_model)) 
     $     call util_file_name('hmmos', .false.,.false., v_high_model)
      if (.not. rtdb_cget(rtdb, 'oniom:v_medium_model', 1, 
     $     v_medium_model)) 
     $     call util_file_name('mmmos', .false.,.false., v_medium_model)
      if (.not. rtdb_cget(rtdb, 'oniom:v_medium_inter', 1, 
     $     v_medium_inter)) 
     $     call util_file_name('mimos', .false.,.false., v_medium_inter)
      if (.not. rtdb_cget(rtdb, 'oniom:v_low_inter', 1, 
     $     v_low_inter)) 
     $     call util_file_name('limos', .false.,.false., v_low_inter)
c
      call util_file_name_resolve(v_low_real, .false.)
      call util_file_name_resolve(v_low_model, .false.)
      call util_file_name_resolve(v_high_model, .false.)
      call util_file_name_resolve(v_medium_model, .false.)
      call util_file_name_resolve(v_medium_inter, .false.)
      call util_file_name_resolve(v_low_inter, .false.)
c     
      if (.not. rtdb_get(rtdb, 'oniom:model:nat', mt_int, 1,
     $     model_nat)) call errquit('oniom: model_nat?',0, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'oniom:inter:nat', mt_int, 1,
     $     inter_nat)) inter_nat = 0
      if (.not. rtdb_get(rtdb, 'oniom:model:charge', mt_dbl, 1, 
     $        model_charge)) call errquit('oniom:model charge',0,
     &       RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'oniom:inter:charge', mt_dbl, 1, 
     $        inter_charge)) inter_charge = 0.0d0
c     
      model_nbonds = 0
      if (rtdb_get(rtdb,'oniom:model:nbonds',mt_int,1,model_nbonds))then
         if (.not. rtdb_get(rtdb, 'oniom:model:bonds', mt_int,
     $        2*model_nbonds, model_bonds)) call errquit
     $        ('oniom: failed to read model-real bonds?',0,
     &       RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'oniom:model:p', mt_dbl,
     $        model_nbonds, model_p)) call errquit
     $        ('oniom: failed to read model-real bonds factors?',0,
     &       RTDB_ERR)
         if (.not. rtdb_cget(rtdb,'oniom:model:tags', 
     $           model_nbonds, model_tags)) call errquit
     $        ('oniom: failed to read model-real bonds tags?',0,
     &       RTDB_ERR)
      end if
c     
      inter_nbonds = 0
      if (rtdb_get(rtdb,'oniom:inter:nbonds',mt_int,1,inter_nbonds))then
         if (.not. rtdb_get(rtdb, 'oniom:inter:bonds', mt_int,
     $        2*inter_nbonds, inter_bonds)) call errquit
     $        ('oniom: failed to read inter-real bonds?',0,
     &       RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'oniom:inter:p', mt_dbl,
     $        inter_nbonds, inter_p)) call errquit
     $        ('oniom: failed to read inter-real bonds factors?',0,
     &       RTDB_ERR)
         if (.not. rtdb_cget(rtdb,'oniom:inter:tags', 
     $           inter_nbonds, inter_tags)) call errquit
     $        ('oniom: failed to read model-real bonds tags?',0,
     &       RTDB_ERR)
      end if
c     
c     Save user indirections for the basis set and geometry
c     
      if (.not. rtdb_cget(rtdb,'geometry',1,gname)) gname = ' '
      if (.not. rtdb_cget(rtdb,'ao basis',1,bname)) bname = ' '
      if (.not. rtdb_cget(rtdb,'ecp basis',1,ename)) ename = ' '
      if (.not. rtdb_get(rtdb, 'charge', mt_dbl, 1, real_charge))
     $     real_charge = 0.0d0
c
      if (model_charge .eq. -999.0) model_charge = real_charge
      if (inter_charge .eq. -999.0) inter_charge = real_charge
c
      nlayer = 2
      if (inter_nat .gt. 0) nlayer = 3
c     
c     Summarize parameters
c
      if (oprint) then
         call util_print_centered(6, 'NWChem ONIOM Module', 40, .true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         end if
c     
         write(6,10)
 10      format(7x,'Theory',/,7x,'------')
         if (high_ecp .eq. ' ') then
            write(6,1) '   high', 
     $           high_theory(1:inp_strlen(high_theory)), 
     $           high_basis(1:inp_strlen(high_basis))
         else
            write(6,1) '   high', 
     $           high_theory(1:inp_strlen(high_theory)), 
     $           high_basis(1:inp_strlen(high_basis)),
     $           high_ecp(1:inp_strlen(high_ecp))
         end if
         if (high_input .ne. ' ') 
     $        write(6,101) high_input(1:inp_strlen(high_input))
         if (medium_theory .ne. ' ') then
            if (medium_ecp .eq. ' ') then
               write(6,1) ' medium', 
     $              medium_theory(1:inp_strlen(medium_theory)), 
     $              medium_basis(1:inp_strlen(medium_basis))
            else
               write(6,1) ' medium', 
     $              medium_theory(1:inp_strlen(medium_theory)), 
     $              medium_basis(1:inp_strlen(medium_basis)),
     $              medium_ecp(1:inp_strlen(medium_ecp))
            end if
            if (medium_input .ne. ' ') 
     $           write(6,101) medium_input(1:inp_strlen(medium_input))
         end if
         if (low_ecp .eq. ' ') then
            write(6,1) '    low', low_theory(1:inp_strlen(low_theory)), 
     $           low_basis(1:inp_strlen(low_basis))
         else
            write(6,1) '    low', low_theory(1:inp_strlen(low_theory)), 
     $           low_basis(1:inp_strlen(low_basis)),
     $           low_ecp(1:inp_strlen(low_ecp))
         end if
         if (low_input .ne. ' ') 
     $        write(6,101) low_input(1:inp_strlen(low_input))
 1       format(6x,a7,2x,'theory=',a,2x,'basis="',a,'"',:,2x,
     $        'ecp="',a,'"')
 101     format(15x,'input string="',a,'"')
c
         write(6,*)
         write(6,15)
 15      format(6x,'Vectors',/,6x,'-------')
         if (nlayer .eq. 2) then
            write(6,21) 'High-model', 
     $           v_high_model(1:inp_strlen(v_high_model))
            write(6,21) 'Low-real', 
     $           v_low_real(1:inp_strlen(v_low_real))
            write(6,21) 'Low-model', 
     $           v_low_model(1:inp_strlen(v_low_model))
         else
            write(6,21) 'High-model', 
     $           v_high_model(1:inp_strlen(v_high_model))
            write(6,21) 'Medium-model', 
     $           v_medium_model(1:inp_strlen(v_medium_model))
            write(6,21) 'Medium-inter', 
     $           v_medium_inter(1:inp_strlen(v_medium_inter))
            write(6,21) 'Low-real', 
     $           v_low_real(1:inp_strlen(v_low_real))
            write(6,21) 'Low-inter', 
     $           v_low_inter(1:inp_strlen(v_low_inter))
         end if
 21      format(8x,a12,2x,a)
c     
         write(6,*)
         write(6,20)
 20      format(5x,'Geometry',/,5x,'--------')
         write(6,2) '   real', real_nat, real_charge
         if (inter_nat .gt. 0) write(6,2) '  inter', inter_nat, 
     $        inter_charge
         if (inter_nbonds .gt. 0) then
            write(6,*) '                 ',
     $           'Bonds between intermediate and real'
            write(6,3) (inter_bonds(1,i),inter_bonds(2,i),inter_p(i),
     $           inter_tags(i)(1:inp_strlen(inter_tags(i))),
     $           i=1,inter_nbonds)
         end if
         write(6,2) '  model', model_nat, model_charge
         if (model_nbonds .gt. 0) then
            write(6,*) '            Bonds between model and',
     $           ' real/intermediate'
            write(6,3) (model_bonds(1,i),model_bonds(2,i),model_p(i),
     $           model_tags(i)(1:inp_strlen(model_tags(i))), 
     $           i=1,model_nbonds)
         end if
 2       format(6x,a7,2x,i4,' atoms with charge',f8.4)
 3       format(20x,2i4,f6.2,2x,'"',a,'"')
         call util_flush(6)
      end if
c     
c     Make the model geometry resetting the symmetry & autoz info
c     
      if (oprint) then
         write(6,*)
         write(6,*) ' Making model geometry '
         write(6,*)
      end if
c     
      if (.not. geom_create(geom,'geometry')) 
     $     call errquit('oniom: failed creating geometry?',0,
     &       GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'oniom_geometry'))
     $     call errquit('oniom: no geometry ', 0,
     &       RTDB_ERR)
      if (.not. geom_cart_get(geom,real_nat,tags,g,charge))
     $     call errquit('oniom: geom_cart_get failed', 0, GEOM_ERR)
      call dcopy(3*real_nat, g, 1, coords, 1)
c     
      nat = model_nat
      do bond = 1, model_nbonds
         i = model_bonds(1,bond)
         j = model_bonds(2,bond)
         p = model_p(bond)
         nat = nat + 1
         tags(nat) = model_tags(bond) !!!!!!!!
         ignore = geom_tag_to_element(tags(nat), symbol, element, atn)
         charge(nat) = atn
         do k = 1, 3
            coords(k,nat) = g(k,i)*(1.0d0-p) + g(k,j)*p
         end do
      end do
c     
      if (.not. geom_cart_set(geom,nat,tags,coords,charge))
     $     call errquit('oniom: geom_cart_get failed', 0, GEOM_ERR)
c     
      nattest = nat             ! sym_nwc corrupts the input #atoms
      call sym_nwc(geom, rtdb, nattest, .false., 1.0d0, 1d-4,
     ,     nops, usesymmol)
      if (.not. geom_ncent(geom,nattest)) call errquit
     $     ('oniom: getting ncent',0, GEOM_ERR)
      if (nat .ne. nattest) call errquit
     $     ('oniom: model does not have correct symmetry', 0, GEOM_ERR)
      ignore = geom_zmtmak(rtdb, geom, .false.)
      if (.not. geom_rtdb_store(rtdb, geom, 'oniom model'))
     $     call errquit('oniom: failed storing geometry',0, RTDB_ERR)
      if (ga_nodeid().eq.0 .and.
     $     util_print('geometry',print_default)) then
         if (.not. geom_print(geom)) 
     $        call errquit('oniom: geom print?',0, GEOM_ERR)
      end if
      if (.not. geom_destroy(geom))
     $     call errquit('oniom: destroying geometry',0, GEOM_ERR)
c     
c     Make the intermediate geometry resetting the symmetry & autoz info
c     
      if (inter_nat .gt. 0) then
         if (oprint) then
            write(6,*)
            write(6,*) ' Making intermediate geometry '
            write(6,*)
         end if
         if (.not. geom_create(geom,'geometry')) 
     $        call errquit('oniom: failed creating geometry?',0,
     &       GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geom, 'oniom_geometry'))
     $        call errquit('oniom: no geometry ', 0, RTDB_ERR)
         if (.not. geom_cart_get(geom,real_nat,tags,g,charge))
     $        call errquit('oniom: geom_cart_get failed', 0, GEOM_ERR)
         call dcopy(3*real_nat, g, 1, coords, 1)
         nat = inter_nat
         do bond = 1, inter_nbonds
            i = inter_bonds(1,bond)
            j = inter_bonds(2,bond)
            p = inter_p(bond)
            nat = nat + 1
            tags(nat) = inter_tags(bond) !!!! 'H oniom'
            charge(nat) = 1.0d0
            do k = 1, 3
               coords(k,nat) = g(k,i)*(1.0d0-p) + g(k,j)*p
            end do
         end do
         if (.not. geom_cart_set(geom,nat,tags,coords,charge))
     $        call errquit('oniom: geom_cart_get failed', 0, GEOM_ERR)
         nattest = nat          ! sym_nwc corrupts the input #atoms
         call sym_nwc(geom, rtdb, nattest, .false., 1.0d0, 1d-4,
     ,        nops, usesymmol)
         if (.not. geom_ncent(geom,nattest)) call errquit
     $        ('oniom: getting ncent',0, GEOM_ERR)
         if (nat .ne. nattest) call errquit
     $        ('oniom: inter does not have correct symmetry', 0,
     &       GEOM_ERR)
         ignore = geom_zmtmak(rtdb, geom, .false.)
         if (.not. geom_rtdb_store(rtdb, geom, 'oniom inter'))
     $        call errquit('oniom: failed storing geometry',0, RTDB_ERR)
         if (ga_nodeid().eq.0 .and.
     $        util_print('geometry',print_default)) then
            if (.not. geom_print(geom)) 
     $           call errquit('oniom: geom print?',0, GEOM_ERR)
         end if
         if (.not. geom_destroy(geom))
     $        call errquit('oniom: destroying geometry',0,
     &       GEOM_ERR)
      end if
c     
c     Two-layer model
c     E(High,Real) = E(High,Model) + E(Low,Real) - E(Low,Model)
c     
c     Three-layer model
c     E(High,Real) = E(High,Model) - E(Medium,Model)
c     .            + E(Medium,Inter) + E(Real,Low) - E(Low,Inter) 
c     
c     To avoid recursive calls thru Fortran, call the task*doit layer.
c     
      oniom = .false.
      call dfill(3*nw_max_atom, 0.0d0, g, 1)
      call dfill(3*nw_max_atom, 0.0d0, g_low_model, 1)
      call dfill(3*nw_max_atom, 0.0d0, g_high_model, 1)
      call dfill(3*nw_max_atom, 0.0d0, g_low_real, 1)
      call dfill(3*nw_max_atom, 0.0d0, g_medium_model, 1)
      call dfill(3*nw_max_atom, 0.0d0, g_medium_inter, 1)
      energy = 0.0d0
      e_low_model = 0.0d0
      e_high_model = 0.0d0
      e_low_real = 0.0d0
      e_medium_model = 0.0d0
      e_medium_inter = 0.0d0
c     
c     Low + Real
c     
      if (oprint_head)
     $     call util_print_centered(6, 'ONIOM LOW+REAL', 20, .true.)
      call oniom_param(rtdb, ' ', low_basis, low_ecp, real_charge,
     $     low_input)
      call oniom_manage_mos(rtdb,low_theory,v_low_real)
      if (task .eq. 'energy') then
         status = task_energy_doit(rtdb, low_theory, e_low_real)
      else
         call oniom_munge_actlist(rtdb,real_nat,0,dummy)
         status = task_gradient_doit(rtdb, low_theory, e_low_real, 
     $        g_low_real)
      end if
      if (.not. status) 
     $     call errquit('oniom: low theory + real geometry failed',0,
     &       GEOM_ERR)
c     
c     Low + Model
c     
      if (oprint_head)
     $     call util_print_centered(6, 'ONIOM LOW+MODEL', 20, .true.)
      call oniom_param(rtdb, 'oniom model',low_basis, low_ecp,
     $     model_charge, low_input)
      call oniom_manage_mos(rtdb,low_theory,v_low_model)
      if (task .eq. 'energy') then
         status = task_energy_doit(rtdb, low_theory, e_low_model)
      else
         call oniom_munge_actlist(rtdb,model_nat,
     $        model_nbonds,model_bonds)
         status = task_gradient_doit(rtdb, low_theory, e_low_model, 
     $        g_low_model)
      end if
      if (.not. status) 
     $     call errquit('oniom:low theory + model geometry failed',0,
     &       GEOM_ERR)
c     
c     High + Model
c     
      if (oprint_head)
     $     call util_print_centered(6, 'ONIOM HIGH+MODEL', 20, .true.)
      call oniom_param(rtdb, 'oniom model', 
     $     high_basis, high_ecp, model_charge, high_input)
      call oniom_manage_mos(rtdb,high_theory, v_high_model)
      if (task .eq. 'energy') then
         status = task_energy_doit(rtdb, high_theory, e_high_model)
      else
         call oniom_munge_actlist(rtdb,model_nat,
     $        model_nbonds,model_bonds)
         status = task_gradient_doit(rtdb, high_theory, e_high_model, 
     $        g_high_model)
      end if
      if (.not. status) 
     $     call errquit('oniom: high theory + model geometry failed',0,
     &       GEOM_ERR)
c
      if (nlayer .eq. 3) then
c     
c     Low + Inter
c     
         if (oprint_head)
     $        call util_print_centered(6, 'ONIOM LOW+INTER', 20, .true.)
         call oniom_param(rtdb, 'oniom inter', 
     $        low_basis, low_ecp, inter_charge, low_input)
         call oniom_manage_mos(rtdb, low_theory, v_low_inter)
         if (task .eq. 'energy') then
            status = task_energy_doit(rtdb, low_theory, e_low_inter)
         else
            call oniom_munge_actlist(rtdb,inter_nat,
     $           inter_nbonds,inter_bonds)
            status = task_gradient_doit(rtdb, low_theory, e_low_inter, 
     $           g_low_inter)
         end if
         if (.not. status) 
     $        call errquit('oniom:low theory + model geometry failed',0,
     &       GEOM_ERR)
c     
c     Medium + Inter
c     
         if (oprint_head) call util_print_centered
     $        (6, 'ONIOM MEDIUM+INTER', 20, .true.)
         call oniom_param(rtdb, 'oniom inter', 
     $        medium_basis, medium_ecp, inter_charge, medium_input)
         call oniom_manage_mos(rtdb,medium_theory,v_medium_inter)
         if (task .eq. 'energy') then
            status = task_energy_doit(rtdb,medium_theory,e_medium_inter)
         else
            call oniom_munge_actlist(rtdb,inter_nat,
     $           inter_nbonds,inter_bonds)
            status = task_gradient_doit(rtdb, medium_theory, 
     $           e_medium_inter, g_medium_inter)
         end if
         if (.not. status) 
     $        call errquit('oniom: medium theory + '//
     $        'inter geometry failed',0, GEOM_ERR)
c     
c     Medium + Model
c     
         if (oprint_head) call util_print_centered
     $        (6, 'ONIOM MEDIUM+MODEL', 20, .true.)
         call oniom_param(rtdb, 'oniom model',
     $        medium_basis, medium_ecp, model_charge, medium_input)
         call oniom_manage_mos(rtdb,medium_theory,v_medium_model)
         if (task .eq. 'energy') then
            status = task_energy_doit(rtdb,medium_theory,e_medium_model)
         else
            call oniom_munge_actlist(rtdb,model_nat,
     $           model_nbonds,model_bonds)
            status = task_gradient_doit(rtdb, medium_theory, 
     $           e_medium_model, g_medium_model)
         end if
         if (.not. status) 
     $        call errquit('oniom: medium theory + '//
     $        'model geometry failed',0, GEOM_ERR)
      end if
c     
c     Assemble the ONIOM energy and gradient
c     
      e2 = e_high_model + e_low_real - e_low_model
      energy = e2
      if (nlayer .eq. 3) then
         e3 = e_high_model - e_medium_model
     $        + e_medium_inter + e_low_real - e_low_inter 
         energy = e3
      end if
c     
      if (.not. rtdb_put(rtdb,'oniom:energy',mt_dbl,1,energy))
     $     call errquit('oniom: failed storing energy', 0, RTDB_ERR)
c     
      if (task .eq. 'gradient') then
         if (nlayer .eq. 2) then
            call dcopy(3*real_nat, g_low_real, 1, g, 1)
            do i = 1, model_nat
               do k = 1, 3
                  g(k,i) = g(k,i) + g_high_model(k,i) - g_low_model(k,i)
               end do
            end do
            nat = model_nat
            do bond = 1, model_nbonds
               i = model_bonds(1,bond)
               j = model_bonds(2,bond)
               p = model_p(bond)
               nat = nat + 1
               do k = 1, 3
                  gdummy = g_high_model(k,nat) - g_low_model(k,nat)
                  g(k,i) = g(k,i) + gdummy*(1.0d0-p)
                  g(k,j) = g(k,j) + gdummy*p
               end do
            end do
         else
            call dcopy(3*real_nat, g_low_real, 1, g, 1)
            do i = 1, model_nat
               do k = 1, 3
                  g(k,i) = g(k,i) + 
     $                 g_high_model(k,i) - g_medium_model(k,i)
               end do
            end do
            nat = model_nat
            do bond = 1, model_nbonds
               i = model_bonds(1,bond)
               j = model_bonds(2,bond)
               p = model_p(bond)
               nat = nat + 1
               do k = 1, 3
                  gdummy = g_high_model(k,nat) - g_medium_model(k,nat)
                  g(k,i) = g(k,i) + gdummy*(1.0d0-p)
                  g(k,j) = g(k,j) + gdummy*p
               end do
            end do
            do i = 1, inter_nat
               do k = 1, 3
                  g(k,i) = g(k,i) + 
     $                 g_medium_inter(k,i) - g_low_inter(k,i)
               end do
            end do
            nat = inter_nat
            do bond = 1, inter_nbonds
               i = inter_bonds(1,bond)
               j = inter_bonds(2,bond)
               p = inter_p(bond)
               nat = nat + 1
               do k = 1, 3
                  gdummy = g_medium_inter(k,nat) - g_low_inter(k,nat)
                  g(k,i) = g(k,i) + gdummy*(1.0d0-p)
                  g(k,j) = g(k,j) + gdummy*p
               end do
            end do
         end if
c
c     Enforce active atom constraints on the gradient and also restore
c     the full list of active atoms from the copy made earlier.
c
         if (rtdb_ma_get(rtdb, 'oniom:actlist', ma_type, 
     $        nactive, l_actlist)) then
            if (.not. ma_get_index(l_actlist, k_actlist))
     $           call errquit('oniom: ma_get_index failed',l_actlist,
     &       MA_ERR)
            do i = 1, real_nat
               do j = 1, nactive
                  if (int_mb(k_actlist+j-1) .eq. i) goto 444
               end do
               do k = 1, 3
                  g(k,i) = 0.0d0
               end do
 444           continue
            end do
            if (.not. rtdb_put(rtdb, 'geometry:actlist', mt_int, 
     $           nactive, int_mb(k_actlist))) call errquit
     $           ('oniom: saving actlist failed',0,
     &       RTDB_ERR)
            if (.not. ma_free_heap(l_actlist)) call errquit
     $           ('oniom: free of actlist?',0, RTDB_ERR)
         end if
c         
         if (.not. rtdb_put(rtdb,'oniom:gradient',mt_dbl,3*real_nat,g))
     $        call errquit('oniom: failed storing gradient', 0,
     &       RTDB_ERR)
      end if
c     
c     Reset the geometry and restore any user indirection of basis
c     or geometry. 
c     
      ignore = rtdb_delete(rtdb, 'geometry')
      ignore = rtdb_delete(rtdb, 'ao basis')
      ignore = rtdb_delete(rtdb, 'ecp basis')
      ignore = rtdb_delete(rtdb, 'charge')
      if (gname .ne. ' ') then 
         if (.not. rtdb_cput(rtdb,'geometry',1,gname)) 
     $        call errquit('oniom:failed resetting user geom name',0,
     &       RTDB_ERR)
      end if
      if (bname .ne. ' ') then 
         if (.not. rtdb_cput(rtdb,'ao basis',1,bname)) call errquit
     $        ('oniom:failed resetting user basis name',0, RTDB_ERR)
      end if
      if (ename .ne. ' ') then 
         if (.not. rtdb_cput(rtdb,'ecp basis',1,ename)) call errquit
     $        ('oniom:failed resetting user ecp name',0,
     &       RTDB_ERR)
      end if
      if (real_charge .ne. -999.0d0) then
         if (.not. rtdb_put(rtdb, 'charge', mt_dbl, 1, real_charge))
     $        call errquit('oniom: resetting real charge?', 0,
     &       RTDB_ERR)
      end if
c
c     Don't leave MO vectors info lying around in the database.
c
      call oniom_unmanage_mos(rtdb)
c     
c     Get the geometry in preparation for printing
c     
      if (.not. geom_create(geom,'geometry')) 
     $     call errquit('oniom: failed creating geometry?',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'oniom_geometry'))
     $     call errquit('oniom: no geometry ', 0, RTDB_ERR)
      if (.not. geom_cart_get(geom,real_nat,tags,coords,charge))
     $     call errquit('oniom: geom_cart_get failed', 0, GEOM_ERR)
      if (.not. geom_destroy(geom))
     $     call errquit('oniom: geom destroy?',0, GEOM_ERR)
c     
      if (oprint) then
         call util_print_centered(6, 'ONIOM summary', 20, .true.)
         write(6,*)
         write(6,5) e_low_real, e_high_model, e_low_model, e2
 5       format(
     $        8x, '+   LOW + REAL  ', f20.12,/,
     $        8x, '+  HIGH + MODEL ', f20.12,/,
     $        8x, '-   LOW + MODEL ', f20.12,/,
     $        8x, '   ------------ ',/,
     $        8x, '=        ONIOM2 ', f20.12,/)
         if (nlayer .eq. 3) then
            write(6,55) e_low_real, e_high_model, e_medium_model, 
     $           e_medium_inter, e_low_inter, e3
 55         format(
     $           8x, '+   LOW + REAL  ', f20.12,/,
     $           8x, '+  HIGH + MODEL ', f20.12,/,
     $           8x, '-   MED + MODEL ', f20.12,/,
     $           8x, '+   MED + INTER ', f20.12,/,
     $           8x, '-   LOW + INTER ', f20.12,/,
     $           8x, '   ------------ ',/,
     $           8x, '=        ONIOM3 ', f20.12,/)
         end if
         if (task .eq. 'gradient') then
            write(6,1000) 'ONIOM',
     $           'x','y','z','x','y','z'
            do 30, i=1, real_nat
               write(6,2000) i, tags(i),(coords(j,i),j=1,3),
     $              (g(j,i),j=1,3)
 30         continue
            write(6,*)
 1000       format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $           'coordinates',
     $           24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000       format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
         end if
         call util_flush(6)
      end if
c     
      if (.not. rtdb_cput(rtdb, 'task:theory', 1, 'oniom'))
     $     call errquit('oniom: re-setting theory?',0,
     &       RTDB_ERR)
c     
      oniom = .true.
      call util_print_pop
      call ecce_print_module_exit('oniom', 'ok')
c
      end
      subroutine oniom_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb
c     
      logical status
      character*32 test
      integer maxbonds
      parameter (maxbonds = 50)
c
      integer model_bonds(2,maxbonds), model_nbonds, model_nat
      double precision model_p(maxbonds), model_charge
      character*32 model_tags(maxbonds)
c
      integer inter_bonds(2,maxbonds), inter_nbonds, inter_nat
      double precision inter_p(50), inter_charge
      character*32 inter_tags(maxbonds)
c
      character*32 high_theory, high_basis, high_ecp
      character*32 medium_theory, medium_basis, medium_ecp
      character*32 low_theory, low_basis, low_ecp
      character*1024 high_input, medium_input, low_input
c
      character*256 v_low_real, v_low_model, v_high_model,
     $     v_medium_model, v_medium_inter, v_low_inter
c
      integer bond, junk, i, j
      character*2 symbol
      character*8 element
      integer atn
c
      high_theory = ' '
      high_basis = ' '
      high_ecp = ' '
      high_input = ' '
      medium_theory = ' '
      medium_basis = ' '
      medium_ecp = ' ' 
      medium_input = ' '
      low_theory = ' '
      low_basis = ' '
      low_ecp = ' '
      low_input = ' '
      model_nat = 0
      model_nbonds = 0
      model_charge = -999.0     ! Used to indicate no charge
      inter_nat = 0
      inter_nbonds = 0
      inter_charge = -999.0
c
      v_low_real = ' '
      v_low_model = ' '
      v_high_model = ' '
      v_medium_model = ' '
      v_medium_inter = ' '
      v_low_inter = ' '
c
c     While input is available
c
 10   if (.not. inp_read()) 
     $     call errquit('oniom_input: premature EOF',0, INPUT_ERR)
      if (.not. inp_a(test))
     $     call errquit('oniom_input: failed to read keyword',0,
     &       INPUT_ERR)
c
      if (inp_compare(.false.,'high',test)) then
c
         call oniom_theory_input(high_theory,high_basis,high_ecp,
     $        high_input)
c     
      else if (inp_compare(.false.,'medium',test)) then
c
         call oniom_theory_input(medium_theory,medium_basis,medium_ecp,
     $        medium_input)
c     
      else if (inp_compare(.false.,'low',test)) then
c
         call oniom_theory_input(low_theory,low_basis,low_ecp,
     $        low_input)
c     
      else if (inp_compare(.false.,'print', test) .or.
     $        inp_compare(.false.,'noprint', test)) then
         call util_print_input(rtdb, 'oniom')
c     
      else if (inp_compare(.false.,'model',test)) then
c
c     model natom [charge q] [i1 j1 p1 [tag1] ...] 
c
         if (.not. inp_i(model_nat))
     $        call errquit('oniom_input:model natoms [i1 j1 p1 ...]',0,
     &       INPUT_ERR)
         if (inp_a(test)) then
            if (inp_compare(.false.,'charge',test)) then
               if (.not. inp_f(model_charge)) call errquit
     $              ('oniom_input: reading charge', 0, INPUT_ERR)
            else
               call inp_prev_field()
            end if
         end if
         model_nbonds = 0
         do bond = 1, maxbonds
            if (.not. inp_i(model_bonds(1,bond))) goto 501
            if (.not. inp_i(model_bonds(2,bond))) call errquit
     $           ('oniom_input: reading second atom for bond',bond,
     &       INPUT_ERR)
            if (.not. inp_f(model_p(bond))) call errquit
     $           ('oniom_input: reading scale factor for bond',bond,
     &       INPUT_ERR)
            model_tags(bond) = 'H oniom'
            if (inp_cur_field() .lt. inp_n_field()) then
               if (inp_i(junk)) then
                  call inp_prev_field()
               else if (.not. inp_a(model_tags(bond))) then
                  call errquit('oniom_input: reading tag for bond',bond,
     &       INPUT_ERR)
               end if
            end if
            if (.not. geom_tag_to_element(model_tags(bond), symbol, 
     $           element, atn)) then
               if (symbol.ne.'bq') call errquit
     $              ('oniom_input: tag must be bq* or element*',0,
     &       INPUT_ERR)
            end if
            model_nbonds = model_nbonds + 1
         end do
 501     continue
      else if (inp_compare(.false.,'inter',test)) then
c
c     inter natom [i1 j1 p1 [tag1] ...]
c
         if (.not. inp_i(inter_nat))
     $        call errquit('oniom_input:inter natoms [i1 j1 p1 ...]',0,
     &       INPUT_ERR)
         if (inp_a(test)) then
            if (inp_compare(.false.,'charge',test)) then
               if (.not. inp_f(inter_charge)) call errquit
     $              ('oniom_input: reading charge', 0, INPUT_ERR)
            else
               call inp_prev_field()
            end if
         end if
         inter_nbonds = 0
         do bond = 1, maxbonds
            if (.not. inp_i(inter_bonds(1,bond))) goto 503
            if (.not. inp_i(inter_bonds(2,bond))) call errquit
     $           ('oniom_input: reading second atom for bond',bond,
     &       INPUT_ERR)
            if (.not. inp_f(inter_p(bond))) call errquit
     $           ('oniom_input: reading scale factor for bond',bond,
     &       INPUT_ERR)
            inter_tags(bond) = 'H oniom'
            if (inp_cur_field() .lt. inp_n_field()) then
               if (inp_i(junk)) then
                  call inp_prev_field()
               else if (.not. inp_a(inter_tags(bond))) then
                  call errquit('oniom_input: reading tag for bond',bond,
     &       INPUT_ERR)
               end if
            end if
            if (.not. geom_tag_to_element(inter_tags(bond), symbol, 
     $           element, atn)) then
               if (symbol.ne.'bq') call errquit
     $              ('oniom_input: tag must be bq* or element*',0,
     &       INPUT_ERR)
            end if
            inter_nbonds = inter_nbonds + 1
         end do
 503     continue
c
      else if (inp_compare(.false.,'vectors',test)) then
 111     if (inp_a(test)) then
            if (inp_compare(.false.,'low-real',test)) then
               if (.not. inp_a(v_low_real)) call errquit
     $              ('oniom_input: error reading low_real vectors', 0,
     &       INPUT_ERR)
            else if (inp_compare(.false.,'low-model',test)) then
               if (.not. inp_a(v_low_model)) call errquit
     $              ('oniom_input: error reading low_model vectors', 0,
     &       INPUT_ERR)
            else if (inp_compare(.false.,'high-model',test)) then
               if (.not. inp_a(v_high_model)) call errquit
     $              ('oniom_input: error reading high_model vectors', 0,
     &       INPUT_ERR)
            else if (inp_compare(.false.,'medium-model',test)) then
               if (.not. inp_a(v_medium_model)) call errquit
     $              ('oniom_input:error reading medium_model vectors',0,
     &       INPUT_ERR)
            else if (inp_compare(.false.,'medium-inter',test)) then
               if (.not. inp_a(v_medium_inter)) call errquit
     $              ('oniom_input:error reading medium_inter vectors',0,
     &       INPUT_ERR)
            else if (inp_compare(.false.,'low-inter',test)) then
               if (.not. inp_a(v_low_inter)) call errquit
     $              ('oniom_input: error reading low_inter vectors', 0,
     &       INPUT_ERR)
            else
               call errquit('oniom_input:unknown keyword for vectors',0,
     &       INPUT_ERR)
            end if
            goto 111
         end if
      else if (.not. inp_compare(.false.,'end',test)) then
         call errquit('oniom_input: unknown directive', 0, INPUT_ERR)
      end if
      if (.not. inp_compare(.false.,'end',test)) goto 10
c
c     Sanity tests
c
      if (high_theory.eq.' ' .or. low_theory.eq.' ') call errquit
     $     ('oniom_input: must specify both high and low theory',0,
     &       INPUT_ERR)
      if (model_nat .le. 0) call errquit
     $     ('oniom_input: must specify no. of model atoms',0, INPUT_ERR)
      if ( (medium_theory.ne.' ' .and. inter_nat.eq.0) .or.
     $     (medium_theory.eq.' ' .and. inter_nat.gt.0) ) call errquit
     $     ('oniom_input: 3-layer model requires both medium theory'//
     $     ' and intermediate layer',0, INPUT_ERR)
c
c     Put the bond info into the right order 1=model, 2=real and
c     also do a little checking.
c
      do bond = 1, model_nbonds
         i = min(model_bonds(1,bond),model_bonds(2,bond))
         j = max(model_bonds(1,bond),model_bonds(2,bond))
         if (i.le.0 .or. i.gt.model_nat .or. j.le.model_nat) 
     $        call errquit
     $        ('oniom_input: bad atom info for model bond number', bond,
     &       INPUT_ERR)
         model_bonds(1,bond) = i
         model_bonds(2,bond) = j
      end do
      do bond = 1, inter_nbonds
         i = min(inter_bonds(1,bond),inter_bonds(2,bond))
         j = max(inter_bonds(1,bond),inter_bonds(2,bond))
         if (i.le.0 .or. i.gt.inter_nat .or. j.le.inter_nat) 
     $        call errquit
     $        ('oniom_input: bad atom info for inter bond number', bond,
     &       INPUT_ERR)
         inter_bonds(1,bond) = i
         inter_bonds(2,bond) = j
      end do
c
c     Store info to the database
c
      status = rtdb_cput(rtdb, 'oniom:high:theory', 1, high_theory)
      status = rtdb_cput(rtdb, 'oniom:low:theory', 1, low_theory)
     $     .and. status
      if (medium_theory .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:medium:theory', 1, 
     $     medium_theory) .and. status
c
      if (high_basis .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:high:basis', 1, high_basis)
     $     .and. status
      if (medium_basis .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:medium:basis', 1, 
     $     medium_basis).and. status
      if (low_basis .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:low:basis', 1, low_basis)
     $     .and. status
c
      if (high_ecp .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:high:ecp', 1, high_ecp)
     $     .and. status
      if (medium_ecp .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:medium:ecp', 1, 
     $     medium_ecp).and. status
      if (low_ecp .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:low:ecp', 1, low_ecp)
     $     .and. status
c
      if (high_input .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:high:input', 1, high_input)
     $     .and. status
      if (medium_input .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:medium:input', 1, 
     $     medium_input).and. status
      if (low_input .ne. ' ') 
     $     status = rtdb_cput(rtdb, 'oniom:low:input', 1, low_input)
     $     .and. status
c
      if (v_low_real .ne. ' ') status = status .and. 
     $     rtdb_cput(rtdb, 'oniom:v_low_real', 1, v_low_real)
      if (v_low_model .ne. ' ')  status = status .and. 
     $     rtdb_cput(rtdb, 'oniom:v_low_model', 1, v_low_model)
      if (v_high_model .ne. ' ')  status = status .and. 
     $     rtdb_cput(rtdb, 'oniom:v_high_model', 1, v_high_model)
      if (v_medium_model .ne. ' ')  status = status .and. 
     $     rtdb_cput(rtdb, 'oniom:v_medium_model', 1, v_medium_model)
      if (v_medium_inter .ne. ' ')  status = status .and. 
     $     rtdb_cput(rtdb, 'oniom:v_medium_inter', 1, v_medium_inter)
      if (v_low_inter .ne. ' ')  status = status .and. 
     $     rtdb_cput(rtdb, 'oniom:v_low_inter', 1, v_low_inter)
c
      status = rtdb_put(rtdb, 'oniom:model:nat', mt_int, 1, model_nat)
     $     .and. status
      status = rtdb_put(rtdb, 'oniom:model:charge', mt_dbl, 1, 
     $     model_charge) .and. status
      if (model_nbonds .gt. 0) then
         status = rtdb_put(rtdb, 'oniom:model:nbonds', mt_int, 1, 
     $        model_nbonds) .and. status
         status = rtdb_put(rtdb, 'oniom:model:bonds', mt_int, 
     $        2*model_nbonds, model_bonds) .and. status
         status = rtdb_put(rtdb, 'oniom:model:p', mt_dbl, 
     $        model_nbonds, model_p) .and. status
         status = rtdb_cput(rtdb,'oniom:model:tags', 
     $        model_nbonds, model_tags) .and. status
      end if
c
      if (inter_nat .gt. 0) then
         status = rtdb_put(rtdb, 'oniom:inter:nat', mt_int, 1, 
     $        inter_nat) .and. status
         status = rtdb_put(rtdb, 'oniom:inter:charge', mt_dbl, 1, 
     $        inter_charge) .and. status
         if (inter_nbonds .gt. 0) then
            status = rtdb_put(rtdb, 'oniom:inter:nbonds', mt_int, 1, 
     $           inter_nbonds) .and. status
            status = rtdb_put(rtdb, 'oniom:inter:bonds', mt_int, 
     $           2*inter_nbonds, inter_bonds) .and. status
            status = rtdb_put(rtdb, 'oniom:inter:p', mt_dbl, 
     $           inter_nbonds, inter_p) .and. status
            status = rtdb_cput(rtdb,'oniom:inter:tags', 
     $           inter_nbonds, inter_tags) .and. status
         end if
      end if
c
      if (.not. status) call errquit
     $     ('oniom_input: error putting data into RTDB',0, RTDB_ERR)
c
      end
      subroutine oniom_unmanage_mos(rtdb)
#include "rtdb.fh"
      integer rtdb
c
c     Delete the vectors input/output keywords from the
c     popular (scf,dft,mcscf) places so that there is
c     no confusion or overwriting down the road
c
      logical ignore
c
      ignore = rtdb_delete(rtdb, 'scf:input vectors')
      ignore = rtdb_delete(rtdb, 'dft:input vectors')
      ignore = rtdb_delete(rtdb, 'pspw:input vectors')
      ignore = rtdb_delete(rtdb, 'band:input vectors')
      ignore = rtdb_delete(rtdb, 'mcscf:input vectors')
      ignore = rtdb_delete(rtdb, 'scf:output vectors')
      ignore = rtdb_delete(rtdb, 'dft:output vectors')
      ignore = rtdb_delete(rtdb, 'pspw:output vectors')
      ignore = rtdb_delete(rtdb, 'band:output vectors')
      ignore = rtdb_delete(rtdb, 'mcscf:output vectors')
c
      end
      subroutine oniom_manage_mos(rtdb,theory,movecs)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
      integer rtdb
      character*(*) theory
      character*(*) movecs
c
c     Defeat the convergence tests on the MO methods which
c     are confused by some ONIOM calculations and tell the
c     MO codes to use the correct MO vectors.
c
      logical ignore
      logical oexists
      logical do_mcscf
      character*32 keyin, keyout
c
      do_mcscf = .false.
      if (ga_nodeid() .eq. 0) then
         keyin  = 'scf:input vectors'
         keyout = 'scf:output vectors'
         if (inp_compare(.false.,'dft',theory)) then
            keyin  = 'dft:input vectors'
            keyout = 'dft:output vectors'
         else if (inp_compare(.false.,'pspw',theory)) then
            keyin  = 'pspw:input vectors'
            keyout = 'pspw:output vectors'
         else if (inp_compare(.false.,'band',theory)) then
            keyin  = 'band:input vectors'
            keyout = 'band:output vectors'
         else if (inp_compare(.false.,'mcscf',theory)) then
            keyin  = 'mcscf:input vectors'
            keyout = 'mcscf:output vectors'
            do_mcscf = .true.
         else if (inp_compare(.false.,'selci',theory)) then
            keyin  = 'mcscf:input vectors'
            keyout = 'mcscf:output vectors'
            do_mcscf = .true.
         end if
         inquire(file=movecs,exist=oexists)
c
         ignore = rtdb_parallel(.false.)
         ignore = rtdb_delete(rtdb, 'dft:converged')  ! DO WE REALLY NEED TO DO THIS?
         ignore = rtdb_delete(rtdb, 'scf:converged')
         ignore = rtdb_delete(rtdb, 'mcscf:converged')
         ignore = rtdb_delete(rtdb, keyin)
         ignore = rtdb_delete(rtdb, keyout)
         if (oexists) then
            if (.not. rtdb_cput(rtdb,keyin,1,movecs)) call errquit
     $           ('oniom_manage_mos: failed settting input vectors',0,
     &       RTDB_ERR)
         else if (do_mcscf) then
c           we need to prevent the MCSCF code picking up random SCF vectors
            if (.not. rtdb_cput(rtdb,keyin,1,'atomic')) call errquit
     $           ('oniom_manage_mos: failed settting input vectors',0,
     &       RTDB_ERR)
         end if
         if (.not. rtdb_cput(rtdb,keyout,1,movecs)) call errquit
     $        ('oniom_manage_mos: failed settting output vectors',0,
     &       RTDB_ERR)
         ignore = rtdb_parallel(.true.)
      end if
c
      call ga_sync()
c
      end
      subroutine oniom_munge_actlist(rtdb,nat,nbonds,bonds)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      integer rtdb
      integer nat, nbonds, bonds(2,*)
c
c     For constrained geometry optimizations and numerical frequencies 
c     the list of active atoms for computing the gradient may be set.
c     For ONIOM subcalculations we need to prune the list of atoms
c     to those actually active in the calculation (otherwise some of
c     the other modules rightfully complain) and if either of the
c     atoms in a bond are active we need to ensure that the link
c     atom is also marked active.
c
      integer l_actlist, k_actlist, nactive, i, j, bond, ma_type
      integer active(nw_max_atom)
c
      if (.not. rtdb_ma_get(rtdb, 'oniom:actlist', ma_type, 
     $     nactive, l_actlist)) return
c
      if (.not. ma_get_index(l_actlist, k_actlist))
     $     call errquit('grad_act_at: ma_get_inded failed',l_actlist,
     &       MA_ERR)
c
c     Prune out atoms not in the subset
c
      do i = 1, nat + nbonds
         active(i) = 0
      end do
      do i = 1, nactive
         j = int_mb(k_actlist+i-1)
         if (j.le.nat) active(j) = 1
      end do
c
c     Ensure link atoms are active if necessary
c
      do bond = 1, nbonds
         i = bonds(1,bond)
         j = bonds(2,bond)
         if (active(i).gt.0 .or. active(j).gt.0) active(nat+bond) = 1
      end do
c
c     Compress the active flags down into a list of atoms
c
      nactive = 0
      do i = 1, nat + nbonds
         if (active(i).gt.0) then
            nactive = nactive + 1
            active(nactive) = i
         end if
      end do
c
      if (nactive .eq. 0) call errquit('oniom: nactive=0?',0,
     &       INPUT_ERR)
c
      if (.not. rtdb_put(rtdb, 'geometry:actlist', mt_int, nactive, 
     $     active)) call errquit('oniom: failed storing actlist',0,
     &       RTDB_ERR)
      if (.not. ma_free_heap(l_actlist)) call errquit
     $     ('oniom: free of actlist?',0, MA_ERR)
c
      end
      subroutine oniom_theory_input(theory,basname,ecpname,input)
      implicit none
#include "errquit.fh"
      character*(*) theory, basname, ecpname, input ! [output]
#include "inp.fh"
c
c     (high||medium||low) theory [basis basname] [ecp ecpname] [input inpstr]
c
c     Read oniom theory input ... defaults assumed set outside
c
c     The theory string should always be converted to lower case as the
c     task layer uses lower case strings only.
c
      character*128 test
c
      if (.not. inp_a(theory)) call errquit
     $     ('oniom_input: failed reading theory',0, INPUT_ERR)
      call inp_lcase(theory)
c
 10   if (.not. inp_a(test)) return
      if (inp_compare(.false.,'basis',test)) then
         if (.not. inp_a(basname)) call errquit
     $        ('oniom_input: failed reading basname',0, INPUT_ERR)
      else if (inp_compare(.false.,'ecp',test)) then
         if (.not. inp_a(ecpname)) call errquit
     $        ('oniom_input: failed reading ecpname',0, INPUT_ERR)
      else if (inp_compare(.false.,'input',test)) then
         if (.not. inp_a(input)) call errquit
     $        ('oniom_input: failed reading input [input]',0, INPUT_ERR)
      else
         call errquit('oniom: unknown keyword for high||medium||low',0,
     &       INPUT_ERR)
      end if
      goto 10
c
      end
      subroutine oniom_param(rtdb, geom, basis, ecp, charge, input)
      implicit none
#include "errquit.fh"
      integer rtdb
      character*(*) geom, basis, ecp, input
      double precision charge
#include "rtdb.fh"
#include "mafdecls.fh"
c
      logical ignore
      logical nw_inp_from_character
      external nw_inp_from_character
c
      if (charge .ne. -999.0d0) then
         if (.not. rtdb_put(rtdb,'charge',mt_dbl,1,charge))
     $        call errquit('oniom:setting charge?',0, RTDB_ERR)
      end if
      if (geom .ne. ' ') then
         if (.not. rtdb_cput(rtdb,'geometry',1,geom))
     $        call errquit('oniom: setting geometry?',0, RTDB_ERR)
      end if
c
      if (basis .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'ao basis', 1, basis)) 
     $        call errquit('oniom: setting basis?',0, RTDB_ERR)
      end if
c
      if (ecp .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'ecp basis', 1, ecp)) 
     $        call errquit('oniom: setting ecp?',0, RTDB_ERR)
      else
         ignore = rtdb_delete(rtdb,'ecp basis')
      end if
c
      if (input .ne. ' ') then
         if (.not. nw_inp_from_character(rtdb,input))
     $        call errquit('oniom: error processing input string',0,
     &       INPUT_ERR)
      endif
c     
      end
