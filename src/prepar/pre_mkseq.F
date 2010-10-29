      logical function pre_mkseq(irtdb,lfnout,source,sysnam,
     + lfnpdb,filpdb,lfnseq,filseq,lfnpar,lfnfrg,lfnsgm,
     + lfntmp,filtmp,lfnmod,filmod,
     + slvnam,slvmdl,maxscf,qscale,altloc,chain,icyren,model,nxlnk,
     + mcount,ncount,icount,scount,mfract,nfract,ifract,iconst,his)
c
c $Id$
c
c     in : integer lfnout    = logical file number output file
c          integer lfnpdb    = logical file number for pdb file
c          char*80 filpdb    = file name of pdb file
c          integer lfnseq    = logical file number for seq file
c          char*80 filseq    = file name of seq file
c          integer lfnpar    = logical file number for typ file
c          integer lfnfrg    = logical file number for frg file
c
      implicit none
c
#include "mafdecls.fh"
#include "util.fh"
#include "geom.fh"
#include "pre_common.fh"
c
      integer irtdb
      logical pre_pdbsiz,pre_pdbseq,pre_rtdbsq,pre_wrtnam
      logical pre_rename,pre_sgmnum,pre_sgmndx,pre_gensgm
      external pre_pdbsiz,pre_pdbseq,pre_rtdbsq,pre_wrtnam
      external pre_rename,pre_sgmnum,pre_sgmndx,pre_gensgm
      logical pre_mkfrg,pre_namiup,pre_links,pre_wrtseq
      external pre_mkfrg,pre_namiup,pre_links,pre_wrtseq
c
      integer matm,mseq,mssb,msgm,mbnd,mang,mdih,mimp,mlnk,mato,nato
      integer natm,nseq,nssb,nsgm,nlnk,maxscf,icyren,model,nxlnk
      integer mcount,ncount,iconst
      integer mfract,nfract,ifract(mfract)
      integer mfrb,i_frb,l_frb
      integer l_lseq,i_lseq,l_cseq,i_cseq
      integer l_lsgm,i_lsgm,l_csgm,i_csgm
      integer l_lssb,i_lssb,l_llnk,i_llnk,l_clnk,i_clnk
      integer l_latm,i_latm,l_catm,i_catm,l_xatm,i_xatm,l_qatm,i_qatm
      integer l_bnd,i_bnd,l_ang,i_ang,l_dih,i_dih,l_imp,i_imp
      integer icount(mcount)
      integer l_lato,i_lato,i_cato,l_cato,l_xato,i_xato,l_qato,i_qato
c
      real*8 qscale
c
      character*255 filpdb,filseq,filtmp,filmod
      character*80 sysnam
      character*3 slvnam
      character*80 source
      character*10 slvmdl
      character*1 altloc,chain,his
      character*4 scount(mcount)
c
      integer lfnpdb,lfnout,lfnfrg,lfnpar,lfnseq,lfnsgm,lfntmp,lfnmod
      integer igeom,numslv,nlnkf
c
      if(util_print('sequence',print_high)) then
      write(lfnout,1000)
 1000 format(' SEQUENCE GENERATION',/)
      endif
c
      if(source(1:3).eq.'pdb') then
      if(.not.pre_pdbsiz(lfnpdb,filpdb,lfnout,model,
     + natm,nseq,nssb,nlnk,slvnam,numslv,altloc,chain))
     + call md_abort('pre_pdbsiz failed',9999)
      else
      if(.not.geom_create(igeom,source))
     + call md_abort('pre_mkseq: geom_create failed',0)
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call md_abort('pre_mkseq: geom_rtdb_load failed',0)
      if(.not.geom_ncent(igeom,natm))
     + call md_abort('pre_mkseq: geom_ncent failed',0)
      nseq=1
      nssb=0
      nlnk=0
      endif
c
c     set dimensions
c     --------------
c
      matm=natm
      mseq=nseq+2
      mssb=nssb
      mlnk=nlnk+5*nseq
c
c     allocate memory for atom list
c     -----------------------------
c
c     integer latm(1,matm) : index in sequence list
c                  2       : atomic number (0=undeterminable)
c                  3       : number of bonds to this atom
c                  4       : atom center type
c                  5       : link number
c
      if(.not.ma_push_get(mt_int,5*matm,'latm',l_latm,i_latm))
     + call md_abort('Memory allocation failed for latm',9999)
c
c     char*6  catm(1,matm) : atom name on pdb
c                  2       : atom name on top
c                  3       : atom type
c
      if(.not.ma_push_get(mt_byte,18*matm,'catm',l_catm,i_catm))
     + call md_abort('Memory allocation failed for catm',9999)
c
c     real*8  xatm(1:3,matm) : atomic coordinates in nm
c
      if(.not.ma_push_get(mt_dbl,3*matm,'xatm',l_xatm,i_xatm))
     + call md_abort('Memory allocation failed for xatm',9999)
c
c     real*8 qatm(matm) : partial atomic charges
c
      if(.not.ma_push_get(mt_dbl,matm,'qatm',l_qatm,i_qatm))
     + call md_abort('Memory allocation failed for qatm',9999)
c
c     allocate memory for secondary segment links
c     -------------------------------------------
c
c     integer lssb(1,mssb) : link segment 1
c                  2       : link segment 2
c
      if(.not.ma_push_get(mt_int,2*mssb,'lssb',l_lssb,i_lssb))
     + call md_abort('Memory allocation failed for lssb',9999)
c
c     integer llnk(1,mlnk) : link segment 1
c                  2       : link segment 2 atom number
c                  3       : link segment 1 global atom number
c                  4       : link segment 2 global atom number
c     char*4  clnk(1,mlnk) : link segment 1 atom name
c                  2       : link segment 2 atom name
c
      if(.not.ma_push_get(mt_int,4*mlnk,'llnk',l_llnk,i_llnk))
     + call md_abort('Memory allocation failed for llnk',9999)
      if(.not.ma_push_get(mt_byte,8*mlnk,'clnk',l_clnk,i_clnk))
     + call md_abort('Memory allocation failed for clnk',9999)
c
c     allocate memory for sequence
c     ----------------------------
c
c     integer lseq(1,mseq) : number of sequence entry
c                  2       : index into unique segment list
c                  3       : index of first atom of sequence entry
c                  4       : molecule number
c
      if(.not.ma_push_get(mt_int,6*mseq,'lseq',l_lseq,i_lseq))
     + call md_abort('Memory allocation failed for lseq',9999)
c
c     char*5  cseq(1,mseq) : name of sequence entry on pdb
c                  2       : name of sequence entry on top
c
      if(.not.ma_push_get(mt_byte,10*mseq,'cseq',l_cseq,i_cseq))
     + call md_abort('Memory allocation failed for cseq',9999)
c
c
c     read the sequence and atom list from the pdb file
c     -------------------------------------------------
c
      if(source(1:3).eq.'pdb') then
      if(.not.pre_pdbseq(lfnpdb,filpdb,lfnout,model,slvnam,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_xatm),matm,natm,
     + int_mb(i_lssb),mssb,nssb,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk,altloc,chain))
     + call md_abort('pre_pdbseq failed',9999)
      if(util_print('files',print_default)) then
      write(lfnout,2000) filpdb(1:index(filpdb,' ')-1)
 2000 format(' PDB geometry',t40,a,/)
      endif
      else
      if(.not.pre_rtdbsq(irtdb,igeom,lfnout,sysnam,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_xatm),matm,natm))
     + call md_abort('pre_rtdbsq failed',9999)
      if(util_print('files',print_default)) then
      write(lfnout,2001)
 2001 format(' RTDB geometry',/)
      endif
      endif
c
c     comply with iupac atom names
c     ----------------------------
c
      if(.not.pre_namiup(lfnout,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_latm),byte_mb(i_catm),matm,natm))
     + call md_abort('pre_namiup failed',9999)
c
c     assign segments if not obtained from pdb
c     ----------------------------------------
c
c     dallocate memory for sequence
c     -----------------------------
c
c      if(.not.ma_pop_stack(l_cseq))
c     + call md_abort('Memory deallocation failed for cseq',9999)
c      if(.not.ma_pop_stack(l_lseq))
c     + call md_abort('Memory deallocation failed for lseq',9999)
c
c     reallocate memory for sequence
c     ------------------------------
c
c     integer lseq(1,mseq) : number of sequence entry
c                  2       : index into unique segment list
c                  3       : index of first atom of sequence entry
c                  4       : molecule number
c
c      if(.not.ma_push_get(mt_int,4*mseq,'lseq',l_lseq,i_lseq))
c     + call md_abort('Memory allocation failed for lseq',9999)
c
c     char*5  cseq(1,mseq) : name of sequence entry on pdb
c                  2       : name of sequence entry on top
c
c      if(.not.ma_push_get(mt_byte,10*mseq,'cseq',l_cseq,i_cseq))
c     + call md_abort('Memory allocation failed for cseq',9999)
c
c     check segment links
c     -------------------
c
      if(.not.pre_links(lfnout,lfnmod,filmod,
     + dbl_mb(i_xatm),byte_mb(i_catm),int_mb(i_latm),
     + matm,int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,int_mb(i_llnk),
     + byte_mb(i_clnk),mlnk,nlnk,nlnkf))
     + call md_abort('pre_links failed',9999)
c
c     comply with force field segment names
c     -------------------------------------
c
      if(.not.pre_rename(lfnout,lfnmod,filmod,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_latm),byte_mb(i_catm),matm,natm,int_mb(i_llnk),
     + byte_mb(i_clnk),mlnk,nlnk,icyren,his))
     + call md_abort('pre_rename failed',9999)
c
c     write atom names to temporary file
c     ----------------------------------
c
      if(.not.pre_wrtnam(lfnout,lfntmp,filtmp,
     + byte_mb(i_catm),matm,natm))
     + call md_abort('pre_wrtnam failed',9999)
c
c     get number of unique segments
c     -----------------------------
c
      if(.not.pre_sgmnum(lfnout,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,nsgm))
     + call md_abort('pre_sgmnum failed',9999)
c
c     set dimension
c     -------------
c
      msgm=nsgm
c
c     allocate memory for the segments
c     --------------------------------
c
c     integer lsgm(1,msgm) : number of segments of this type
c                  2       : source, 0,1,2,3,4=?,s,x,u,t
c                  3       : number of atoms in segment
c
      if(.not.ma_push_get(mt_int,3*msgm,'lsgm',l_lsgm,i_lsgm))
     + call md_abort('Memory allocation failed for lsgm',9999)
c
c     char*5  csgm(msgm)   : segment names
c
      if(.not.ma_push_get(mt_byte,5*msgm,'csgm',l_csgm,i_csgm))
     + call md_abort('Memory allocation failed for csgm',9999)
c
c     scan the database directories for segment files
c     -----------------------------------------------
c
      if(.not.pre_sgmndx(lfnfrg,lfnout,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_lsgm),byte_mb(i_csgm),msgm,nsgm,mato)) then
c
c     set dimensions for bonds, angles, torsions and out-of-plane
c     -----------------------------------------------------------
c
      mato=max(4*mato,6)
      mbnd=2*max(matm,mato)
      mang=2*max(matm,mato)
      mdih=3*max(matm,mato)
      mimp=max(matm,mato)
      mfrb=2*mato
c
c     allocate memory for segment atom lists
c     --------------------------------------
c
      if(.not.ma_push_get(mt_int,5*mato,'lato',l_lato,i_lato))
     + call md_abort('Memory allocation failed for lato',9999)
c
      if(.not.ma_push_get(mt_byte,18*mato,'cato',l_cato,i_cato))
     + call md_abort('Memory allocation failed for cato',9999)
c
      if(.not.ma_push_get(mt_dbl,3*mato,'xato',l_xato,i_xato))
     + call md_abort('Memory allocation failed for xato',9999)
c
      if(.not.ma_push_get(mt_dbl,mato,'qato',l_qato,i_qato))
     + call md_abort('Memory allocation failed for qato',9999)
c
      if(.not.ma_push_get(mt_int,2*mfrb,'frb',l_frb,i_frb))
     + call md_abort('mkseq: Memory allocation failed for frb',9999)
c
c     allocate memory for the bonded lists
c     ------------------------------------
c
c     integer lbnd(1:2,mbnd) : bond indices
c
      if(.not.ma_push_get(mt_int,2*mbnd,'bnd',l_bnd,i_bnd))
     + call md_abort('mkseq: Memory allocation failed for bnd',9999)
c
c     integer lang(1:3,mang) : angle indices
c
      if(.not.ma_push_get(mt_int,3*mang,'ang',l_ang,i_ang))
     + call md_abort('mkseq: Memory allocation failed for ang',9999)
c
c     integer ldih(1:4,mdih) : torsion indices
c
      if(.not.ma_push_get(mt_int,4*mdih,'dih',l_dih,i_dih))
     + call md_abort('mkseq: Memory allocation failed for dih',9999)
c
c     integer limp(1:4,mimp) : improper torsion indices
c
      if(.not.ma_push_get(mt_int,4*mimp,'imp',l_imp,i_imp))
     + call md_abort('mkseq: Memory allocation failed for imp',9999)
c
c     generate fragment files for missing segments in dir_t
c
      if(.not.pre_mkfrg(irtdb,source,lfnfrg,lfnpar,
     + lfnout,iconst,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_lsgm),byte_mb(i_csgm),msgm,nsgm,int_mb(i_latm),
     + byte_mb(i_catm),dbl_mb(i_xatm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_lato),byte_mb(i_cato),dbl_mb(i_xato),dbl_mb(i_qato),
     + mato,nato,int_mb(i_frb),mfrb,
     + int_mb(i_bnd),mbnd,int_mb(i_ang),mang,int_mb(i_dih),mdih,
     + int_mb(i_imp),mimp,maxscf,qscale,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk))
     + call md_abort('pre_mkfrg failed',9999)
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_imp))
     + call md_abort('mkseq: Memory deallocation failed for imp',9999)
      if(.not.ma_pop_stack(l_dih))
     + call md_abort('mkseq: Memory deallocation failed for dih',9999)
      if(.not.ma_pop_stack(l_ang))
     + call md_abort('mkseq: Memory deallocation failed for ang',9999)
      if(.not.ma_pop_stack(l_bnd))
     + call md_abort('mkseq: Memory deallocation failed for bnd',9999)
      if(.not.ma_pop_stack(l_frb))
     + call md_abort('Memory deallocation failed for frb',9999)
      if(.not.ma_pop_stack(l_qato))
     + call md_abort('Memory deallocation failed for qato',9999)
      if(.not.ma_pop_stack(l_xato))
     + call md_abort('Memory deallocation failed for xato',9999)
      if(.not.ma_pop_stack(l_cato))
     + call md_abort('Memory deallocation failed for cato',9999)
      if(.not.ma_pop_stack(l_lato))
     + call md_abort('Memory deallocation failed for lato',9999)
      endif
c
c     generate segment files
c     ----------------------
c
      if(.not.pre_gensgm(lfnout,lfnfrg,lfnsgm,
     + int_mb(i_lsgm),byte_mb(i_csgm),msgm,nsgm))
     + call md_abort('pre_gensgm failed',9999)
c
      if(util_print('where',print_debug)) then
      write(lfnout,1001)
 1001 format(' Segment files generated',/)
      endif
c
c     write the sequence file
c     -----------------------
c
      if(.not.pre_wrtseq(lfnseq,filseq,lfnout,
     + byte_mb(i_cseq),int_mb(i_lseq),mseq,nseq,
     + int_mb(i_lssb),mssb,nssb,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk,nlnkf,slvmdl,
     + mcount,ncount,icount,scount,mfract,nfract,ifract))
     + call md_abort('pre_wrtseq failed',9999)
c
      if(util_print('where',print_debug)) then
      write(lfnout,1002)
 1002 format(' Sequence file generated',/)
      endif
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_csgm))
     + call md_abort('Memory deallocation failed for csgm',9999)
      if(.not.ma_pop_stack(l_lsgm))
     + call md_abort('Memory deallocation failed for lsgm',9999)
      if(.not.ma_pop_stack(l_cseq))
     + call md_abort('Memory deallocation failed for cseq',9999)
      if(.not.ma_pop_stack(l_lseq))
     + call md_abort('Memory deallocation failed for lseq',9999)
      if(.not.ma_pop_stack(l_clnk))
     + call md_abort('Memory deallocation failed for clnk',9999)
      if(.not.ma_pop_stack(l_llnk))
     + call md_abort('Memory deallocation failed for llnk',9999)
      if(.not.ma_pop_stack(l_lssb))
     + call md_abort('Memory deallocation failed for lssb',9999)
      if(.not.ma_pop_stack(l_qatm))
     + call md_abort('Memory deallocation failed for qatm',9999)
      if(.not.ma_pop_stack(l_xatm))
     + call md_abort('Memory deallocation failed for xatm',9999)
      if(.not.ma_pop_stack(l_catm))
     + call md_abort('Memory deallocation failed for catm',9999)
      if(.not.ma_pop_stack(l_latm))
     + call md_abort('Memory deallocation failed for latm',9999)
c
      if(source(1:3).ne.'pdb') then
      if(.not.geom_destroy(igeom))
     + call md_abort('pre_mkseq: geom_destroy failed',9999)
      endif
c
      pre_mkseq=.true.
      return
      end
