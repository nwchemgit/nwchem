      logical function pre_mkseq(irtdb,lfnout,lprint,source,ffield,
     + lfnpdb,filpdb,lfnseq,filseq,lfnpar,filpar,lfnfrg,lfnsgm,
     + lfntmp,filtmp,dir_s,dir_x,dir_u,dir_t,slvnam,slvmdl)
c
c     in : integer lfnout    = logical file number output file
c          char*80 ffield    = force field from [amber]
c          integer lfnpdb    = logical file number for pdb file
c          char*80 filpdb    = file name of pdb file
c          integer lfnseq    = logical file number for seq file
c          char*80 filseq    = file name of seq file
c          integer lfnpar    = logical file number for typ file
c          char*80 filpar    = file name of typ file
c          integer lfnfrg    = logical file number for frg file
c          char*80 dir_s     = directory standard ffield files
c          char*80 dir_x     = directory extensions ffield files
c          char*80 dir_u     = directory user preferences ffield files
c          char*80 dir_t     = directory temporary ffield files
c
      implicit none
c
#include "mafdecls.fh"
#include "geom.fh"
c
      integer irtdb
      logical pre_pdbsiz,pre_pdbseq,pre_rtdbsq,pre_wrtnam
      logical pre_namamb,pre_sgmnum,pre_sgmndx,pre_gensgm
      external pre_pdbsiz,pre_pdbseq,pre_rtdbsq,pre_wrtnam
      external pre_namamb,pre_sgmnum,pre_sgmndx,pre_gensgm
      logical pre_mkfrg,pre_namiup,pre_links,pre_wrtseq
      external pre_mkfrg,pre_namiup,pre_links,pre_wrtseq
c
      integer matm,mseq,mssb,msgm,mbnd,mang,mdih,mimp,mlnk
      integer natm,nseq,nssb,nsgm,nlnk
      integer lprint
      integer l_lseq,i_lseq,l_cseq,i_cseq
      integer l_lsgm,i_lsgm,l_csgm,i_csgm
      integer l_lssb,i_lssb,l_llnk,i_llnk,l_clnk,i_clnk
      integer l_latm,i_latm,l_catm,i_catm,l_xatm,i_xatm,l_qatm,i_qatm
      integer l_bnd,i_bnd,l_ang,i_ang,l_dih,i_dih,l_imp,i_imp
c
      character*255 filpdb,filpar,filseq,filtmp
      character*255 dir_s,dir_x,dir_u,dir_t
      character*80 ffield
      character*3 slvnam
      character*80 source
      character*10 slvmdl
c
      integer lfnpdb,lfnout,lfnfrg,lfnpar,lfnseq,lfnsgm,lfntmp
      integer igeom,numslv
c
      if(lprint.ge.5) write(lfnout,1000)
 1000 format('SEQUENCE GENERATION')
c
      if(source(1:3).eq.'pdb') then
      if(.not.pre_pdbsiz(lfnpdb,filpdb,lfnout,lprint,0,
     + natm,nseq,nssb,nlnk,slvnam,numslv))
     + call errquit('pre_pdbsiz failed',9999)
      else
      if(.not.geom_create(igeom,source))
     + call errquit('pre_mkseq: geom_create failed',0)
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('pre_mkseq: geom_rtdb_load failed',0)
      if(.not.geom_ncent(igeom,natm))
     + call errquit('pre_mkseq: geom_ncent failed',0)
      nseq=1
      nssb=0
      nlnk=0
      endif
c
c     set dimensions
c     --------------
c
      matm=natm
      mseq=nseq+1
      mssb=nssb
      mlnk=nlnk+100
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
      if(.not.ma_push_get(mt_int,5*matm,"latm",l_latm,i_latm))
     + call errquit('Memory allocation failed for latm',9999)
c
c     char*6  catm(1,matm) : atom name on pdb
c                  2       : atom name on top
c                  3       : atom type
c
      if(.not.ma_push_get(mt_byte,18*matm,"catm",l_catm,i_catm))
     + call errquit('Memory allocation failed for catm',9999)
c
c     real*8  xatm(1:3,matm) : atomic coordinates in nm
c
      if(.not.ma_push_get(mt_dbl,3*matm,"xatm",l_xatm,i_xatm))
     + call errquit('Memory allocation failed for xatm',9999)
c
c     real*8 qatm(matm) : partial atomic charges
c
      if(.not.ma_push_get(mt_dbl,matm,"qatm",l_qatm,i_qatm))
     + call errquit('Memory allocation failed for qatm',9999)
c
c     allocate memory for secondary segment links
c     -------------------------------------------
c
c     integer lssb(1,mssb) : link segment 1
c                  2       : link segment 2
c
      if(.not.ma_push_get(mt_int,2*mssb,"lssb",l_lssb,i_lssb))
     + call errquit('Memory allocation failed for lssb',9999)
c
c     integer llnk(1,mlnk) : link segment 1
c                  2       : link segment 1 atom number
c                  3       : link segment 2
c                  4       : link segment 2 atom number
c     char*4  clnk(1,mlnk) : link segment 1 atom name
c                  2       : link segment 2 atom name
c
      if(.not.ma_push_get(mt_int,2*mlnk,"llnk",l_llnk,i_llnk))
     + call errquit('Memory allocation failed for llnk',9999)
      if(.not.ma_push_get(mt_byte,8*mlnk,"clnk",l_clnk,i_clnk))
     + call errquit('Memory allocation failed for clnk',9999)
c
c     allocate memory for sequence
c     ----------------------------
c
c     integer lseq(1,mseq) : number of sequence entry
c                  2       : index into unique segment list
c                  3       : index of first atom of sequence entry
c                  4       : molecule number
c
      if(.not.ma_push_get(mt_int,4*mseq,"lseq",l_lseq,i_lseq))
     + call errquit('Memory allocation failed for lseq',9999)
c
c     char*5  cseq(1,mseq) : name of sequence entry on pdb
c                  2       : name of sequence entry on top
c
      if(.not.ma_push_get(mt_byte,10*mseq,"cseq",l_cseq,i_cseq))
     + call errquit('Memory allocation failed for cseq',9999)
c
c
c     read the sequence and atom list from the pdb file
c     -------------------------------------------------
c
      if(source(1:3).eq.'pdb') then
      if(.not.pre_pdbseq(lfnpdb,filpdb,lfnout,lprint,0,slvnam,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_xatm),matm,natm,
     + int_mb(i_lssb),mssb,nssb,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk))
     + call errquit('pre_pdbseq failed',9999)
      else
      if(.not.pre_rtdbsq(irtdb,igeom,lfnout,lprint,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_xatm),matm,natm))
     + call errquit('pre_rtdbsq failed',9999)
      endif
c
c     comply with iupac atom names
c     ----------------------------
c
      if(.not.pre_namiup(lfnout,lprint,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_latm),byte_mb(i_catm),matm,natm))
     + call errquit('pre_namiup failed',9999)
c
c     assign segments if not obtained from pdb
c     ----------------------------------------
c
c     dallocate memory for sequence
c     -----------------------------
c
c      if(.not.ma_pop_stack(l_cseq))
c     + call errquit('Memory deallocation failed for cseq',9999)
c      if(.not.ma_pop_stack(l_lseq))
c     + call errquit('Memory deallocation failed for lseq',9999)
c
c     reallocate memory for sequence
c     ------------------------------
c
c     integer lseq(1,mseq) : number of sequence entry
c                  2       : index into unique segment list
c                  3       : index of first atom of sequence entry
c                  4       : molecule number
c
c      if(.not.ma_push_get(mt_int,4*mseq,"lseq",l_lseq,i_lseq))
c     + call errquit('Memory allocation failed for lseq',9999)
c
c     char*5  cseq(1,mseq) : name of sequence entry on pdb
c                  2       : name of sequence entry on top
c
c      if(.not.ma_push_get(mt_byte,10*mseq,"cseq",l_cseq,i_cseq))
c     + call errquit('Memory allocation failed for cseq',9999)
c
c     check segment links
c     -------------------
c
      if(.not.pre_links(lfnout,lprint,
     + dbl_mb(i_xatm),byte_mb(i_catm),int_mb(i_latm),
     + matm,int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,int_mb(i_llnk),
     + byte_mb(i_clnk),mlnk,nlnk))
     + call errquit('pre_links failed',9999)
c
c     comply with force field segment names
c     -------------------------------------
c
      if(ffield(1:5).eq.'amber') then
      if(.not.pre_namamb(lfnout,lprint,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_latm),byte_mb(i_catm),matm,natm))
     + call errquit('pre_namamb failed',9999)
      endif
c
c     write atom names to temporary file
c     ----------------------------------
c
      if(.not.pre_wrtnam(lfntmp,filtmp,byte_mb(i_catm),matm,natm))
     + call errquit('pre_wrtnam failed',9999)
c
c     get number of unique segments
c     -----------------------------
c
      if(.not.pre_sgmnum(lfnout,lprint,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,nsgm))
     + call errquit('pre_sgmnum failed',9999)
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
      if(.not.ma_push_get(mt_int,3*msgm,"lsgm",l_lsgm,i_lsgm))
     + call errquit('Memory allocation failed for lsgm',9999)
c
c     char*5  csgm(msgm)   : segment names
c
      if(.not.ma_push_get(mt_byte,5*msgm,"csgm",l_csgm,i_csgm))
     + call errquit('Memory allocation failed for csgm',9999)
c
c     scan the database directories for segment files
c     -----------------------------------------------
c
      if(.not.pre_sgmndx(lfnfrg,dir_s,dir_x,dir_u,dir_t,lfnout,lprint,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_lsgm),byte_mb(i_csgm),msgm,nsgm)) then
c
c     set dimensions for bonds, angles, torsions and out-of-plane
c     -----------------------------------------------------------
c
      mbnd=2*matm
      mang=2*matm
      mdih=3*matm
      mimp=matm
c
c     allocate memory for the bonded lists
c     ------------------------------------
c
c     integer lbnd(1:2,mbnd) : bond indices
c
      if(.not.ma_push_get(mt_int,2*mbnd,"bnd",l_bnd,i_bnd))
     + call errquit('Memory allocation failed for bnd',9999)
c
c     integer lang(1:3,mang) : angle indices
c
      if(.not.ma_push_get(mt_int,3*mang,"ang",l_ang,i_ang))
     + call errquit('Memory allocation failed for ang',9999)
c
c     integer ldih(1:4,mdih) : torsion indices
c
      if(.not.ma_push_get(mt_int,4*mdih,"dih",l_dih,i_dih))
     + call errquit('Memory allocation failed for dih',9999)
c
c     integer limp(1:4,mimp) : improper torsion indices
c
      if(.not.ma_push_get(mt_int,4*mimp,"imp",l_imp,i_imp))
     + call errquit('Memory allocation failed for imp',9999)
c
c     generate fragment files for missing segments in dir_t
c
      if(.not.pre_mkfrg(irtdb,source,lfnfrg,lfnpar,filpar,
     + dir_s,dir_x,dir_u,dir_t,lfnout,lprint,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_lsgm),byte_mb(i_csgm),msgm,nsgm,int_mb(i_latm),
     + byte_mb(i_catm),dbl_mb(i_xatm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_bnd),mbnd,int_mb(i_ang),mang,int_mb(i_dih),mdih,
     + int_mb(i_imp),mimp))
     + call errquit('pre_mkfrg failed',9999)
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_imp))
     + call errquit('Memory deallocation failed for imp',9999)
      if(.not.ma_pop_stack(l_dih))
     + call errquit('Memory deallocation failed for dih',9999)
      if(.not.ma_pop_stack(l_ang))
     + call errquit('Memory deallocation failed for ang',9999)
      if(.not.ma_pop_stack(l_bnd))
     + call errquit('Memory deallocation failed for bnd',9999)
      endif
c
c     generate segment files
c     ----------------------
c
      if(.not.pre_gensgm(lfnout,lprint,lfnfrg,lfnsgm,
     + ffield,dir_s,dir_x,dir_u,dir_t,
     + int_mb(i_lsgm),byte_mb(i_csgm),msgm,nsgm))
     + call errquit('pre_gensgm failed',9999)
c
      if(lprint.ge.5) write(lfnout,1001)
 1001 format('SEGMENTS FILES GENERATED')
c
c     write the sequence file
c     -----------------------
c
      if(.not.pre_wrtseq(lfnseq,filseq,
     + byte_mb(i_cseq),int_mb(i_lseq),mseq,nseq,
     + int_mb(i_lssb),mssb,nssb,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk,slvmdl))
     + call errquit('pre_wrtseq failed',9999)
c
      if(lprint.ge.5) write(lfnout,1002)
 1002 format('SEQUENCE FILE GENERATEd')
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_csgm))
     + call errquit('Memory deallocation failed for csgm',9999)
      if(.not.ma_pop_stack(l_lsgm))
     + call errquit('Memory deallocation failed for lsgm',9999)
      if(.not.ma_pop_stack(l_cseq))
     + call errquit('Memory deallocation failed for cseq',9999)
      if(.not.ma_pop_stack(l_lseq))
     + call errquit('Memory deallocation failed for lseq',9999)
      if(.not.ma_pop_stack(l_clnk))
     + call errquit('Memory deallocation failed for clnk',9999)
      if(.not.ma_pop_stack(l_llnk))
     + call errquit('Memory deallocation failed for llnk',9999)
      if(.not.ma_pop_stack(l_lssb))
     + call errquit('Memory deallocation failed for lssb',9999)
      if(.not.ma_pop_stack(l_qatm))
     + call errquit('Memory deallocation failed for qatm',9999)
      if(.not.ma_pop_stack(l_xatm))
     + call errquit('Memory deallocation failed for xatm',9999)
      if(.not.ma_pop_stack(l_catm))
     + call errquit('Memory deallocation failed for catm',9999)
      if(.not.ma_pop_stack(l_latm))
     + call errquit('Memory deallocation failed for latm',9999)
c
      pre_mkseq=.true.
      return
      end
