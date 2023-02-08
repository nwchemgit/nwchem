      subroutine  nwpwxc_vdw3_dftd3()
      call errquit(' dft-d3 functionality missing',0,0)
      end subroutine
#if 0      
      subroutine nwpwxc_printoptions
      end subroutine nwpwxc_printoptions
      subroutine nwpwxc_setfuncpar(func,version,TZ,s6,rs6,s18,rs18,alp)
      end subroutine nwpwxc_setfuncpar
      subroutine nwpwxc_rdpar(dtmp,version,s6,s18,rs6,rs18,alp)
      end subroutine nwpwxc_rdpar
      subroutine nwpwxc_edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
      end subroutine nwpwxc_edisp
      subroutine nwpwxc_adisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
      end subroutine nwpwxc_adisp
c      subroutine nwpwxc_gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
      subroutine nwpwxc_gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
      end subroutine nwpwxc_gdisp
      subroutine nwpwxc_get_dC6_dCNij(maxc,max_elem,c6ab,mxci,mxcj,
      end subroutine nwpwxc_get_dC6_dCNij
      subroutine nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
      end subroutine nwpwxc_getc6
      subroutine nwpwxc_ncoord(natoms,rcov,iz,xyz,cn,cn_thr)
      end subroutine nwpwxc_ncoord
      subroutine nwpwxc_loadc6(fname,maxc,max_elem,c6ab,maxci)
      end subroutine nwpwxc_loadc6
      subroutine nwpwxc_limit(iat,jat,iadr,jadr)
      end subroutine nwpwxc_limit
#endif
