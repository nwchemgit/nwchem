      character*8 luciasymname
      character*255 luciascfmo
      integer luciageom, luciabasis
      integer ipnt(mxpirr),ipnt2(mxpirr)
      integer ioff(mxpirr),icount(mxpirr)
      integer kirs(mxporb),reord(mxporb)
      double precision nwchem_sc(mxorb,mxpirr)
      integer luciaao,luciamo(2)
      COMMON/NWCHEM_DATAI/luciageom,luciabasis,
     &                    ipnt,ipnt2,ioff,icouno,kirs,
     &                    reord,luciaao,luciamo
      COMMON/NWCHEM_DATAD/nwchem_sc
      COMMON/NWCHEM_DATAC/luciascfmo,luciasymname
