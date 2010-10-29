      character*2 function pre_atnam(number)
c
c $Id$
c
      implicit none
      integer number
c
      character*2 aname(0:105)
c
      data aname / '  ',
     + ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     + 'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     + 'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     + 'Ga','Ge','As','Se','Br',' R','Rb','Sr',' Y','Zr',
     + 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     + 'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     + 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     + 'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     + 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     + 'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     + 'Md','No','Lr','Rf','Ha' /
c
      pre_atnam=aname(number)
c
      return
      end
