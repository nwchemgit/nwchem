      subroutine argos_cafe_num2nam(i,a)
c
      implicit none
c
      integer i
      character*16 a
c
      character*2 symbol(103)
c
      data symbol/
     + 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     + 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     + 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     + 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     + 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     + 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     + 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     + 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     + 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     + 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     + 'Md','No','Lr'/
c
      a='     Bq         '
      if(i.gt.0.and.i.lt.104) a(6:7)=symbol(i)
      if(a(7:7).eq.' ') a(6:7)=' '//a(6:6)
c
      return
      end
c $Id$
