      subroutine aor_write_roadata
     &   (luout, alfa, beta, gp, dipquad, dqpol)
      
      integer luout
      double precision alfa(3,3), beta(3,3), gp(3,3) 
      double precision dipquad(3,6), dqpol(3,3,3)

      write (luout,'(/1x,a)') 'roa begin'
      write (luout,'(/1x,a)') 'alpha:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', alfa(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', alfa(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', alfa(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', alfa(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', alfa(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', alfa(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', alfa(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', alfa(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', alfa(3,3)
      
      write (luout,'(/1x,a)') 'beta:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', beta(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', beta(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', beta(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', beta(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', beta(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', beta(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', beta(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', beta(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', beta(3,3)

      write (luout,'(/1x,a)') 'G-prime = -omega beta:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', gp(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', gp(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', gp(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', gp(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', gp(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', gp(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', gp(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', gp(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', gp(3,3)
      
c     second index of dipquadre:
c     XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6
      write (luout,'(/1x,a)') 'dipole-quadrupole (Cartesian):'
      write (luout,'(1x,a,t10,e20.13)') 'XXX', dipquad(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XXY', dipquad(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XXZ', dipquad(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'XYY', dipquad(1,4)
      write (luout,'(1x,a,t10,e20.13)') 'XYZ', dipquad(1,5)
      write (luout,'(1x,a,t10,e20.13)') 'XZZ', dipquad(1,6)

      write (luout,'(1x,a,t10,e20.13)') 'YXX', dipquad(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YXY', dipquad(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YXZ', dipquad(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'YYY', dipquad(2,4)
      write (luout,'(1x,a,t10,e20.13)') 'YYZ', dipquad(2,5)
      write (luout,'(1x,a,t10,e20.13)') 'YZZ', dipquad(2,6)

      write (luout,'(1x,a,t10,e20.13)') 'ZXX', dipquad(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZXY', dipquad(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZXZ', dipquad(3,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZYY', dipquad(3,4)
      write (luout,'(1x,a,t10,e20.13)') 'ZYZ', dipquad(3,5)
      write (luout,'(1x,a,t10,e20.13)') 'ZZZ', dipquad(3,6)

      write (luout,'(/1x,a)') 'dipole-quadrupole (traceless):'
      write (luout,'(1x,a,t10,e20.13)') 'XXX', dqpol(1,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XXY', dqpol(1,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XXZ', dqpol(1,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'XYY', dqpol(1,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'XYZ', dqpol(1,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'XZZ', dqpol(1,3,3)

      write (luout,'(1x,a,t10,e20.13)') 'YXX', dqpol(2,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'YXY', dqpol(2,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'YXZ', dqpol(2,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YYY', dqpol(2,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YYZ', dqpol(2,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'YZZ', dqpol(2,3,3)

      write (luout,'(1x,a,t10,e20.13)') 'ZXX', dqpol(3,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZXY', dqpol(3,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZXZ', dqpol(3,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZYY', dqpol(3,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZYZ', dqpol(3,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZZZ', dqpol(3,3,3)

      write (luout,'(/1x,a/)') 'roa end'

      return
      end


c F. Krausbeck, January 2016, ETH Zurich:
c new write function for Re+Im Parts for Resonance VROA

      subroutine aor_write_roadata2
     &   (luout, alfare, alfaim, betare, gpre, 
     &     betaim, gpim, dipquadre, dipquadim, dqpol)
      
      integer luout
      double precision alfare(3,3), betare(3,3) 
      double precision gpre(3,3), gpim(3,3) 
      double precision alfaim(3,3), betaim(3,3) 
      double precision dipquadre(3,6), dqpol(3,3,3)
      double precision dipquadim(3,6)

      write (luout,'(/1x,a)') 'roa begin'
      write (luout,'(/1x,a)') 'alpha real:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', alfare(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', alfare(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', alfare(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', alfare(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', alfare(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', alfare(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', alfare(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', alfare(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', alfare(3,3)
      
      write (luout,'(/1x,a)') 'alpha im:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', alfaim(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', alfaim(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', alfaim(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', alfaim(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', alfaim(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', alfaim(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', alfaim(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', alfaim(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', alfaim(3,3)

      write (luout,'(/1x,a)') 'beta real:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', betare(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', betare(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', betare(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', betare(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', betare(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', betare(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', betare(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', betare(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', betare(3,3)

      write (luout,'(/1x,a)') 'G-prime = -omega beta(real):'
      write (luout,'(1x,a,t10,e20.13)') 'XX', gpre(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', gpre(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', gpre(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', gpre(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', gpre(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', gpre(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', gpre(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', gpre(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', gpre(3,3)

      write (luout,'(/1x,a)') 'beta im:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', betaim(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', betaim(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', betaim(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', betaim(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', betaim(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', betaim(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', betaim(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', betaim(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', betaim(3,3)

      write (luout,'(/1x,a)') 'G-prime = -omega beta(imag):'
      write (luout,'(1x,a,t10,e20.13)') 'XX', gpim(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', gpim(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', gpim(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', gpim(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', gpim(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', gpim(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', gpim(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', gpim(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', gpim(3,3)
      
c     second index of dipquadre:
c     XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6
      write (luout,'(/1x,a)') 'dipole-quadrupole real (Cartesian):'
      write (luout,'(1x,a,t10,e20.13)') 'XXX', dipquadre(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XXY', dipquadre(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XXZ', dipquadre(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'XYY', dipquadre(1,4)
      write (luout,'(1x,a,t10,e20.13)') 'XYZ', dipquadre(1,5)
      write (luout,'(1x,a,t10,e20.13)') 'XZZ', dipquadre(1,6)

      write (luout,'(1x,a,t10,e20.13)') 'YXX', dipquadre(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YXY', dipquadre(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YXZ', dipquadre(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'YYY', dipquadre(2,4)
      write (luout,'(1x,a,t10,e20.13)') 'YYZ', dipquadre(2,5)
      write (luout,'(1x,a,t10,e20.13)') 'YZZ', dipquadre(2,6)

      write (luout,'(1x,a,t10,e20.13)') 'ZXX', dipquadre(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZXY', dipquadre(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZXZ', dipquadre(3,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZYY', dipquadre(3,4)
      write (luout,'(1x,a,t10,e20.13)') 'ZYZ', dipquadre(3,5)
      write (luout,'(1x,a,t10,e20.13)') 'ZZZ', dipquadre(3,6)

c     second index of dipquadim:
c     XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6
      write (luout,'(/1x,a)') 'dipole-quadrupole imag (Cartesian):'
      write (luout,'(1x,a,t10,e20.13)') 'XXX', dipquadim(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XXY', dipquadim(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XXZ', dipquadim(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'XYY', dipquadim(1,4)
      write (luout,'(1x,a,t10,e20.13)') 'XYZ', dipquadim(1,5)
      write (luout,'(1x,a,t10,e20.13)') 'XZZ', dipquadim(1,6)

      write (luout,'(1x,a,t10,e20.13)') 'YXX', dipquadim(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YXY', dipquadim(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YXZ', dipquadim(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'YYY', dipquadim(2,4)
      write (luout,'(1x,a,t10,e20.13)') 'YYZ', dipquadim(2,5)
      write (luout,'(1x,a,t10,e20.13)') 'YZZ', dipquadim(2,6)

      write (luout,'(1x,a,t10,e20.13)') 'ZXX', dipquadim(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZXY', dipquadim(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZXZ', dipquadim(3,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZYY', dipquadim(3,4)
      write (luout,'(1x,a,t10,e20.13)') 'ZYZ', dipquadim(3,5)
      write (luout,'(1x,a,t10,e20.13)') 'ZZZ', dipquadim(3,6)

      write (luout,'(/1x,a)') 'dipole-quadrupole (traceless):'
      write (luout,'(1x,a,t10,e20.13)') 'XXX', dqpol(1,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XXY', dqpol(1,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XXZ', dqpol(1,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'XYY', dqpol(1,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'XYZ', dqpol(1,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'XZZ', dqpol(1,3,3)

      write (luout,'(1x,a,t10,e20.13)') 'YXX', dqpol(2,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'YXY', dqpol(2,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'YXZ', dqpol(2,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YYY', dqpol(2,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YYZ', dqpol(2,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'YZZ', dqpol(2,3,3)

      write (luout,'(1x,a,t10,e20.13)') 'ZXX', dqpol(3,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZXY', dqpol(3,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZXZ', dqpol(3,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZYY', dqpol(3,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZYZ', dqpol(3,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZZZ', dqpol(3,3,3)

      write (luout,'(/1x,a/)') 'roa end'

      return
      end
c $Id$
