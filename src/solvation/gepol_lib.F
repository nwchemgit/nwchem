c       =============================================================
c       GEPOL Library - Modified Version
c       Rafael de P. Soares - https://ufrgs.br/lvpp
c       March 2025
c
c       Modification Guidelines:
c       - All modifications are in **lowercase**, while original code
c         remains in UPPERCASE.
c       - The same convention applies to comments: original comments 
c         start with "C", while modified ones start with "c".
c
c       Summary of Modifications:
c       - Added a subroutine to compute the cavity and return results 
c         via arguments, eliminating the need for file I/O.
c       - Added a flag to disable the regrouping of secondary 
c         triangles, allowing for finer tesserae output.
c
c       Further cleanup and improvements may be possible.
c       Original source available at:
c       https://ftp.uv.es/pub/src/gepol/gepol93/
c
c       =============================================================
C       =============================================================
C       GEPOL93 (GEometria POLihedro,1993)
C       Version 8
C       Burjassot September 17, 1993.
C       =============================================================
C*** Written by
C
C       J.L. Pascual-Ahuir, E. Silla and I. Tunon
C       Departamento de Quimica Fisica.
C       Facultad de Quimica.
C       Universidad de Valencia.
C       C/Dr.Moliner 50.
C       Burjassot (Valencia) 46100
C       SPAIN
C       
C       Phone number: (34)-6-3864332
C       FAX   number: (34)-6-3864564
C       EMAIL   PASCUAL@EVALUN11
C                 SILLA@EVALUN11
C
C*** Grants
C
C       A grant has been provided by the Spanish Direccion General de 
C       Investigacion Cientifica Y Tecnica (DGICYT), project PS 90-0264.
C
C*** OLD VERSIONS
C
C     - GEPOL87 by J.L. Pascual-Ahuir, E. Silla, J. Tomasi and R. Bonacorsi
C     - GEPOL92 by J.L. Pascual-Ahuir, E. Silla and I. Tunon.
C 
C*** Modifications with respect to GEPOL92
C     The modifications affect mainly the computation of the Solvent-
C     Excluding Surface (ESURF). There are important changes in the 
C     algorithm of the computation of the ESURF. The set of parameters
C     used for that computation is different from the old versions.
C     The computation of the surfaces of large systems as proteins,
C     with the old versions, takes a lot of time. With GEPOL93 this 
C     problem has been solved. For a clear description of the new
C     algorithm see reference D.
C
C*** References 
C   (Basic references)
C
C    A- J.L.Pascual-Ahuir, E.Silla, J.Tomasi y R.Bonacorsi.
C       Electrostatic Interaction of Solute with a Continuum.
C       Improved Description of the Cavity and of the Surface Cavity
C       Bound Charge Distribution.
C       J. Comput. Chem.,8(1987)778-787
C
C    B- J. L. Pascual-Ahuir and E. Silla
C       GEPOL: An improved description of molecular surfaces. I. Building
C       the spherical surface set.
C       J. Comput. Chem.,11(1990)1047-1060
C
C    C- E. Silla, I. Tunon and J. L. Pascual-Ahuir
C       GEPOL: An improved description of molecular surfaces. II. Computing
C       the molecular area and volume.
C       J. Comput. Chem.,12(1991)1077-1088
C
C    D- J. L. Pascual-Ahuir, I Tunon and E. Silla
C       GEPOL: An improved description of molecular surfaces. III. A New
C       algorithm for the computation of the Solvent-Excluding Surface.
C       (explain the algorithm used in GEPOL93)
C       To be submited to J. Comput. Chem. during 1993
C
C   (References with examples of use)
C
C    
C    E- J. L. Pascual-Ahuir and E. Silla
C       GEPOL: A method to calculate the envelope surface. Computation
C       of changes in conformational area and volume of n-octanol
C       Quantum Chemistry - Basics Aspects, Actual Trends. Ed. R. Carbo.
C       Studies in Physical and Theoretical Chemistry,62(1989)597-603
C       Elsevier Science Publishers B. V., Amsterdam 1989. p
C 
C    F- J. L. Pascual-Ahuir, J. Andres and E. Silla.
C       Calculations of the relative basicities of methylamines in 
C       solution
C       Chem. Phys. Lett., 169(1990)297-300
C
C    G- E. Silla, F. Villar, O. Nilsson, J.L. Pascual-Ahuir and O. Tapia.
C       Molecular volumes and surfaces of biomacromolecules via GEPOL: A
C       fast and efficient algorithm
C       J. Mol. Graphics,8(1990)168-172
C
C    H- F. M. Floris, J. Tomasi, J.L. Pascual-Ahuir
C       Dispersion and repulsion contributions to the solvation energy:
C       Refinements to a simple computational model in the continuum
C       Approximation.
C       J. Comput. Chem.,12(1991)784-791
C
C    I- I. Tunon, E. Silla and J. L. Pascual-Ahuir
C       Molecular surface area and hydrophobic effect
C       Protein  Eng.,5(1992)715-716
C
C    J- I. Tunon, E. Silla and J. L. Pascual-Ahuir
C       Continuum-uniform approach calculations of the solubility
C       of hydrocarbons in water
C       Chem. Phys. Lett., 203(1993)289-294
C
C    K- I. Tunon, E. Silla and J. L. Pascual-Ahuir
C       Theoretical study of the inversion of the alcohol acidity 
C       scale in aqueous solution. Toward an interpretation of the
C       acid-base behavior of organic compounds in solution.
C       J. Am. Chem. Soc. 115(1993)2226-2230 
C
C*** Aim
C
C       This program computes the surface for a molecule as a distribution
C       of points and calculates its area and the volume enclosed. Each
C       point represents a piece of the surface named a tesserae. 
C
C       Three kinds of envelope surfaces can be computed:
C
C         - THE VAN DER WAALS MOLECULAR SURFACE. This is the envelope
C       surface of a set of intersecting spheres with given atomic radii 
C       centered on the nuclei of selected atoms of the molecule.
C
C         - THE ACCESSIBLE MOLECULAR SURFACE. This was defined by
C       B.Lee and F.M.Richards (J.MOL.BIOL.55(1971)379-400). It is the
C       surface defined by the center of the solvent, considered as a 
C       rigid sphere (probe sphere), when it rolls around the van der 
C       Waals surface. 
C
C         - THE SOLVENT-EXCLUDING SURFACE. This is the surface envelope 
C       of the volume excluded to the solvent, considered as a rigid 
C       sphere (probe sphere), when it rolls around the van der Waals 
C       surface. It was defined initially by F. M. Richards and named 
C       as the MOLECULAR SURFACE (Ann. Rev. Biophys. Bioeng.,6 (1977)
C       151-176) 
C       GEPOL constructs the SOLVENT-EXCLUDING SURFACE creating a 
C       set of new spheres located among the original spheres defined in 
C       the input geometry.
C        
C        
C*** Computacional especifications
C
C       This program is written in FORTRAN 77 and has been run on 
C        VAX (VMS)
C        SILICON GRAPHICS (UNIX)
C        IBM RISC (AIX)       
C
C*** Input
C       
C       GENERAL INPUT (FOR005)
C       ----------------------
C       This file will have the commands and parameters that control 
C       the program.
C       Every record is a command and it is divided in two fields:The first
C       five characters of the field form a KEYWORD.The second field is 
C       reserved for parameters.
C       You do not need to put the commands in a given order.
C       There are two types of KEYWORD depending on whether a parameter 
C       is needed in the second field of the record or not. Those that need
C       a parameter have the symbol equal(=) included in the KEYWORD.
C       Every command has a default defined within the program. The default
C       values for the parameters are selected to achieve a good compromise
C       between time and accuracy.
C 
C       Each KEYWORD is explained, as follows:
C 
C  
C  TITL=     This Keyword denotes that the title of the calculation goes 
C            into the second field of the record. You may use up to 20 
C            records as titles, but each record should start with the 
C            keyword TITL=.If you do not use this command the program 
C            will provide a title.
C
C  WSURF     The van der Waals molecular surface is calculated.
C
C  ASURF     The accessible molecular surface is calculated.
C
C  ESURF     The Solvent-Excluding Surface (Molecular Surface) is 
C            calculated.
C
C            ( if you do not put in any of the latter keywords
C              the program will calculate the van der Waals surface)
C
C  NDIV=     In the second field goes an integer that can take values 
C            between 1 and 5. It specifies the division level for the 
C            triangles on the surface. The accuracy of the calculation
C            improves as NDIV rises. The default for this command is 3.
C
C  OFAC=     This parameter will be used only if ESURF is computed.
C            In the second field goes a real number that can take values 
C            between 0.0 and 1.0. This parameter is the Overlapping FACtor.
C            The accuracy improves as the OFAC value increases.The default
C            value is 0.8
C
C  RMIN=     This parameter will be used only if ESURF is computed.
C            In the second field goes a real number that can take values 
C            larger than 0.0. This parameter is the radius of the smallest
C            sphere that can be created.The accuracy improves as the RMIN
C            value decreases.The default value is 0.50.
C
C            (OFAC,RMIN  are the parameters that control the creation
C             of new spheres)
C 
C  RSOL=     This parameter will be used only if ASURF or ESURF is computed.
C            In the second field goes a real number. It is the probe or 
C            solvent radius. The default value is 1.4
C
C  COOF=     In the second field goes the name of the file with the 
C            coordinates and radii.
C
C  SPHF=     The area of every sphere will be printed in the file
C            indicated in the second field. As default, the area of every
C            sphere will not be printed. For more information see later.
C            
C  VECF=     The coordinates of the points, the components of the vectors 
C            perpendicular to the surface and the areas of every tesserae, 
C            will be printed in the file indicated in the second field. 
C            As default the information will not be printed.The file
C            is written in binary format.
C  
C  DVEC=     Only if you have selected VECF=. The number in the second 
C            field gives the direction of the vectors and their module. 
C            If the number is positive they will go outward and if it 
C            is negative inward. The module of the vectors will be equal 
C            to the absolutevalue of the number. As default the value
C            is +1.
C
C  DISF=     This option prints, in the file indicated in the second field,
C            information needed for a program, developed in our laboratory,
C            named HELIOS, that displays the surfaces. As default the 
C            information is not printed. The file is written in binary 
C            format.
C
C  REDUC     Only if  SPHF= has been selected. The use of this command
C            will REDUCe the printing, giving only information about the
C            initial spheres. This option only has sense when the ESURF 
C            is computed.
C
C  ASSG1     Only if ESURF has been selected. This command assigns 
C            the tesserae of the new spheres to the initial ones.
C            Thus the value of the area of each initial sphere will be 
C            the sum of the area of the tesserae of itself plus the area
C            of the tesserae, from the new spheres, assigned to it. 
C            A tesserae is assigned to the sphere whose surface is closest
C            to the tesserae.
C            As DEFAULT no assignation is made.
C
C  LPRIN     This command produce a longer printing in the general output.
C            As default a shorter printing is made
C
C
C       COORDINATES FILE   
C       ----------------
C
C       The name of this file has been given in the general input after
C       the keyword COOF=. The format of the file is given next:
C
C     -- TITLE
C        You can put in as many records of titles as you want, but each 
C        record should start with the symbol *.
C    
C     -- NATOM (I8)
C        Number of atoms
C     
C     -- One record per atom (4F10.5,A33)
C        The Variables in each record are:
C      
C           XE        Coordinate X of the atom.
C           YE        Coordinate Y of the atom.
C           ZE        Coordinate Z of the atom.
C           RE        Radius of the sphere centered on the atom.
C           LABEL     Atom number, atom name, residue number... This 
C                     variable is not used in GEPOL93, it is only read to
C                     be written in the file SPHF= to identify each atom.
C         
C         OBSERVATIONS!!.
C         - With GEPOL93 we provide a program named PGEPOL that prepares
C           this file from others formats (PDB,GROMOS,CHARMM...)
C         - If you give a value to RE of 0.00000 this atom will not 
C           participate at all in the calculation.
C         - If you give a negative value to RE this atom will have a
C           GHOST sphere. Sometimes we are interested only in computing
C           the surface of part of a molecule, but taking into account 
C           the presence of the rest of the molecule. Thus, those atoms
C           that belong to the part of the molecule, which interests us,
C           will have positive RE and the rest of the atoms will have 
C           their proper radii but negative. 
C         - We are developing a set of programs to analyze the output of
C           GEPOL that makes use of the variable LABEL with a specific 
C           format and information. If you are interested in using these
C           programs the format and information of LABEL should be the
C           following:
C            (I8,1x,A4,I7,1x,A4,I3,1x,A4)
C            IAT  Atom number(order)
C            ATN  Atom name
C            IRE  Residue number(order)
C            REN  Residue name
C            ISE  Molecule number
C            SEN  Molecule name
C                     
C
C*** Output files
C  - The general output is connected to the fortran unit FOR006
C
C  - Information about the spheres
C    -----------------------------
C    The name of this file is given after the keyword SPHF=. It contains
C    information about every sphere.
C
C     -- TITLE (A80)
C          The same lines that were read in the coordinates file.
C 
C     -- NATOM (I8)
C        Number of atoms
C     
C     -- One record by atom or sphere (4F10.5,A33,F10.5)
C        The Variables in each record are:
C      
C           XE        Coordinate X of the atom(center).
C           YE        Coordinate Y of the atom.
C           ZE        Coordinate Z of the atom.
C           RE        Radius of the sphere centred on the atom.
C           LABEL     Atom number, atom name, residue number... 
C           AE        Surface Area of this atom or sphere
C    
C     OBSERVATIONS!!
C       - If you have used REDUC you will get only information about the 
C         original spheres.
C
C  -  Information about the tesserae forming the surface.
C    ---------------------------------------------------
C    The name of this file is given after the keyword VECF=. It contains
C    the coordinates of the center of each tesserae, the vectors 
C    perpendicular to the surface at this point, the area of each tesserae.
C    It is written in binary format. The records have the following
C    information:
C      Record 1.     LT      (Integer*4) is the number of title record
C      Next LT rec.  TITLE   (charac.*80) title
C      Rec. LT+2     NP      (Integer*4) Number of tesserae or points.
C      Next NP rec.          (3integer*4,7real*4) each record has the 
C                            following variables:
C
C                    ITO  each sphere can have up to 60 tesserae, ITO gives
C                         the order number of this tesserae among the 60.
C                    ISO  gives the sphere number where the point lies.
C                    ISA  gives the sphere number assigned to this point
C                    XP   coordinate x of the point(center of the tesserae)
C                    YP   coordinate y of the point
C                    ZP   coordinate z of the point
C                    AP   area of the tesserae
C                    XVEC component X of the vector
C                    YVEC component Y of the vector
C                    ZVEC component Z of the vector
C
C*** UNITS
C    The units of the area and volume will depend on the units of your
C    coordinates. Remember that XE,YE,ZE,RE,RSOL and RMIN should have 
C    the same units. Thus if you use Angstroms the area will appear 
C    in Angstroms**2 and the volume in Angstrom**3
C
C*** Other information
C    This information may be useful if you want to add something to the
C    program. 
C    MC = Maximum number of Centers (atoms, spheres)
C    MV = Maximum number of Vectors (points/tesserae)
C    IUSE(I) indicates the type of the center I
C          IUSE=1 Initial center with radius=0  (Set in Sub. READCOOR) 
C          IUSE=2 Sphere engulfed by another    (Set in Sub. CREA and BULK)
C          IUSE=3 Ghost sphere                  (Set in Sub. READCOOR)
C          IUSE=4 Sphere with final area 0      (Set in Sub. MZERO5  )
C          IUSE=5 Semi ghost sphere             (Set in Sub. SHELL   )
C          IUSE=6 Real sphere                   (Set in Sub. READCOOR)
C***Problems
C   If you have any problems or suggestions, please do not hesitate to 
C   contact us.
C***************************************************************************
      subroutine gepol(KSURF,REDU,PVEC,PSPH,PDIS,TIT,LT
     &                  ,RMIN,VECF,SPHF,DISF,OFAC,RD,NDIV,DVEC
     &                  ,ASS1,LPR,
c additional inputs here
     &                  natom,ratm,xyzatm,mxseg,
c outputs from here on
     &                  np,xyzseg,areaseg,atomseg,stot,vol,rgrp)
      IMPLICIT NONE
      
      LOGICAL ASS1
      LOGICAL GHOST
      LOGICAL LPR
      LOGICAL PDIS,PSPH,PVEC
      LOGICAL REDU

      INTEGER*2 IUSE

      INTEGER*4 ISA,ISO,ITO
      INTEGER*4 LT
      INTEGER*4 MC,MV
      INTEGER*4 NATOM,NCOR,NDIV,NP

      integer mxseg
      double precision radius(natom)
      double precision xyzatm(3,natom)
      double precision ratm(natom)
      double precision xyzseg(3,mxseg)
      double precision areaseg(mxseg)
      integer atomseg(mxseg)
      double precision dx,dy,dz
      double precision dseg, dsegmin
      integer i,j
      logical rgrp
       
      REAL*4 AP
      REAL*4 DVEC
      REAL*4 OFAC
      REAL*4 RD,RE,RMIN
      REAL*4 XC1,XE,XP
      REAL*4 YC1,YE,YP
      REAL*4 ZC1,ZE,ZP

      REAL*8 CV
      REAL*8 STOT,VOL

      CHARACTER*80 DISF,COOF,SPHF,TIT,VECF
      CHARACTER*33 LABEL
      CHARACTER*5  KSURF
      
      PARAMETER (MC=100000,MV=100000)
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)
    
      DIMENSION LABEL(MC),TIT(20)

c      WRITE(6,'(A )') ' ============================='
c      WRITE(6,'(A )') ' ====  GEPOL93 Version 8  ===='
c      WRITE(6,'(A/)') ' ============================='

C*****Read input****

C Set defaults and Read General Input
c      CALL READIN (KSURF,REDU,PVEC,PSPH,PDIS,TIT,COOF,LT
c     &            ,RMIN,VECF,SPHF,DISF,OFAC,RD,NDIV,DVEC
c     &            ,ASS1,LPR)

C Write general output
c      CALL GWRITE(KSURF,REDU,PVEC,PSPH,PDIS,TIT,COOF,LT
c     &           ,RMIN,VECF,SPHF,DISF,OFAC,RD,NDIV,DVEC
c     &           ,ASS1,LPR)
      
C 
C Read COOR File
c      CALL READCOOR(COOF,NATOM,LABEL,GHOST,LPR)
c Atom coordinates from arguments and not from file
       do i=1,natom
            xe(i)=xyzatm(1,i)
            ye(i)=xyzatm(2,i)
            ze(i)=xyzatm(3,i)

            re(i)=ratm(i)
            iuse(i)=6
       end do
c

      NCOR=NATOM
      
      IF(LPR)CALL PCOUNT(NATOM,NCOR,KSURF,LPR)

C*****Make tesselation for sphere of radius 1.0*********

      CALL TES(LPR)

      CALL DIVIDE(NDIV,LPR)

C*****Beginning  the calculations ******

C Create the new set of spheres if we are interested in the 
C Solvent-excluding Surface

      IF(KSURF.EQ.'ESURF')THEN

       IF(GHOST)CALL GSHELL(NCOR,LPR,RD)
      
       CALL BULK(NATOM,NCOR,NDIV,OFAC,RD,LPR)

       IF(LPR)CALL PCOUNT(NATOM,NCOR,KSURF,LPR)

       CALL CLEAN5(NATOM,NCOR,NDIV,LPR)

       IF(LPR)CALL PCOUNT(NATOM,NCOR,KSURF,LPR)

       CALL CREA(NATOM,NCOR,RMIN,OFAC,RD,LPR,NDIV) 

       IF(LPR)CALL PCOUNT(NATOM,NCOR,KSURF,LPR)

       CALL CLEAN5(NATOM,NCOR,NDIV,LPR)

c       WRITE(15,'(/A)')' About final set of coordinates'
c       WRITE(15,'( A)')' ------------------------------'
       CALL PCOUNT(NATOM,NCOR,KSURF,LPR)

      END IF 

C
C  Compute the surface

      IF(KSURF.EQ.'ASURF')CALL SUM(NCOR,RD,'SUMA',LPR)

      CALL GEOCAV(NCOR,NP,NDIV,GHOST,LPR,rgrp)

      IF(KSURF.EQ.'ESURF')THEN

         IF(ASS1)CALL ASSIGN1(NP,NATOM,GHOST,LPR)

      END IF

C
C Compute the area and volume
      CALL VOLARE(NP,LPR,GHOST,STOT,VOL)

C
C Print several files
      IF(PSPH)CALL PRISPH(NATOM,NCOR,NP,REDU,SPHF,LABEL,TIT,LT,LPR)

      IF(PVEC)CALL PRIVEC(NP,DVEC,TIT,LT,VECF,LPR)

      IF(PDIS)CALL PRIDIS(NP,DISF,TIT,LT,LPR)
     

c Copy back segment information
      do i=1,np
        xyzseg(1,i) = xp(i)
        xyzseg(2,i) = yp(i)
        xyzseg(3,i) = zp(i)

        areaseg(i) = ap(i)
        atomseg(i) = isa(i)
      end do

c      STOP
      return
      END
C

      SUBROUTINE PCOUNT(NATOM,NCOR,KSURF,LPR)
C     -------------------------------------------------------------------
C     This prints general counters
C     -------------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL FIRST
      LOGICAL LPR

      INTEGER*2 IUSE

      INTEGER*4 I,IUC
      INTEGER*4 J
      INTEGER*4 MC,MI
      INTEGER*4 NATOM,NCOR,NEWS

      REAL*4 RE
      REAL*4 XE
      REAL*4 YE
      REAL*4 ZE

      CHARACTER*5 KSURF

      PARAMETER (MC=100000,MI=6)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)

      DIMENSION IUC(MI)

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine Pcount '

      NEWS=NCOR-NATOM

c      IF(KSURF.EQ.'ESURF')THEN
c         WRITE(6,'(3(A,I7/))')
c     & ' Number of INITIAL coordinates          ',NATOM,
c     & ' Number of NEW     coordinates          ',NEWS,
c     & ' Number of TOTAL   coordinates          ',NCOR
c      ELSE
c         WRITE(6,'(A,I7/)')
c     & ' Number of TOTAL coordinates            ',NCOR
c      END IF      

      IF(LPR)THEN
      
        DO J=1,MI
        IUC(J)=0
        END DO

        DO I=1,NATOM
        J=IUSE(I)
        IUC(J)=IUC(J)+1
        END DO
     
        FIRST=.TRUE.
        DO J=1,MI
        IF(IUC(J).NE.0)THEN

          IF(FIRST)THEN
          FIRST=.FALSE.
          WRITE(6,'(A)')
     &          ' MORE INFORMATION ABOUT INITIAL SET OF COORDINATES'
          END IF

          WRITE(6,'(A,I2,A,I7)')
     &    ' Number of coordinates with IUSE ',J,'=',IUC(J)
        END IF
        END DO

        IF(KSURF.EQ.'ESURF')THEN

          DO J=1,MI
          IUC(J)=0
          END DO

          DO I=NATOM+1,NCOR
          J=IUSE(I)
          IUC(J)=IUC(J)+1
          END DO
     
          FIRST=.TRUE.
          DO J=1,MI
          IF(IUC(J).NE.0)THEN

              IF(FIRST)THEN
                FIRST=.FALSE.
                WRITE(6,'(A)')
     &        ' MORE INFORMATION ABOUT NEW SET OF COORDINATES'
              END IF

            WRITE(6,'(A,I2,A,I7)')
     &    ' Number of coordinates with IUSE ',J,'=',IUC(J)
          END IF
          END DO

        END IF

      IF(LPR)CALL GSTAT(NCOR)

      END IF
     
      RETURN
      END

      BLOCK DATA
C     -----------------------------------------------------------------
C     This has the information about the vertices
C     -----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 JVT1,JVT2
      COMMON/PENTA/JVT1(3,60),JVT2(3,4)
      DATA JVT1 /  1,   6,   2,  1,   2,   3,  1,   3,   4,
     &             1,   4,   5,  1,   5,   6,  7,   2,   6,
     &             8,   3,   2,  9,   4,   3, 10,   5,   4,
     &            11,   6,   5,  8,   2,  12,  9,   3,  13,
     &            10,   4,  14, 11,   5,  15,  7,   6,  16,
     &             7,  12,   2,  8,  13,   3,  9,  14,   4,
     &            10,  15,   5, 11,  16,   6,  8,  12,  18,
     &             9,  13,  19, 10,  14,  20, 11,  15,  21,
     &             7,  16,  17,  7,  17,  12,  8,  18,  13,
     &             9,  19,  14, 10,  20,  15, 11,  21,  16,
     &            22,  12,  17, 23,  13,  18, 24,  14,  19,
     &            25,  15,  20, 26,  16,  21, 22,  18,  12,
     &            23,  19,  13, 24,  20,  14, 25,  21,  15,
     &            26,  17,  16, 22,  17,  27, 23,  18,  28,
     &            24,  19,  29, 25,  20,  30, 26,  21,  31,
     &            22,  28,  18, 23,  29,  19, 24,  30,  20,
     &            25,  31,  21, 26,  27,  17, 22,  27,  28,
     &            23,  28,  29, 24,  29,  30, 25,  30,  31,
     &            26,  31,  27, 32,  28,  27, 32,  29,  28,
     &            32,  30,  29, 32,  31,  30, 32,  27,  31 /
      DATA JVT2 /  1,   5,   4,
     &             5,   2,   6,
     &             4,   6,   3,
     &             6,   4,   5 /
      END
C 
      SUBROUTINE TES(LPR)
C     --------------------------------------------------------------------
C     This computes the triangle vertex coordinates for a sphere of radius
C     one, projecting the pentakisdodecahedro onto it.
C     --------------------------------------------------------------------
      IMPLICIT NONE
      
      LOGICAL LPR

      INTEGER*4 I,II
      INTEGER*4 J

      REAL*4 XC1,YC1,ZC1

      REAL*8 CTH,CV
      REAL*8 FI,FIR,FIV
      REAL*8 STH
      REAL*8 TH,THEV
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      
      DIMENSION THEV(6),FIV(6)
      
      DATA THEV/0.6523581397843682D0,1.1071487177940905D0,
     $          1.3820857960113345D0,1.7595068575784587D0,
     $          2.0344439357957027D0,2.4892345138054251D0/
      DATA FIV/0.6283185307179586D0,0.0 D0              ,
     $         0.6283185307179586D0,0.0 D0              ,
     $         0.6283185307179586D0,0.0 D0              /
      DATA FIR/1.2566370614359173 D0/

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine TES '

      CV(1,1)=0.D0
      CV(1,2)=0.D0
      CV(1,3)=1.D0
      CV(32,1)=0.D0
      CV(32,2)=0.D0
      CV(32,3)=-1.D0
      II=1
      DO 520 I=1,6
      TH=THEV(I)
      FI=FIV(I)
      CTH=DCOS(TH)
      STH=DSIN(TH)
      DO 521 J=1,5
      FI=FI+FIR
      IF(J.EQ.1) FI=FIV(I)
      II=II+1
      CV(II,1)=STH*DCOS(FI)
      CV(II,2)=STH*DSIN(FI)
      CV(II,3)=CTH
  521 CONTINUE
  520 CONTINUE
      RETURN
      END
C
      SUBROUTINE DIVIDE(NDIV,LPR)
C     ---------------------------------------------------------------
C     This divides the initial 60 spherical triangles to the level
C     indicated by NDIV
C     ---------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL LPR

      INTEGER*4 IJ
      INTEGER*4 J,J2,J3,J4,J5,JVT1,JVT2
      INTEGER*4 NDIV,NV1,NV2,NV21,NV22,NV23,NV3,NV31,NV32,NV33
      INTEGER*4 NV41,NV42,NV43,NV51,NV52,NV53
     
      REAL*4 XC1,YC1,ZC1
   
      REAL*8 CC,CV,CVN2,CVN3,CVN4,CVN5
      REAL*8 FOUR     
      REAL*8 PI
      REAL*8 XV1,XV2,XV3      
      REAL*8 YV1,YV2,YV3     
      REAL*8 ZERO,ZV1,ZV2,ZV3     

      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/PENTA/JVT1(3,60),JVT2(3,4)

      DIMENSION CVN2(6,3),CVN3(6,3),CVN4(6,3),CVN5(6,3),CC(3)
 
      DATA ZERO/0.0D0/
      DATA PI/3.1415926535897932D0/
      
      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine DIVIDE '

      IJ=0
C*****Level 1****************************
      DO 10 J=1,60
      NV1=JVT1(1,J)
      NV2=JVT1(2,J)
      NV3=JVT1(3,J)
      XV1=CV(NV1,1)
      YV1=CV(NV1,2)
      ZV1=CV(NV1,3)
      XV2=CV(NV2,1)
      YV2=CV(NV2,2)
      ZV2=CV(NV2,3)
      XV3=CV(NV3,1)
      YV3=CV(NV3,2)
      ZV3=CV(NV3,3)
      IF(NDIV.GT.1) GO TO 20
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 10
C*****Level 2**********************
   20 CONTINUE
      CVN2(1,1)=XV1
      CVN2(1,2)=YV1
      CVN2(1,3)=ZV1
      CVN2(2,1)=XV2
      CVN2(2,2)=YV2
      CVN2(2,3)=ZV2
      CVN2(3,1)=XV3
      CVN2(3,2)=YV3
      CVN2(3,3)=ZV3
      CALL CALVER(CVN2)
      DO 21 J2=1,4
      NV21=JVT2(1,J2)
      NV22=JVT2(2,J2)
      NV23=JVT2(3,J2)
      XV1=CVN2(NV21,1)
      YV1=CVN2(NV21,2)
      ZV1=CVN2(NV21,3)
      XV2=CVN2(NV22,1)
      YV2=CVN2(NV22,2)
      ZV2=CVN2(NV22,3)
      XV3=CVN2(NV23,1)
      YV3=CVN2(NV23,2)
      ZV3=CVN2(NV23,3)
      IF(NDIV.GT.2) GO TO 30
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 21
C*****Level 3**********************************
   30 CONTINUE
      CVN3(1,1)=XV1
      CVN3(1,2)=YV1
      CVN3(1,3)=ZV1
      CVN3(2,1)=XV2
      CVN3(2,2)=YV2
      CVN3(2,3)=ZV2
      CVN3(3,1)=XV3
      CVN3(3,2)=YV3
      CVN3(3,3)=ZV3
      CALL CALVER(CVN3)
      DO 31 J3=1,4
      NV31=JVT2(1,J3)
      NV32=JVT2(2,J3)
      NV33=JVT2(3,J3)
      XV1=CVN3(NV31,1)
      YV1=CVN3(NV31,2)
      ZV1=CVN3(NV31,3)
      XV2=CVN3(NV32,1)
      YV2=CVN3(NV32,2)
      ZV2=CVN3(NV32,3)
      XV3=CVN3(NV33,1)
      YV3=CVN3(NV33,2)
      ZV3=CVN3(NV33,3)
      IF(NDIV.GT.3) GO TO 40
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 31
C*****Level 4******************************
   40 CONTINUE
      CVN4(1,1)=XV1
      CVN4(1,2)=YV1
      CVN4(1,3)=ZV1
      CVN4(2,1)=XV2
      CVN4(2,2)=YV2
      CVN4(2,3)=ZV2
      CVN4(3,1)=XV3
      CVN4(3,2)=YV3
      CVN4(3,3)=ZV3
      CALL CALVER(CVN4)
      DO 41 J4=1,4
      NV41=JVT2(1,J4)
      NV42=JVT2(2,J4)
      NV43=JVT2(3,J4)
      XV1=CVN4(NV41,1)
      YV1=CVN4(NV41,2)
      ZV1=CVN4(NV41,3)
      XV2=CVN4(NV42,1)
      YV2=CVN4(NV42,2)
      ZV2=CVN4(NV42,3)
      XV3=CVN4(NV43,1)
      YV3=CVN4(NV43,2)
      ZV3=CVN4(NV43,3)
      IF(NDIV.GT.4) GO TO 50
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 41
C*****Level 5*************************************
   50 CONTINUE
      CVN5(1,1)=XV1
      CVN5(1,2)=YV1
      CVN5(1,3)=ZV1
      CVN5(2,1)=XV2
      CVN5(2,2)=YV2
      CVN5(2,3)=ZV2
      CVN5(3,1)=XV3
      CVN5(3,2)=YV3
      CVN5(3,3)=ZV3
      CALL CALVER(CVN5)
      DO 51 J5=1,4
      NV51=JVT2(1,J5)
      NV52=JVT2(2,J5)
      NV53=JVT2(3,J5)
      XV1=CVN5(NV51,1)
      YV1=CVN5(NV51,2)
      ZV1=CVN5(NV51,3)
      XV2=CVN5(NV52,1)
      YV2=CVN5(NV52,2)
      ZV2=CVN5(NV52,3)
      XV3=CVN5(NV53,1)
      YV3=CVN5(NV53,2)
      ZV3=CVN5(NV53,3)
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
   51 CONTINUE
   41 CONTINUE
   31 CONTINUE
   21 CONTINUE
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALVER(CVN)
C     ---------------------------------------------------------------------
C     This divides one triangle into four.
C     ---------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 CVN,XXX,YYY,ZZZ,RRR,FC
      INTEGER*4 N,N1,N2
      DIMENSION CVN(6,3)
      DO 7 N=1,3
      N2=N+3
      N1=N-1
      IF(N.EQ.1)N1=3
      XXX=(CVN(N,1)+CVN(N1,1))/2.0D0
      YYY=(CVN(N,2)+CVN(N1,2))/2.0D0
      ZZZ=(CVN(N,3)+CVN(N1,3))/2.0D0
      RRR=SQRT(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      FC=1.0D0/RRR
      CVN(N2,1)=XXX*FC
      CVN(N2,2)=YYY*FC
      CVN(N2,3)=ZZZ*FC
    7 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
C     ---------------------------------------------------------------------
C     This computes the center of a spherical triangle.
C     ---------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC,XXX,YYY,ZZZ,RRR,FC
      DIMENSION CC(3)
      XXX=(XV1+XV2+XV3)/3.0D0
      YYY=(YV1+YV2+YV3)/3.0D0
      ZZZ=(ZV1+ZV2+ZV3)/3.0D0
      RRR=SQRT(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      FC=1.0D0/RRR
      CC(1)=XXX*FC
      CC(2)=YYY*FC
      CC(3)=ZZZ*FC
      RETURN
      END

C
      SUBROUTINE GEOCAV(NCOR,NP,NDIV,GHOST,LPR,rgrp)
C     ------------------------------------------------------------------
C     This computes the surface.
C     ------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL GHOST
      LOGICAL LPR
      logical rgrp

      INTEGER*2 IUSE

      INTEGER*4 I,IIJ,IJE,ITO,ISA,ISO
      INTEGER*4 J
      INTEGER*4 K 
      INTEGER*4 L1
      INTEGER*4 MC,MV
      INTEGER*4 N3,N4,NCOR,NDIV,NEJCI,NINF,NP,NSUP,NTRIAN,NTS
      INTEGER*4 UN3

      REAL*4 AP,ATP,ATS
      REAL*4 DD,DIJ2
      REAL*4 FC,FNDIV
      REAL*4 PI
      REAL*4 RE,REI,RREJ,RRR
      REAL*4 SRE2
      REAL*4 XC1,XE,XEI,XP,XPL,XSL,XSM
      REAL*4 YC1,YE,YEI,YP,YPL,YSL,YSM
      REAL*4 ZC1,ZE,ZEI,ZP,ZPL,ZSL,ZSM

      REAL*8 CV

      PARAMETER (MC=100000,MV=100000)
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

      DIMENSION IJE(MC)

      DATA PI/3.141593E0/

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine GEOCAV '

      NP=0

C begin
      NTRIAN=4**(NDIV-1)
      FNDIV=60*NTRIAN

C It selects one sphere
      DO 1 I=1,NCOR
      
      IF(IUSE(I).EQ.6) THEN
      REI=RE(I)
      XEI=XE(I)
      YEI=YE(I)
      ZEI=ZE(I)
      ATS=4.0E0*PI*REI*REI/FNDIV

      IIJ=0
C
C  It determines which spheres are linked to sphere I
      DO J=1,NCOR

      IF(IUSE(J).GE.3) THEN

      IF(I.NE.J)THEN

      DIJ2=(XEI-XE(J))*(XEI-XE(J))+
     &     (YEI-YE(J))*(YEI-YE(J))+
     &     (ZEI-ZE(J))*(ZEI-ZE(J))
      SRE2=(REI+RE(J))*(REI+RE(J))

          IF(DIJ2.LT.SRE2) THEN
            IIJ=IIJ+1
            IJE(IIJ)=J
            NEJCI=IIJ
          END IF

      END IF
      END IF
      END DO

C 
C It selects one main triangle.
      NSUP=0
      UN3=1
      DO 2 J=1,60
      XPL=0.0E0
      YPL=0.0E0
      ZPL=0.0E0
      NTS=0
      NINF=NSUP+1
      NSUP=NINF+NTRIAN-1
C
C It selects one secondary triangle.
      DO 3 K=NINF,NSUP
      XSL=XC1(K)*REI
      YSL=YC1(K)*REI
      ZSL=ZC1(K)*REI
      XSM=XSL+XEI
      YSM=YSL+YEI
      ZSM=ZSL+ZEI
C
C It fixes if the secondary triangle is inside or outside.
      L1=UN3
      DO N3=L1,NEJCI
      UN3=N3
      N4=IJE(N3)
      DD=(XSM-XE(N4))*(XSM-XE(N4))+
     &   (YSM-YE(N4))*(YSM-YE(N4))+
     &   (ZSM-ZE(N4))*(ZSM-ZE(N4))
      RREJ=RE(N4)*RE(N4)
      IF(DD.LT.RREJ) GO TO 3
      END DO

      DO N3=1,L1-1
      UN3=N3
      N4=IJE(N3)
      DD=(XSM-XE(N4))*(XSM-XE(N4))+
     &   (YSM-YE(N4))*(YSM-YE(N4))+
     &   (ZSM-ZE(N4))*(ZSM-ZE(N4))
      RREJ=RE(N4)*RE(N4)
      IF(DD.LT.RREJ) GO TO 3
      END DO
C
C It prepares the coordinates for the main triangle
      XPL=XPL+XSL
      YPL=YPL+YSL
      ZPL=ZPL+ZSL
      NTS=NTS+1
c
c In case we want all secondary triangles
      if(.not.rgrp) then
        NP=NP+1
        XP(NP)=XSL+XEI
        YP(NP)=YSL+YEI
        ZP(NP)=ZSL+ZEI
        AP(NP)=ATS
        ITO(NP)=J 
        ISO(NP)=I
        ISA(NP)=I
      end if
c
    3 CONTINUE
C
C It reduces the secondary triangles to the main triangle.
      IF(NTS.EQ.0)GO TO 2

c If the regroup of secondary triangles is on
      if(rgrp) then
        ATP=ATS*NTS
        XPL=XPL/NTS
        YPL=YPL/NTS
        ZPL=ZPL/NTS
        RRR=SQRT(XPL*XPL+YPL*YPL+ZPL*ZPL)
        FC=REI/RRR

        NP=NP+1
        XP(NP)=XPL*FC+XEI
        YP(NP)=YPL*FC+YEI
        ZP(NP)=ZPL*FC+ZEI
        AP(NP)=ATP
        ITO(NP)=J 
        ISO(NP)=I
        ISA(NP)=I
      end if

      
    2 CONTINUE

      END IF

    1 CONTINUE

      RETURN
      END
C
      
C 
      SUBROUTINE PRIVEC(NP,DVEC,TIT,LT,VECF,LPR)
C     --------------------------------------------------------------
C     This computes the vectors over the surface with module DVEC. 
C     Prints the coordinates of the centers of the tesserae, the 
C     area of tesserae and the components of the perpendicular vectors 
C     to the surface at each center of the tesserae.
C     --------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL LPR

      INTEGER*2 IUSE
   
      INTEGER*4 I,ISCH,ISA,ISO,ITO
      INTEGER*4 J
      INTEGER*4 LT
      INTEGER*4 MC,MV
      INTEGER*4 NP

      REAL*4 AP
      REAL*4 DVEC
      REAL*4 FDR
      REAL*4 RE
      REAL*4 XE,XP,XVEC
      REAL*4 YE,YP,YVEC
      REAL*4 ZE,ZP,ZVEC

      CHARACTER*80 TIT,VECF

      PARAMETER (MC=100000,MV=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)
     
      DIMENSION TIT(20)

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine PRIVEC '
      WRITE(6,'(A/1X,A)')' Printing vectors in file:',VECF

      OPEN(UNIT=8,FILE=VECF,FORM='UNFORMATTED',STATUS='UNKNOWN')

      WRITE(8)LT
      DO I=1,LT
      WRITE(8)TIT(I)
      END DO

      WRITE(8)NP

C   
      ISCH=0
      DO I=1,NP

      J=ISO(I)

      IF(ISCH.NE.J)THEN

         FDR=DVEC/RE(J)
         ISCH=J

      END IF
      
      XVEC=(XP(I)-XE(J))*FDR
      YVEC=(YP(I)-YE(J))*FDR
      ZVEC=(ZP(I)-ZE(J))*FDR

      WRITE(8)ITO(I),J,ISA(I),XP(I),YP(I),ZP(I),AP(I),XVEC,YVEC,ZVEC
     
      END DO

      CLOSE(8)
      RETURN
      END

C
      SUBROUTINE PRISPH(NATOM,NCOR,NP,REDU,SPHF,LABEL,TIT,LT,LPR)
C     ---------------------------------------------------------------
C     This prints the coordinates, radius, label and area of each 
C     sphere
C     ---------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL LPR
      LOGICAL REDU
   
      INTEGER*2 IUSE

      INTEGER*4 I,ISA,ISO,ITO
      INTEGER*4 J
      INTEGER*4 L,LT
      INTEGER*4 MC,MV
      INTEGER*4 NATOM,NCOR,NP

      REAL*4 AE,AP
      REAL*4 RE
      REAL*4 XE,XP
      REAL*4 YE,YP
      REAL*4 ZE,ZP

      CHARACTER*4 LNEW
      CHARACTER*33 LABEL
      CHARACTER*80 TIT,SPHF

      PARAMETER (MC=100000,MV=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

      DIMENSION AE(MC),LABEL(MC),TIT(20)

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine PRISPH '
      WRITE(6,'(A/1X,A)')' Printing Information of Spheres in file:',
     &                     SPHF

      DO J=1,NCOR
      AE(J)=0.0E0
      END DO

      DO I=1,NP
      J=ISA(I)
      AE(J)=AE(J)+AP(I)
      END DO
      
      OPEN(UNIT=7,FILE=SPHF,FORM='FORMATTED',STATUS='UNKNOWN')

      DO I=1,LT
      WRITE(7,'(2A)')'*',TIT(I)
      END DO

      IF(REDU)THEN
      WRITE(7,'(I8)')NATOM
      ELSE
      WRITE(7,'(I8)')NCOR
      END IF

      DO I=1,NATOM
      WRITE(7,'(4F10.5,A,F10.5)') 
     & XE(I),YE(I),ZE(I),RE(I),LABEL(I),AE(I)
      END DO

      IF(.NOT.REDU)THEN
         LNEW='NEW '
         L=0
         DO  I=NATOM+1,NCOR
         WRITE(7,'(4F10.5,I8,1X,A,I7,1X,A,I3,1X,A,F10.5)')
     &   XE(I),YE(I),ZE(I),RE(I),I,LNEW,L,LNEW,L,LNEW,AE(I)
         END DO
      END IF

      CLOSE(7)


      RETURN
      END

C  
      SUBROUTINE PRIDIS(NP,DISF,TIT,LT,LPR)
C     -------------------------------------------------------------
C     This prints the information needed by the Helios program to 
C     display the surface.
C     -------------------------------------------------------------

      IMPLICIT NONE
      
      LOGICAL LPR

      INTEGER*2 IUSE

      INTEGER*4 FFF
      INTEGER*4 I,ISA,ISO,ITO
      INTEGER*4 J,JV,JVT1,JVT2
      INTEGER*4 LE,LT
      INTEGER*4 MC,MV
      INTEGER*4 NP,NVX1,NVX2,NVX3
   
      REAL*4 AP
      REAL*4 RE,REI
      REAL*4 XP,XC1,XE,XEI
      REAL*4 YP,YC1,YE,YEI
      REAL*4 ZP,ZC1,ZE,ZEI
      REAL*4 VERTEX

      REAL*8 CV,PLUS1

      CHARACTER*80 TIT,DISF
      
      PARAMETER (MC=100000,MV=100000)
      
      COMMON/PENTA/JVT1(3,60),JVT2(3,4)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

      DIMENSION VERTEX(6,3),PLUS1(6,3),TIT(20)
      
      IF(LPR)WRITE(6,'(/A)')' ==> Start Subroutine PRIDIS '
      WRITE(6,'(A/1x,A)')' Printing Display information in file:',DISF

      OPEN(UNIT=15,FILE=DISF,FORM='UNFORMATTED',STATUS='UNKNOWN')

      WRITE(15)LT

      DO I=1,LT
      WRITE(15)TIT(I)
      END DO
      
      WRITE(15)NP

C  Here we calculate the vertices of the main triangles of each sphere 
C  which belongs to the molecular surface.

      DO I=1,NP

      J=ITO(I)
      LE=ISO(I)
      REI=RE(LE)
      XEI=XE(LE)
      YEI=YE(LE)
      ZEI=ZE(LE)

      NVX1=JVT1(1,J)
      NVX2=JVT1(2,J)
      NVX3=JVT1(3,J)

      DO FFF=1,3
          PLUS1(1,FFF)=CV(NVX1,FFF)
          PLUS1(2,FFF)=CV(NVX2,FFF)
          PLUS1(3,FFF)=CV(NVX3,FFF)
      END DO

      CALL CALVER(PLUS1)

      DO JV=1,6
           VERTEX(JV,1)=PLUS1(JV,1)*REI+XEI
           VERTEX(JV,2)=PLUS1(JV,2)*REI+YEI
           VERTEX(JV,3)=PLUS1(JV,3)*REI+ZEI
      END DO

      WRITE(15)ISO(I),ISA(I),((VERTEX(JV,FFF),FFF=1,3),JV=1,6)
      
      END DO

      CLOSE(15)
      RETURN
      END
C

C 
      SUBROUTINE ASSIGN1(NP,NATOM,GHOST,LPR)
C     --------------------------------------------------------------------
C     This assigns the tesserae of the new spheres to the initial ones.
C     --------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL GHOST
      LOGICAL LPR

      INTEGER*2 IU,IUSE

      INTEGER*4 I,IASS,ISA,ISO,ITO
      INTEGER*4 J
      INTEGER*4 MC,MV
      INTEGER*4 NATOM,NI,NP

      REAL*4 AP
      REAL*4 DIS,DMI
      REAL*4 RE
      REAL*4 XE,XP
      REAL*4 YE,YP
      REAL*4 ZE,ZP

      PARAMETER (MC=100000,MV=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)
      
      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine ASSIGN1 '
      DO I=1,NP
 
           IF(ISO(I).GT.NATOM)THEN
                  DMI=99999.9
                  IASS=0

                  DO J=1,NATOM

                   IF(IUSE(J).NE.1)THEN
                   DIS= (XE(J)-XP(I))*(XE(J)-XP(I))+
     &                  (YE(J)-YP(I))*(YE(J)-YP(I))+
     &                  (ZE(J)-ZP(I))*(ZE(J)-ZP(I))
                   DIS=SQRT(DIS)-RE(J)

                     IF(DIS.LT.DMI)THEN
                     IASS=J
                     DMI=DIS
                     END IF

                   END IF

                  END DO

           ISA(I)=IASS
           END IF
      END DO
     
      IF(GHOST)THEN
        J=0
        DO I=1,NP
        NI=I-J
        ITO(NI)=ITO(I)
        ISO(NI)=ISO(I)
        ISA(NI)=ISA(I)
        XP(NI)=XP(I)
        YP(NI)=YP(I)
        ZP(NI)=ZP(I)
        AP(NI)=AP(I)
        IU=IUSE(ISA(I))        

         IF(IU.EQ.2)THEN
            IUSE(ISA(I))=6        
         ELSE IF(IU.EQ.3)THEN 
            J=J+1
         ELSE IF(IU.EQ.4)THEN 
            IUSE(ISA(I))=6        
         ELSE IF(IU.EQ.5)THEN 
            J=J+1
         END IF

        END DO
        NP=NP-J
      END IF

      RETURN
      END

C
      SUBROUTINE GSHELL(NCOR,LPR,RD)
C     --------------------------------------------------------------------
C     This determines which ghost spheres are around the real ones.
C     --------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL LPR

      INTEGER*2 IUSE

      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 MC
      INTEGER*4 NCOR

      REAL*4 DD
      REAL*4 RD,RE
      REAL*4 TEST
      REAL*4 XE
      REAL*4 YE
      REAL*4 ZE

      PARAMETER (MC=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine SHELL '


      DO I=1,NCOR-1
      DO J=I+1,NCOR

        IF((IUSE(I).EQ.3).AND.(IUSE(J).EQ.6))THEN

              DD = (XE(I)-XE(J)) * (XE(I)-XE(J)) +
     &             (YE(I)-YE(J)) * (YE(I)-YE(J)) +
     &             (ZE(I)-ZE(J)) * (ZE(I)-ZE(J))
              TEST=(RE(I)+RE(J)+2*RD)*(RE(I)+RE(J)+2*RD)
              IF(DD.LT.TEST)THEN
              IUSE(I)=5
              END IF

        ELSE IF((IUSE(I).EQ.6).AND.(IUSE(J).EQ.3))THEN

              DD = (XE(I)-XE(J)) * (XE(I)-XE(J)) +
     &             (YE(I)-YE(J)) * (YE(I)-YE(J)) +
     &             (ZE(I)-ZE(J)) * (ZE(I)-ZE(J))
              TEST=(RE(I)+RE(J)+2*RD)*(RE(I)+RE(J)+2*RD)

              IF(DD.LT.TEST)THEN
              IUSE(J)=5
              END IF

        END IF

      END DO
      END DO

      RETURN 
      END
C
      SUBROUTINE VOLARE(NP,LPR,GHOST,STOT,VOL)
C     ---------------------------------------------------------------------
C     This calculates total area and volume.
C     ---------------------------------------------------------------------

      IMPLICIT NONE    

      LOGICAL LPR,GHOST

      INTEGER*2 IUSE

      INTEGER*4 I,ISA,ISO,ITO
      INTEGER*4 J,JU
      INTEGER*4 NP
      INTEGER*4 MV,MC

      REAL*4 AP
      REAL*4 DP
      REAL*4 RE,REJ
      REAL*4 VD,VN
      REAL*4 XE,XP,XEJ
      REAL*4 YE,YP,YEJ
      REAL*4 ZE,ZP,ZEJ

      REAL*8 STOT,VOL

      PARAMETER (MC=100000,MV=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine VOLARE '

      STOT=0.0E0
      VOL=0.0E0

      DO I=1,NP
      STOT=STOT+AP(I)
      END DO


      IF(.NOT.GHOST)THEN
      
        JU=0
        DO I=1,NP
        J=ISO(I)

           IF(JU.NE.J)THEN
            JU=J
            XEJ=XE(J)
            YEJ=YE(J)
            ZEJ=ZE(J)
            REJ=RE(J)
           END IF

        VOL=VOL+((XP(I)-XEJ)*XP(I)+
     &           (YP(I)-YEJ)*YP(I)+
     &           (ZP(I)-ZEJ)*ZP(I)) *AP(I)/(3.0E0*REJ)


        END DO

      END IF


c      WRITE(6,'(//A)')' ----------------- RESULTS -----------------'
c      WRITE(6,'(A/)') ' -------------------------------------------'
c      WRITE(6,'(A)')  ' *******************************************'
c      WRITE(6,'(A)')  ' *******************************************'
c      WRITE(6,'(A,F17.3,A)')' ** Area             =',STOT,'   **'
c      IF(.NOT.GHOST)THEN
c      WRITE(6,'(A,F17.3,A)')' ** Volume           =',VOL,'   **'
c      END IF
c      WRITE(6,'(A,I11,A)')' ** Number of Points =  ',NP,'       **'
c      WRITE(6,'(A)')' *******************************************'
c      WRITE(6,'(A)')' *******************************************'

      RETURN 
      END
C


      SUBROUTINE GSTAT(NCOR)
C     ------------------------------------------------------------- 
C     This prepares some statistics about the set of spheres
C     ------------------------------------------------------------- 
      IMPLICIT NONE

      INTEGER*2 IUSE

      INTEGER*4 CT6,CT4
      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 MAXI,MC
      INTEGER*4 NCOR,NINT
      INTEGER*4 TMIN,TMAX

      REAL*4 C
      REAL*4 R1,R2,RE
      REAL*4 VINT
      REAL*4 XE
      REAL*4 YE
      REAL*4 ZE

      PARAMETER (MC=100000,MAXI=100)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)

      DIMENSION CT4(MAXI),CT6(MAXI)

      WRITE(6,'(/A)')' ==> Start subroutine STAT'
      
      VINT=0.1
      TMIN=9000
      TMAX=0

      DO J=1,MAXI
      CT4(J)=0      
      CT6(J)=0      
      END DO

      DO I=1,NCOR
      C=RE(I)/VINT
      J=INT(C)+1
      IF(IUSE(I).EQ.6)THEN

        IF(J.LT.MAXI)THEN
          CT6(J)=CT6(J)+1
          TMAX=MAX(TMAX,J)
          TMIN=MIN(TMIN,J)
        END IF

      ELSE IF(IUSE(I).EQ.4)THEN 

        IF(J.LT.MAXI)THEN
          CT4(J)=CT4(J)+1
          TMAX=MAX(TMAX,J)
          TMIN=MIN(TMIN,J)
        END IF

      END IF
      END DO

      WRITE(6,'(A)')' RADII .GE. and .LT.    TYPE 4    TYPE 6'

      DO J=TMIN,TMAX
      R1=(J-1)*VINT
      R2=J*VINT
      WRITE(6,'(2F10.5,2I10)')R1,R2,CT4(J),CT6(J)
      END DO
 
      RETURN
      END
C
      SUBROUTINE SUM(NCOR,RD,OP,LPR)
C     --------------------------------------------------------------------
C     Add the solvent radius to every sphere raddi
C     Or subtract the solvent radius.
C     --------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL LPR

      INTEGER*2 IUSE

      INTEGER*4 I
      INTEGER*4 MC
      INTEGER*4 NCOR

      REAL*4 RD,RE
      REAL*4 XE
      REAL*4 YE
      REAL*4 ZE
      
      CHARACTER*4 OP
 
      PARAMETER (MC=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)     

      IF(LPR)WRITE(6,'(/A)')' ==> Start subroutine SUM'

      IF(OP.EQ.'SUMA')THEN

        IF(LPR)WRITE(6,'(A)')' Adding the radius of the solvent'
        DO I=1,NCOR
         IF(IUSE(I).NE.1)RE(I)=RE(I)+RD
        END DO

      ELSE IF(OP.EQ.'REST')THEN

        IF(LPR)WRITE(6,'(A)')' Subtracting the radius of the solvent'
        DO I=1,NCOR
         IF(IUSE(I).NE.1)RE(I)=RE(I)-RD
        END DO

      ELSE

        WRITE(6,'(A)')' Variable OP badly defined'
        STOP

      END IF
 
      RETURN
      END

C
      SUBROUTINE CREA(NATOM,NCOR,RMIN,OFAC,RD,LPR,NDIV)
C     ----------------------------------------------------------------
C     This subroutine creates the new spheres
C     ----------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL LPR
      INTEGER*2 IUSE
  
      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 K,KK,KG,KP,KENG
      INTEGER*4 LK
      INTEGER*4 MC
      INTEGER*4 NATOM,NCOR,NEK,NGE,NL1,NL2,NDIV,NENG
      
      REAL*4 DI
      REAL*4 FC,FC1
      REAL*4 OFAC,OFACT
      REAL*4 RD,RE,REK,REN,REG,REG2,REGD2,REND2A,REND2C
      REAL*4 REP,REP2,REPD2,RMIN,RMID2,RGN
      REAL*4 RIJ,RIJ2,RIK,RIN,RNK2
      REAL*4 TEST
      REAL*4 XE,XEN
      REAL*4 YE,YEN
      REAL*4 ZE,ZEN

      PARAMETER (MC=100000)
      
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
 
      DIMENSION DI(MC),KENG(MC)
      
      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine CREA '

      OFACT=1.0E0-OFAC*2.0E0

      RMID2=(RMIN+RD)*(RMIN+RD)

      NGE=0

      NL1=2   
      NL2=NCOR

  600 CONTINUE

C Loop to select the first sphere of the pair.
      DO 602 I=NL1,NL2
      IF(IUSE(I).LE.4) GO TO 602

      DO  J=1,NCOR
          DI(J)= SQRT( (XE(I)-XE(J)) * (XE(I)-XE(J)) +
     &                 (YE(I)-YE(J)) * (YE(I)-YE(J)) +
     &                 (ZE(I)-ZE(J)) * (ZE(I)-ZE(J)) )
      END DO

C
C Loop to select the second sphere of the pair.

      DO 603 J=1,I-1
      IF(IUSE(J).LE.4)GO TO 603

      RIJ= DI(J)

C If the solvent can pass through the pair this pair is discarded
      TEST=RD+RD+RE(I)+RE(J)
      IF(RIJ.GE.TEST) GO TO 603


C Determine which is the largest and smallest sphere
      IF(RE(I).GT.RE(J))THEN
        REG=RE(I)
        REP=RE(J)
        KG=I
        KP=J
      ELSE
        REG=RE(J)
        REP=RE(I)
        KG=J
        KP=I
      END IF

      REG2=REG*REG
      REP2=REP*REP
      REGD2=(REG+RD)*(REG+RD)

C Determine whether the spheres are overlapped
      IF(RIJ.LE.(REP+REG))THEN

C Test of overlapping
         TEST=REG+REP*OFACT
         IF(RIJ.LT.TEST)GO TO 603

C Test of the small sphere
         RIJ2=RIJ*RIJ
         REPD2=(REP+RD)*(REP+RD)
         RGN=(RIJ-REP+REG)*0.5E0
         REND2A=REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ)
         IF(REND2A.LE.RMID2) GO TO 603

C Spheres A
         FC=(RIJ-REP+REG)/(RIJ+REP-REG)
         FC1=FC+1.0E0
         XEN=(XE(KG)+FC*XE(KP))/FC1
         YEN=(YE(KG)+FC*YE(KP))/FC1
         ZEN=(ZE(KG)+FC*ZE(KP))/FC1
         REN=SQRT(REND2A)-RD
         RIN=(RIJ-RE(J)+RE(I))*0.5E0

      ELSE

C Test of the small sphere
         RIJ2=RIJ*RIJ
         REPD2=(REP+RD)*(REP+RD)
         REND2C=REGD2+REG2-(REG/RIJ)*(REGD2+RIJ2-REPD2)
         IF(REND2C.LE.RMID2) GO TO 603

C Calculate radius for sphere of kind B and separate B and C
         RGN=(RIJ-REP+REG)*0.5E0
         REND2A=REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ)

         IF(REND2A.GT.RMID2)THEN

C Spheres B
            FC=(RIJ-REP+REG)/(RIJ+REP-REG)
            FC1=FC+1.0E0
            XEN=(XE(KG)+FC*XE(KP))/FC1
            YEN=(YE(KG)+FC*YE(KP))/FC1
            ZEN=(ZE(KG)+FC*ZE(KP))/FC1
            REN=SQRT(REND2A)-RD
            RIN=(RIJ-RE(J)+RE(I))*0.5E0

         ELSE

C Spheres C
            FC=REG/(RIJ-REG)
            FC1=FC+1.0E0
            XEN=(XE(KG)+FC*XE(KP))/FC1
            YEN=(YE(KG)+FC*YE(KP))/FC1
            ZEN=(ZE(KG)+FC*ZE(KP))/FC1
            REN=SQRT(REND2C)-RD
            IF(KG.EQ.I)THEN
               RIN=REG
            ELSE
               RIN=(RIJ-REG)
            END IF

         END IF

      END IF

C
C Test of overlapping for the new sphere
      NENG=0
      DO 604 K=1,NCOR

      RIK=DI(K)

      IF(RIK.GE.(RIN+REN+RE(K))) GO TO 604

      IF(IUSE(K).LE.3) GO TO 604

      RNK2 = (XEN-XE(K)) * (XEN-XE(K)) +
     &       (YEN-YE(K)) * (YEN-YE(K)) +
     &       (ZEN-ZE(K)) * (ZEN-ZE(K))
      
      REK=RE(K)

      IF(RNK2.GE.((REK+REN)*(REK+REN))) GO TO 604
      
      TEST=(REN-REK)*(REN-REK)
      IF(RNK2.LE.TEST)THEN

           IF(REN.GT.REK)THEN
              IF(TEST.LT.4.0E-4) GO TO 603
              NENG=NENG+1
              KENG(NENG)=K
              GO TO 604
           ELSE
              GO TO 603
           END IF

      END IF

      TEST=REK+REN*OFACT
      IF(TEST.LT.0.0E0) GO TO 604
      TEST=TEST*TEST
      IF(RNK2.LE.TEST)GO TO 603

  604 CONTINUE

C
C Mark spheres engulfed by the new sphere
C
         DO LK=1,NENG
         K=KENG(LK)
         IUSE(K)=2
         END DO
       
C
C Creates the new spheres
C
         NCOR=NCOR+1
         XE(NCOR)=XEN
         YE(NCOR)=YEN
         ZE(NCOR)=ZEN
         RE(NCOR)=REN
         IUSE(NCOR)=6
         DI(NCOR)=RIN

  603 CONTINUE

  602 CONTINUE


      IF(NCOR.GT.MC)THEN
        WRITE(6,'(A)')' %-ERROR-% DIMENSION. TOO MANY SPHERES CREATED'
        STOP
      END IF

      NGE=NGE+1

      IF(NCOR.NE.NL2) THEN 
C Mark spheres with area zero
      CALL MZERO5(NATOM,NCOR,NL2,NDIV,LPR)
        IF(NCOR.NE.NL2) THEN 
           NL1=NL2+1
           NL2=NCOR
           GO TO 600
        END IF
      END IF
      
      IF(LPR)THEN
       WRITE(6,'(A,I10)')' Number of generations     =',NGE
      END IF

      RETURN
      END
C
      SUBROUTINE BULK(NATOM,NCOR,NDIV,OFAC,RD,LPR)
C     ----------------------------------------------------------------
C     This subroutine creates new spheres 
C     ----------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL LPR
      INTEGER*2 IUSE
  
      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 K,KK,KENG
      INTEGER*4 L,L1,LK,LL
      INTEGER*4 MC
      INTEGER*4 NATOM,NCOR,NEK,NGE,NDIV,NL1,NL2,NTRI,NENG
      INTEGER*4 SC
      INTEGER*4 ULL
      
      REAL*4 DD,DI
      REAL*4 FC,FC1
      REAL*4 OFAC,OFACT
      REAL*4 RD,RE,REI,REJ,REK,REN,RNK2
      REAL*4 RIJ,RIK,RIN
      REAL*4 TEST
      REAL*4 XC1,XE,XEN,XPN
      REAL*4 YC1,YE,YEN,YPN
      REAL*4 ZC1,ZE,ZEN,ZPN

      REAL*8 CV

      PARAMETER (MC=100000)
      
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
 
      DIMENSION DI(MC),SC(MC),KENG(MC)
      

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine BULK '
     
      OFACT=1.0E0-OFAC*2.0E0

      NTRI=60*4**(NDIV-1)
      NGE=0

      NL1=2   
      NL2=NCOR

  600 CONTINUE
C
C Loop to select the first sphere of the pair.
      DO 602 I=NL1,NL2
      IF(IUSE(I).LE.4) GO TO 602
      REI=RE(I)

      DO  J=1,NCOR
          DI(J)=SQRT( (XE(I)-XE(J)) * (XE(I)-XE(J)) +
     &                (YE(I)-YE(J)) * (YE(I)-YE(J)) +
     &                (ZE(I)-ZE(J)) * (ZE(I)-ZE(J)) )
      END DO
C
C Loop to select the second sphere of the pair.

      DO 603 J=1,I-1
      IF(IUSE(J).LE.4)GO TO 603
      REJ=RE(J)
      RIJ= DI(J)
C
C If the solvent can pass through the pair this pair is discarded

      IF(RIJ.GE.(RD+RD+REI+REJ)) GO TO 603
C
C Test of overlapping and assignment of the radius of the new sphere

      IF(REI.GT.REJ)THEN
          TEST=REI+REJ*OFACT
          IF(RIJ.LT.TEST)GO TO 603
          REN=REI
      ELSE
          TEST=REJ+REI*OFACT
          IF(RIJ.LT.TEST)GO TO 603
          REN=REJ
      END IF
C
C Computes coordinates of the new sphere 

      FC=(RIJ-REJ+REI)/(RIJ+REJ-REI)
      FC1=FC+1.0E0
      XEN=(XE(I)+FC*XE(J))/FC1
      YEN=(YE(I)+FC*YE(J))/FC1
      ZEN=(ZE(I)+FC*ZE(J))/FC1
      RIN=(RIJ-REJ+REI)*0.5E0

C
C Test of overlapping for the new sphere

      NENG=0
      DO 604 K=1,NCOR

      RIK=DI(K)

      IF(RIK.GE.(RIN+REN+RE(K))) GO TO 604

      IF(IUSE(K).LE.3) GO TO 604

      RNK2 = (XEN-XE(K)) * (XEN-XE(K)) +
     &       (YEN-YE(K)) * (YEN-YE(K)) +
     &       (ZEN-ZE(K)) * (ZEN-ZE(K))
      
      REK=RE(K)

      IF(RNK2.GE.((REK+REN)*(REK+REN))) GO TO 604
      
      TEST=(REN-REK)*(REN-REK)
      IF(RNK2.LE.TEST)THEN

           IF(REN.GT.REK)THEN
              IF(TEST.LT.4.0E-4) GO TO 603
              NENG=NENG+1
              KENG(NENG)=K
              GO TO 604
           ELSE
              GO TO 603
           END IF

      END IF

      TEST=REK+REN*OFACT
      IF(TEST.LT.0.0E0) GO TO 604
      TEST=TEST*TEST
      IF(RNK2.LE.TEST)GO TO 603

  604 CONTINUE
C
C Find spheres that are overlapped with the new sphere
C Use radius of the sphere plus RD
C
      NEK=0
      DO K=1,NATOM

         IF(IUSE(K).GE.2)THEN
            DD= (XEN-XE(K))*(XEN-XE(K)) +
     &          (YEN-YE(K))*(YEN-YE(K)) + 
     &          (ZEN-ZE(K))*(ZEN-ZE(K))

            TEST=(REN+RD+RD+RE(K))*(REN+RD+RD+RE(K))
               
                IF(DD.LT.TEST)THEN
                  NEK=NEK+1
                  SC(NEK)=K
                END IF

         END IF

      END DO


C
C Determine if the new sphere has accessible surface area zero.

      ULL=1
      DO 3 L=1,NTRI
      XPN=XC1(L)*(REN+RD)+XEN
      YPN=YC1(L)*(REN+RD)+YEN
      ZPN=ZC1(L)*(REN+RD)+ZEN

            L1=ULL
            DO LL=L1,NEK
            ULL=LL
            K=SC(LL)
            DD=(XPN-XE(K))*(XPN-XE(K)) +
     &         (YPN-YE(K))*(YPN-YE(K)) + 
     &         (ZPN-ZE(K))*(ZPN-ZE(K)) 
            TEST=(RE(K)+RD)*(RE(K)+RD)
            IF(DD.LT.TEST)GO TO 3
            END DO

            DO LL=1,L1-1
            ULL=LL
            K=SC(LL)
            DD=(XPN-XE(K))*(XPN-XE(K)) +
     &         (YPN-YE(K))*(YPN-YE(K)) + 
     &         (ZPN-ZE(K))*(ZPN-ZE(K)) 
            TEST=(RE(K)+RD)*(RE(K)+RD)
            IF(DD.LT.TEST)GO TO 3
            END DO
            
      GO TO 603
      
    3 CONTINUE

C
C Marks spheres that are engulfed by the new sphere
C
      DO LK=1,NENG
      K=KENG(LK)
      IUSE(K)=2
      END DO

C
C Save information about the new sphere
C
      NCOR=NCOR+1
      XE(NCOR)=XEN
      YE(NCOR)=YEN
      ZE(NCOR)=ZEN
      RE(NCOR)=REN
      IUSE(NCOR)=6
      DI(NCOR)=RIN

  603 CONTINUE

  602 CONTINUE

      IF(NCOR.GT.MC)THEN
        WRITE(6,'(A)')' %-ERROR-% DIMENSION. TOO MANY SPHERES CREATED'
        STOP
      END IF


C Compute number of generations
      NGE=NGE+1

C Check if there are new spheres
      IF(NCOR.NE.NL2) THEN 
C Mark spheres with area zero
        CALL MZERO5(NATOM,NCOR,NL2,NDIV,LPR)
        IF(NCOR.NE.NL2) THEN 
           NL1=NL2+1
           NL2=NCOR
           GO TO 600
        END IF
      END IF

C Write
C Number of generations
      IF(LPR)THEN
      WRITE(6,'(A,I10)')' Number of generations       =',NGE
      END IF


      RETURN
      END
C

C
      SUBROUTINE MZERO5(NATOM,NCOR,NL2,NDIV,LPR)
C     ----------------------------------------------------------------
C     This discards spheres(IUSE=2) engulfed by another
C     Marks (IUSE=4) the spheres with total area zero.
C     ----------------------------------------------------------------
      IMPLICIT NONE
    
      LOGICAL LPR,CHECK

      INTEGER*2 IUSE
  
      INTEGER*4 I,J,K
      INTEGER*4 MC
      INTEGER*4 L,LL,L1
      INTEGER*4 NATOM,NCOR,NDIV,NEJ,NEW,NI,NL2,NTRI,NZERO
      INTEGER*4 SC
      INTEGER*4 ULL

      REAL*4 DD
      REAL*4 RE,REI
      REAL*4 TEST
      REAL*4 XE,XEI,XC1,XP
      REAL*4 YE,YEI,YC1,YP
      REAL*4 ZE,ZEI,ZC1,ZP

      REAL*8 CV

      PARAMETER (MC=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)


      DIMENSION SC(MC),CHECK(MC)

      IF(LPR)WRITE(6,'(/A)')' ==> Start subroutine MZERO5'

C
C Discard Spheres that have been engulfed by another (iuse=2)
      J=0
      K=0
      DO I=NATOM+1,NCOR
      NI=I-J
      XE(NI)=XE(I)
      YE(NI)=YE(I)
      ZE(NI)=ZE(I)
      RE(NI)=RE(I)
      IUSE(NI)=IUSE(I)

        IF(IUSE(I).EQ.2) THEN
           J=J+1
           IF(I.LE.NL2)THEN
              K=K+1
           END IF
        END IF

      END DO

      NEW=NCOR-NL2
      NCOR=NCOR-J
      NL2=NL2-K

C Print information about spheres discarded
      IF(LPR)THEN
      WRITE(6,'(A,I7)')' New spheres in last generation       =',NEW
      NEW=NEW-(J-K)
      WRITE(6,'(A,I7)')' Final new spheres in last generation =',NEW
      WRITE(6,'(A,I7)')' Spheres discarded (engulfed)         =',J
      END IF

C Find spheres with zero area.
      NZERO=0
      NTRI=60*4**(NDIV-1)

      DO I=1,NL2
      CHECK(I)=.FALSE.
      END DO

      DO I=NL2+1,NCOR
      CHECK(I)=.FALSE.
      IF(IUSE(I).GT.4)CHECK(I)=.TRUE.
      END DO

C
C Start
C
      DO 5 I=NCOR,1,-1
      IF(CHECK(I))THEN

      XEI=XE(I)
      YEI=YE(I)
      ZEI=ZE(I)
      REI=RE(I)

C Find spheres that overlap sphere I.
C  
          NEJ=0
          DO 7 J=1,NCOR
          IF(IUSE(J).GE.3)THEN

              DD=(XEI-XE(J))*(XEI-XE(J)) +
     &           (YEI-YE(J))*(YEI-YE(J)) +
     &           (ZEI-ZE(J))*(ZEI-ZE(J)) 
              TEST=(REI+RE(J))*(REI+RE(J))

              IF(DD.LT.TEST)THEN
              IF(I.NE.J)THEN

                  NEJ=NEJ+1
                  SC(NEJ)=J

              END IF
              END IF

          END IF
     
    7     CONTINUE
C          
C Mark spheres that should be checked
C
          IF(I.GT.NL2)THEN
          DO LL=1,NEJ
          J=SC(LL)
          IF(IUSE(J).GT.4)CHECK(J)=.TRUE.
          END DO
          END IF

C
C Determine if the sphere I has area zero
          ULL=1
          DO 4 L=1,NTRI
          XP=XC1(L)*REI+XEI
          YP=YC1(L)*REI+YEI
          ZP=ZC1(L)*REI+ZEI

             L1=ULL
             DO LL=L1,NEJ 
             ULL=LL
             J=SC(LL)
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)
             IF(DD.LT.TEST) GO TO 4
             END DO

             DO LL=1,L1-1 
             ULL=LL
             J=SC(LL)
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)
             IF(DD.LT.TEST) GO TO 4
             END DO
          
             GO TO 5

    4     CONTINUE

      IUSE(I)=4
      NZERO=NZERO+1

      END IF
    5 CONTINUE



      IF(LPR)THEN
      WRITE(6,'(A,I7)')' Spheres marked with zero area        =',NZERO
      END IF

      RETURN
      END
C

C
      SUBROUTINE CLEAN5(NATOM,NCOR,NDIV,LPR)
C     ----------------------------------------------------------------
C     This discard the spheres of area zero that are not needed for
C     the computation of the surface.
C     ----------------------------------------------------------------
      IMPLICIT NONE
    
      LOGICAL LPR,USEFUL

      INTEGER*2 IUSE
  
      INTEGER*4 I,J
      INTEGER*4 MC
      INTEGER*4 L,LJ,L1
      INTEGER*4 NATOM,NCOR,NDIV,NEJ,NI,NTRI
      INTEGER*4 SC
      INTEGER*4 ULJ

      REAL*4 DD
      REAL*4 RE,REI
      REAL*4 TEST
      REAL*4 XE,XEI,XC1,XP
      REAL*4 YE,YEI,YC1,YP
      REAL*4 ZE,ZEI,ZC1,ZP

      REAL*8 CV

      PARAMETER (MC=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)

      DIMENSION SC(MC),USEFUL(MC)

      IF(LPR)WRITE(6,'(/A)')' ==> Start subroutine CLEAN5  '

      NTRI=60*4**(NDIV-1)


      DO I=1,NATOM
         USEFUL(I)=.TRUE.   
      END DO

      DO I=NATOM+1,NCOR
          IF(IUSE(I).NE.4)THEN
            USEFUL(I)=.TRUE.   
          ELSE
            USEFUL(I)=.FALSE.   
          END IF
      END DO

      DO I=1,NCOR
      IF(IUSE(I).GE.4)THEN      
      XEI=XE(I)
      YEI=YE(I)
      ZEI=ZE(I)
      REI=RE(I)

C Find spheres that overlap sphere I.
          NEJ=0
          DO J=1,NCOR
          IF(IUSE(J).GE.3)THEN

              DD=(XEI-XE(J))*(XEI-XE(J)) +
     &           (YEI-YE(J))*(YEI-YE(J)) +
     &           (ZEI-ZE(J))*(ZEI-ZE(J)) 
              TEST=(REI+RE(J))*(REI+RE(J))

              IF(DD.LT.TEST)THEN
              IF(I.NE.J)THEN
                  NEJ=NEJ+1
                  SC(NEJ)=J
              END IF
              END IF

          END IF
          END DO

C Find spheres that are needed to discard the triangle L
          ULJ=1
          DO 40 L=1,NTRI
          XP=XC1(L)*REI+XEI
          YP=YC1(L)*REI+YEI
          ZP=ZC1(L)*REI+ZEI

C Among the useful spheres
             L1=ULJ
             DO LJ=L1,NEJ
             ULJ=LJ
             J=SC(LJ)
             IF(USEFUL(J))THEN
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)

                IF(DD.LT.TEST) GO TO 40

             END IF
             END DO

             DO LJ=1,L1-1
             ULJ=LJ
             J=SC(LJ)
             IF(USEFUL(J))THEN
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)

                IF(DD.LT.TEST) GO TO 40

             END IF
             END DO

C Among the not useful
             DO LJ=1,NEJ
             ULJ=LJ
             J=SC(LJ)
             IF(.NOT.USEFUL(J))THEN
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)

                IF(DD.LT.TEST)THEN
                USEFUL(J)=.TRUE.
                GO TO 40
                END IF

             END IF
             END DO



   40     CONTINUE
      END IF
      END DO

C Discard spheres totally inside of others
      J=0
      DO I=NATOM+1,NCOR
      NI=I-J
      XE(NI)=XE(I)
      YE(NI)=YE(I)
      ZE(NI)=ZE(I)
      RE(NI)=RE(I)
      IUSE(NI)=IUSE(I)

        IF(.NOT.USEFUL(I)) THEN
           J=J+1
        END IF

      END DO

      NCOR=NCOR-J
    
      IF(LPR)THEN
      WRITE(6,'(A,I7)')' Spheres discarded (not useful)       =',J
      END IF

      RETURN
      END
