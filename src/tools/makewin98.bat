mkdir lib
cd lib
mkdir win32
cd ..

mkdir bin
cd bin
mkdir win32
cd ..

cd ma
nmake
cd ..

cd tcgmsg-mpi
nmake
cd ..

cd armci\src
nmake
cd ..\..

cd global\src
nmake
cd ..\..

cd pario

cd elio
nmake
cd ..

cd dra
nmake
cd ..

cd sf
nmake
cd ..

cd eaf
nmake
cd ..

nmake
cd ..

