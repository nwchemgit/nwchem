mkdir -p nwchem-6.5
cd nwchem-6.5
svn checkout --non-recursive https://svn.pnl.gov/svn/nwchem/branches/release-6-5-patches/src
cd src
svn update \
 property \
 config \
 rimp2 \
 moints \
 include \
 blas \
 lapack \
 tools \
 peigs \
 perfm \
 optim \
 basis \
 geom \
 inp \
 input \
 pstat \
 rtdb \
 task \
 symmetry \
 bq \
 cons \
 NWints \
 atomscf \
 gradients \
 nwdft \
 nwxc \
 stepper \
 driver \
 hessian \
 vib \
 mp2_grad \
 util \
 ddscf \
 tce 
cd tools
./get-tools
cd ..
cp  build4xeonphi.sh nwchem-6.5
tar --exclude=.svn -cf nwchem-6.5.xeonphi.tar nwchem-6.5/*
gzip nwchem-6.5.xeonphi.tar

