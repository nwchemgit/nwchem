How to fix bizarre check-ins
====

* Case when committer got in bizarre state on this git repository and
the changes have been pushed already to github.com

* git log showed

$ git log  
commit d890621e5d1b664ebb94683e804e859f527a41aa  
Merge: 86c8a48 a7055a2  
Author: Niri Govind <niri.govind@pnnl.gov>  
Date:   Mon Nov 27 15:07:35 2017 -0800  

Merge branch 'master' of https://github.com/nwchemgit/nwchem

commit a7055a28bc6efa8e6ceff79ae541d5e0b8d8a8b8    
Author: edoapra <edoardo.apra@gmail.com>  
Date:   Mon Nov 27 14:14:54 2017 -0800  

svn replaced by git
    
*  let's go back to the previos commit a7055a28bc6efa8e6ceff79ae541d5e0b8d8a8b8

$ git reset --hard a7055a28bc6efa8e6ceff79ae541d5e0b8d8a8b8

* Check that everything is good with gitg

* Push everything back to github.com

$ git push -f
   

* Full log was ...  
Branch: refs/heads/master  
  Home:   https://github.com/nwchemgit/nwchem  
  Commit: 68ac6c4ac3f284a8c1fa7dad82f41b22c039927e  
      https://github.com/nwchemgit/nwchem/commit/68ac6c4ac3f284a8c1fa7dad82f41b22c039927e  
  Author: Niri Govind <niri.govind@pnnl.gov>  
  Date:   2017-11-20 (Mon, 20 Nov 2017)  

  Changed paths:    
src/nwdft/rt_tddft/rtutils/rt_tddft_movecs_import.F   

  Log Message:   
  minor edit to reflect the initial density matrix  


  Commit: 2447af6b3899c46408f645dc32d5928036177987  
      https://github.com/nwchemgit/nwchem/commit/2447af6b3899c46408f645dc32d5928036177987  
  Author: Niri Govind <niri.govind@pnnl.gov>  
  Date:   2017-11-20 (Mon, 20 Nov 2017)  

  Changed paths:   
A docs/Capabilities  
A docs/_config.yml  
A docs/_layouts/default.html  
A docs/assets/css/style.scss  
A docs/assets/images/checker.png  
A docs/assets/js/scale.fix.js  
A docs/index.html  
A docs/index2.html  
M src/nwpw/nwpw_input.F  

  Log Message:  
  Merge branch 'master' of https://github.com/nwchemgit/nwchem  


  Commit: 86c8a484492300231fb8d18666e17bfa77f1f1c9  
  https://github.com/nwchemgit/nwchem/commit/86c8a484492300231fb8d18666e17bfa77f1f1c9  
  Author: Niri Govind <niri.govind@pnnl.gov>  
  Date:   2017-11-20 (Mon, 20 Nov 2017)  

  Changed paths:  
  R docs/Capabilities  
  R docs/_config.yml  
  R docs/_layouts/default.html  
  R docs/assets/css/style.scss  
  R docs/assets/images/checker.png  
  R docs/assets/js/scale.fix.js  
  R docs/index.html  
  R docs/index2.html  

  Log Message:  
  Merge branch 'master' of https://github.com/nwchemgit/nwchem  


  Commit: d890621e5d1b664ebb94683e804e859f527a41aa  
      https://github.com/nwchemgit/nwchem/commit/d890621e5d1b664ebb94683e804e859f527a41aa  
  Author: Niri Govind <niri.govind@pnnl.gov>  
  Date:   2017-11-27 (Mon, 27 Nov 2017)  

  Changed paths:  
    A .travis.yml  
    M README.md  
    M src/util/util_nwchem_version.bash  
    A travis/build_openblas.sh  
    A travis/build_scalapack.sh  
    A travis/sleep_loop.sh  

  Log Message:  
  Merge branch 'master' of https://github.com/nwchemgit/nwchem  


Compare: https://github.com/nwchemgit/nwchem/compare/a7055a28bc6e...d890621e5d1b
