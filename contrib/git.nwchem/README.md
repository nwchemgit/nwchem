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
   
