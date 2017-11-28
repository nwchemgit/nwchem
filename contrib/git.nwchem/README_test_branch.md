How to work on a test branch
=====





How to re-sync a forked repository
====

Adapted from https://help.github.com/articles/syncing-a-fork/
  git remote -v
  git remote add upstream https://github.com/nwchemgit/nwchem.git
  git remote -v
  git fetch upstream
  git checkout  master
  git merge  upstream/master
  git push
(or for new branch: git push origin mynewbranch)
  git pull

