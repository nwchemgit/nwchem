
Some git config options:

git config --global user.email "you@example.com"  
git config --global user.name "Your Name"  
git config --global core.editor "emacs"  (default on most OS is vi)  
git config --global merge.tool "meld"  

Using --local instead of --global above allows you to set these preferences on a per repository basis. --global sets "system wide preferences"  

Branching Notes:

Assuming we are currently on *master*

- git checkout -b test  (Creates new branch test based on master branch)
 
We are now on branch test unless specified otherwise
- git add f1 f2 ... (add some files)
- git commit -m "First commit to test branch"

Occasionally, we might want to get updates from master into branch *test* using:
 - git pull origin master

- Few more commits to *test* branch

When ready to merge into master:

    git checkout master (we are now on master branch)  
    git pull (make sure master is upto date)  
    git merge test  
    If there are any conflicts, fix those files manually, commit the fixed files and then: (git merge abort if something seems wrong to abort the merge process)  
    git push  

------------------------------

- OPTIONAL  
For some reason, if you want to push your local branch *test* to the github repo,  
  - When you push *test* branch for 1st time, use: git push --set-upstream origin test  
  - For Subsequent pushes to *test* branch "git push" will do.  

To delete a branch that has been pushed to the github repo:  
- git checkout master (switch to master branch, since you cannot delete a branch you are currently on)  
**We are now on master branch**  
- git push --delete origin test --> Will delete the *test* branch from remote(github) repo  

At this point, your local copy of branch *test* still exists - To remove it use:  
- git branch -d test (delete a local branch *test* that has been merged with master)  

Caveat: Usually you delete a branch after merging it back into master as we did above. For some reason, If *test* branch is not merged into master and you want to delete it, use:  
- git branch -D test  (same cmd as above but using -D instead of -d )  


How to re-sync a forked repository
====

Adapted from https://help.github.com/articles/syncing-a-fork/   
  
Your fork of nwchem -> git@github.com:USERNAME/nwchem.git

  git clone git@github.com:USERNAME/nwchem.git  
  cd nwchem  

  git remote -v  
  git remote add upstream https://github.com/nwchemgit/nwchem.git  
  git remote -v  
  git fetch upstream  
  git checkout  master  
  git merge  upstream/master (OR git pull upstream master)  
  git push  
(or for new branch: git push origin mynewbranch)  
  git pull 
  

How to merge a Pull
=====

1.  git clone https://github.com/nwchemgit/nwchem.git
    cd nwchem/
2.  git checkout -b mjw99-master master
    git branch
3.  git pull https://github.com/mjw99/nwchem-1.git master
    gitk
4.  git checkout master
    gitk
5.  git merge --no-ff mjw99-master
    gitk
6.  git push origin master
  
  or
  
5.  git rebase -p mjw99-master 
    gitk
6.  git push origin master
		

More commands
====

* Revision information

$ git describe

 nwchem_on_git-14-g2b641c8
 
(this needs tags to be present)

$ git describe --all

 heads/master


* What changes are scheduled to be sent by git push?

$ git diff --stat --cached remotes/origin/release-6-8

or

$ git diff --stat --cached origin/master

git diff branch1 branch2 --name-only (shows diff b/w 2 branches)

How to git clone without git

 curl -LJO https://github.com/nwchemgit/nwchem/tarball/master
 curl -LJO https://github.com/nwchemgit/nwchem/tarball/hotfix/release-6-8