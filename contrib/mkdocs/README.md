# Updating the website https://nwchemgit.github.io using mkdocs

The wiki source hosted on the repository https://github.com/nwchemgit/nwchem-wiki (and mirrored in https://github.com/nwchemgit/nwchem/wiki) is used to generated the content of the https://nwchemgit.github.io using [mkdocs](https://mkdocs.readthedocs.io/). 

This directory contains scripts to generate the mkdocs content:  
 - install_mkdocs.sh  
 - build_website.sh  
 - commit_changes.sh  
 
 This scripts are used by the Github Action https://github.com/nwchemgit/nwchem-wiki/blob/master/.github/workflows/website_update.yml so that every time a new commit is pushed to the NWChem wiki, the content of the website  https://nwchemgit.github.io is updated.

`IMPORTANT`: The Information below is for reference ONLY, since the wiki source is now automatically uploaded to the https://nwchemgit.github.io by means of the Github Action https://github.com/nwchemgit/nwchem-wiki/blob/master/.github/workflows/website_update.yml

## ~~Install mkdocs https://www.mkdocs.org/~~ (for reference ONLY)

```
 python3 -m pip install --user mkdocs
 python3 -m pip install --user mkdocs-material
``` 

## ~~Install the MkDocs plugin that simplifies internal links creation~~ (for reference ONLY)

```
 git clone https://github.com/cmitu/mkdocs-altlink-plugin
 cd mkdocs-altlink-plugin
 python3 -m pip install --user -e .
```

## ~~Install the  python-markdown-math plugin for MathJax~~ (for reference ONLY)

```
python3 -m pip  install --user python-markdown-math
```
# ~~Edit the configuration~~ (for reference ONLY)
```
mkdocs.yml
```
https://www.mkdocs.org/#theming-our-documentation  
https://www.mkdocs.org/user-guide/plugins/  

## ~~Checkout the Markdown source of the Wiki pages~~ (for reference ONLY)
```
git clone https://github.com/nwchemgit/nwchem-wiki docs
```

## ~~Checkout the html files of the archive forum~~ (for reference ONLY)
```
git clone https://github.com/nwchemgit/archivedforum
rsync -av archivedforum/Special_AWCforum docs/.
```

## ~~Test the changes~~ (for reference ONLY)

The file `index.md` is the main file

```
mkdocs serve
```
Point your browser to  http://127.0.0.1:8000



## ~~push files to the nwchemgit.github.io: method #1~~ (for reference ONLY)

### Build the html source
```
mkdocs build
```
### push the html source to the nwchemgit.github.io repository
```
git clone https://github.com/nwchemgit.github.io
cd nwchemgit.github.io
rsync -av /path/to/mkdocs/site/* .
 git add -A 
git commit -m
git push
```

## ~~push files to the nwchemgit.github.io: method #2~~ (for reference ONLY)

```
git clone https://github.com/nwchemgit/nwchemgit.github.io
cd nwchemgit.github.io
mkdocs -v gh-deploy --config-file /path/to/mkdocs/mkdocs.yml --remote-branch master
```

Browse new web pages at
https://nwchemgit.github.io/

## tools for rendering latex math equations (for reference ONLY)

* python-markdown-math
https://github.com/mitya57/python-markdown-math
https://pypi.org/project/python-markdown-math/
```
python3 -m pip install python-markdown-math
```
https://github.com/mkdocs/mkdocs/issues/253
* arithmatex from pymdown-extensions
```
python3 -m pip install pymdown-extensions
```
https://facelessuser.github.io/pymdown-extensions  

Caveat: gets confused by the existing latex2svg embedded images

* katex
https://pypi.org/project/markdown-katex/
https://gitlab.com/mbarkhau/markdown-katex
```
python3 -m pip install markdown-katex
```
Not tried yet

* more in https://github.com/Python-Markdown/markdown/wiki/Third-Party-Extensions

