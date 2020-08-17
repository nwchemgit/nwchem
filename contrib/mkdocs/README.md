## Install mkdocs https://www.mkdocs.org/

```
 python3 -m pip install --user mkdocs
 python3 -m pip install --user mkdocs-material
``` 

## Install the MkDocs plugin that simplifies internal links creation

```
 git clone https://github.com/cmitu/mkdocs-altlink-plugin
 cd mkdocs-altlink-plugin
 python3 -m pip install --user -e .
```

## Install the  python-markdown-math plugin for MathJax (work in progress)

```
python3 -m pip  install --user python-markdown-math
```
# Edit the configuration
```
mkdocs.yml
```
https://www.mkdocs.org/#theming-our-documentation  
https://www.mkdocs.org/user-guide/plugins/  

## Checkout the Markdown source of the Wiki pages
```
git clone https://github.com/nwchemgit/nwchem-wiki docs
```

## Checkout the html files of the archive forum
```
git clone https://github.com/nwchemgit/archivedforum
rsync -av archivedforum/Special_AWCforum docs/.
```

## Test the changes 

The file `index.md` is the main file

```
mkdocs serve
```
Point your browser to  http://127.0.0.1:8000



## push files to the nwchemgit.github.io: method #1

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

## push files to the nwchemgit.github.io: method #2

```
git clone https://github.com/nwchemgit/nwchemgit.github.io
cd nwchemgit.github.io
mkdocs -v gh-deploy --config-file /path/to/mkdocs/mkdocs.yml --remote-branch master
```

Browse new web pages at
https://nwchemgit.github.io/

## tools for rendering latex math equations

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
