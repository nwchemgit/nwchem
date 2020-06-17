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
git clone https://github.com/nwchemgit.github.io
cd nwchemgit.github.io
mkdocs -v gh-deploy --config-file /path/to/mkdocs/mkdocs.yml --remote-branch master
```

Browse new web pages at
https://nwchemgit.github.io/
