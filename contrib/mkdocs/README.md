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

## Build the html
```
mkdocs build
```

## Create gh-pages branch (if needed)

**Need to do it only once!
The second time will wipe out the whole existing branch ...**
```
cd /path/to/repo-gh-pages
git symbolic-ref HEAD refs/heads/gh-pages
rm .git/index
git clean -fdx
echo "My GitHub Page" > index.html
git add .
git commit -a -m "First pages commit"
git push origin gh-pages
```

## push files to gh-pages branch

```
cd /path/to/repo-gh-pages
rsync -av /path/to/mkdocs/site/* .
git commit -m
git push
```

# Browse new web pages
https://edoapra.github.io/nwchem/
or in future
https://nwchemgit.github.io/nwchem/
