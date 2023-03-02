#!/bin/bash
#if [[ $(command -v pandoc 2 > /dev/null) ]]; then
#    echo need to install pandoc
#    sudo apt-get install -y pandoc
#    sudo apt-get install -y pandoc-citeproc
#fi
#pandoc --version
#if [[ "$?" != 0 ]]; then
#    echo broken pandoc installation
#    exit 1
#fi
if [[ -z "$VIRTUAL_ENV" ]]; then
    USEROPT=--user
else
    USEROPT=
fi
python3 -m pip install --upgrade pip --quiet  $USEROPT
python3 -m pip install  MarkupSafe --use-pep517 --quiet $USEROPT
#python3 -m pip install mkdocs==1.3.1 --quiet $USEROPT
python3 -m pip install mkdocs==1.4.0 --quiet $USEROPT
python3 -m pip install mkdocs-material --quiet $USEROPT
python3 -m pip install git+https://github.com/cmitu/mkdocs-altlink-plugin/ --use-pep517 --quiet $USEROPT
python3 -m pip install python-markdown-math --quiet $USEROPT
python3 -m pip install install pymdown-extensions --quiet $USEROPT
#python3 -m pip install install mkdocs-with-pdf --quiet $USEROPT
python3 -m pip install install mkdocs-print-site-plugin --quiet $USEROPT
python3 -m pip install  mkdocs-bibtex --use-pep517 --quiet $USEROPT
python3 -m pip install  pandoc --use-pep517 --quiet $USEROPT
python3 -m pip install  pypandoc_binary --use-pep517 --quiet $USEROPT
python3 -m pip install install git+https://github.com/flywire/caption --use-pep517 --quiet $USEROPT
