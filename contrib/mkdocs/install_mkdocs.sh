#!/bin/bash
if [[ -z "$VIRTUAL_ENV" ]]; then
    USEROPT=--user
else
    USEROPT=
fi
python3 -m pip install --upgrade pip --quiet  $USEROPT
python3 -m pip install  MarkupSafe --use-pep517 --quiet $USEROPT
python3 -m pip install mkdocs --quiet $USEROPT
python3 -m pip install mkdocs-material --quiet $USEROPT
python3 -m pip install git+https://github.com/cmitu/mkdocs-altlink-plugin/ --use-pep517 --quiet $USEROPT
python3 -m pip install python-markdown-math --quiet $USEROPT
python3 -m pip install install pymdown-extensions --quiet $USEROPT
#python3 -m pip install install mkdocs-with-pdf --quiet $USEROPT
python3 -m pip install install mkdocs-print-site-plugin --quiet $USEROPT
