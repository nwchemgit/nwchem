#!/bin/bash
python3 -m pip install --upgrade pip --quiet
python3 -m pip install --user  MarkupSafe --use-pep517 --quiet
python3 -m pip install --user mkdocs --quiet
python3 -m pip install --user mkdocs-material --quiet
python3 -m pip install --user git+https://github.com/cmitu/mkdocs-altlink-plugin/ --use-pep517 --quiet
python3 -m pip install --user python-markdown-math --quiet
python3 -m pip install --user install pymdown-extensions --quiet
