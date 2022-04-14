#!/bin/bash
python3 -m pip install --user mkdocs
python3 -m pip install --user mkdocs-material
git clone https://github.com/cmitu/mkdocs-altlink-plugin
cd mkdocs-altlink-plugin
python3 -m pip install --user -e .
cd ..
python3 -m pip install --user python-markdown-math
