#!/bin/bash -f
#git log v6.8...v6.8-release --oneline --reverse
#args are tag1 and tag2
git log $1...$2 --oneline --reverse
