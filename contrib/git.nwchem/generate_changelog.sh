#!/usr/bin/env bash 
#git log v6.8...v6.8-release --oneline --reverse
#args are tag1 and tag2
#or branches, eg:
#$1 origin/hotfix/release-7-0-0
#$2  origin/release-7-2-0
git log $1..$2 --oneline --reverse --no-merges --pretty="- %s (%as:%an)" | grep -v "ci skip" | grep -v "travis skip" | uniq|more
