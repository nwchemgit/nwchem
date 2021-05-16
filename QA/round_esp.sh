#!/bin/bash
input_file=$1
read -r line < "$input_file"
echo "$line"
while IFS= read -r line || [[ -n "$line" ]]
do
    read  -n 2 -r elem
    if [ ! -z "$elem" ]; then
	printf "%2s\t " $elem
	for i in {1..8}; do
	    read  -n 12 -r x ; printf "%9.4f"  $x
	done
	printf "\n"
    fi
done < "$input_file"
