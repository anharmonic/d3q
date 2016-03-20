#!/bin/bash

line=$1
shift

col=$1
shift

echo "Getting column $col, line $line"

files="$@"
echo "Files: $files"

nwords=$(grep -v "^ *#" $1|head -n1 |wc -w)
nat3=$[($nwords-5)/3]
nat=$[$nat3/3]

echo "$nwords, $nat3, $nat"
echo "idx: 1"
echo "path: 2"
echo "q_x..q_z: 3,4,5"
echo "omega: $[5+1]..$[5+$nat3]"
echo "width: $[5+$nat3+1]..$[5+2*$nat3]"
echo "shift: $[5+2*$nat3+1]..$[5+3*$nat3]"

#+1 because of the file name when grepping multiple files
colomega=$[1+5+$col]
colwidth=$[1+5+$nat3+$col]
colshift=$[1+5+2*$nat3+$col]

echo "Columns to get (+1): $colomega, $colwidth, $colshift"

grep "^ *$line "  $files /dev/null |\
	awk "{print \$1,\$$colomega,\$$colwidth,\$$colshift}"|\
	sed -e "s/^.*_T//" -e "s/_s/  /" -e "s/[.]out:/ /"|\
	sed -re "s/ +/\t/g"|\
	sort -n
