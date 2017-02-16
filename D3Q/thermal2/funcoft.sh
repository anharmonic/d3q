#!/bin/bash

error_and_quit(){
echo
echo "ERROR: $1"
echo
cat << EOF
This script collects the linewidth/shifts of a certain point q
among all the files give in input and print on stdout a list
of linewidths/shifts as function of T and sigma.
Syntax:

 $0 point mode file [file2 [file3...]] 

 point:	 index of the point to collect among the files
 mode:	 number of the phonon mode 1...3*nat
 file1..n: list of files produced by lw.x
 The output will contain 5 columns: 
 1. temperature
 2. smearing
 3. omega
 4. linewidth (gamma, HWHM)
 5. shifted frequency (if available in the file)
EOF
exit $2

}

line=$1
shift

col=$1
shift

if [ $# -lt 1 ]; then
 error_and_quit "Not enough arguments" 255
fi

echo "Getting mode $col, line $line" >&2

files="$@"
echo "Files: $files" >&2

for i in $files; do
  test -f "$i" || error_and_quit "File not found: '$i'" 256
done

nwords=$(grep "^ *1 " $1|wc -w)
if grep -q lineshift $1; then
  nat3=$[($nwords-5)/3]
else
  nat3=$[($nwords-5)/2]
fi  
nat=$[$nat3/3]

echo "$nwords, $nat3, $nat" >&2
#echo "idx: 1"
#echo "path: 2"
#echo "q_x..q_z: 3,4,5"
#echo "omega: $[5+1]..$[5+$nat3]"
#echo "width: $[5+$nat3+1]..$[5+2*$nat3]"
#echo "shift: $[5+2*$nat3+1]..$[5+3*$nat3]"

#+1 because of the file name when grepping multiple files
colomega=$[1+5+$col]
colwidth=$[1+5+$nat3+$col]
colshift=$[1+5+2*$nat3+$col]

echo "Columns to get (+1): $colomega, $colwidth, $colshift" >&2

grep "^ *$line "  $files /dev/null |\
	awk "{print \$1,\$$colomega,\$$colwidth,\$$colshift}"|\
	sed -e "s/^.*_T//" -e "s/_s/  /" -e "s/[.]out:/ /"|\
	sed -re "s/ +/\t/g"|\
	sort -k 1n -k 2n 


