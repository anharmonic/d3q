#!/bin/bash
input="fc"
output="mat2R_fc"
sum_rule="no"

rand=$RANDOM$RANDOM
header="header-${rand}"
dyn="dyn-${rand}_"
q2r="q2r-${rand}"
tmpdir="tmp"

print_help(){
        echo "Take a force constant file from q2r (or import_phonopy.py) and"
        echo "convert it to d3_q2r.x format, doing a back and forward FFT"
        echo
	echo "Syntax:"
	echo "  $0 [-i FC.input] [-o FC.out] [-a ASR_TYPE]"
        echo
        echo "FC.in:  input force constants (normal q2r.x format) default: 'fc'"
        echo "FC.out: output force constants (optimized d3_q2r.x format) default: 'mat2R_fc'"
        echo "ASR_TYPE: sum rule method to apply default: no"
        exit $1
}

while getopts "h?a:i:o:" opt; do
    case "$opt" in
    h|\?)
        print_help 0
        ;;
    a)  sum_rule=$OPTARG
        ;;
    i)  input=$OPTARG
        ;;
    o)  output=$OPTARG
        ;;
    esac
done
shift $((OPTIND-1))


mkdir -p ${tmpdir}


nqline=$(grep -E "^ *[1-9]+[0-9]* +[1-9]+[0-9]* +[1-9]+[0-9]* *$" ${input}|head -n 1)

n1=$(echo $nqline|awk '{print $1}')
n2=$(echo $nqline|awk '{print $2}')
n3=$(echo $nqline|awk '{print $3}')

nn=$[n1*n2*n3]
echo "Number of files: $nn, grid:  $n1, $n2, $n3"

function n2q {
  q1=$(echo $1./$4.|/usr/bin/bc -l)
  q2=$(echo $2./$5.|/usr/bin/bc -l)
  q3=$(echo $3./$6.|/usr/bin/bc -l)
  echo $q1 $q2 $q3
}

#prepare ${header} file
grep -B100000 "^ *${n1}  *${n2}  *${n3} *$" ${input} > ${tmpdir}/${header}-0
nl=0

echo "Artificial dyn mat file " > ${tmpdir}/${header}
echo " " >> ${tmpdir}/${header}
  head -n1 "${tmpdir}/${header}-0" >> ${tmpdir}/${header}
# add a damn line between celldm and lattice vectors if ibrav=0
if grep -qE "^ +[1-9]+ +[1-9]+ +0 +([0-9.]+ *){6,6} *$" ${tmpdir}/${header}-0; then
  echo "lattice vectors " >> ${tmpdir}/${header} 
fi
  nl=$(wc -l ${tmpdir}/${header}-0|awk '{print $1}')
  tail -n $[$nl-1] ${tmpdir}/${header}-0 |head -n $[$nl-3]>> ${tmpdir}/${header}


cat << eof > ${tmpdir}/${q2r}
 &input
   fildyn  =  ''
   flfrc   =  '${output}'
 /
$n1 $n2 $n3
$nn
eof

# generate, one by one, the dynamical matrices files
# for the list of q-points in the dyn0 file
n=0
for i1 in $(seq 0 $[n1-1]); do
for i2 in $(seq 0 $[n2-1]); do
for i3 in $(seq 0 $[n3-1]); do
  let n++
  q=$(n2q $i1 $i2 $i3 $n1 $n2 $n3)
  echo $q
  echo ${tmpdir}/${dyn}${n} >> ${tmpdir}/${q2r}

cat << eof | matdyn.x > ${tmpdir}/out-${rand}.${n}
&input
  flfrc='${input}'
  asr='${sum_rule}' 
  fldyn='${tmpdir}/${dyn}${n}'
  flfrq='${tmpdir}/freq-${rand}'
  flvec='${tmpdir}/modes-${rand}'
  q_in_cryst_coord = .true.
/
1
${q}
eof

cat ${tmpdir}/${header}  >  ${tmpdir}/ref_${dyn}${n}
cat ${tmpdir}/${dyn}${n} >> ${tmpdir}/ref_${dyn}${n}
mv ${tmpdir}/ref_${dyn}${n}  ${tmpdir}/${dyn}${n}

done
done
done

d3_q2r.x < ${tmpdir}/${q2r} > ${tmpdir}/out-${rand}.d3_q2r

rm -f ${tmpdir}/*${rand}*
ls ${tmpdir}/* >/dev/null 2>&1 || rm -r ${tmpdir}

