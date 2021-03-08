#!/bin/bash

# Starting from dynamical matrix file dyn* and force constants in "fc", 
# create a set of files with the sum rules already applied to them
# It uses a modified matdyn.x that also calls write_epsilon_and_zeu when q=0
# to the fldyn file. Change matdyn.f90 like this:
#        if(iout_dyn.ne.0) then
#          call write_dyn_on_file(q(1,n),dyn,nat, iout_dyn)
#          if(sum(abs(q(:,n)))==0._dp) call  write_epsilon_and_zeu (zeu, epsil, nat, iout_dyn)
#        endif
# (e.g. add the 3rd line) and recompile.

input="dyn"
output="asr_dyn"
sum_rule="crystal"

rand=$RANDOM$RANDOM
header="header-${rand}"
postd="postd-${rand}"
fcf="fc-${rand}"
tmpdir="tmp"

mkdir -p ${tmpdir}

print_help(){
        echo "Take a series of dynamical matrix files, apply sum rule to them"
        echo "and save the resulting dynamical matrix to new files"
        echo
	echo "Syntax:"
	echo "  $0 [-i IN] [-o OUT] [-a ASR_TYPE]"
        echo
        echo "IN:       input fildyn file prefix (default: 'dyn')"
        echo "OUT:      output fildyn file prefix (default: 'asr_dyn')"
        echo "ASR_TYPE: type of sum rule to use (default: crystal)"
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

ls ${input}* >/dev/null 2>&1 || print_help 1

echo "Starting to apply sum rule '${sum_rule}' on:"
ls ${input}*
echo "files with sum rule applied will be saved as '${output}'"
echo

n=0
IFS="
"

#rm -f *${postd}*
#rm -f ref_*
#rm -f rot_*
#rm -f out.*
#rm -f ${header}
#rm -f fc.*


# prepare force constants from initial matrices
cat << eof |q2r.x > ${tmpdir}/out-${rand}.0
 &input
   fildyn  =  '${input}'
   zasr    =  '${sum_rule}', 
   flfrc   =  '${tmpdir}/${fcf}'
 /
eof

if [ ! -s "${tmpdir}/${fcf}" ]; then
 echo -e "\n\nError: Could not generate force constants file" >&2
 exit 1
fi

#prepare ${header} file
nhead=$(grep -B100 "Dynamical  Matrix in cartesian axes" ${input}1|wc -l)
head -n $[$nhead-2] ${input}1 > ${tmpdir}/${header}


# generate, one by one, the dynamical matrices files
# for the list of q-points in the dyn0 file
n=0
for q in $(tail -n +3 ${input}0); do

let n++
echo "== $n: $q =="

# generate the dynamical matrix (without ${header})
cat << eof | matdyn.x > ${tmpdir}/out-${rand}.${n}
&input
  flfrc='${tmpdir}/${fcf}'
  asr='${sum_rule}' 
  fldyn='${tmpdir}/${postd}${n}'
  flfrq='${tmpdir}/freq-${rand}'
  flvec='${tmpdir}/modes-${rand}'
/
1
${q}
eof

# Add ${header} 
cat ${tmpdir}/${header}    >  ${tmpdir}/ref_${postd}${n}
cat ${tmpdir}/${postd}${n} >> ${tmpdir}/ref_${postd}${n}

# Generate the star of q
q2qstar.x ${tmpdir}/ref_${postd}${n} ${tmpdir}/rot_${postd}${n} >> ${tmpdir}/out-${rand}.${n}

# Add effective charges for dyn1
if [ ${n} -eq 1 ] && grep -q "Dielectric" ${input}1; then
 grep -B10000 "^     Diagonalizing the dynamical matrix" ${tmpdir}/rot_${postd}1 |grep -v Diagonalizing > ${tmpdir}/tmp-${rand}
 grep -A10000 Dielectric ${tmpdir}/${postd}1 >> ${tmpdir}/tmp-${rand}
 grep -A10000 "^     Diagonalizing the dynamical matrix" ${tmpdir}/rot_${postd}1 >> ${tmpdir}/tmp-${rand}
 mv ${tmpdir}/tmp-${rand} ${tmpdir}/rot_${postd}1
 #exit 0
fi

mv ${tmpdir}/rot_${postd}${n} ${output}${n}

done

rm -f ${tmpdir}/*${rand}*
ls ${tmpdir}/* >/dev/null 2>&1 || rm -r ${tmpdir}

# Copy the list of points
cp -f ${input}0 ${output}0




