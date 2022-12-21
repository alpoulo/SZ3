#!usr/bin/env bash

version="v5"

qbins="b16"
errmode="ABS"
prediction="lorenzo"

appidx=0
file_name="vy"
fname=${file_name}
#fname=""
prefix1="/home/alpoulo/datasets"
prefix2="/zfs/fthpc/alpoulo/datasets"

apps=("HACC" "CESM2D" "exaalt-copper" "ISABEL" "CESM3D" "NYX" "MIRANDA")
dim=("-1 1073726487" "-2 3600 1800" "-2 5243 3137" "-3 500 500 100" "-3 3600 1800 26" "-3 512 512 512" "-3 384 384 256")
exts=(".f32" "_1_1800_3600.f32" ".f32.dat" ".bin.f32" "_1_26_1800_3600.f32" ".dat" ".d64")
dtypes=("-f" "-f" "-f" "-f" "-f" "-f" "-d")
prefixes=("${prefix2}" "${prefix1}" "${prefix1}" "${prefix1}" "${prefix2}" "${prefix2}" "${prefix1}")

app=${apps[$appidx]}
dims=${dim[$appidx]}
ext=${exts[$appidx]}
prefix=${prefixes[$appidx]}
infile="${prefix}/${app}/${file_name}${ext}"

baselinefile="${fname}_baseline"
outfile="${baselinefile}_${errmode}_${prediction}_${qbins}.out"

scratch="/scratch1/alpoulo/adaptive_quantization/output/${app}/tests/${errmode}/${prediction}/${qbins}"
mkdir -p $scratch

entropydir="${scratch}/entropy/baseline"
mkdir -p $entropydir

histdir="${scratch}/histogram/baseline"
mkdir -p $histdir

cmpdir="${scratch}/cmp/baseline"
mkdir -p $cmpdir

for i in {1..5}
do
    histfile="${histdir}/${baselinefile}_1e-${i}_hist.out"
    #histfile="${baselinefile}_1e-${i}_hist.out"
    export QUANTHISTOGRAM="${histfile}"
    quantentropyfile="${entropydir}/${baselinefile}_1e-${i}_quant.out"
    export QUANTENTROPY=${quantentropyfile}
    
    cmpfile="${cmpdir}/${baselinefile}_1e-${i}.dat.sz.out"
    
    #./sz3_$version -f -i $infile -o $cmpfile $dims -c sz3.config -M $errmode $i -a >> $outfile
    ./sz3_${version} ${dtype} -i ${infile} -o ${cmpfile} ${dims} -c sz3.config -M ${errmode} 1e-${i} -a >> ${outfile}
    printProperty -f $quantentropyfile 3 >> $outfile
    rm $cmpfile
    rm $quantentropyfile
    printf "\n" >> $outfile
done


