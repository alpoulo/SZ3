#!usr/bin/env bash

appidx=0
srcversion="v5"
errmode="REL"
prediction="lorenzo"
qbins="b16"
walltime="4:00:00"
mem="64gb"

configfile="sz3_${qbins}.config"

prefix1="/home/alpoulo/datasets"
prefix2="/zfs/fthpc/alpoulo/datasets"

apps=("HACC" "CESM2D" "exaalt-copper" "ISABEL" "CESM3D" "NYX" "Miranda" "exaalt-copper")
dim=("-1 1073726487" "-2 3600 1800" "-2 5243 3137" "-3 500 500 100" "-3 3600 1800 26" "-3 512 512 512" "-3 384 384 256" "-2 83 1077290")
exts=(".f32" "_1_1800_3600.f32" ".f32.dat" ".bin.f32" "_1_26_1800_3600.f32" ".dat" ".d64" ".f32.dat")
dtypes=("-f" "-f" "-f" "-f" "-f" "-f" "-d" "-f")
prefixes=("${prefix2}" "${prefix1}" "${prefix1}" "${prefix1}" "${prefix2}" "${prefix2}" "${prefix1}" "${prefix1}")

app=${apps[$appidx]} dims=${dim[$appidx]}
dtype=${dtypes[$appidx]}
prefix=${prefixes[$appidx]}
ext=${exts[$appidx]}
inputdir="${prefix}/${app}"

srcdir="/home/alpoulo/repositories/compression/SZ3/release/examples"
outdir="/scratch1/alpoulo/adaptive_quantization/output/${app}/${errmode}/${prediction}/${qbins}"

cmp="${outdir}/cmp/baseline"
mkdir -p $cmp
baselinedir="${outdir}/baseline"
mkdir -p $baselinedir
entropydir="${outdir}/entropy/baseline"
mkdir -p $entropydir
histdir="${outdir}/histogram/baseline"
mkdir -p $histdir
pbsdir="${outdir}/pbs"
mkdir -p $pbsdir

#inputdir="$inputdir/dataset1-5423x3137"
#inputdir="${inputdir}/dataset2-83x1077290"
#for infile in "$inputdir".*"${ext}"
for infile in "$inputdir"/*"${ext}"
do
    fname="${infile##*/}"
    fname="${fname%.bin.f32}"
    fname="${fname%.f32}"
    fname="${fname%%_1*}"
    fname="${fname%.dat}"
    pbsout="${pbsdir}/${fname}_baseline.pbs"
    echo "#!/bin/bash" > $pbsout
    echo "#PBS -l select=1:ncpus=8:mem=${mem},walltime=${walltime}" >> $pbsout
    echo "#PBS -m a" >> $pbsout
    echo "#PBS -N ${fname}_baseline" >> $pbsout 
    echo "#PBS -e stderr_${fname}.txt" >> $pbsout
    printf "\n" >> $pbsout
    echo "cd ${srcdir}" >> $pbsout
    printf "\n" >> $pbsout
    
    for eb in {1..5}
    do
        baselinefile="${baselinedir}/${fname}_${errmode}_1e-${eb}_baseline.out"
        entropyfile="${entropydir}/${fname}_${errmode}_1e-${eb}_baseline.out"
        histfile="${histdir}/${fname}_${errmode}_1e-${eb}_baseline.out"
        echo "export QUANTENTROPY=${entropyfile}" >> ${pbsout}
        echo "export QUANTHISTOGRAM=${histfile}" >> ${pbsout}
        printf "\n" >> $pbsout
        cmpfile="${cmp}/${fname}_baseline_1e-${eb}.dat.sz.out"
        echo "./sz3_${srcversion} ${dtype} -i ${infile} -o ${cmpfile} ${dims} -c ${configfile} -M ${errmode} 1e-${eb} -a > ${baselinefile}" >> ${pbsout}
        echo "printProperty -f ${entropyfile} 3 >> ${baselinefile}" >> ${pbsout}
        echo "rm ${entropyfile}" >> ${pbsout}
        #echo "rm ${cmpfile}" >> ${pbsout}
    done
    qsub $pbsout
done

