#!usr/bin/env bash

appidx=7
srcversion="v5"
errmode="ABS"
prediction="lorenzo"
qbins="b10"
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

app=${apps[$appidx]}
dims=${dim[$appidx]}
dtype="${dtypes[$appidx]}"
prefix=${prefixes[$appidx]}
ext=${exts[$appidx]}
inputdir="${prefix}/${app}"

srcdir="/home/alpoulo/repositories/compression/SZ3/release/examples"

outdir="/scratch1/alpoulo/adaptive_quantization/output/${app}/${errmode}/${prediction}/${qbins}"

cmp="${outdir}/cmp/${srcversion}"
mkdir -p $cmp

adaptivedir="${outdir}/${srcversion}"
mkdir -p $adaptivedir

entropydir="${outdir}/entropy/${srcversion}"
mkdir -p $entropydir

histogramdir="${outdir}/histogram/${srcversion}"
mkdir -p $histogramdir

pbsdir="${outdir}/pbs"
mkdir -p $pbsdir

################# FOR exaalt-copper data #################
#inputdir="${inputdir}/dataset1-5423x3137"
inputdir="${inputdir}/dataset2-83x1077290"
for infile in "$inputdir".*"${ext}"
##########################################################
#for infile in "$inputdir"/*"${ext}"
do
    fname="${infile##*/}"
    fname="${fname%.bin.f32}"
    fname="${fname%.f32}"
    fname="${fname%.d64}"
    fname="${fname%%_1*}"
    fname="${fname%.dat}"
    pbsout="${pbsdir}/${fname}_adaptive_${srcversion}.pbs"
    echo "#!/bin/bash" > $pbsout
    echo "#PBS -l select=1:ncpus=8:mem=${mem},walltime=${walltime}" >> $pbsout
    echo "#PBS -m a" >> $pbsout
    echo "#PBS -N ${fname}_adaptive" >> $pbsout
    echo "#PBS -e stderr_${fname}.txt" >> $pbsout
    printf "\n" >> $pbsout
    echo "cd ${srcdir}" >> $pbsout
    printf "\n" >> $pbsout

#    for eb in "${errbounds[@]}"
    for eb in {1..5}
    do
        for ab in {1..4}
        do
            let adapt_eb=2**${ab}
            adaptivefile="${adaptivedir}/${fname}_${errmode}_1e-${eb}_q${ab}_adaptive.out"
            quantentropyfile="${entropydir}/${fname}_${errmode}_1e-${eb}_q${ab}_quant.out"
            adaptentropyfile="${entropydir}/${fname}_${errmode}_1e-${eb}_q${ab}_adapt.out"
            quanthistfile="${histogramdir}/${fname}_${errmode}_1e-${eb}_q${ab}_quant_histogram.out"
            adapthistfile="${histogramdir}/${fname}_${errmode}_1e-${eb}_q${ab}_adapt_histogram.out"
            echo "export QUANTENTROPY=${quantentropyfile}" >> ${pbsout}
            echo "export ADAPTENTROPY=${adaptentropyfile}" >> ${pbsout}
            echo "export QUANTHISTOGRAM=${quanthistfile}" >> ${pbsout}
            echo "export ADAPTHISTOGRAM=${adapthistfile}" >> ${pbsout}
            printf "\n" >> $pbsout
            cmpfile="${cmp}/${fname}_1e-${eb}_q${ab}_adaptive.dat.sz.out"
            echo "./sz3_${srcversion} ${dtype} -i ${infile} -o ${cmpfile} ${dims} -c ${configfile} -M ${errmode} ${adapt_eb}e-${eb} -a -q ${ab} > ${adaptivefile}" >> ${pbsout}
            echo "echo QUANTENTROPY >> ${adaptivefile}" >> ${pbsout}
            echo "printProperty -f ${quantentropyfile} 3 >> ${adaptivefile}" >> ${pbsout}
            echo "echo ADAPTENTROPY >> ${adaptivefile}" >> ${pbsout}
            echo "printProperty -f ${adaptentropyfile} 3 >> ${adaptivefile}" >> ${pbsout}
            echo "rm ${quantentropyfile}" >> ${pbsout}
            echo "rm ${adaptentropyfile}" >> ${pbsout}
            echo "rm ${cmpfile}" >> ${pbsout}
        done
    done
    qsub $pbsout
done

