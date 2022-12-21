#!usr/bin/env bash

prefix1="/home/alpoulo/datasets"
prefix2="/zfs/fthpc/alpoulo/datasets"

#apps=("HACC" "CESM2D" "exaalt-copper" "ISABEL" "CESM3D" "NYX" "MIRANDA")
apps=("NYX")
#exts=(".f32" "_1_1800_3600.f32" ".f32.dat" ".bin.f32" "_1_26_1800_3600.f32" ".dat" ".d64")
exts=(".dat")
prefixes=("${prefix2}" "${prefix1}" "${prefix1}" "${prefix1}" "${prefix2}" "${prefix2}" "${prefix1}")

#for infile in "$inputdir"/*.f32
for appidx in {0..0}
do
    app=${apps[$appidx]}
    ext=${exts[$appidx]}
    prefix=${prefixes[$appidx]}
    inputdir="${prefix}/${app}"

    outdir="/scratch1/alpoulo/adaptive_quantization/output/${app}"
    entropydir="${outdir}/original/entropy"
    mkdir -p $entropydir
    echo ${app}
    
    for infile in "$inputdir"/*"$ext"
    do
        fname="${infile##*/}"
        fname="${fname%.bin.f32}"
        fname="${fname%.f32}"
        fname="${fname%.d64}"
        fname="${fname%%_1*}"
        fname="${fname%.dat}"
        echo ${fname}
        outfile="${entropydir}/${fname}.out"
        printProperty -f $infile 3 > $outfile 
    done
done
